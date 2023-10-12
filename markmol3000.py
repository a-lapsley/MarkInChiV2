from rdkit import Chem
from rdkit.Chem.rdchem import EditableMol
from rdkit.Chem.rdchem import Atom
from rdkit.Chem.Draw import ShowMol, MolsToImage #only for debugging
import os

FILENAME = "molfiles\\exx12.1.mol"

class MarkMol3000(object):

    def __init__(self) -> None:
        self.molfile_lines = [] #Lines of v3000 .mol file that will be converted

    def generate_markinchi(self):
        #Generates the MarkInChI for the loaded v3000 .mol file

        core, rgroups = self.read_molfile_lines(self.molfile_lines)
        

        for rgroup in rgroups:
            rgroup.add_xe()
            rgroup.generate_component_inchis()
            rgroup.sort_components()
            rgroup.generate_inchi()

        rgroups, rgroup_mapping = self.sort_rgroups(rgroups)

        core = self.relabel_pseudoatoms(core, rgroup_mapping)

        for rgroup in rgroups:
            rgroup.generate_links(core)

        core, pseudoatom_indices = self.pseudoatoms_to_xe(core)
        core_inchi, core_aux = Chem.MolToInchiAndAuxInfo(core)
        #ShowMol(core)

        mapping = get_canonical_map(core_aux)
    
        core, rgroup_indices = self.set_atom_mapping(core, mapping)

        for rgroup in rgroups:
            rgroup.remap_attachments(mapping)
            rgroup.set_index_from_map(rgroup_indices)

        new_pseudoatom_indices = self.canonical_pseudoatom_indices(
            pseudoatom_indices, mapping
            )
        core_inchi = self.remove_pseudoatom_isotopes(
            core_inchi, new_pseudoatom_indices
            )
        core_inchi = core_inchi.replace("Xe", "Zz")

        rgroups = sorted(rgroups, key=lambda item: item.get_index())

        final_inchi = ""
        final_inchi += core_inchi
        

        for rgroup in rgroups:
            final_inchi += rgroup.get_final_inchi()

        if final_inchi.find("<M>") != -1:
            # Check it's not just a normal InChI
            final_inchi = final_inchi.replace("InChI=1S", "MarkInChI=1B")

        return final_inchi
        
        #ShowMols([core, relabelled_core])

    def load_from_file(self, file_path):
        # Loads V3000 molfile with path 'file_path'
        # Stores contents of the molfile in self.molfile_lines
        # Returns 1 if file is valid and is loaded succesfully
        # Returns -1 otherwise

        # Check file is a .mol file
        if file_path.endswith(".mol") != True:
            print("Invalid file format. Please use a .mol file")
            return -1
        try:
            with open(file_path) as file:
                lines = file.readlines()

                # Check file uses the V3000 format
                if lines[3].find("V3000") == -1:
                    print("Invalid file provided. Please ensure the file is "
                          "using the V3000 molfile format."
                          )
                    return -1
                else:
                    self.molfile_lines = lines 
                    return 1
        except:
            print("Unable to open file \'%s\'" % file_path)
            return -1

    def read_molfile_lines(self, molfile_lines):
        # Reads the contents of molfile_lines to generate:
        # - An RDKit Mol object for the core of the molecule, core_mol
        # - A list of RGroup objects for the R groups, rgroups

        writing_core = 0
        #status of writing the core molecule string 
        #0=core not yet written
        #1=writing core
        #2=core writing finished
        writing_rgroups = False

        core_ctab = ""
        rgroups = []

        for line in molfile_lines:

            #Write the CTAB for the core of the molecule
            if writing_core == 0:
                if line.find("BEGIN CTAB") != -1:
                    writing_core = 1

            if writing_core == 1:
                if line.find("END CTAB") != -1:
                    writing_core = 2

                core_ctab += line

            #Get the CTABS for each component of each R group
            #Also store attachment information
            if line.find("BEGIN RGROUP") != -1:
                rgroup_number = int(line.split()[4])
                rgroup = RGroup(rgroup_number)
                writing_rgroups = True
                rgroup_ctab = ""

            if line.find("END RGROUP") != -1:
                rgroups.append(rgroup)

            if writing_rgroups:

                if line.find("BEGIN CTAB") != -1:
                    rgroup_ctab = ""
                    rgroup_attachments = []#[(attchpt, atom_no),..]

                #Read attachment information
                #Strip attachment information from the line as RDKit doesn't
                #like atoms with multiple attchpts
                if line.find("ATTCHPT") != -1:
                    new_line = ""
                    line_parts = line.split()
                    atom_idx = int(line_parts[2])
                    for line_part in line_parts:
                        if line_part.find("ATTCHPT") != -1:
                            attchpt = int(line_part.split("=")[1])
                            #-1 means attached to 1 and 2
                            if attchpt == -1:
                                rgroup_attachments.append(
                                    (1, atom_idx)
                                )
                                rgroup_attachments.append(
                                    (2, atom_idx)
                                )
                                
                            else:
                                rgroup_attachments.append(
                                    (attchpt, atom_idx)
                                )

                        else:
                            if line_part == "M":
                                line_part += " "
                            new_line += line_part
                            new_line += " "
        
                    line = new_line + "\n"

                rgroup_ctab += line

                if line.find("END CTAB") != -1:
                    rgroup.add_component(rgroup_ctab, rgroup_attachments)

                

                    

        core_molblock = ctab_to_molblock(core_ctab)
        core_mol = Chem.MolFromMolBlock(core_molblock)
        
        return core_mol, rgroups
    
    def pseudoatoms_to_xe(self, mol):
        """
        Changes the Rn pseudoatoms in mol to Xe atoms

        Isotopically labels the Xe atoms according to the R group number
        For Rn, Z = 232 - n

        Xe standard isotope weight: 131
        InChI supports isotopes from -100 to +100 relative to the 
        standard isotope weight, ie. 31 to 231
        We want R1 to have the highest priority, so highest isotopic 
        weight, ie. 231, all the way down to a possible R201 with an 
        isotopic weight of 31

        Also generate a list of the indices of these atoms.
        """
        pseudoatom_indices = []
        for atom in mol.GetAtoms():
            if atom.HasProp("_MolFileRLabel"):
                rlabel = atom.GetIntProp("_MolFileRLabel")
                atom.SetAtomicNum(54) #Change atom to Xe
                atom.SetIsotope(232-rlabel)
                idx = atom.GetIdx()
                pseudoatom_indices.append(idx + 1)
        return mol, pseudoatom_indices
        
    def set_atom_mapping(self, mol, mapping):
        # Sets the atom map number for each atom in mol to the corresponding
        # canonical index according to mapping
        # This allows finding atoms based on their canonical InChI index
        # Also generates a map between the number of an R group and its index
        # (rlabel, index)
        rgroup_indices = {}
        for atom in mol.GetAtoms():
            original_idx = atom.GetIdx() + 1
            canonical_idx = mapping[original_idx]
            atom.SetAtomMapNum(canonical_idx)
            if atom.HasProp("_MolFileRLabel"):
                rlabel = atom.GetIntProp("_MolFileRLabel")
                rgroup_indices[rlabel] = canonical_idx

        return mol, rgroup_indices

    def canonical_pseudoatom_indices(self, pseudoatom_indices, mapping):
        # Returns list of the canonical indices for the pseudoatoms
        new_indices = []
        for idx in pseudoatom_indices:
            new_indices.append(mapping[idx])
        return new_indices

    def remove_pseudoatom_isotopes(self, inchi, indices):
        # Removes isotope labels for the pseudoatoms from inchi
        # indices gives the canonical index of the pseudoatoms

        inchi_parts = inchi.split("/")
        new_inchi = ""
        for inchi_part in inchi_parts:
            if inchi_part[0] != "i":
                new_inchi += inchi_part
                new_inchi += "/"
            else:
                isotope_layer = "i"
                fragments = inchi_part.split("i")[1].split(",")
                for fragment in fragments:
                    if fragment.find("+") != -1:
                        idx = fragment.split("+")[0]
                        if int(idx) not in indices:
                            isotope_layer += fragment
                            isotope_layer += ","
                    else:
                        idx = fragment.split("-")[0]
                        if int(idx) not in indices:
                            isotope_layer += fragment
                            isotope_layer += ","
                if isotope_layer != "i":
                    isotope_layer = isotope_layer[:len(isotope_layer)-1]
                    new_inchi += isotope_layer
                    new_inchi += "/"
        new_inchi = new_inchi[:len(new_inchi) - 1]
        return new_inchi

    def sort_rgroups(self, rgroups):
        rgroup_mapping = {}
        rgroups = sorted(rgroups, key=lambda item: item.get_inchi())
        for i, rgroup in enumerate(rgroups):
            old_id = rgroup.get_id()
            new_id = i + 1
            rgroup_mapping[old_id] = new_id
            rgroup.set_id(new_id)
        return rgroups, rgroup_mapping

    def relabel_pseudoatoms(self, mol, mapping):
        for atom in mol.GetAtoms():
            if atom.HasProp("_MolFileRLabel"):
                old_label = atom.GetIntProp("_MolFileRLabel")
                new_label = mapping[old_label]
                atom.SetIntProp("_MolFileRLabel", new_label)
        return mol
    

class RGroup():

    def __init__(self, id) -> None:
        self.id = id
        self.components = []
        self.inchi = ""
        self.index = None
     
    def get_id(self):
        return self.id

    def set_id(self, id):
        self.id = id
        return None

    def add_component(self, ctab, attachments):
        molblock = ctab_to_molblock(ctab)
        mol = Chem.MolFromMolBlock(molblock)
        component = {}
        component["mol"] = mol
        component["attachments"] = attachments
        self.components.append(component)
        return None
    
    def add_xe(self):
        for component in self.components:
            mol = component["mol"]
            attachments = component["attachments"]
            atoms = mol.GetAtoms()
            if len(atoms) > 1:
                # Single atom R groups are treated differently
                for atom in atoms:
                    atom_idx = atom.GetIdx()
                    # RDKit labels start at 0 but molfile labels start at 1
                    n = atom_idx + 1
                    for attachment in attachments:
                        if attachment[1] == n:
                            edit_mol = EditableMol(mol)
                            xe_atom = Atom(54)
                            xe_idx = edit_mol.AddAtom(xe_atom)
                            edit_mol.AddBond(
                                atom_idx,
                                xe_idx, 
                                order=Chem.rdchem.BondType.SINGLE
                                )
                            mol = edit_mol.GetMol()
                component["mol"] = mol
                
    def generate_component_inchis(self):
        for component in self.components:
            mol = component["mol"]
            if len(mol.GetAtoms()) == 1:
                atom = mol.GetAtomWithIdx(0)
                inchi = atom.GetSymbol()
                component["inchi"] = inchi
                component["mapping"] = None
            else:
                inchi, aux = Chem.MolToInchiAndAuxInfo(mol)
                inchi = inchi.replace("InChI=1S/","")
                inchi = inchi.replace("Xe", "Zz")
                component["inchi"] = inchi
                mapping = get_canonical_map(aux)
                component["mapping"] = mapping
                #print(inchi)
                #print(mapping)

    def sort_components(self):
        sorted_list = sorted(self.components, key=lambda item: item["inchi"])
        self.components = sorted_list 

    def generate_inchi(self):
        inchi = "<M>"
        for component in self.components:
            inchi += component["inchi"]
            inchi += "!"
        inchi = inchi[:len(inchi)-1]
        self.inchi = inchi
        return inchi

    def get_inchi(self):
        return self.inchi

    def generate_links(self, mol):
        # Get the original indices of the atoms in the core that this
        # pseudoatom is bonded to
        linked_atoms = []
        for atom in mol.GetAtoms():
            if atom.HasProp("_MolFileRLabel"):
                rlabel = atom.GetIntProp("_MolFileRLabel")
                if rlabel == self.get_id():
                    bonded_atoms = atom.GetNeighbors()
                    #print(rlabel)

        for atom in bonded_atoms:
            linked_atoms.append(atom.GetIdx() + 1)
            linked_atoms = sorted(linked_atoms)

        for component in self.components:
            # convert ranked atom indices from the v3000 file format to the 
            # absolute atom indices
            new_attachments = [] 
            for attachment in component["attachments"]:
                # attachments are of the form (attchpt, atom_idx)
                # Difficult to explain in words but for example:
                # If component is bonded to atom numbers 1, 8, 9 on the core 
                # structure, 
                # attchpt 1 = atom 1
                # attchpt 2 = atom 8
                # attchpt 3 = atom 9
                # Represents position in ordered list of atom numbers
                i = attachment[0] - 1
                core_atom = linked_atoms[i]
                component_atom = attachment[1]
                new_attachments.append((core_atom, component_atom))
            component["attachments"] = new_attachments
                   
    def remap_attachments(self, core_mapping):
        # relabel the links defined in attachments
        # original core index -> canonical core index defined by core_mapping
        # original component index -> canonical component index defined by
        # component["mapping"], which was generated when the InChI for the 
        # component was generated
        for component in self.components:
            new_attachments = []
            #print(component["attachments"])
            for attachment in component["attachments"]:
                old_core_idx = attachment[0]
                new_core_idx = core_mapping[old_core_idx]

                old_component_idx = attachment[1]
                if component["mapping"] != None:
                    new_component_idx = component["mapping"][old_component_idx]
                else: 
                    new_component_idx = old_component_idx
                new_attachment = (new_core_idx, new_component_idx)
                new_attachments.append(new_attachment)
            component["attachments"] = new_attachments

    def get_final_inchi(self):
        inchi = self.inchi
        return inchi
    
    def set_index_from_map(self, rgroup_mapping):
        index = rgroup_mapping[self.id]
        self.index = index
    
    def get_index(self):
        return self.index
        

def ctab_to_molblock(ctab):
        #Adds necessary beginning and end to a CTAB to allow RDKit to read it
        #as a molblock
        molblock = "\n\n\n  0  0  0     0  0            999 V3000\n"
        molblock += ctab
        molblock += "M  END"
        return molblock

def get_canonical_map(auxinfo):
    #Generates a map between the original atom indices and the canonical
    #indices found by the InChI algorithm
    #These are extracted from the aux info generated by the InChI algorithm
    #Stored in a dictionary with entries original_label:new_label
    mapping_string = auxinfo.split("/")[2]
    mapping_string = mapping_string[2:]
    original_labels = mapping_string.split(",")
    canonical_map = {}
    for i, original_label in enumerate(original_labels):
        new_label = i + 1
        original_label = int(original_label)
        canonical_map[original_label] = new_label
    return canonical_map

def ShowMols(mols, subImgSize=(200, 200), title='RDKit Molecule',
            stayInFront=True, **kwargs):
  """
  Generates a picture of molecules and displays it in a Tkinter window.

  This function is a copy of the ShowMol function from the RDKit source code
  that displays multiple Mols in the same Tkinter window.

  It is only used for debugging purposes.
  """
  import tkinter

  from PIL import ImageTk

  img = MolsToImage(mols, subImgSize, **kwargs)

  tkRoot = tkinter.Tk()
  tkRoot.title(title)
  tkPI = ImageTk.PhotoImage(img)
  tkLabel = tkinter.Label(tkRoot, image=tkPI)
  tkLabel.place(x=0, y=0, width=img.size[0], height=img.size[1])
  tkRoot.geometry('%dx%d' % (img.size))
  tkRoot.lift()
  if stayInFront:
    tkRoot.attributes('-topmost', True)
  tkRoot.mainloop()

if __name__ == "__main__":
    #If file run independently
    
    filedir = os.path.join(os.getcwd(), FILENAME)
    markinchi = MarkMol3000()
    markinchi.load_from_file(filedir)
    inchi = markinchi.generate_markinchi()
    print(inchi)

    