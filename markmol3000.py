from rdkit import Chem
from rdkit.Chem.rdchem import EditableMol
from rdkit.Chem.rdchem import Atom
from rdkit.Chem.Draw import ShowMol, MolsToImage #only for debugging
from copy import deepcopy
import os

FILENAME = "molfiles\\test6a.mol"
DEBUG = True

class MarkMol3000(object):

    def __init__(self):
        self.molfile_lines = [] #Lines of v3000 .mol file that will be converted

    def generate_markinchi(self):
        #Generates the MarkInChI for the loaded v3000 .mol file

        #Parse the molfile to generate:
        #core - Mol for the core of the molecule
        #rgroups - dictionary of RGroup objects, key = R Group number
        #varattachs - list of VarAttach objects
        read_output = self.read_molfile_lines(self.molfile_lines)

        core = read_output[0]
        rgroups = read_output[1]
        varattachs = read_output[2]

        markinchi = MarkInChI(
            core, rgroups=rgroups, varattachs=varattachs, final=True
            )

        return markinchi.get_markinchi()

        for varattach in varattachs:
            varattach.set_endpts()
            varattach.set_atom_indices()
            core = varattach.split(core)

        for varattach in varattachs:
            varattach.generate_struct_inchi()
            varattach.set_ranker(mapping)

        varattachs = sorted(varattachs, key=lambda item: item.get_ranker())
        
        for varattach in varattachs:
            final_inchi += varattach.get_final_inchi(mapping)

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
        rgroups = {}
        varattchs = []

        for line in molfile_lines:

            #Write the CTAB for the core of the molecule
            if writing_core == 0:
                if line.find("BEGIN CTAB") != -1:
                    writing_core = 1

            if writing_core == 1:
                if line.find("END CTAB") != -1:
                    writing_core = 2

                if line.find("ENDPTS") != -1:
                    varattach = VarAttach(line)
                    varattchs.append(varattach)
                core_ctab += line

            #Get the CTABS for each component of each R group
            #Also store attachment information
            if line.find("BEGIN RGROUP") != -1:
                rgroup_number = int(line.split()[4])
                rgroup = RGroup(rgroup_number)
                writing_rgroups = True
                rgroup_ctab = ""

            if line.find("END RGROUP") != -1:
                rgroups[rgroup_number] = rgroup

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
        return core_mol, rgroups, varattchs


    


class MarkInChI():
    def __init__(
            self, mol, rgroups={}, varattachs=[], final=False, 
            parent_linker_indices = []
            ):
        self.mol = mol
        self.rgroups = deepcopy(rgroups)
        #Use a deepcopy to stop any changes to the ordering of R groups within
        #this fragment affecting its parent
        self.varattachs = varattachs
        self.final = final
        self.parent_linker_indices = parent_linker_indices

    def get_markinchi(self):
        
        child_rgroups = self.get_child_rgroups(self.mol)

        # Generate MarkInChIs for each child R group, sort them, and relabel
        # the pseudoatoms in this molecule accordingly

        if len(child_rgroups) != 0:
            for rlabel in child_rgroups:
                rgroup = self.rgroups[rlabel]
                rgroup.add_xe()
                rgroup.generate_component_inchis(self.rgroups)
                rgroup.sort_components()
                rgroup.generate_inchi()
            

            self.rgroups, rgroup_mapping = self.sort_rgroups(self.rgroups)

            self.mol = self.relabel_pseudoatoms(
                self.mol, rgroup_mapping
            )

            
        # Convert pseudoatoms to Xe and isotopically label according to R label
        # rlabel_map is a dictionary that gives the R group label of the 
        # pseudoatom at a given index
        self.mol, rlabel_map = self.rgroup_pseudoatoms_to_xe(self.mol)

        # Canonicalise the atom indices of self.mol by converting to InChI and 
        # back again
        core_inchi, aux = Chem.MolToInchiAndAuxInfo(self.mol)
        self.mol = Chem.MolFromInchi(core_inchi)

        # Re-map any lists of indices to the new canonical labels
        mapping = self.get_canonical_map(aux)

        rlabel_map = self.remap_rlabel_map(rlabel_map, mapping)

        self.parent_linker_indices = self.remap(
            self.parent_linker_indices, mapping
        )


        # --- Formatting the final string correctly ---
        core_inchi = core_inchi.replace("Xe", "Zz")

        # Obtain the InChI string before the isotope layer, the isotope layer 
        # itself, and any additional layers after the isotope layer
        parts = core_inchi.split("/i")
        final_inchi = parts[0]
        markush_strings = ""

        if len(parts) == 2:
            sub_parts = parts[1].split("/")
            isotope_layer = sub_parts[0]
            if len(sub_parts) == 2:
                additional_layers = sub_parts[1]
            else:
                additional_layers = ""
        else:
            isotope_layer = ""
            additional_layers = ""

        # Iterate through the pseudoatoms in the Mol and add the appropriate 
        # Markush information.
        # As the atom indices are their canonical indices, these parts will be
        # added in the correct order

        xe_atom_count = 0
        for atom in self.mol.GetAtoms():
            if atom.GetAtomicNum() == 54:
                isotope = atom.GetIsotope()
                idx = atom.GetIdx() + 1
                if isotope == 31:
                    markush_strings += "<M></M>"
                    replace_string = "%i-100" % idx
                else:
                    rlabel = isotope - 31
                    markush_strings += self.rgroups[rlabel].get_final_inchi()
                    replace_string = str(idx)
                    if isotope < 131:
                        replace_string += str(isotope-131)
                    elif isotope > 131:
                        replace_string += "+"
                        replace_string += str(isotope-131)
                # This seemed like the simplest way to cover the cases that the 
                # isotope is in the middle of the list or at the end 
                isotope_layer = isotope_layer.replace(replace_string + ",", "")
                isotope_layer = isotope_layer.replace(replace_string, "")

                xe_atom_count += 1

        # If there is still relevant isotopic information left, format it
        if isotope_layer != "":
            if isotope_layer.strip()[len(isotope_layer)-1] == ",":
                isotope_layer = isotope_layer[:len(isotope_layer)-1]
            isotope_layer = "/i" + isotope_layer

        # Add the isotope, additional, and Markush layers to the final string
        final_inchi = final_inchi + isotope_layer
        if additional_layers != "":
            final_inchi += "/"
            final_inchi += additional_layers
        final_inchi += markush_strings

        # MarkInChI strings that are sub parts of a greater MarkInChI should be
        # formatted differently - no InChI prefix, single atoms are just the 
        # symbol
        if not self.final:
            final_inchi = final_inchi.replace("InChI=1S/", "")
            if xe_atom_count == 0 and len(self.mol.GetAtoms()) == 1:
                atom = self.mol.GetAtomWithIdx(0)
                #This tidies up the isotope information for single atoms
                if isotope_layer != "":
                    isotope_layer = isotope_layer.replace("/i1","")
                    isotope_layer = isotope_layer.replace("/", "")
                final_inchi = atom.GetSymbol() + isotope_layer

        # If final string has Markush information, label it as a MarkInChI 
        if final_inchi.find("<M>") != 0:
            final_inchi = final_inchi.replace("InChI=1S/", "MarkInChI=1B/")

        #Show(self.mol, indices=True)
        return final_inchi


    def sort_rgroups(self, rgroups):
        rgroup_mapping = {}
        rgroups = dict(sorted(
            rgroups.items(), key=lambda item: item[1].get_inchi()
            ))
        for i in rgroups.keys():
            
            rgroup = rgroups[i]
            old_id = rgroup.get_id()
            new_id = i
            rgroup_mapping[old_id] = new_id
            rgroup.set_id(new_id)
        return rgroups, rgroup_mapping
            
    def get_child_rgroups(self, mol):
        #Gets a list of the R group labels that are immediately referenced
        #by this structure
        
        child_rgroups = []
        for atom in mol.GetAtoms():
            if atom.HasProp("_MolFileRLabel"):
                rlabel = atom.GetIntProp("_MolFileRLabel")
                if rlabel not in child_rgroups:
                    child_rgroups.append(rlabel)

        return child_rgroups
            
    def relabel_pseudoatoms(self, mol, rgroup_mapping):
        for atom in mol.GetAtoms():
            if atom.HasProp("_MolFileRLabel"):
                old_label = atom.GetIntProp("_MolFileRLabel")
                new_label = rgroup_mapping[old_label]
                atom.SetIntProp("_MolFileRLabel", new_label)
        return mol
    
    def rgroup_pseudoatoms_to_xe(self, mol):
        """
        Changes the Rn pseudoatoms in mol to Xe atoms

        Isotopically labels the Xe atoms according to the R group number
        For Rn, Z = 31 + n

        Xe standard isotope weight: 131
        InChI supports isotopes from -100 to +100 relative to the 
        standard isotope weight, ie. 31 to 231

        Xe-31 we reserve for representing links back to parent fragments.

        Also generate a map (atom_idx: rlabel).
        """
        index_to_rlabel_map = {}

        for atom in mol.GetAtoms():
            if atom.HasProp("_MolFileRLabel"):
                rlabel = atom.GetIntProp("_MolFileRLabel")
                atom.SetAtomicNum(54) #Change atom to Xe
                atom.SetIsotope(31+rlabel)
                idx = atom.GetIdx()
                index_to_rlabel_map[idx+1] = rlabel
        return mol, index_to_rlabel_map

    def get_canonical_map(self, auxinfo):
        #Generates a map between the original atom indices and the canonical
        #indices found by the InChI algorithm
        #These are extracted from the aux info generated by the InChI algorithm
        #Stored in a dictionary with entries original_label:new_label
        if auxinfo == None:
            return {1: 1}
        
        mapping_string = auxinfo.split("/")[2]
        mapping_string = mapping_string[2:]
        original_labels = mapping_string.split(",")
        canonical_map = {}
        for i, original_label in enumerate(original_labels):
            new_label = i + 1
            original_label = int(original_label)
            canonical_map[original_label] = new_label
        return canonical_map
    
    def remove_pseudoatom_isotopes(self, inchi, rgroup_indices):

        # Removes isotope labels for the pseudoatoms from inchi

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
                        if int(idx) not in rgroup_indices:
                            isotope_layer += fragment
                            isotope_layer += ","
                    else:
                        idx = fragment.split("-")[0]
                        if int(idx) not in rgroup_indices:
                            isotope_layer += fragment
                            isotope_layer += ","
                if isotope_layer != "i":
                    isotope_layer = isotope_layer[:len(isotope_layer)-1]
                    new_inchi += isotope_layer
                    new_inchi += "/"
        new_inchi = new_inchi[:len(new_inchi) - 1]
        return new_inchi

    def remap(self, index_array, mapping):
        #Simple function to remap indices in a list according to a mapping
        new_index_array = []
        for old_idx in index_array:
            new_idx = mapping[old_idx]
            new_index_array.append(new_idx)
        return new_index_array
    
    def remap_rlabel_map(self, rlabel_map, mapping):
        new_rlabel_map = {}
        for old_idx in rlabel_map.keys():
            new_idx = mapping[old_idx]
            new_rlabel_map[new_idx] = rlabel_map[old_idx]
        return new_rlabel_map
    
class RGroup():

    def __init__(self, id) -> None:
        self.id = id
        self.components = []
        self.inchi = ""
        self.index = None
        self.child_rgroups = []
     
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
        component["linker_indices"] = []
        self.components.append(component)
        return None
    
    def add_xe(self):
        # For each component, 
        # adds pseudoatoms to the Mol structure to represent the linkage points
        # back to the parent structure
        # Also sets list of indices for these atoms. 
        parent_linker_indices = []
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
                            xe_atom.SetIsotope(31)
                            xe_idx = edit_mol.AddAtom(xe_atom)
                            edit_mol.AddBond(
                                atom_idx,
                                xe_idx, 
                                order=Chem.rdchem.BondType.SINGLE
                                )
                            mol = edit_mol.GetMol()
                            parent_linker_indices.append(xe_idx + 1)
                component["linker_indices"] = parent_linker_indices
                component["mol"] = mol
                
    def generate_component_inchis(self, rgroups):

        for component in self.components:
            mol = component["mol"]
            linker_indices = component["linker_indices"]
            markinchi = MarkInChI(
                mol, rgroups=rgroups, parent_linker_indices=linker_indices)
            inchi = markinchi.get_markinchi()
            component["inchi"] = inchi

    def sort_components(self):
        sorted_list = sorted(self.components, key=lambda item: item["inchi"])
        self.components = sorted_list 

    def generate_inchi(self):
        inchi = "<M>"
        for component in self.components:
            inchi += component["inchi"]
            inchi += "!"
        inchi = inchi[:len(inchi)-1]
        inchi += "</M>"
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
        

class VarAttach():

    def __init__(self, line):
        self.line = line
        self.endpts = []
        self.atom_indices = []
        self.struct = None #Mol for the variable attachment
        self.struct_inchi = ""

    def set_endpts(self):
        #Gets the indices of the atoms on the core structure this variable
        #attachment can attach to and stores them in self.endpts
        endpts = []
        line = self.line.split("ENDPTS")[1]
        line = line.split("(")[1]
        line = line.split(")")[0]
        endpt_strings = line.split()
        endpt_strings = endpt_strings[1:]
        for endpt_string in endpt_strings:
            endpts.append(int(endpt_string))
        self.endpts = sorted(endpts)

    def set_atom_indices(self):
        #Gets the indices of the placeholder atom for this variable attahcment
        #and the atom on the variable attachment it is connected to. It is 
        #not clear at this point which index is which - we work that out later
        line_parts = self.line.split()
        atom1_idx = int(line_parts[4])
        atom2_idx = int(line_parts[5])
        self.atom_indices = [atom1_idx, atom2_idx]

    def split(self, core):
        #Splits the molecule core into the core of the molecule and the 
        #variable attachment. Returns the stripped core and stores the 
        #attachment as a Mol in self.attachment


        frag_indices = []
        frags = Chem.GetMolFrags(
            core, asMols=True, fragsMolAtomMapping = frag_indices
            )
        frags = list(zip(frags, frag_indices))
        for frag in frags:
            #Work out which fragment is the attachment by checking to see
            #whether it contains one of the indices stored in self.atom_indices
            is_attachment = False
            for atom_idx in self.atom_indices:
                if atom_idx - 1 in frag[1]:
                    is_attachment = True
            
            #If it is the attachment, remove the placeholder atom and add a
            #xenon atom placeholder instead - this will be turned into a Zz
            #later. 
            #This is maybe not the most elegant way of doing this but other 
            #methods seem to cause issues as it seems like RDKit is still 
            #trying to store variable attachment data for the molecule which
            #causes issues if we try to use ReplaceAtom() instead of doing it
            #manually. I'm not sure why that is but this method seems to work.
            if is_attachment:
                attachment = frag[0]
                for i, atom in enumerate(attachment.GetAtoms()):
                    if atom.GetAtomicNum() == 0:
                        neighbor_indices = []
                        for neighbor in atom.GetNeighbors():
                            neighbor_indices.append(neighbor.GetIdx())
                        
                        edit_mol = EditableMol(attachment)
                        xe_atom = Atom(54)
                        xe_idx = edit_mol.AddAtom(xe_atom)
                        for neighbor_idx in neighbor_indices:
                            edit_mol.AddBond(
                                neighbor_idx,
                                xe_idx, 
                                order=Chem.rdchem.BondType.SINGLE
                                )
                        edit_mol.RemoveAtom(i)
                        attachment = edit_mol.GetMol()
                        
                

            else:
                core = frag[0]
         
        self.struct = attachment
        return core
        
    def generate_struct_inchi(self):
        #Gets the InChI for the variable attachment structure
        atoms = self.struct.GetAtoms()
        if len(atoms) == 2:
            for atom in atoms:
                if atom.GetAtomicNum() != 54:
                    self.struct_inchi = atom.GetSymbol()
        else:
            self.struct_inchi = Chem.MolToInchi(self.struct)
        self.struct_inchi = self.struct_inchi.replace("InChI=1S/", "")
        self.struct_inchi = self.struct_inchi.replace("Xe", "Zz")
        
    def get_final_inchi(self, mapping):
        #Constructs the full MarkInChI label including the connections
        final_inchi = "<M>"
        canonical_endpts = []
        for endpt in self.endpts:
            canoncical_endpt = mapping[endpt]
            canonical_endpts.append(canoncical_endpt)

        canonical_endpts = sorted(canonical_endpts)
        for endpt in canonical_endpts:
            final_inchi += "%iH," % endpt
        final_inchi = final_inchi[:len(final_inchi)-1]
        final_inchi += "-"
        final_inchi += self.struct_inchi
        return final_inchi

    def set_ranker(self, mapping):
        #Used for determining order of variable attachments in the final
        #MarkInChI string
        sum = 0
        for endpt in self.endpts:
            sum += mapping[endpt]
        ranker = str(sum)
        ranker += self.struct_inchi
        self.ranker = ranker

    def get_ranker(self):
        return self.ranker



def ctab_to_molblock(ctab):
        #Adds necessary beginning and end to a CTAB to allow RDKit to read it
        #as a molblock
        molblock = "\n\n\n  0  0  0     0  0            999 V3000\n"
        molblock += ctab
        molblock += "M  END"
        return molblock





def Show(mols, subImgSize=(200, 200), title='RDKit Molecule',
            stayInFront=True, indices=False, **kwargs):
  """
  Generates a picture of molecule(s) and displays it in a Tkinter window.
  Only if in debug mode, otherwise does nothing. 

  This function is a copy of the ShowMol function from the RDKit source code
  that displays multiple Mols in the same Tkinter window.

  It is only used for debugging purposes.
  """
  if DEBUG:
    mols = deepcopy(mols)
    import tkinter

    from PIL import ImageTk

    if type(mols) is list:
        
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
    else:
        if indices:
            for atom in mols.GetAtoms():
                atom.SetProp("molAtomMapNumber", str(atom.GetIdx() + 1))
        ShowMol(mols)
    

if __name__ == "__main__":
    #If file run independently
    
    filedir = os.path.join(os.getcwd(), FILENAME)
    markinchi = MarkMol3000()
    markinchi.load_from_file(filedir)
    inchi = markinchi.generate_markinchi()
    print(inchi)

    