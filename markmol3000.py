from rdkit import Chem
from rdkit.Chem.rdchem import EditableMol
from rdkit.Chem.rdchem import Atom
from rdkit.Chem.Draw import ShowMol, MolsToImage #only for debugging
from copy import deepcopy
import os
import sys, getopt

#Default settings
debug = False


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


        markinchi = MarkInChI(
            core, rgroups=rgroups, final=True
            )

        return markinchi.get_markinchi()

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
        return core_mol, rgroups


    


class MarkInChI():
    def __init__(
            self, mol, rgroups={}, final=False
            ):
        self.mol = mol
        self.rgroups = deepcopy(rgroups)
        #Use a deepcopy to stop any changes to the ordering of R groups within
        #this fragment affecting its parent
        self.final = final


    def get_markinchi(self):
        
        Show(self.mol)
        self.mol, varattachs = self.get_varattachs(self.mol)
        
        
        

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
        self.mol = self.rgroup_pseudoatoms_to_xe(self.mol)

        # Canonicalise the atom indices of self.mol by converting to InChI and 
        # back again
        core_inchi, aux = Chem.MolToInchiAndAuxInfo(self.mol)
        self.mol = Chem.MolFromInchi(core_inchi)

        # Re-map any indices to the new canonical labels
        mapping = self.get_canonical_map(aux)


        for varattach in varattachs:
            varattach.generate_inchi(self.rgroups)
            varattach.remap_endpts(mapping)
        
        varattachs = self.sort_varattachs(varattachs)

        final_inchi = self.finalise_markinchi(core_inchi, varattachs)

        

        #print(final_inchi)
        
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
        """
        for atom in mol.GetAtoms():
            if atom.HasProp("_MolFileRLabel"):
                rlabel = atom.GetIntProp("_MolFileRLabel")
                atom.SetAtomicNum(54) #Change atom to Xe
                atom.SetIsotope(31+rlabel)
        return mol

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

    def finalise_markinchi(self, inchi, varattachs):
        # This function constructs the final MarkInChI, making sure everything
        # is formatted correctly etc.

        core_inchi = inchi.replace("Xe", "Zz")

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
        
        # Add Markush layer for each variable attachment

        for varattach in varattachs:
            markush_strings += varattach.get_final_inchi()

        if xe_atom_count == 1:
            markush_strings = markush_strings.replace("<M></M>","")

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
        # formatted differently - no InChI prefix
        if not self.final:
            final_inchi = final_inchi.replace("InChI=1S/", "")

        # If final string has Markush information, label it as a MarkInChI 
        if final_inchi.find("<M>") != 0:
            final_inchi = final_inchi.replace("InChI=1S/", "MarkInChI=1B/")

        return final_inchi

    def get_varattachs(self, mol):

        # Breaks the mol up into the core and the variable attachments
        # Returns the core and a list of the variable attachments that directly
        # join onto the core

        # Label each atom to track its original index once the mol is split
        for atom in mol.GetAtoms():
            atom.SetProp("molAtomMapNumber", str(atom.GetIdx() + 1))

        # Get the variable attachment information for this molecule
        # A variable attachment bond is represented as a bond between the
        # attachment and a placeholder atom with Z=0, and the bond has a
        # special property '_MolFileBondEndPts' which is a list of the possible
        # indices the attachment could be bonded to (the endpoints).
        # We find the bonds that have this property, and for those bonds get
        # the beginning and end atoms to find the placeholder atom, and also
        # parse the value of the property to get the index of the endpoints.
        # Store these in a dictionary {placeholder_idx:[endpt1, endpt2, ...], }
        attachments = {}
        for bond in mol.GetBonds():
            if bond.HasProp("_MolFileBondEndPts"):
                endpts_string = bond.GetProp("_MolFileBondEndPts")
                endpts_string = endpts_string[1:len(endpts_string)-1]
                endpts = endpts_string.split()
                endpts = endpts[1:]
                for i, endpt in enumerate(endpts):
                    endpts[i] = int(endpt)

                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                if begin_atom.GetAtomicNum() == 0:
                    placeholder_idx = begin_atom.GetIdx() + 1
                else:
                    placeholder_idx = end_atom.GetIdx() + 1                
                attachments[placeholder_idx] = endpts
        
        # Split the mol based on the fragments that are not covalently bonded to
        # each other (the variable attachments are not considered bonded to the
        # core in the RDKit Mol)
        # frag_indices stores the original indices of the atoms in each fragment
        frag_indices = []
        frags = Chem.GetMolFrags(
            mol, asMols=True, fragsMolAtomMapping = frag_indices
            )
        frags = list(zip(frags, frag_indices))

        all_varattachs = []

        # Work out which fragment is the core and which are variable attachments
        # by checking to see if each fragment contains any variable attachment
        # placeholder atoms. 
        for frag in frags:

            is_core = True

            for i in frag[1]:
                if i + 1 in attachments.keys():
                    is_core = False
                    varattach = VarAttach(frag[0], frag[1], attachments[i+1])
                    all_varattachs.append(varattach)
            
            if is_core:
                core = frag[0]
        
        # Check whether a variable attachment is nested in another variable 
        # attachment - if it is, recombine this attachment into the Mol for its
        # parent attachment - this will get dealt with later when the MarkInChI
        # for the parent attachment is generated. 
        for a in all_varattachs:
            endpts = a.get_endpts()

            for b in all_varattachs:
                indices = b.get_original_indices()

                if all(endpt in indices for endpt in endpts):
                    a.set_nested()
                    b.add_nested_attachment(a)

        # Make a list of the variable attachments that link directly to the core
        non_nested_varattachs = []
        for v in all_varattachs:
            if not v.is_nested():
                v.reform_structure()
                non_nested_varattachs.append(v)
        

        return core, non_nested_varattachs

    def sort_varattachs(self, varattachs):

        # Sort the variable attachments. They are sorted first by the sum of the
        # indices of the endpoints, and then alphabetically
        # (in the code we sort alphabetically first to break ties, then the
        # endpoint sum takes precedence)

        varattachs = sorted(varattachs, key=lambda v: v.get_inchi())
        varattachs = sorted(varattachs, key=lambda v: v.get_endpt_sum())

        return varattachs

class RGroup():

    def __init__(self, id):
        self.id = id
        self.components = []
        self.inchi = ""
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
        self.components.append(component)
        return None
    
    def add_xe(self):
        # For each component, 
        # adds pseudoatoms to the Mol structure to represent the linkage points
        # back to the parent structure
        # Also sets list of indices for these atoms. 
        for component in self.components:
            mol = component["mol"]
            attachments = component["attachments"]
            atoms = mol.GetAtoms()
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
                component["mol"] = mol

    def generate_component_inchis(self, rgroups):

        for component in self.components:
            mol = component["mol"]
            markinchi = MarkInChI(
                mol, rgroups=rgroups)
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

        
class VarAttach():

    def __init__(self, mol, original_indices, endpts):
        self.mol = mol
        self.original_indices = original_indices
        self.endpts = endpts
        self.nested = False
        self.inchi = ""
        self.ranker = ""       

        
    def get_mol(self):
        return self.mol
    
    def get_endpts(self):
        return self.endpts
    
    def get_original_indices(self):
        return self.original_indices
        
    def set_nested(self):
        self.nested = True

    def is_nested(self):
        return self.nested
    
    def get_inchi(self):
        return self.inchi

    def reform_structure(self):
        # For any nested variable attachments, reforms the variable attachment
        # bond, now referring to the new indices on this fragment
        # Also turn the linker back to the parent into a labelled Xe atom
        for bond in self.mol.GetBonds():
            if bond.HasProp("_MolFileBondEndPts"):
                
                endpts_string = bond.GetProp("_MolFileBondEndPts")
                endpts_string = endpts_string[1:len(endpts_string)-1]
                endpts = endpts_string.split()
                endpts = endpts[1:]
                for i, endpt in enumerate(endpts):
                    endpts[i] = int(endpt)

                if endpts == self.endpts:
                    bond.ClearProp("_MolFileBondEndPts")
                    for atom in (bond.GetBeginAtom(), bond.GetEndAtom()):
                        if atom.GetSymbol() == "*":
                            atom.SetAtomicNum(54)
                            atom.SetIsotope(31)

                else:
                    new_endpts = []
                    for atom in self.mol.GetAtoms():
                        if int(atom.GetProp("molAtomMapNumber")) in endpts:
                            new_endpts.append(str(atom.GetIdx() + 1))
                    new_string = "("
                    new_string += str(len(new_endpts)) + " "
                    for new_endpt in new_endpts:
                        new_string += str(new_endpt) + " "
                    new_string = new_string[:len(new_string)-1] + ")"
                    bond.SetProp("_MolFileBondEndPts", new_string)

    def add_nested_attachment(self, attachment):
        self.mol = Chem.CombineMols(self.mol, attachment.get_mol())
        
    def generate_inchi(self, rgroups):

        markinchi = MarkInChI(
            self.mol, rgroups
        )

        inchi = markinchi.get_markinchi()

        self.inchi = inchi
        return inchi

    def remap_endpts(self, mapping):

        new_endpts = []
        for endpt in self.endpts:
            new_endpt = mapping[endpt]
            new_endpts.append(new_endpt)
        
        self.endpts = sorted(new_endpts)

    def get_endpt_sum(self):

        # Used for canonically sorting variable attachments

        endpt_sum = 0

        for endpt in self.endpts:
            endpt_sum += endpt

        return endpt_sum

    def get_final_inchi(self):

        final_inchi = "<M>"

        for endpt in self.endpts:
            final_inchi += "%iH," % endpt
        
        final_inchi = final_inchi[:len(final_inchi) - 1] + "-"
        final_inchi += self.inchi
        final_inchi += "</M>"

        return final_inchi
    
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
  if debug:
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

        
    argv = sys.argv[1:]
    opts, args = getopt.getopt(argv, "i:d")
    for opt, arg in opts:
        if opt == "-d":
            debug = True
        elif opt == "-i":
            filename = arg
    
    if len(argv) == 0:
        filename = "molfiles\\structures_for_testing\\ext9.mol"
        debug = True

    filedir = os.path.join(os.getcwd(), filename)
    markinchi = MarkMol3000()
    markinchi.load_from_file(filedir)
    markinchi_string = markinchi.generate_markinchi()
    print(markinchi_string)

    