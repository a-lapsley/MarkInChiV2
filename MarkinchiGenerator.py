from rdkit import Chem
from rdkit.Chem.rdchem import EditableMol, Atom, Mol
from copy import deepcopy
from MarkinchiUtils import show, ctab_to_molblock
import os
import sys
import getopt

# Default settings
debug = False



class MarkinchiGenerator(object):

    def __init__(self) -> None:
        self.core = None  # Core of the molecule
        self.rgroups = None  # List of RGroups

    def generate_markinchi(self) -> str:
        # Generates the MarkInChI for the core and rgroups
        if self.core != None:
            markinchi = MarkInChI(
                self.core, rgroups=self.rgroups, final=True
            )
            return markinchi.get_markinchi()
        else:
            return "No or invalid file provided"

    def read_molfile_lines(self, molfile_lines: list) -> None:
        # Reads the contents of the lines to set:
        # self.core
        # self.rgroups

        writing_core = 0
        # status of writing the core molecule string
        # 0=core not yet written
        # 1=writing core
        # 2=core writing finished
        writing_rgroups = False

        core_ctab = ""
        rgroups = {}

        for line in molfile_lines:

            # Write the CTAB for the core of the molecule
            if writing_core == 0:
                if line.find("BEGIN CTAB") != -1:
                    writing_core = 1

            if writing_core == 1:
                if line.find("END CTAB") != -1:
                    writing_core = 2

                core_ctab += line

            # Get the CTABS for each component of each R group
            # Also store attachment information
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
                    rgroup_attachments = []  # [(attchpt, atom_no),..]

                # Read attachment information
                # Strip attachment information from the line as RDKit doesn't
                # like atoms with multiple attchpts
                if line.find("ATTCHPT") != -1:
                    new_line = ""
                    line_parts = line.split()
                    atom_idx = int(line_parts[2])
                    for line_part in line_parts:
                        if line_part.find("ATTCHPT") != -1:
                            attchpt = int(line_part.split("=")[1])
                            # -1 means attached to 1 and 2
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

        self.core = core_mol
        self.rgroups = rgroups

    def get_from_molfile(self, file_path: str) -> None:
        # Loads V3000 molfile with path 'file_path'
        # Splits into individual lines
        # If file is valid, parses the molfile lines to get the core and the 
        # R groups

        # Check file is a .mol file
        if file_path.endswith(".mol") != True:
            print("Invalid file format. Please use a .mol file")
            return None
        try:
            with open(file_path) as file:
                lines = file.readlines()

                # Check file uses the V3000 format
                if lines[3].find("V3000") == -1:
                    print("Invalid file provided. Please ensure the file is "
                          "using the V3000 molfile format."
                          )
                    return None
                else:
                    self.read_molfile_lines(lines)
                    return None
        except:
            print("Unable to open file \'%s\'" % file_path)
            return None

    def get_from_molblock(self, molblock: str) -> None:

        # Splits molblock into individual lines
        # If molblock is valid, parses the molfile lines to get the core and the 
        # R groups

        lines = molblock.split("\n")
        new_lines = []
        for line in lines:
            line = line + "\n"
            new_lines.append(line)
        
        try:
            self.read_molfile_lines(new_lines)
        except:
            pass

class MarkInChI():

    # Class for generating a MarkInChI string from a fragment

    # Core functions

    def __init__(self,
                 mol: Mol,
                 rgroups: dict = {},
                 final: bool = False) -> None:

        self.mol = mol  # The molecule to get the MarkInChI of
        self.rgroups = deepcopy(rgroups)
        # Use a deepcopy to stop any changes to the ordering of R groups within
        # this fragment affecting its parent
        self.final = final  # boolean, if False, MarkInChI is nested in another
        self.child_rgroups = {} 

    def get_markinchi(self) -> str:

        # This is the main algorithm for generating the MarkInChI for a fragment
        # or a molecule

        # Find any variable attachments, atom lists and child R groups on this
        # fragment
        self.mol, varattachs = self.get_varattachs(self.mol)

        self.mol, listatoms = self.get_listatoms(self.mol)

        child_rgroup_labels = self.get_child_rgroups(self.mol)

        # Generate MarkInChIs for each child R group, sort them, and relabel
        # the pseudoatoms in this molecule accordingly

        

        if len(child_rgroup_labels) != 0:
            for rlabel in child_rgroup_labels:
                rgroup = self.rgroups[rlabel]
                if not rgroup.is_processed():
                    rgroup.add_xe()
                    rgroup.generate_component_inchis(self.rgroups)
                    rgroup.sort_components()
                    rgroup.generate_inchi()
                    rgroup.set_processed(True)

                self.child_rgroups[rlabel] = rgroup

            self.child_rgroups, rgroup_mapping = self.sort_rgroups(
                self.child_rgroups)
            
            self.mol = self.relabel_pseudoatoms(
                self.mol, rgroup_mapping
            )
            
        # Generate MarkInChIs for each Variable Attachment, and sort them
        # alphabetically
        if len(varattachs) != 0:
            for varattach in varattachs:
                varattach.generate_inchi(self.rgroups)

            varattachs = sorted(varattachs, key=lambda v: v.get_inchi())

        # Sort the listatoms according to their atomic numbers (see function
        # for detailed explanation)
        listatoms = self.sort_listatoms_by_atomic_nums(listatoms)

        # Convert pseudoatoms to Xe and isotopically label according to R label
        self.mol = self.rgroup_pseudoatoms_to_xe(self.mol)

        # Canonicalize the molecule indices according to the InChI algorithm
        # This is to ensure any variable attachments and listatoms with multiple
        # equivalent indices are canonicalized
        # Represent variable attachment points with chains of Rn atoms
        # Represent list atoms with chains of Ne atoms
        self.mol, varattachs, listatoms = self.canonize_inchi(
            self.mol, varattachs, listatoms
        )

        # Canonicalise the atom indices of self.mol by converting to InChI and
        # back again

        core_inchi, aux = Chem.MolToInchiAndAuxInfo(self.mol)
        self.mol = Chem.MolFromInchi(core_inchi)

        # Re-map any indices to the new canonical labels
        mapping = self.get_canonical_map(aux)

        listatoms = self.remap_listatoms(listatoms, mapping)
        listatoms = self.get_atom_symbols(listatoms)

        # Deal with variable attachments
        for varattach in varattachs:
            varattach.generate_inchi(self.rgroups)
            varattach.remap_endpts(mapping)

        varattachs = self.sort_varattachs(varattachs)

        # Construct the final markinchi string
        final_inchi = self.finalise_markinchi(
            core_inchi, varattachs, listatoms
        )
        
        if debug:
            show(self.mol, indices=True)


        return final_inchi

    def finalise_markinchi(self,
                           inchi: str,
                           varattachs: list,
                           listatoms: list) -> str:
        # This function constructs the final MarkInChI, making sure everything
        # is formatted correctly etc.

        core_inchi = inchi.replace("Xe", "Zz")
        # If the fragment is just a single R group, we don't need the core bit
        # (But only if it's not final, otherwise it's not a proper MarkInChI)
        # Also check the Xe doesn't have a mass of 31, which would mean it is a
        # linker to parent, and so the fragment is an H atom, so we do need the
        # core
        mol_is_rgroup = False
        if len(self.mol.GetAtoms()) == 1 and not self.final:
            if self.mol.GetAtomWithIdx(0).GetAtomicNum() == 54:
                if self.mol.GetAtomWithIdx(0).GetIsotope() != 31:
                    core_inchi = ""
                    mol_is_rgroup = True

        # Obtain the InChI string before the isotope layer, the isotope layer
        # itself, and any additional layers after the isotope layer
        parts = core_inchi.split("/i")
        final_inchi = parts[0]
        markush_strings = ""

        if len(parts) == 2:
            sub_parts = parts[1].split("/")
            isotope_layer = sub_parts[0]

            if len(sub_parts) > 1:

                additional_layers = ""
                for sub_part in sub_parts[1:]:
                    additional_layers += "/" + sub_part
                
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
                xe_atom_count += 1

        for atom in self.mol.GetAtoms():
            if atom.GetAtomicNum() == 54:
                isotope = atom.GetIsotope()
                idx = atom.GetIdx() + 1
                if isotope == 31 and xe_atom_count > 1:
                    markush_strings += "<M></M>"
                    replace_string = "%i-100" % idx
                elif isotope == 31 and xe_atom_count == 1:
                    replace_string = "%i-100" % idx
                else:
                    rlabel = isotope - 31
                    string = self.child_rgroups[rlabel].get_final_inchi()
                    # If the mol is just an R group, strip the <M> markers as
                    # they are not needed
                    if mol_is_rgroup:
                        string = string[3:]
                        string = string[:len(string)-4]
                    markush_strings += string
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

        # Add Markush layer for the atom lists

        for listatom in listatoms:
            markush_strings += self.get_listatom_string(listatom)

        # Add Markush layer for each variable attachment

        for varattach in varattachs:
            markush_strings += varattach.get_final_inchi()

        # If there is still relevant isotopic information left, format it
        if isotope_layer != "":
            if isotope_layer.strip()[len(isotope_layer)-1] == ",":
                isotope_layer = isotope_layer[:len(isotope_layer)-1]
            isotope_layer = "/i" + isotope_layer

        # Add the isotope, additional, and Markush layers to the final string
        final_inchi = final_inchi + isotope_layer
        final_inchi += additional_layers
        final_inchi += markush_strings

        # MarkInChI strings that are sub parts of a greater MarkInChI should be
        # formatted differently - no InChI prefix
        if not self.final:
            final_inchi = final_inchi.replace("InChI=1S/", "")

        # If final string has Markush information, label it as a MarkInChI
        if final_inchi.find("<M>") != -1:
            final_inchi = final_inchi.replace("InChI=1S/", "MarkInChI=1B/")

        return final_inchi

    def canonize_rdkit(self,
                       mol: Mol,
                       varattachs, listatoms) -> tuple[Mol, list, list]:
        #
        #   THIS FUNCTION IS NOT CALLED
        #   IT IS HERE FOR TESTING PURPOSES ONLY
        #
        # Canonicalizes the indices of the molecule according to the RDKit
        # canonicalization algorithm.

        # For each Varattach, add an Og atom, isotopically labelled according
        # to the order of the Varattach in the list
        # We reverse the order of the list, so that the highest priority group
        # (the first one alphabetically) has the highest mass number, which will
        # give the connection points the lower indices after canonization

        for i, varattach in enumerate(varattachs[::-1]):
            for endpt in varattach.get_endpts():
                mol.UpdatePropertyCache()

                # Remove an implicit hydrogen
                core_atom = mol.GetAtomWithIdx(endpt - 1)
                hydrogens = core_atom.GetTotalNumHs()

                if hydrogens > 0:
                    core_atom.SetNumExplicitHs(hydrogens - 1)

                # Check whether the atom has already been marked by a labelled
                # Og. If it has, change the isotope to be the larger of the
                # 2, and otherwise, add a new Og.

                atom_already_marked = False

                for neighbor in core_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 118:
                        if neighbor.GetIsotope() < (i + 1):
                            neighbor.SetIsotope(i + 1)
                        atom_already_marked = True

                if not atom_already_marked:
                    edit_mol = EditableMol(mol)
                    atom = Atom(118)
                    atom.SetIsotope(i + 1)
                    atom_idx = edit_mol.AddAtom(atom)
                    edit_mol.AddBond(
                        endpt - 1, atom_idx, order=Chem.rdchem.BondType.SINGLE
                    )
                    mol = edit_mol.GetMol()

        # For each listatom, turn the atom into an isotopically labelled Tn
        # atom
        for i, listatom in enumerate(listatoms[::-1]):
            core_atom = mol.GetAtomWithIdx(listatom["idx"]-1)

            core_atom.SetAtomicNum(117)
            core_atom.SetIsotope(i + 1)

        # Generate the RDKit canonical ranking and renumber the atoms
        mol.UpdatePropertyCache()
        ranking = Chem.CanonicalRankAtoms(mol)
        new_order = zip(*sorted([(j, i) for i, j in enumerate(ranking)]))

        new_order = tuple(new_order)[1]
        new_order = new_order[::-1]

        mol = Chem.RenumberAtoms(mol, new_order)
        # Remap the varattach endpoints
        for varattach in varattachs:
            endpts = varattach.get_endpts()
            new_endpts = list(
                map(lambda i: new_order.index(i - 1) + 1, endpts))
            new_endpts = sorted(new_endpts)
            varattach.set_endpts(new_endpts)

        # Remap the listatom indices
        for listatom in listatoms:
            old_index = listatom["idx"]
            new_index = new_order.index(old_index - 1) + 1
            listatom["idx"] = new_index

        # Turn the placeholder atoms back into Hs or original atoms
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 118:
                atom.SetAtomicNum(1)
                atom.SetIsotope(0)
            if atom.GetAtomicNum() == 117:
                i = atom.GetIsotope() - 1
                atom.SetIsotope(0)
                atomic_nums = listatoms[::-1][i]["atomic_nums"]
                if 6 in atomic_nums:
                    atom.SetAtomicNum(6)
                elif 1 in atomic_nums:
                    atom.SetAtomicNum(atomic_nums[1])
                else:
                    atom.SetAtomicNum(atomic_nums[0])

        return mol, varattachs, listatoms

    def canonize_inchi(self, 
                       mol: Mol, 
                       varattachs: list, 
                       listatoms: list) -> tuple[Mol, list, list]:
        
        
        # Canonicalizes the indices of the molecule using the InChI algorithm,
        # labelling the variable attachments and listatoms to break any
        # symmetry to ensure the end result is canonical.

        # For each variable attachment, add a chain of Rn atoms to each endpoint
        # The chain length corresponds to the highest priority variable
        # attachment that has this endpoint
        #
        # Iterate through attachments in alphabetical order (they have already
        # been sorted) - this is the order of priority
        for i, varattach in enumerate(varattachs):
            for endpt in varattach.get_endpts():
                mol.UpdatePropertyCache(strict=False)
                mol = Chem.rdmolops.AddHs(mol)

                core_atom = mol.GetAtomWithIdx(endpt - 1)

                # If the endpoint has already been labelled by a higher priority
                # group, don't do anything further
                already_marked = False
                for neighbor in core_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 86:
                        already_marked = True

                if not already_marked:

                    idx_a = None
                    # Replace one of the hydrogens on this endpoint with a Rn
                    for neighbor in core_atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 1 and idx_a == None:
                            neighbor.SetAtomicNum(86)
                            neighbor.SetProp("firstInChain", "True")
                            idx_a = neighbor.GetIdx()

                    # Extend the chain of Rn atoms according to the priority of
                    # the group
                    # E.g. if there are 3 groups, highest priority group has a
                    # chain length of 3, so add 2 more Rn (we have already
                    # added the first one in the previous step)

                    for j in range(len(varattachs) - i - 1):
                        edit_mol = EditableMol(mol)
                        idx_b = edit_mol.AddAtom(Chem.Atom(86))
                        edit_mol.AddBond(
                            idx_a,
                            idx_b,
                            order=Chem.rdchem.BondType.SINGLE
                        )
                        mol = edit_mol.GetMol()
                        idx_a = idx_b

                mol = Chem.rdmolops.RemoveHs(mol, sanitize=False)

        # Add chains of Ne atoms for each listatom

        for i, listatom in enumerate(listatoms):

            mol.UpdatePropertyCache(strict=False)
            mol = Chem.rdmolops.AddHs(mol)
            core_atom = mol.GetAtomWithIdx(listatom["idx"] - 1)

            idx_a = None
            for neighbor in core_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 1 and idx_a == None:
                    neighbor.SetAtomicNum(10)
                    neighbor.SetProp("firstInChain", "True")
                    idx_a = neighbor.GetIdx()

            if idx_a == None:
                edit_mol = EditableMol(mol)
                new_atom = Chem.Atom(10)
                new_atom.SetProp("firstInChain", "True")
                idx_b = edit_mol.AddAtom(new_atom)
                edit_mol.AddBond(
                    core_atom.GetIdx(),
                    idx_b,
                    order=Chem.rdchem.BondType.SINGLE
                )
                mol = edit_mol.GetMol()
                idx_a = idx_b

            for j in range(len(listatoms) - i - 1):
                edit_mol = EditableMol(mol)
                new_atom = Chem.Atom(10)
                idx_b = edit_mol.AddAtom(new_atom)
                edit_mol.AddBond(
                    idx_a,
                    idx_b,
                    order=Chem.rdchem.BondType.SINGLE
                )
                mol = edit_mol.GetMol()
                idx_a = idx_b
            
            mol = Chem.rdmolops.RemoveHs(mol, sanitize=False)

        # Generate the InChI for the molecule to get the AuxInfo, and get the
        # mapping from the original indices to the canonical indices
        #
        # With this mapping, the element of index i in the list is the index
        # of the atom in the original molecule that gets mapped to the canonical
        # index i
        aux = Chem.MolToInchiAndAuxInfo(mol)[1]
        aux = aux.split("/N:")[1]
        aux = aux.split("/")[0]
        new_indices = []
        for idx in aux.split(","):
            new_indices.append(int(idx)-1)

        # Reverse the list - arbitrary but means the highest priority group now
        # has the lowest index
        if new_indices == [0]:
            # Case of a single Xe atom
            new_indices = new_indices
        elif new_indices == [1]:
            # Case of a single H atom
            new_indices.append(0)
        else:
            new_indices = new_indices[::-1]

        new_indices = tuple(new_indices)
        mol = Chem.RenumberAtoms(mol, new_indices)
        
        # Label the atoms to keep track of their new index when we remove the
        # Rn atoms
        for atom in mol.GetAtoms():
            atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))

            if atom.HasProp("firstInChain"):
                atom.SetAtomicNum(1)

        # Remove any Rn and Ne
        # Remove any H atoms that aren't labelled (i.e. don't remove if the
        # fragment is just an H atom)
        edit_mol = EditableMol(mol)
        edit_mol.BeginBatchEdit()
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 86:
                edit_mol.RemoveAtom(atom.GetIdx())

            if (atom.GetAtomicNum() == 1 and
                    not atom.HasProp("molAtomMapNumber")):
                edit_mol.RemoveAtom(atom.GetIdx())

            if atom.GetAtomicNum() == 10:
                edit_mol.RemoveAtom(atom.GetIdx())

        edit_mol.CommitBatchEdit()
        mol = edit_mol.GetMol()
        if len(mol.GetAtoms()) > 2:
            mol = Chem.rdmolops.RemoveHs(mol, sanitize=False)

        # Get a mapping from the new index from the canonicalization to the
        # index of the atom in the stripped molecule without the Rn atoms
        new_idx_to_stripped_idx_map = {}
        for atom in mol.GetAtoms():
            new_idx = int(atom.GetProp("molAtomMapNumber")) + 1
            stripped_idx = atom.GetIdx() + 1
            new_idx_to_stripped_idx_map[new_idx] = stripped_idx

        # Generate a final map from the original atom indices to the indices in
        # the stripped molecule
        final_map = []
        for i in range(len(mol.GetAtoms())):
            new_idx = new_indices.index(i)
            final_idx = new_idx_to_stripped_idx_map[new_idx + 1]
            final_map.append(final_idx)

        # Remap the varattach endpoints using this map
        for varattach in varattachs:
            endpts = varattach.get_endpts()
            new_endpts = list(map(lambda i: final_map[i-1], endpts))
            new_endpts = sorted(new_endpts)
            varattach.set_endpts(new_endpts)

        # Remap the listatoms
        for listatom in listatoms:
            old_idx = listatom["idx"]
            new_idx = final_map[old_idx - 1]
            listatom["idx"] = new_idx

        return mol, varattachs, listatoms

    def get_canonical_map(self, auxinfo: str) -> dict:
        # Generates a map between the original atom indices and the canonical
        # indices found by the InChI algorithm
        # These are extracted from the aux info generated by the InChI algorithm
        # Stored in a dictionary with entries original_label:new_label.

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

    # R Group handling

    def sort_rgroups(self, rgroups: dict) -> tuple[dict, dict]:
        # Sort and relabel the rgroups according to the MarkInChI for the group
        # Also return a mapping from the old R label to the new R label

        sorted_rgroups = sorted(
            rgroups.items(), key=lambda item: item[1].get_inchi()
        )

        rgroup_mapping = {}
        rgroups = {}

        for i, rgroup in enumerate(sorted_rgroups):
            rgroup = rgroup[1]
            old_id = rgroup.get_id()
            new_id = i + 1
            rgroup_mapping[old_id] = new_id
            rgroup.set_id(new_id)
            rgroups[new_id] = rgroup

        return rgroups, rgroup_mapping

    def get_child_rgroups(self, mol: Mol) -> list:
        # Gets a list of the R group labels that are immediately referenced
        # by this structure

        child_rgroups = []
        for atom in mol.GetAtoms():
            if atom.HasProp("_MolFileRLabel"):
                rlabel = atom.GetIntProp("_MolFileRLabel")
                if rlabel not in child_rgroups:
                    child_rgroups.append(rlabel)

        return child_rgroups

    def relabel_pseudoatoms(self, mol: Mol, rgroup_mapping: dict) -> Mol:
        # Relabel the R group pseudoatoms in mol according to the mapping

        for atom in mol.GetAtoms():
            if atom.HasProp("_MolFileRLabel"):
                old_label = atom.GetIntProp("_MolFileRLabel")
                new_label = rgroup_mapping[old_label]
                atom.SetIntProp("_MolFileRLabel", new_label)

        return mol

    def rgroup_pseudoatoms_to_xe(self, mol: Mol) -> Mol:
        """
        Changes the Rn pseudoatoms in mol to Xe atoms

        Isotopically labels the Xe atoms according to the R group number
        For Rn, Z = 31 + n

        Xe standard isotope weight: 131
        InChI supports isotopes from -100 to +100 relative to the 
        standard isotope weight, ie. 31 to 231

        Xe-31 we reserve for representing links back to parent fragments.

        So Xe-32 is R1, Xe-33 is R2, etc.
        """
        for atom in mol.GetAtoms():
            if atom.HasProp("_MolFileRLabel"):
                rlabel = atom.GetIntProp("_MolFileRLabel")
                atom.SetAtomicNum(54)  # Change atom to Xe
                atom.SetIsotope(31+rlabel)

        return mol

    def remove_pseudoatom_isotopes(self, 
                                   inchi: str, 
                                   rgroup_indices: list) -> str:

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

    # Variable attachment handling

    def get_varattachs(self, mol: Mol) -> tuple[Mol, list]:

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
            mol, asMols=True, fragsMolAtomMapping=frag_indices
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

    def sort_varattachs(self, varattachs: list) -> list:

        # Sort the variable attachments. They are sorted first by the sum of the
        # indices of the endpoints, and then alphabetically
        # (in the code we sort alphabetically first to break ties, then the
        # endpoint sum takes precedence)

        varattachs = sorted(varattachs, key=lambda v: v.get_inchi())
        varattachs = sorted(varattachs, key=lambda v: v.get_endpt_sum())

        return varattachs

    # Listatom handling

    def get_listatoms(self, mol: Mol) -> tuple[Mol, list]:

        # Gets the list of atomic numbers for the index of each listatom
        # Removes the listatom and replaces it with a normal atom according to
        # the following rules:
        # If carbon present in the list, replace with C
        # Otherwise, replace with lowest atomic number element other than H

        listatoms = []
        for atom in mol.GetAtoms():
            atom.ClearProp("molAtomMapNumber")
            if atom.HasQuery():
                smarts = atom.GetSmarts()
                if smarts != "*":

                    smarts = smarts.replace("[", "")
                    smarts = smarts.replace("]", "")
                    smarts = smarts.replace("#", "")
                    smarts_parts = smarts.split(",")

                    atomic_nums = []

                    for part in smarts_parts:
                        atomic_nums.append(int(part))

                    atomic_nums = sorted(atomic_nums)

                    idx = atom.GetIdx() + 1

                    edit_mol = EditableMol(mol)

                    if 6 in atomic_nums:
                        new_atom = Atom(6)
                    elif 1 in atomic_nums:
                        new_atom = Atom(atomic_nums[1])
                    else:
                        new_atom = Atom(atomic_nums[0])

                    # Arbitrary property so we can track this atom
                    new_atom.SetProp("isList", "1")
                    new_atom_idx = edit_mol.AddAtom(new_atom)

                    for bond in atom.GetBonds():

                        if bond.GetBeginAtomIdx() == idx - 1:
                            new_bond_end = bond.GetEndAtomIdx()
                        else:
                            new_bond_end = bond.GetBeginAtomIdx()

                        bond_type = bond.GetBondType()

                        edit_mol.AddBond(
                            new_atom_idx, new_bond_end, order=bond_type
                        )

                    edit_mol.RemoveAtom(idx - 1)

                    mol = edit_mol.GetMol()

                    new_indices = []

                    for i in range(len(mol.GetAtoms())):
                        if i == idx - 1:
                            new_indices.append(new_atom_idx - 1)
                        elif i < idx:
                            new_indices.append(i)
                        elif i >= idx:
                            new_indices.append(i - 1)

                    mol = Chem.RenumberAtoms(mol, tuple(new_indices))

                    listatom = {}
                    for atom in mol.GetAtoms():
                        if atom.HasProp("isList"):
                            listatom["idx"] = atom.GetIdx() + 1
                            atom.ClearProp("isList")

                    listatom["atomic_nums"] = atomic_nums
                    listatoms.append(listatom)

        # I don't know why Atom Map Numbers get added here but we need to remove
        # them as they cause issues with the canonicalization later

        for atom in mol.GetAtoms():
            atom.ClearProp("molAtomMapNumber")

        return mol, listatoms

    def remap_listatoms(self, listatoms: list, mapping: dict) -> list:
        # Remaps the listatoms to the new indices of the molecule, and sort
        # by this new index

        new_listatoms = []

        for listatom in listatoms:
            new_listatom = {}
            new_listatom["idx"] = mapping[listatom["idx"]]
            new_listatom["atomic_nums"] = listatom["atomic_nums"]
            new_listatoms.append(new_listatom)

        new_listatoms = sorted(new_listatoms, key=lambda item: item["idx"])

        return new_listatoms

    def get_atom_symbols(self, listatoms: list) -> list:

        # Adds a list of the atomic symbols for each atom in each list, in the
        # correct order (C, H, then ascending atomic)

        periodic_table = Chem.rdchem.GetPeriodicTable()

        for listatom in listatoms:
            atom_symbols = []
            if 6 in listatom["atomic_nums"]:
                atom_symbols.append("C")
            if 1 in listatom["atomic_nums"]:
                atom_symbols.append("H")
            for i in listatom["atomic_nums"]:
                if i not in (1, 6):
                    atom_symbols.append(periodic_table.GetElementSymbol(i))
            listatom["atom_symbols"] = atom_symbols

        return listatoms

    def get_listatom_string(self, listatom: dict) -> str:
        # Constructs the Markush string for a listatom

        markush_string = "<M>%i-" % listatom["idx"]

        for symbol in listatom["atom_symbols"]:
            markush_string += "%s!" % symbol
        markush_string = markush_string[:len(markush_string) - 1]
        markush_string += "</M>"

        return markush_string

    def sort_listatoms_by_atomic_nums(self, listatoms: list) -> list:
        # Sorts listatoms according to the atomic numbers it refers to.
        # We want highest priority to go to the list containing the highest
        # atomic number, and if two lists have the same highest, go to the next
        # highest, and so on
        # This can be easily done by summing over 2^(n-1) for n in each list,
        #  which represents each possible list by a unique integer

        for listatom in listatoms:
            sum = 0
            for n in listatom["atomic_nums"]:
                sum += 2**(n-1)
            listatom["sum"] = sum

        listatoms = sorted(listatoms, key=lambda item: item["sum"])
        return listatoms

class RGroup():

    # Class for representing an R group

    def __init__(self, id: int) -> None:
        self.id = id  # The R label of the R group
        self.components = []  # List of each possible component in the R group
        self.inchi = ""  # MarkInChI string for the R group
        self.child_rgroups = []  # R groups directly referenced by this one
        self.processed = False # To ensure it is only dealt with once

    def get_id(self) -> int:
        return self.id

    def set_id(self, id: int) -> None:
        self.id = id

    def is_processed(self) -> bool:
        return self.processed
    
    def set_processed(self, processed: bool) -> None:
        self.processed = processed

    def add_component(self, ctab: str, attachments: list) -> None:
        # Adds a component to the R group from a CTAB
        # attachments is not currently used but can be used in future for
        # multiply attached R groups

        molblock = ctab_to_molblock(ctab)
        mol = Chem.MolFromMolBlock(molblock)
        component = {}
        component["mol"] = mol
        component["attachments"] = attachments
        self.components.append(component)

    def add_xe(self) -> None:
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

    def generate_component_inchis(self, rgroups: dict) -> None:

        # Generates the MarkInChI for each component

        for component in self.components:
            mol = component["mol"]
            markinchi = MarkInChI(
                mol, rgroups=rgroups)
            inchi = markinchi.get_markinchi()
            component["inchi"] = inchi

    def sort_components(self) -> None:
        # Sorts components alphabetically by their MarkInChI string
        sorted_list = sorted(self.components, key=lambda item: item["inchi"])
        self.components = sorted_list

    def generate_inchi(self) -> str:
        # This MarkInChI is used used for sorting the R groups in this R group's
        # parent fragment

        inchi = "<M>"
        for component in self.components:
            inchi += component["inchi"]
            inchi += "!"
        inchi = inchi[:len(inchi)-1]
        inchi += "</M>"
        self.inchi = inchi
        return inchi

    def get_inchi(self) -> str:
        return self.inchi

    def generate_links(self, mol: Mol) -> None:
        # Get the original indices of the atoms in the core that this
        # pseudoatom is bonded to
        # Not currently used but may be used for multiply attached R groups
        linked_atoms = []
        for atom in mol.GetAtoms():
            if atom.HasProp("_MolFileRLabel"):
                rlabel = atom.GetIntProp("_MolFileRLabel")
                if rlabel == self.get_id():
                    bonded_atoms = atom.GetNeighbors()
                    # print(rlabel)

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

    def remap_attachments(self, core_mapping: dict) -> None:
        # relabel the links defined in attachments
        # original core index -> canonical core index defined by core_mapping
        # original component index -> canonical component index defined by
        # component["mapping"], which was generated when the InChI for the
        # component was generated
        for component in self.components:
            new_attachments = []
            # print(component["attachments"])
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

    def get_final_inchi(self) -> str:
        # Gets the final MarkInChI string for the R group
        inchi = self.inchi
        return inchi

class VarAttach():

    # Class for representing a variable attachment

    def __init__(self, 
                 mol: Mol, 
                 original_indices: list, 
                 endpts: list) -> None:
        
        self.mol = mol  # The Mol object for the attachment
        self.original_indices = original_indices  # What the atoms in this Mol
        # were labelled in the parent fragment
        self.endpts = endpts  # Indices on the parent fragment of the connection
        # points of this attachment
        self.nested = False  # If true, this is nested within another VarAttach
        self.inchi = ""  # MarkInChI string for the attachment
        self.ranker = ""  # String for canonicalising variable attachments

    def get_mol(self) -> Mol:
        return self.mol

    def get_endpts(self) -> list:
        return self.endpts

    def set_endpts(self, endpts: list) -> None:
        self.endpts = endpts

    def get_original_indices(self) -> list:
        return self.original_indices

    def set_nested(self) -> None:
        self.nested = True

    def is_nested(self) -> bool:
        return self.nested

    def get_inchi(self) -> str:
        return self.inchi

    def reform_structure(self) -> None:
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
                        if (atom.GetAtomicNum() == 0 and
                            not atom.HasProp("_MolFileRLabel")):
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

        # If variable attachment is just an R group, treat it differently
        atoms = self.mol.GetAtoms()
        if len(atoms) == 2:
            for atom in atoms:
                if atom.HasProp("_MolFileRLabel"):
                    new_mol = Chem.Mol()
                    new_mol = EditableMol(new_mol)
                    new_mol.AddAtom(atom)
                    new_mol = new_mol.GetMol()
                    self.mol = new_mol

    def add_nested_attachment(self, attachment: 'VarAttach') -> None:
        # Adds a nested attachment to the Mol for this VarAttach
        self.mol = Chem.CombineMols(self.mol, attachment.get_mol())

    def generate_inchi(self, rgroups: dict) -> str:

        # Generates the MarkInChI string for this VarAttach

        markinchi = MarkInChI(
            self.mol, rgroups
        )

        inchi = markinchi.get_markinchi()

        self.inchi = inchi
        return inchi

    def remap_endpts(self, mapping: dict) -> None:

        # Relabel this attachment's endpoints using the mapping, and sort them

        new_endpts = []
        for endpt in self.endpts:
            new_endpt = mapping[endpt]
            new_endpts.append(new_endpt)

        self.endpts = sorted(new_endpts)

    def get_endpt_sum(self) -> int:

        # Used for canonically sorting variable attachments

        endpt_sum = 0

        for endpt in self.endpts:
            endpt_sum += endpt

        return endpt_sum

    def get_final_inchi(self) -> str:

        # Constructs the final MarkInChI string for this attachment

        final_inchi = "<M>"

        for endpt in self.endpts:
            final_inchi += "%iH," % endpt

        final_inchi = final_inchi[:len(final_inchi) - 1] + "-"
        final_inchi += self.inchi
        final_inchi += "</M>"

        return final_inchi

# If this script is called directly
if __name__ == "__main__":

    # Parse script arguments
    argv = sys.argv[1:]

    opts, args = getopt.getopt(argv, "d")
    
    for opt, arg in opts:
        if opt == "-d":
            debug = True

    # This is just for testing purposes (e.g. when this script is run directly
    # from an IDE)
    if len(args) == 0:
            filename = "molfiles\\test36 .mol"
            debug = False
    else:
        filename = args[0]
    

    
    

    # Generate and print the MarkInChI
    filedir = os.path.join(os.getcwd(), filename)
    markinchi_generator = MarkinchiGenerator()
    markinchi_generator.get_from_molfile(filedir)
    markinchi_string = markinchi_generator.generate_markinchi()
    print(markinchi_string)
