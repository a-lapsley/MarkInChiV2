from rdkit import Chem
from rdkit.Chem.rdchem import EditableMol, Mol
from rdkit.Chem.MolStandardize import rdMolStandardize
from MarkinchiUtils import show

import os
import sys
import getopt

debug = False

class MarkinchiParser(object):

    # Object for converting a MarkInChI string into an RDKit molecule and R
    # groups

    def __init__(self, markinchi: str = "") -> None:
        self.markinchi = markinchi
        self.rgroups = []
        self.rgroups_dict = {}
        self.nested_rgroups = []
        self.core_mol = None

    def set_markinchi(self, markinchi: str) -> None:
        self.markinchi = markinchi

    def parse_markinchi(self) -> tuple[Mol, dict]:
        
        # Main script for parsing the MarkInChI

        parts = self.markinchi.split("<")

        core_mol = self.get_core(parts[0])
        
        markush_strings = self.get_markush_strings(parts[1:])
        string_lists = self.sort_markush_strings(markush_strings)
        markush_stereo = string_lists[0]
        rgroup_strings = string_lists[1]
        varattach_strings = string_lists[2]
        listatom_strings = string_lists[3]

        # Re-add Markush stereochemistry

        if markush_stereo != "":

            core_mol = self.add_markush_stereo(
                core_mol,
                markush_stereo,
                rgroup_strings,
                varattach_strings,
                listatom_strings
            )
        



        core_mol, rgroup_strings = self.xe_to_rgroups(core_mol, rgroup_strings)

        for rgroup_string in rgroup_strings:
            self.add_rgroup_from_string(rgroup_string)

        for varattach_string in varattach_strings:
            core_mol = self.add_varattach(core_mol, varattach_string)

        for rgroup, nested_rgroups in zip(self.rgroups, self.nested_rgroups):
            for component in rgroup:
                component = self.update_rlabels(component)
            
            nested_rgroups = self.update_nested_rlabels(nested_rgroups)

            self.rgroups += nested_rgroups

        for i, rgroup in enumerate(self.rgroups):
            self.rgroups_dict[i + 1] = rgroup


        for listatom_string in listatom_strings:
            core_mol = self.add_listatom(core_mol, listatom_string)
        
        if debug:
            show(core_mol, indices=True)

        self.core_mol = core_mol

        return self.core_mol, self.rgroups_dict

    def get_core(self, markinchi_part: str) -> Mol:

        # Generates the core of the molecule from the first part of the
        # MarkInChI

        if markinchi_part.find("MarkInChI=1B") != -1:
            inchi = markinchi_part.replace("MarkInChI=1B","InChI=1S")
        elif markinchi_part.find("InChI=1S") != -1:
            inchi = markinchi_part
        else:
            inchi = "InChI=1S/" + markinchi_part

        inchi = inchi.replace("Zz","Xe")
        
        core_mol = Chem.MolFromInchi(inchi)

        return core_mol

    def add_markush_stereo(
            self, 
            mol: Mol, 
            markush_stereo: str, 
            rgroup_strings: list, 
            varattach_strings: list,
            listatom_strings: list) -> Mol:
        

        xe_counter = 0
        
        mol.UpdatePropertyCache(strict=False)

        edit_mol = EditableMol(mol)

        rgroup_count = 0

        for rgroup_string in rgroup_strings:
            if rgroup_string != "":
                rgroup_count += 1
        
        dummy_atom_indices = []

        if rgroup_strings != []:
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 54:
                    if rgroup_strings[xe_counter] != "":
                        idx_a = atom.GetIdx()
                        for i in range(rgroup_count - xe_counter):
                            new_atom = Chem.Atom(54)
                            idx_b = edit_mol.AddAtom(new_atom)
                            edit_mol.AddBond(
                                idx_a,
                                idx_b,
                                order=Chem.rdchem.BondType.SINGLE
                            )
                            idx_a = idx_b
                            dummy_atom_indices.append(idx_b)
                        xe_counter += 1

        mol = edit_mol.GetMol()

        varattach_count = len(varattach_strings)

        for i, varattach_string in enumerate(varattach_strings):
            endpts_str = varattach_string.split("-", 1)[0]
        
            endpts = []
            for endpt in endpts_str.split(","):
                endpt = endpt.replace("H","")
                endpts.append(int(endpt))

            for endpt in endpts:
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
                    for j in range(varattach_count - i - 1):
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

        first_in_chain_indices = []    
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() in (86, 10):
                if atom.HasProp("firstInChain"):
                    first_in_chain_indices.append(atom.GetIdx())
                else:
                    dummy_atom_indices.append(atom.GetIdx())
        
        new_stereo = {
            "main": {"b":"", "t": "", "m": "", "s": ""},
            "isotope": {"b":"", "t": "", "m": "", "s": ""}
        }

        markush_stereo_parts = markush_stereo.split("/")

        layer = "main"
        for part in markush_stereo_parts:
            key = part[0]
            if key == "i":
                layer = "isotope"
            if key in ("b", "t", "m", "s"):
                new_stereo[layer][key] = part

        inchi, aux = Chem.MolToInchiAndAuxInfo(mol)
        aux = aux.split("/N:")[1]
        aux = aux.split("/")[0]
        
        new_indices = []
        for idx in aux.split(","):
            new_indices.append(int(idx)-1)
        
        reverse_mapping = []
        for i in range(len(new_indices)):
            reverse_mapping.append(new_indices.index(i))
        reverse_mapping = tuple(reverse_mapping)

        new_inchi_parts = {"main": "", "isotope": ""}
        inchi_parts = inchi.split("/")

        layer = "main"
        for part in inchi_parts:
            key = part[0]
            if key == "i":
                layer = "isotope"
            if key in ("b", "t", "m", "s"):
                if new_stereo[layer][key] != "":
                    part = new_stereo[layer][key]

            new_inchi_parts[layer] += part + "/"

        for layer in ("main", "isotope"):
            for key in ("b", "t", "m", "s"):
                if new_inchi_parts[layer].find("/%s" % key) == -1:
                    if new_stereo[layer][key] != "":
                        new_inchi_parts[layer] += new_stereo[layer][key]
                        new_inchi_parts[layer] += "/"
        
        
        new_inchi = new_inchi_parts["main"] + new_inchi_parts["isotope"]
        new_inchi = new_inchi[:-1]

        mol = Chem.MolFromInchi(new_inchi)
        mol = Chem.RenumberAtoms(mol, reverse_mapping)
        mol.UpdatePropertyCache()


        for idx in first_in_chain_indices:
            atom = mol.GetAtomWithIdx(idx)
            atom.SetAtomicNum(1)

        edit_mol = EditableMol(mol)
        edit_mol.BeginBatchEdit()
        for dummy_idx in dummy_atom_indices:
            edit_mol.RemoveAtom(dummy_idx)

        

        edit_mol.CommitBatchEdit()
        mol = edit_mol.GetMol()
        mol = Chem.rdmolops.RemoveHs(mol)

        return mol

    def xe_to_rgroups(self, mol: Mol, rgroup_strings: list) -> tuple[Mol, list]:
        
        # Converts Xenon atoms to placeholder R group atoms
        # Removes any empty R group strings from the list, as these represent a
        # link back to a parent molecule, and replaces these with a pseudoatom
        # labelled *

        # Find the position of the parent linker in the list, and for all other
        # R groups add to a new list without the linker

        new_rgroup_strings = []
        if len(rgroup_strings) == 0:
            parent_linker_rank = 0
        else:
            parent_linker_rank = -1
            for i, rgroup_string in enumerate(rgroup_strings):
                if rgroup_string == "":
                    parent_linker_rank = i
                else:
                    new_rgroup_strings.append(rgroup_string)

        # For each Xe, convert to an R group with the appropriate R label, 
        # except for the linker back to the parent
        
        pseudoatoms_counted = 0
        rlabel = 1

        

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 54:
                if pseudoatoms_counted == parent_linker_rank:
                    atom.SetAtomicNum(0)
                    atom.SetIsotope(0)
                    atom.SetProp("dummyLabel", "*")
                    atom.SetProp("isParentLinker", "True")
                    atom.SetNumRadicalElectrons(0)
                    atom.SetNoImplicit(True)
                else:
                    atom.SetAtomicNum(0)
                    atom.SetIsotope(0)
                    atom.SetProp("dummyLabel", "R%i" % rlabel)
                    atom.SetIntProp("_MolFileRLabel", rlabel)
                    atom.SetNumRadicalElectrons(0)
                    atom.SetNoImplicit(True)
                    rlabel += 1
                pseudoatoms_counted += 1
        
        if len(mol.GetAtoms()) == 1 and parent_linker_rank == 0:
            mol = Chem.MolFromSmiles("[*][H]")
            mol.GetAtomWithIdx(0).SetProp("isParentLinker", "True")
            
        return mol, new_rgroup_strings

    def get_markush_strings(self, parts: list) -> list:

        # Gets a list of the direct Markush substructures of this mol
        # This involves parsing the brackets to keep further layers of nested
        # substructures intact

        depth = 0
        markush_components = []
        for part in parts:
            if part[:2] == "M>":

                if depth == 0:
                    component = ""

                depth += 1
        
            component += "<" + part

            if part[:3] == "/M>":
                depth -= 1
                if depth == 0:
                    component = component[3:len(component)-4]
                    markush_components.append(component)

        return markush_components

    def sort_markush_strings(self, strings: list) -> tuple[list, list, list]:
        # Sorts the Markush strings by the type of Markush structure, and 
        # returns a list of each (R groups, variable attachments, atom lists)
        # Also return the Markush stereo layer

        stereo_layer = ""
        rgroups = []
        varattachs = []
        listatoms = []

        for string in strings:
            # Get any text before the molecula formula, as this tells us what
            # type of Markush structure it is.

            markush_type = None
            core = string.split("<")[0]
            formula = core.split("/")[0]
            if formula != "" and formula[0] in ("t", "b", "m", "s"):
                stereo_layer = string
            elif formula.find("-") != -1: 
                pre_formula = formula.split("-")[0]
                if pre_formula.find("H") != -1:
                    varattachs.append(string)
                else:
                    listatoms.append(string)
            else:
                rgroups.append(string)

        return stereo_layer, rgroups, varattachs, listatoms

    def add_varattach(self, mol: Mol, varattach: str) -> Mol:
        # Adds a variable attachment to mol according to the string varattach

        # Extract the endpoints and the core InChI from the string
        # Construct the endpoint label used by RDKit and in molfiles to
        # represent the variable attachment endpoints

        endpts_str, core_str = varattach.split("-", 1)
        
        bond_label_str = "("
        bond_label_str += str(len(endpts_str.split(","))) + " "
        endpt_indices = []
        for endpt in endpts_str.split(","):
            endpt = endpt.replace("H","")
            endpt_indices.append(int(endpt) - 1)
            bond_label_str += endpt + " "
        
        mol.UpdatePropertyCache()
        
        # If there is no H on the endpoint atom, this means we don't have the 
        # correct tautomer, so we need to find all tautomers for this molecule
        # and find one that works
        needs_tautomerising = False
        for endpt_idx in endpt_indices:
            atom = mol.GetAtomWithIdx(endpt_idx)
            if atom.GetTotalNumHs() == 0:
                needs_tautomerising = True
        
        if needs_tautomerising:
            tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
            
            # Isotopically label each atom so they are all unique, and store the
            # original isotopes so we can restore them afterwards
            isotopes = []
            for i, atom in enumerate(mol.GetAtoms()):
                isotopes.append(atom.GetIsotope())
                atom.SetIsotope(1000 + i)           
           
            # Find all tautomers for the molecule and make a list of those where
            # the endpoint atom has an H 
            tautomers = tautomer_enumerator.Enumerate(mol)
            valid_tautomers = []
            for tautomer in tautomers:
                atom = tautomer.GetAtomWithIdx(endpt_idx)
                if atom.GetTotalNumHs() > 0:
                    valid_tautomers.append(tautomer)
            
            # Choose the first of these (arbitrary, but shouldn't matter)
            # Restore the isotopic information
            mol = valid_tautomers[0]
            for atom, isotope in zip(mol.GetAtoms(), isotopes):
                atom.SetIsotope(isotope)


        bond_label_str = bond_label_str[:len(bond_label_str) - 1] + ")"

        # Check if Variable attachment is just an R group, and if so deal with
        # separately
        depth = 0
        is_rgroup = False
        for i in range(len(core_str)):
            if core_str[i:i+3] == "<M>":
                depth += 1
            if core_str[i:i+4] == "</M>":
                depth -= 1
            if core_str[i] == "!" and depth == 0:
                is_rgroup = True
        
        if is_rgroup:
            self.add_rgroup_from_string(core_str)
            
            

            varattach_mol = Chem.MolFromSmiles("ClCl")
            atom1 = varattach_mol.GetAtomWithIdx(0)
            atom2 = varattach_mol.GetAtomWithIdx(1)
            atom2.SetProp("dummyLabel","R")
            atom2.SetAtomicNum(0)
            atom2.SetIntProp("_MolFileRLabel", 0)
            atom1.SetAtomicNum(0)
            atom1.SetProp("isParentLinker", "True")
            varattach_mol = self.update_rlabels(varattach_mol)

        else:

            # Get the mol from the MarkInChI of the core
            parser = MarkinchiParser(core_str)
            varattach_mol, rgroups = parser.parse_markinchi()

            varattach_mol = self.update_rlabels(varattach_mol)
            for i in rgroups.keys():
                self.rgroups.append(rgroups[i])

        index_offset = len(mol.GetAtoms())

        # Combine the core mol with the variable attachment
        mol = Chem.CombineMols(mol, varattach_mol)
        #show(mol)

        # Update the indices of any nested variable attachments on the variable
        # attachment we have just added
        for bond in mol.GetBonds():
            if bond.HasProp("_MolFileBondEndPts"):
                old_endpts = bond.GetProp("_MolFileBondEndPts")
                
                old_endpts = old_endpts.split()

                if (
                    bond.GetBeginAtomIdx() >= index_offset or
                    bond.GetEndAtomIdx() >= index_offset
                ):
                    new_endpts = old_endpts[0] + " " 
                    for old_endpt in old_endpts[1:]:
                        old_endpt = old_endpt.replace(")","")
                        old_endpt = int(old_endpt)
                        new_endpt = old_endpt + index_offset
                        new_endpts += str(new_endpt) + " "
                
                    new_endpts = new_endpts[:len(new_endpts) - 1] + ")"
                    bond.SetProp("_MolFileBondEndPts", new_endpts)

        # Add the endpoint information to the variable attachment bond
        parent_linker_indices = []
        nested_linker_indices = []
        for atom in mol.GetAtoms():
            if (atom.HasProp("isParentLinker") and
                atom.GetIdx() >= index_offset):
                atom.ClearProp("_MolFileRLabel")
                atom.SetProp("dummyLabel", "*")
                nested_linker_indices.append(atom.GetIdx())
                for bond in atom.GetBonds():
                        bond.SetProp("_MolFileBondType", str(1))
                        bond.SetProp("_MolFileBondEndPts", bond_label_str)
                        bond.SetProp("_MolFileBondAttach", "ANY")
            if (atom.HasProp("isParentLinker") and
                atom.GetIdx() < index_offset):
                parent_linker_indices.append(atom.GetIdx())

        molblock = Chem.MolToV3KMolBlock(mol)
        mol = Chem.MolFromMolBlock(molblock, sanitize=False)
        Chem.rdmolops.SanitizeMol(mol)

        for idx in parent_linker_indices:
            atom = mol.GetAtomWithIdx(idx)
            atom.SetProp("dummyLabel", "*")
            atom.SetProp("isParentLinker", "True")
        
        for idx in nested_linker_indices:
            atom = mol.GetAtomWithIdx(idx)
            atom.SetProp("dummyLabel", "*")

        return mol

    def add_rgroup_from_string(self, rgroup_string: str) -> None:
        # Parses a Markinchi substring for an R group to generate a list of Mols
        # Adds this R group to the list of R groups, as well as the list of any
        # nested R groups within this one.

        components = []
        component = ""
        depth = 0
        
        for i in range(len(rgroup_string)):

            if rgroup_string[i] == "!" and depth == 0:
                components.append(component)
                component = ""
            else:
                component += rgroup_string[i]
            
            if rgroup_string[i:i+3] == "<M>":
                depth += 1

            if rgroup_string[i:i+4] == "</M>":
                depth -= 1

        components.append(component)

        rgroup = []
        nested_rgroups = []
        for component in components:
            if component != "":
                parser = MarkinchiParser(component)
                mol, rgroups = parser.parse_markinchi()
                rgroup.append(mol)
                for i in rgroups.keys():
                    nested_rgroups.append(rgroups[i])

        self.rgroups.append(rgroup)
        self.nested_rgroups.append(nested_rgroups)

    def update_rlabels(self, mol: Mol) -> Mol:
        
        # Increase the R label of any nested R groups in this fragment to 
        # avoid conflicts with the R groups in the parent

        for atom in mol.GetAtoms():
            if atom.HasProp("_MolFileRLabel"):
                
                old_label = atom.GetIntProp("_MolFileRLabel")
                new_label = old_label + len(self.rgroups)
                atom.SetProp("dummyLabel", "R%i" % new_label)
                atom.SetIntProp("_MolFileRLabel", new_label)
        
        return mol

    def update_nested_rlabels(self, rgroups: list) -> list:

        updated_rgroups = []
        for rgroup in rgroups:

            updated_rgroup = []
            for component in rgroup:

                updated_component = self.update_rlabels(component)
                updated_rgroup.append(updated_component)

            updated_rgroups.append(updated_rgroup)
        
        return updated_rgroups

    def add_listatom(self, mol: Mol, listatom: str) -> Mol:

        # For each listatom, turns the atom in mol into a query atom with a 
        # list of elements 

        index = int(listatom.split("-")[0])
        elements = listatom.split("-")[1].split("!")

        molfile_string = "["
        for element in elements:
            molfile_string += element
            molfile_string += ","

        atom_symbol = mol.GetAtomWithIdx(index - 1).GetSymbol()
        
        molfile_string = molfile_string[:len(molfile_string) - 1] + "]"
        replace_string = "M  V30 %i %s" % (index, atom_symbol)
        new_string = "M  V30 %i %s" % (index, molfile_string)

        parent_linker_indices = []
        for atom in mol.GetAtoms():
            if atom.HasProp("isParentLinker"):
                parent_linker_indices.append(atom.GetIdx())

        molblock = Chem.MolToV3KMolBlock(mol)
        molblock = molblock.replace(replace_string, new_string)
        mol = Chem.MolFromMolBlock(molblock)

        for idx in parent_linker_indices:
            atom = mol.GetAtomWithIdx(idx)
            atom.SetProp("isParentLinker", "True")
            atom.SetProp("dummyLabel", "*")
            atom.ClearProp("_MolFileRLabel")


        return mol

    def get_molblock(self) -> str:
        
        core_block = Chem.MolToV3KMolBlock(self.core_mol)
        core_block = core_block.replace(" VAL=1","")
        core_block = core_block.replace("\nM  END", "")
        
        for i, rgroup in enumerate(self.rgroups):
            rgroup_block = "M  V30 BEGIN RGROUP %i" % (i + 1)
            rgroup_block += "\nM  V30 RLOGIC 0 0 \"\"\n"
            for component in rgroup:
                for atom in component.GetAtoms():
                    if atom.HasProp("isParentLinker"):
                        for neighbor in atom.GetNeighbors():
                            neighbor.SetProp("isAttachPoint", "True")
                        
                        linker_idx = atom.GetIdx()
                        editmol = EditableMol(component)
                        editmol.RemoveAtom(linker_idx)
                        new_component = editmol.GetMol()

                if len(new_component.GetAtoms()) == 0:
                    new_component = Chem.MolFromSmiles("[H]")
                    atom = new_component.GetAtomWithIdx(0)
                    atom.SetProp("isAttachPoint", "True")
                
                for atom in new_component.GetAtoms():
                    if atom.HasProp("isAttachPoint"):
                        idx = atom.GetIdx() + 1
                        element = atom.GetSymbol()
                        
                if (len(new_component.GetAtoms()) == 1 and
                    new_component.GetAtomWithIdx(0).GetAtomicNum() == 1):
                    new_component = Chem.AddHs(new_component)
      
                component_block = Chem.MolToV3KMolBlock(new_component)
                new_component_block = ""
                lines = component_block.split("\n")
                lines = lines[:len(lines) - 1]
                for i, line in enumerate(lines):
                    if line.find("M  V30 %i %s" % (idx, element)) != -1:
                        line = line + " ATTCHPT=1"
                        if line.find("VAL") != -1:
                            new_line = ""
                            lineparts = line.split()
                            for part in lineparts:
                                if part.find("VAL") == -1:
                                    new_line += part
                                    new_line += " "
                            line = new_line.replace("M ", "M  ")
                    if i > 3:
                        if line.find("M  END") == -1:
                            new_component_block += line + "\n"        
                
                rgroup_block += new_component_block
            
            rgroup_block += "M  V30 END RGROUP\n"
            core_block += rgroup_block

        core_block += "M  END"

        return core_block

if __name__ == "__main__":

    from MarkinchiGenerator import MarkinchiGenerator

    debug = False

    filename = "molfiles\\structures_for_testing\\exx2.mol"
    filedir = os.path.join(os.getcwd(), filename)
    markinchi_generator = MarkinchiGenerator()

    markinchi_generator.get_from_molfile(filedir)
    
    markinchi = markinchi_generator.generate_markinchi()

    print(markinchi)

    parser = MarkinchiParser(markinchi)
    mol, rgroups = parser.parse_markinchi()
    molblock = parser.get_molblock()
    for i in rgroups.keys():
        print("R%i" % (i))
        for component in rgroups[i]:
            show(component, indices=True)

    markinchi_generator = MarkinchiGenerator()
    markinchi_generator.get_from_molblock(molblock)
    new_markinchi = markinchi_generator.generate_markinchi()
    print(new_markinchi)
    


    
    
    import MarkinchiUtils as MIUtils

    mol_list = MIUtils.enumerate_markush_mol(mol, rgroups)
    new_inchi_list = list(
    set(MIUtils.inchis_from_mol_list(mol_list)))
    new_inchi_list = sorted(new_inchi_list)


    ref_mol, ref_rgroups = MIUtils.parse_molfile(filename)
    ref_list = MIUtils.enumerate_markush_mol(ref_mol, ref_rgroups)
    ref_inchi_list = MIUtils.inchis_from_mol_list(ref_list)
    ref_inchi_list = list(set(MIUtils.inchis_from_mol_list(ref_list)))
    ref_inchi_list = sorted(ref_inchi_list)
    print(ref_inchi_list)
    print(new_inchi_list)
    
    print(new_inchi_list == ref_inchi_list)


    for inchi in new_inchi_list:
        show(Chem.MolFromInchi(inchi))


    