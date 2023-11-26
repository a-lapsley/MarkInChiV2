from rdkit import Chem
from rdkit.Chem.rdchem import EditableMol, Mol
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
        self.nested_rgroups = []
        self.core_mol = None

    def set_markinchi(self, markinchi: str) -> None:
        self.markinchi = markinchi

    def parse_markinchi(self) -> tuple[Mol, list]:
        
        # Main script for parsing the MarkInChI

        parts = self.markinchi.split("<")

        core_mol = self.get_core(parts[0])
        
        markush_strings = self.get_markush_strings(parts[1:])
        string_lists = self.sort_markush_strings(markush_strings)
        rgroup_strings = string_lists[0]
        varattach_strings = string_lists[1]
        listatom_strings = string_lists[2]

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

        for listatom_string in listatom_strings:
            core_mol = self.add_listatom(core_mol, listatom_string)
        
        if debug:
            show(core_mol, indices=True)

        self.core_mol = core_mol

        return self.core_mol, self.rgroups

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


        if inchi.find("Xe") != -1:
            xe_count = inchi.split("/")[1].split("Xe")[1]
            if xe_count == "":
                xe_count = 1
            else:
                xe_count = int(xe_count)
        else:
            xe_count = 0

        stripped_inchi = ""

        for i, part in enumerate(inchi.split("/")):
            if i < 5:
                stripped_inchi += part + "/"
        stripped_inchi = stripped_inchi[:-1]

        atom_count = len(Chem.MolFromInchi(stripped_inchi).GetAtoms())
        print(inchi)
        inchi_parts = inchi.split("/i")
        new_inchi = inchi_parts[0]

        if len(inchi_parts) == 1:
            isotope_layer = "/i"
            isotope_parts = []
        else:
            isotope_parts = inchi_parts[1].split("/")
            isotope_layer = "/i" + isotope_parts[0]

        if len(isotope_parts) > 1:
            stereo_layer = ""
            for part in isotope_parts[1:]:
                stereo_layer += "/" + part
        else:
            stereo_layer = ""

        for i in range(atom_count - xe_count, atom_count):
            idx = i + 1
            if i <= 100:
                isotope_layer += ",%i-%i" % (
                    idx, 100 - (i - (atom_count - xe_count)) )
            if i >= 100:
                isotope_layer += ",%i+%i" % (
                    i + 1, atom_count - xe_count + idx)

        isotope_layer = isotope_layer.replace("i,", "i")       
        new_inchi += isotope_layer
        new_inchi += stereo_layer
        print(new_inchi)
        core_mol = Chem.MolFromInchi(new_inchi)
        return core_mol

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
                else:
                    atom.SetAtomicNum(0)
                    atom.SetIsotope(0)
                    atom.SetProp("dummyLabel", "R%i" % rlabel)
                    atom.SetIntProp("_MolFileRLabel", rlabel)
                    atom.SetNumRadicalElectrons(0)
                    rlabel += 1
                pseudoatoms_counted += 1
        
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

        rgroups = []
        varattachs = []
        listatoms = []

        for string in strings:
            # Get any text before the molecula formula, as this tells us what
            # type of Markush structure it is.

            markush_type = None

            core = string.split("<")[0]
            formula = core.split("/")[0]
            if formula.find("-") != -1: 
                pre_formula = formula.split("-")[0]
                if pre_formula.find("H") != -1:
                    varattachs.append(string)
                else:
                    listatoms.append(string)
            else:
                rgroups.append(string)

        return rgroups, varattachs, listatoms

    def add_varattach(self, mol: Mol, varattach: str) -> Mol:
        # Adds a variable attachment to mol according to the string varattach

        # Extract the endpoints and the core InChI from the string
        # Construct the endpoint label used by RDKit and in molfiles to
        # represent the variable attachment endpoints

        endpts_str, core_str = varattach.split("-", 1)
        
        bond_label_str = "("
        bond_label_str += str(len(endpts_str.split(","))) + " "
        for endpt in endpts_str.split(","):
            endpt = endpt.replace("H","")
            bond_label_str += str(int(endpt)) + " "
        
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
            self.rgroups += rgroups

        index_offset = len(mol.GetAtoms())

        # Combine the core mol with the variable attachment
        mol = Chem.CombineMols(mol, varattach_mol)

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
        mol = Chem.MolFromMolBlock(molblock)

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
                nested_rgroups += rgroups

        self.rgroups.append(rgroup)
        self.nested_rgroups.append(rgroups)

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

    filename = "molfiles\\ext92.mol"
    filedir = os.path.join(os.getcwd(), filename)
    markinchi_generator = MarkinchiGenerator()

    markinchi_generator.get_from_molfile(filedir)
    
    markinchi = markinchi_generator.generate_markinchi()

    print(markinchi)

    parser = MarkinchiParser(markinchi)
    mol, rgroups = parser.parse_markinchi()
    show(mol)
    for i, rgroup in enumerate(rgroups):
        print("R%i" % (i + 1))
        for component in rgroup:
            show(component)
    
    from MarkinchiUtils import enumerate_markush_mol

    mol_list = enumerate_markush_mol(mol, rgroups)

    for mol in mol_list:
        show(mol)

    molblock = parser.get_molblock()

    #print(molblock)
    #show(mol, indices=True)
    markinchi_generator = MarkinchiGenerator()
    markinchi_generator.get_from_molblock(molblock)
    print(markinchi_generator.generate_markinchi())

    