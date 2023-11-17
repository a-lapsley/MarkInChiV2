from rdkit import Chem
from rdkit.Chem.rdchem import EditableMol

from MarkinchiGenerator import Show

import os
import sys, getopt

from typing import TypeAlias

Mol: TypeAlias = Chem.rdchem.Mol


class MarkinchiParser(object):

    # Object for converting a MarkInChI string into an RDKit molecule and R
    # groups

    def __init__(self, markinchi: str = "") -> None:
        self.markinchi = markinchi
        self.mol = None
        self.rgroups = None

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
        
        for varattach_string in varattach_strings:
            core_mol = self.add_varattach(core_mol, varattach_string)

        core_mol = self.xe_to_rgroups(core_mol)

        rgroups = []

        return core_mol, rgroups

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

    def xe_to_rgroups(self, mol: Mol) -> Mol:
        
        # Converts Xenon atoms to placeholder R group atoms

        rlabel = 1
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 54:
                atom.SetAtomicNum(0)
                atom.SetProp("dummyLabel", "R%i" % rlabel)
                atom.SetIntProp("_MolFileRLabel", rlabel)
                atom.SetNumRadicalElectrons(0)
        
        return mol

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

        # Get the mol from the MarkInChI of the core
        parser = MarkinchiParser(core_str)
        varattach_mol, rgroups = parser.parse_markinchi()

        # If variable attachment has no child R groups, the only pseudoatom is 
        # the linker of the variable attachment to the parent

        if len(rgroups) == 0:
            parent_linker_rank = 0

        # Change the appropriate pseudoatom to the linker back to the parent        
        pseudoatoms_counted = 0

        for atom in varattach_mol.GetAtoms():
            if atom.HasProp("_MolFileRLabel"):
                if pseudoatoms_counted == parent_linker_rank:
                    atom.ClearProp("_MolFileRLabel")
                    atom.SetProp("dummyLabel", "*")
                    atom.SetProp("isParentLinker", "True")

                pseudoatoms_counted += 1


        index_offset = len(mol.GetAtoms())

        # Combine the core mol with the variable attachment
        mol = Chem.CombineMols(mol, varattach_mol)
        mol.UpdatePropertyCache()

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
        for atom in mol.GetAtoms():
            if atom.HasProp("isParentLinker"):
                atom.ClearProp("_MolFileRLabel")
                atom.SetProp("dummyLabel", "*")
                for bond in atom.GetBonds():
                        bond.SetProp("_MolFileBondType", str(1))
                        bond.SetProp("_MolFileBondEndPts", bond_label_str)
                        bond.SetProp("_MolFileBondAttach", "ANY")

        molblock = Chem.MolToV3KMolBlock(mol)
        mol = Chem.MolFromMolBlock(molblock)

        return mol

if __name__ == "__main__":

    from MarkinchiGenerator import MarkinchiGenerator

    filename = "molfiles\\test10.mol"
    filedir = os.path.join(os.getcwd(), filename)
    markinchi_generator = MarkinchiGenerator()

    markinchi_generator.get_from_molfile(filedir)
    
    markinchi = markinchi_generator.generate_markinchi()

    print(markinchi)

    parser = MarkinchiParser(markinchi)
    mol, rgroups = parser.parse_markinchi()

    
    molblock = Chem.MolToV3KMolBlock(mol)
    
    markinchi_generator = MarkinchiGenerator()
    markinchi_generator.get_from_molblock(molblock)
    print(markinchi_generator.generate_markinchi())

    Show(mol, overrideDebug=True, indices=True)
    