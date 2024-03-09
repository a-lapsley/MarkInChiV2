from rdkit.Chem.Draw import ShowMol, MolsToImage
from rdkit.Chem.rdchem import Mol
from rdkit import Chem
from copy import deepcopy
import json
import os

ABBREVIATION_FILE = "abbreviations.json"

def show(mols: list, 
         subImgSize: tuple = (200, 200), 
         title: str = 'RDKit Molecule',
         stayInFront: bool = True, 
         indices: bool = False, 
         tidyCoords: bool = True,
         **kwargs) -> None:
    """
    Generates a picture of molecule(s) and displays it in a Tkinter window.

    This function is a copy of the ShowMol function from the RDKit source code
    that displays multiple Mols in the same Tkinter window.

    It is only used for debugging purposes.
    """
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
        # If specified, label each atom with it's index
        if indices:
            for atom in mols.GetAtoms():
                atom.SetProp("molAtomMapNumber", str(atom.GetIdx() + 1))
        # If specified, automatically compute coordinates for atoms to make the
        # image more legible

        if tidyCoords:
            Chem.rdDepictor.Compute2DCoords(mols)
        ShowMol(mols)

def ctab_to_molblock(ctab: str) -> str:

    # Adds necessary beginning and end to a CTAB to allow RDKit to read it
    # as a molblock
    molblock = "\n\n\n  0  0  0     0  0            999 V3000\n"
    molblock += ctab
    molblock += "M  END"

    return molblock

def parse_molfile(filename: str) -> tuple[Mol, dict]:

    with open(filename) as file:
                molfile_lines = file.readlines()

    writing_core = 0
    # status of writing the core molecule string
    # 0=core not yet written
    # 1=writing core
    # 2=core writing finished
    writing_rgroups = False

    core_ctab = ""
    rgroups = {}
    abbr_rgroups = {}

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
            rgroup = []
            writing_rgroups = True
            rgroup_ctab = ""

        if line.find("END RGROUP") != -1:
            rgroups[rgroup_number] = rgroup
            pass

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
                component_molblock = ctab_to_molblock(rgroup_ctab)
                component, nested_abbrs = parse_molblock(component_molblock)
                abbr_rgroups = {**nested_abbrs, **abbr_rgroups}

                component = Chem.rdmolops.AddHs(component)
                attach_idx = rgroup_attachments[0][1] - 1
                atom = component.GetAtomWithIdx(attach_idx)

                linker_added = False
                for neighbor in atom.GetNeighbors():
                    if (neighbor.GetAtomicNum() == 1 and not 
                        linker_added):
                        neighbor.SetAtomicNum(0)
                        neighbor.SetProp("dummyLabel","*")
                        neighbor.SetProp("isParentLinker", "True")
                        linker_added = True
                
                if atom.GetAtomicNum() == 0:
                    linker_atom = Chem.Atom(0)
                    linker_atom.SetProp("dummyLabel", "*")
                    linker_atom.SetProp("isParentLinker", "True")
                    edit_mol = Chem.EditableMol(component)
                    linker_idx = edit_mol.AddAtom(linker_atom)
                    edit_mol.AddBond(
                        atom.GetIdx(),
                        linker_idx,
                        Chem.rdchem.BondType.SINGLE
                    )
                    component = edit_mol.GetMol()
                    linker_added = True

                if not linker_added:
                    raise Exception("Invalid R group connection point")
                
                atom.SetChiralTag(
                    Chem.rdchem.ChiralType.CHI_UNSPECIFIED
                )

                component = Chem.rdmolops.RemoveHs(component)
                

                # Set the geometry of all bonds to the attachment point
                # to be undefined
                # We need to do this as the mol file format is ambiguous
                # with its definition of R group connections
                new_atom = component.GetAtomWithIdx(atom_idx)
                for bond in new_atom.GetBonds():
                    if(
                        bond.GetBondType() == 
                        Chem.rdchem.BondType.DOUBLE
                        ):
                        
                        bond.SetStereo(
                            Chem.rdchem.BondStereo.STEREOANY
                        )
                
                component.UpdatePropertyCache()
                Chem.rdmolops.SanitizeMol(component)

                rgroup.append(component)

                pass

    core_molblock = ctab_to_molblock(core_ctab)
    core_mol, core_abbrs = parse_molblock(core_molblock)

    abbr_rgroups = {**abbr_rgroups, **core_abbrs}

    for id in abbr_rgroups.keys():
        rgroup = []
        for component in abbr_rgroups[id]:
            for atom in component.GetAtoms():
                if atom.GetAtomicNum() == 54:
                    atom.SetAtomicNum(0)
                    atom.SetProp("dummyLabel","*")
                    atom.SetProp("isParentLinker", "True")
            rgroup.append(component)
        rgroups[id] = rgroup

    
    return core_mol, rgroups

def parse_molblock(molblock: str) -> tuple[Mol, dict]:
     
    molblock, rgroups = abbreviations_to_rgroups(molblock)
    mol = Chem.MolFromMolBlock(molblock, sanitize=False)
    Chem.rdmolops.SanitizeMol(mol)
    Chem.rdmolops.AssignStereochemistry(mol)

    return mol, rgroups

def abbreviations_to_rgroups(molblock: str) -> tuple[str, dict]:

    abbr_dict = get_abbr_dict()

    new_molblock = ""
    present_abbrs = []
    lines = molblock.split("\n")
    
    atom_block = False
    for line in lines:

        if atom_block:
            element = line.split()[3]
            if element.lower() in abbr_dict.keys():
                id = abbr_dict[element.lower()]["id"]
                line = line.replace(element, "R%i" % id)
                line += " RGROUPS=(1 %i)" % id
                present_abbrs.append(element.lower())

        if line.find("BEGIN ATOM") != -1:
            atom_block = True
        if line.find("END ATOM") != -1:
            atom_block = False

        
        new_molblock += line + "\n"

    rgroups = {}
    for abbr in present_abbrs:
        rgroup = []
        for option in abbr_dict[abbr]["options"]:
            mol = Chem.MolFromInchi(option)
            rgroup.append(mol)
        id = abbr_dict[abbr]["id"]
        rgroups[id] = rgroup
    
    return new_molblock, rgroups

def get_abbr_dict() -> dict:

    file = os.path.join(os.getcwd(), ABBREVIATION_FILE)
    
    with open(file) as f:
        definitions = json.load(f)

    abbr_dict = {}    
    keys = sorted(definitions.keys())
    
    for i, key in enumerate(keys):
        abbr_dict[key] = {}
        abbr_dict[key]["id"] = 10000 + i
        abbr_dict[key]["options"] = enumerate_options(key, definitions)

    return abbr_dict

def enumerate_options(key: str, definitions: dict) -> list:
    
    options = []

    for option in definitions[key]:
        if option in definitions.keys():
            options += enumerate_options(option, definitions)
        else:
            options += [option]
    
    return options
    
