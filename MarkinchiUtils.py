from rdkit.Chem.Draw import ShowMol, MolsToImage
from rdkit.Chem.rdchem import Mol
from rdkit import Chem
from copy import deepcopy

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
                component = Chem.MolFromMolBlock(
                    component_molblock, sanitize=False
                    )
                
                Chem.rdmolops.SanitizeMol(component)
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
    core_mol = Chem.MolFromMolBlock(core_molblock, sanitize=False)
    Chem.rdmolops.SanitizeMol(core_mol)
    
    return core_mol, rgroups

    