from rdkit.Chem.Draw import ShowMol, MolsToImage
from rdkit.Chem.rdchem import Mol, Atom, EditableMol
from rdkit.Chem import rdMolEnumerator, rdDistGeom, AllChem
from rdkit import Chem
from copy import deepcopy

from molzip import molzip


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

def enumerate_markush_mol(mol: Mol, rgroups: dict) -> list:
    # Gets a list of all possible Mols for a Markush Mol and list of R groups
    
    # Get any R groups directly referenced by this Mol
    child_rgroups = {}
    rgroup_count = 0
    for atom in mol.GetAtoms():
        if atom.HasProp("_MolFileRLabel"):
            rlabel = atom.GetProp("_MolFileRLabel")
            child_rgroups[rlabel] = rgroups[int(rlabel)]
            rgroup_count += 1
    
    # For each child R group, get a list of all possible structures it could be
    enumerated_rgroups = {}
    for r in child_rgroups.keys():

        enumerated_components = []
        rgroup = child_rgroups[r]

        for component in rgroup:
            enumerated_components += enumerate_markush_mol(component, rgroups)

        enumerated_rgroups[r] = enumerated_components
    

    #Enumerate the varattachs using the built in RDKit functionality
    new_mol_list = enumerate_varattachs(mol)

    mol_list = new_mol_list

    #Enumerate all R groups
    for rlabel in child_rgroups.keys():
        rgroup = enumerated_rgroups[rlabel]
        new_mol_list = []
        for listmol in mol_list:
            new_mol_list += enumerate_rgroups(listmol, rlabel, rgroup)
        mol_list = new_mol_list

    # Enumerate all list atoms
    new_mol_list = []
    for mol in mol_list:
        enumerated_list = enumerate_listatoms(mol)
        new_mol_list += enumerated_list

    mol_list = new_mol_list

    # Sanitize molecules to remove any invalid molecules
    new_mol_list = []

    for mol in mol_list:
        try:
            Chem.rdmolops.SanitizeMol(mol)
            new_mol_list.append(mol)
        except:
            #show(mol)
            print("Skipping an invalid molecule")
    
    mol_list = new_mol_list

    return mol_list

def enumerate_listatoms(mol: Mol) -> list:
    #Enumerates all possibilities of listatoms in a Mol
    #Find the atoms with an atom list query
    for atom in mol.GetAtoms():
        atom.ClearProp("molAtomMapNumber")
        if atom.HasQuery():
            #Parse the atom SMARTS to get a list of elements it could be
            smarts = atom.GetSmarts()
            print(smarts)
            if smarts != "*":

                smarts = smarts.replace("[", "")
                smarts = smarts.replace("]", "")
                smarts = smarts.replace("#", "")
                smarts_parts = smarts.split(",")

                atomic_nums = []

                for part in smarts_parts:
                    atomic_nums.append(int(part))

                listatom_idx = atom.GetIdx()

                new_list = []

                #For each element, replace the list atom with an atom of that
                #element
                for atomic_num in atomic_nums:

                    new_atom = Atom(atomic_num)
                    edit_mol = EditableMol(mol)
                    new_atom_idx = edit_mol.AddAtom(new_atom)

                    for bond in atom.GetBonds():

                        if bond.GetBeginAtomIdx() == listatom_idx:
                            new_bond_end = bond.GetEndAtomIdx()
                        else:
                            new_bond_end = bond.GetBeginAtomIdx()

                        bond_type = bond.GetBondType()

                        edit_mol.AddBond(
                            new_atom_idx, new_bond_end, order=bond_type
                        )

                    edit_mol.RemoveAtom(listatom_idx)

                    new_mol = edit_mol.GetMol()

                    # For this new Mol with the replaced listatom, recursively
                    # replace further listatoms to ensure all listatoms are 
                    # dealt with
                    new_list += enumerate_listatoms(new_mol)
                return new_list

    # If there are no listatoms left in this Mol, return it as a single element
    # list            
    return [mol]
                
def enumerate_rgroups(mol: Mol, rlabel: int, rgroup: list) -> list:
    # Enumerates possibilities for an R group in a Mol 
    
    # Find the R group pseudoatoms
    for atom in mol.GetAtoms():
        if atom.HasProp("_MolFileRLabel"):
            atom_rlabel = atom.GetProp("_MolFileRLabel")
            # Replace the pseudoatom with an Xe atom as these are easier to work
            # with later
            if atom_rlabel == rlabel:
                atom.SetProp("isLinker", "True")
                atom.SetProp("dummyLabel", "*")

                mol = mark_stereo(mol, atom.GetIdx())
                        
                new_list = []


                # For each component the R group could be, add the component to 
                # the core and add this new Mol to the list of possibilities
                for component in rgroup:

                    # Find the pseudoatom that represents the connection point

                    for comp_atom in component.GetAtoms():
                        if comp_atom.HasProp("isParentLinker"):
                            # Mark the atom which connects to the core
                            comp_atom.SetIntProp("molAtomMapNumber", 1)

                            component = mark_stereo(
                                component, comp_atom.GetIdx())

                    # Use RDKit replace substructs to replace the Xe atom with 
                    # the component
                    # We need to turn the molecule into a 3D conformation and
                    # then recompute the stereochemistry afterwards to ensure 
                    # it is correctly conserved
                    try:  
                        mol_copy = Chem.Mol(mol)

                        for atom_copy in mol_copy.GetAtoms():
                            if atom_copy.HasProp("isLinker"):
                                atom_copy.SetIntProp("molAtomMapNumber", 1)

                        new_mol = molzip(mol_copy, component)
                        new_mol = new_mol.GetMol()
                        new_mol = update_stereo(new_mol)
                        new_mol = cleanup_flags(new_mol)
                        
                        new_list.append(new_mol)
                    except:
                        print("WARNING: Skipping invalid structure")
    
    # Remove any markers from the atoms in the Mols in the list we have
    # generated to prevent it causing further issues 
    for new_mol in new_list:
        for atom in new_mol.GetAtoms():
            atom.ClearProp("coreLinkPoint")
            atom.ClearProp("rgroupLinkPoint")

        new_mol.UpdatePropertyCache()
        
    return new_list

def cleanup_flags(mol: Mol) -> Mol:

    for atom in mol.GetAtoms():
        atom.ClearProp("molAtomMapNumber")
        atom.ClearProp("bond_me")
        atom.ClearProp("delete_me")
        atom.ClearProp("fixed")

    return mol

def mark_stereo(mol: Mol, idx: int) -> Mol:

    # Label any bonds whose stereochemistry depends on atom with index idx
    # This allows us to correctly relabel them if this atom is replaced

    for bond in mol.GetBonds():
        stereo_atoms = bond.GetStereoAtoms()
        if len(stereo_atoms) != 0:
            if idx in stereo_atoms:
                bond.SetProp("updateStereo", "True")

                if idx == stereo_atoms[0]:
                    fixed_atom = mol.GetAtomWithIdx(stereo_atoms[1])
                else:
                    fixed_atom = mol.GetAtomWithIdx(stereo_atoms[0])

                fixed_atom.SetProp("fixed", "True")

    return mol

def update_stereo(mol: Mol) -> Mol:

    # Update the stereo labels for any bonds which were affected by combining
    # two fragments

    for bond in mol.GetBonds():
        if bond.HasProp("updateStereo"):
            bgn = bond.GetBeginAtom()
            end = bond.GetEndAtom()

            
            # atom_1 is the one which has been reconnected
            # atom_2 is the one that has not been changed
            if bgn.HasProp("bond_me"):    
                atom_1 = bgn
                atom_2 = end
            else:
                atom_1 = end
                atom_2 = bgn
               
            for neighbor in atom_1.GetNeighbors():
                if neighbor.HasProp("bond_me"):
                    marker_1 = neighbor.GetIdx()
            
            for neighbor in atom_2.GetNeighbors():
                if neighbor.HasProp("fixed"):
                    marker_2 = neighbor.GetIdx()

            if atom_1.GetIdx() == bgn.GetIdx():
                bond.SetStereoAtoms(marker_1, marker_2)
            else:
                bond.SetStereoAtoms(marker_2, marker_1)
            
        bond.ClearProp("updateStereo")
    
    return mol

def enumerate_varattachs(mol: Mol) -> list:

    # Enumerates all possibilities for variable attachments

    new_mol_list = []
    for bond in mol.GetBonds():
        if bond.HasProp("_MolFileBondEndPts"):

            atom_a = bond.GetBeginAtom()
            atom_b = bond.GetEndAtom()

            if (atom_a.GetAtomicNum() == 0 and not
                atom_a.HasProp("_MolFileRLabel")):
                delete_idx = atom_a.GetIdx()
                attach_idx = atom_b.GetIdx()
            elif (atom_b.GetAtomicNum() == 0 and not
                atom_b.HasProp("_MolFileRLabel")):
                delete_idx = atom_b.GetIdx()
                attach_idx = atom_a.GetIdx()
            
            endpt_str = bond.GetProp("_MolFileBondEndPts")[1:-1]
            endpts = endpt_str.split()[1:]


            for endpt in endpts:
                edit_mol = EditableMol(mol)
                edit_mol.BeginBatchEdit()
                end_idx = int(endpt) - 1

                end_atom = mol.GetAtomWithIdx(end_idx)
                if end_atom.GetTotalNumHs() > 0:
                    edit_mol.AddBond(
                                attach_idx,
                                end_idx,
                                order=Chem.rdchem.BondType.SINGLE
                            )
                    edit_mol.RemoveAtom(delete_idx)
                    edit_mol.CommitBatchEdit()
                    new_mol = edit_mol.GetMol()
                    new_mol.UpdatePropertyCache()
                    new_mol_list += enumerate_varattachs(new_mol)

            return new_mol_list
        
    return [mol]
    
def inchis_from_mol_list(mol_list: list) -> list:

    # Gets the InChI for each mol in a list of Mols
    # If the structure is invalid, the InChI is given as None

    inchi_list = []

    for mol in mol_list:
        try:
            inchi = Chem.MolToInchi(mol)
            inchi_list.append(inchi)
        except:
            inchi = None
            inchi_list.append(inchi)

    return inchi_list

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

if __name__ == "__main__":

    ref_mol, ref_rgroups = parse_molfile("D:\\alexl\\Documents\\MarkInChiV2\\molfiles\\test62.mol")
    ref_list = enumerate_markush_mol(ref_mol, ref_rgroups)

    for mol in ref_list:
        show(mol)
    ref_inchi_list = sorted(inchis_from_mol_list(ref_list))
    