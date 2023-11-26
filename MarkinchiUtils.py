from rdkit.Chem.Draw import ShowMol, MolsToImage
from rdkit.Chem.rdchem import Mol, Atom, EditableMol
from rdkit.Chem import rdMolEnumerator, rdDistGeom
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

def enumerate_markush_mol(mol: Mol, rgroups: list) -> list:
    # Gets a list of all possible Mols for a Markush Mol and list of R groups
    
    # Get any R groups directly referenced by this Mol
    child_rgroups = {}
    rgroup_count = 0
    for atom in mol.GetAtoms():
        if atom.HasProp("_MolFileRLabel"):
            rlabel = atom.GetProp("_MolFileRLabel")
            child_rgroups[rlabel] = rgroups[int(rlabel)-1]
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

    return mol_list

def enumerate_listatoms(mol: Mol) -> list:
    #Enumerates all possibilities of listatoms in a Mol
    #Find the atoms with an atom list query
    for atom in mol.GetAtoms():
        atom.ClearProp("molAtomMapNumber")
        if atom.HasQuery():
            #Parse the atom SMARTS to get a list of elements it could be
            smarts = atom.GetSmarts()
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
                atom.SetAtomicNum(54)
                new_list = []

                # For each component the R group could be, add the component to 
                # the core and add this new Mol to the list of possibilities
                for component in rgroup:
                    # Find the pseudoatom that represents the connection point
                    link_idx = None
                    for comp_atom in component.GetAtoms():
                        if comp_atom.HasProp("isParentLinker"):
                            # Mark the atom which connects to the core
                            for comp_neighbor in comp_atom.GetNeighbors():
                                comp_neighbor.SetProp("rgroupLinkPoint", "True")
                                
                            # Remove the pseudoatom
                            edit_mol = EditableMol(component)
                            edit_mol.RemoveAtom(comp_atom.GetIdx())
                            component = edit_mol.GetMol()
                    
                    # Get the index of the connecting atom on the component
                    for comp_atom in component.GetAtoms():
                        if comp_atom.HasProp("rgroupLinkPoint"):
                            link_idx = comp_atom.GetIdx()
                    
                    # If the component is empty, it is just an H atom
                    if link_idx == None:
                        component = Chem.MolFromSmiles("[H]")
                        link_idx = 0

                    # Use RDKit replace substructs to replace the Xe atom with 
                    # the component
                    # We need to turn the molecule into a 3D conformation and
                    # then recompute the stereochemistry afterwards to ensure 
                    # it is correctly conserved
                    # try:
                    if True:
                        replace_struct = Chem.MolFromSmiles("[Xe]")
                        mol_copy = Chem.rdmolops.AddHs(mol, addCoords=True)
                        #rdDistGeom.EmbedMolecule(mol_copy, maxAttempts=20)
                        new_mol = Chem.ReplaceSubstructs(
                            mol_copy,
                            replace_struct,
                            component,
                            replaceAll = True,
                            replacementConnectionPoint = link_idx
                        )[0]

                        Chem.rdmolops.AssignStereochemistryFrom3D(new_mol)
                        new_mol = Chem.rdmolops.RemoveHs(new_mol)
                        new_mol.UpdatePropertyCache()
                        new_list.append(new_mol)
                    #except:
                    else:
                        print("WARNING: Skipping invalid structure")
    
    # Remove any markers from the atoms in the Mols in the list we have
    # generated to prevent it causing further issues 
    for new_mol in new_list:
        for atom in new_mol.GetAtoms():
            atom.ClearProp("coreLinkPoint")
            atom.ClearProp("rgroupLinkPoint")

        new_mol.UpdatePropertyCache()
        
    return new_list
        
def enumerate_varattachs(mol: Mol) -> list:

    # Enumerates all possibilities for variable attachments
    # This is straightforwardly done using RDKit's in built method
    bundle = rdMolEnumerator.Enumerate(mol)

    if len(bundle) == 0:
        return [mol]
    else:
        return bundle
    
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