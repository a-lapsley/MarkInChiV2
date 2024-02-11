from rdkit.Chem.rdchem import Mol, Atom, EditableMol
from rdkit import Chem
from molzip import molzip

class MolEnumerator():

    def __init__(self, mol: Mol, rgroups: dict) -> None:
        
        self.mol = mol
        self.rgroups = rgroups
        self.mol_list = self.enumerate_markush_mol(self.mol, self.rgroups)

    def enumerate_markush_mol(self, mol: Mol, rgroups: dict) -> list:
        # Gets a list of all possible Mols for a Markush Mol and list of Rgroups
        
        # Get any R groups directly referenced by this Mol
        child_rgroups = {}
        rgroup_count = 0
        for atom in mol.GetAtoms():
            if atom.HasProp("_MolFileRLabel"):
                rlabel = atom.GetProp("_MolFileRLabel")
                child_rgroups[rlabel] = rgroups[int(rlabel)]
                rgroup_count += 1
        
        # For each child R group, get a list of all possible structures it could
        # be
        enumerated_rgroups = {}
        for r in child_rgroups.keys():

            enumerated_components = []
            rgroup = child_rgroups[r]

            for component in rgroup:
                enumerated_components += self.enumerate_markush_mol(
                    component, rgroups)

            enumerated_rgroups[r] = enumerated_components
        

        #Enumerate the varattachs using the built in RDKit functionality
        new_mol_list = self.enumerate_varattachs(mol)

        mol_list = new_mol_list

        #Enumerate all R groups
        for rlabel in child_rgroups.keys():
            rgroup = enumerated_rgroups[rlabel]
            new_mol_list = []
            for listmol in mol_list:
                new_mol_list += self.enumerate_rgroups(listmol, rlabel, rgroup)
            mol_list = new_mol_list

        # Enumerate all list atoms
        new_mol_list = []
        for mol in mol_list:
            enumerated_list = self.enumerate_listatoms(mol)
            new_mol_list += enumerated_list

        mol_list = new_mol_list

        # Sanitize molecules to remove any invalid molecules
        new_mol_list = []

        for mol in mol_list:
            try:
                Chem.rdmolops.SanitizeMol(mol)
                new_mol_list.append(mol)
            except:
                print("Skipping an invalid molecule")
        
        mol_list = new_mol_list

        return mol_list

    def enumerate_listatoms(self, mol: Mol) -> list:
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

                    
                    #For each element, replace the list atom with an atom of
                    # that element
                    for atomic_num in atomic_nums:
                        new_mol = Chem.Mol(mol)
                        Chem.rdmolops.Kekulize(new_mol, clearAromaticFlags=True)
                        listatom = new_mol.GetAtomWithIdx(listatom_idx)
                        new_atom = Atom(atomic_num)
                        edit_mol = EditableMol(new_mol)
                        new_atom_idx = edit_mol.AddAtom(new_atom)

                        for bond in listatom.GetBonds():

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
                        new_mol.UpdatePropertyCache()
                        Chem.rdmolops.SanitizeMol(new_mol)
                        # For this new Mol with the replaced listatom, 
                        # recursively replace further listatoms to ensure all
                        # listatoms are dealt with
                        new_list += self.enumerate_listatoms(new_mol)
                    return new_list

        # If there are no listatoms left in this Mol, return it as a single
        # element list           
        return [mol]
    
    def enumerate_rgroups(self, mol: Mol, rlabel: int, rgroup: list) -> list:
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

                    mol = self.mark_stereo(mol, atom.GetIdx())
                            
                    new_list = []


                    # For each component the R group could be, add the component to 
                    # the core and add this new Mol to the list of possibilities
                    for component in rgroup:

                        # Find the pseudoatom that represents the connection point

                        for comp_atom in component.GetAtoms():
                            if comp_atom.HasProp("isParentLinker"):
                                # Mark the atom which connects to the core
                                comp_atom.SetIntProp("molAtomMapNumber", 1)

                                component = self.mark_stereo(
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
                            new_mol = self.update_stereo(new_mol)
                            new_mol = self.cleanup_flags(new_mol)
                            
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

    def enumerate_varattachs(self, mol: Mol) -> list:

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
                    mol_copy = Chem.rdmolops.AddHs(mol)
                    edit_mol = EditableMol(mol_copy)
                    edit_mol.BeginBatchEdit()
                    end_idx = int(endpt) - 1

                    end_atom = mol_copy.GetAtomWithIdx(end_idx)

                    h_idx = None
                    for neighbor in end_atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 1 and h_idx == None:
                            h_idx = neighbor.GetIdx()

                    if h_idx != None:
                        edit_mol.AddBond(
                                    attach_idx,
                                    end_idx,
                                    order=Chem.rdchem.BondType.SINGLE
                                )
                        edit_mol.RemoveAtom(delete_idx)
                        edit_mol.RemoveAtom(h_idx)
                        edit_mol.CommitBatchEdit()
                        new_mol = edit_mol.GetMol()
                        new_mol = Chem.rdmolops.RemoveHs(new_mol)
                        new_mol_list += self.enumerate_varattachs(new_mol)

                return new_mol_list
            
        return [mol]
    
    def update_stereo(self, mol: Mol) -> Mol:

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
    
    def cleanup_flags(self, mol: Mol) -> Mol:

        for atom in mol.GetAtoms():
            atom.ClearProp("molAtomMapNumber")
            atom.ClearProp("bond_me")
            atom.ClearProp("delete_me")
            atom.ClearProp("fixed")

        return mol
    
    def mark_stereo(self, mol: Mol, idx: int) -> Mol:

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
    
    def get_inchi_list(self, tidy=True) -> list:

        # Gets the InChI for each mol in list of Mols
        # If the structure is invalid, the InChI is given as None

        inchi_list = []

        for mol in self.mol_list:
            try:
                inchi = Chem.MolToInchi(mol)
                inchi_list.append(inchi)
            except:
                inchi = None
                inchi_list.append(inchi)

        # Remove duplicates and sort list
        if tidy:
            inchi_list = list(set(inchi_list))
            inchi_list = sorted(inchi_list)

        return inchi_list
    
    def get_mol_list(self) -> list:

        return self.mol_list