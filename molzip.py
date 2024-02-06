# Copyright (C) 2019 Brian P Kelley
#
# @@ All Rights Reserved @@
#
#  The contents are covered by the terms of the 3-Clause BSD license
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met: 
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.  
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution. 
#
# 3. Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
# AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
# WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE. 

"""Zip multiple molecules together.

Mols are in Matched Pair or RGroupDecomposition  form


  a= '[*:1]CCC'
  b= '[:1]CCC'

Stereo around chiral stereo centers is preserved

     >>> smi1 = '[C@H](Br)([*:1])F'
     >>> smi2 = '[*:1]N'
     >>> a = Chem.MolFromSmiles(smi1)
     >>> b = Chem.MolFromSmiles(smi2)
     >>> Chem.MolToSmiles(molzip(a, b, preserveChirality=True))
     N[C@@H](F)Br

If preserveChirality is False, then the bond order in the reassembled
atom will be used, sometimes with incorrect results

     >>> Chem.MolToSmiles(molzip(a, b, preserveChirality=False))
     N[C@H](F)Br
  
With thanks to Andrew Dalke and Brian Cole for discussions.

Andrew's blog post here:

 http://www.dalkescientific.com/writings/diary/
"""
from collections import defaultdict
from rdkit import Chem

def get_mappings(a):
    """Return a mapping of {atom_map_num: atom} for the molecule"""
    maps = {}
    for atom in a.GetAtoms():
        mapno = atom.GetAtomMapNum()
        if mapno and atom.GetAtomicNum() == 0:
            if mapno in maps:
                raise ValueError("duplicate atommaps for %s", mapno)
            maps[mapno] = atom
    return maps

def get_other_atom(atom):
    """Get the other atoms attached to this atom, but fail if there are more than one"""
    bonds = list(atom.GetBonds())
    if len(bonds) > 1:
        raise ValueError("Multiple bonds to mapped atom, can't currently handle fusions")
    bond = bonds[0]
    a = bond.GetBeginAtom()
    if a.GetIdx() == atom.GetIdx():
        return bond.GetEndAtom()
    return a

def countSwapsToInterConvert(orders1, orders2):
    nSwaps = 0
    ref_idx = 0
    probe_idx = 0
    while ref_idx < len(orders1):
        if orders1[ref_idx] != orders2[ref_idx]:
            try:
                idx = orders2.index(orders1[ref_idx])
            except:
                raise ValueError("orders1 and orders2 don't contain the same elements")
            orders2[probe_idx], orders2[idx] = orders2[idx], orders2[probe_idx]
            nSwaps += 1
        probe_idx += 1
        ref_idx += 1
    return nSwaps

def needs_inversion(orders1, orders2):
    return countSwapsToInterConvert(orders1, orders2) % 2 == 1

def bookmark_chirality(a, b, remap_chiral, deletions):      
    """Bookmark the chirality around atom a
    
     :param a: chiral atom
     :param b: future atom to be bonded
     :param remap_chiral: dict of the chiral atom: atom, orders to be preserved
     :params deletions: if we know the atom that is being deleted
                          and that atom that is being bonded, mark
                          the 'new' atom with the chiral order of the one
                          being deleted
    """
    tag = a.GetProp("chiral_tag")
    assert tag
    orders = get_orders(a)
    remap_chiral[a.GetIdx()] = [a, orders]
    # which atom is getting deleted
    for bond in a.GetBonds():
        oatom = bond.GetOtherAtom(a)
        if oatom.GetIdx() in deletions:
            order = oatom.GetIntProp(tag)
            b.SetIntProp(tag, order)

def mark_chiral_orders(atom):
    """Mark the chiral orders around this atom"""
    chiral_tag = "chiral_order_a_%s"%atom.GetIdx()
    atom.SetProp("chiral_tag", chiral_tag)
    for i,bond in enumerate(atom.GetBonds()):
        bond.GetOtherAtom(atom).SetIntProp(chiral_tag, i)

def get_orders(atom):
    """Return the marked chiralorders for this atom"""
    tag = atom.GetProp("chiral_tag")
    orders = []
    for bond in atom.GetBonds():
        oatom = bond.GetOtherAtom(atom)
        orders.append(oatom.GetIntProp(tag))
    return orders
            
def restore_chirality(atom, orders1):
    """Restore the chirality given the original marked orders"""
    orders2 = get_orders(atom)
    if needs_inversion(orders1, orders2):#not_rotation(orders, orders2):
        atom.InvertChirality()    
    # cleanup tags?
    
def molzip(a, b, preserveChirality=True):
    """zip two molecules together using their atom mappings and preserve any original chirality
    if possible.
    
    Molecule attachments must be using atom map numbers and dummy atoms.
    The atoms attached to these dummy atoms will be linked and the stereo in relationship to the dummy
    atom will be preserved.
    
    Example:
    
     >>> smi1 = '[C@H](Br)([*:1])F'
     >>> smi2 = '[*:1]N'
     >>> Chem.MolToSmiles(molzip(Chem.MolFromSmiles(smi1), Chem.MolFromSmiles(smi2), preserveChirality=False))
     N[C@H](F)Br
     >>> Chem.MolToSmiles(molzip(Chem.MolFromSmiles(smi1), Chem.MolFromSmiles(smi2), preserveChirality=True))
     N[C@@H](F)Br
    
    Now, is this right?
    
     >>> smi_manual = '[C@H](Br)(N)F'
     >>> Chem.MolToSmiles(Chem.MolFromSmiles(smi_manual))
     N[C@@H](F)Br
    
    Looks like it!
    """
    mapsa = get_mappings(a)
    mapsb = get_mappings(b)

    count = 0
    for mapno, atom_a in mapsa.items():
        if mapno in mapsb:
            atom_b = mapsb[mapno]

            other_a = get_other_atom(atom_a)
            other_b = get_other_atom(atom_b)
               
            if preserveChirality:
                if other_a.GetChiralTag():
                    mark_chiral_orders(other_a)

                if other_b.GetChiralTag():
                    mark_chiral_orders(other_b)

            # resolve bond order
            count += 1
            order = 1 # only single bonds for now, perhaps take max of existing bonds
            
            other_a.SetIntProp("bond_me", count)
            other_b.SetIntProp("bond_me", count)
            atom_a.SetIntProp("delete_me", count)
            atom_b.SetIntProp("delete_me", count)

    actions = Chem.RWMol(Chem.CombineMols(a, b))       
    if count:
        bonds = defaultdict(list)
        deletions = {}
        for atom in actions.GetAtoms():
            if atom.HasProp("bond_me"):
                bonds[atom.GetIntProp("bond_me")].append(atom)
            elif atom.HasProp("delete_me"):
                deletions[atom.GetIdx()] = atom

        remap_chiral = {}
        for count, (a,b) in bonds.items():
            orders = None
            if a.HasProp("chiral_tag"):
                bookmark_chirality(a, b, remap_chiral, deletions)
            if b.HasProp("chiral_tag"):
                bookmark_chirality(b, a, remap_chiral, deletions)
                    
            actions.AddBond(a.GetIdx(), b.GetIdx(), Chem.BondType.SINGLE)

            
        for atom in deletions.values():
            actions.RemoveAtom(atom.GetIdx())

        for idx, (atom, orders) in remap_chiral.items():
            restore_chirality(atom, orders)                
        actions.UpdatePropertyCache()
    return actions