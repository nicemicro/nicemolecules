# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 15:52:19 2022

@author: nicemicro
"""

from math import sqrt 

class Cov_bond:
    def __init__(self, atoms, electrons):
        self.atoms = atoms
        self.electrons = electrons
    
    def describe(self):
        desc  = f"Bond order: {self.order()}\n"
        for atom, electron in zip(self.atoms, self.electrons):
            desc += f"  Atom symbol: {atom.symbol}"
            desc += f" at ({atom.coord_x}, {atom.coord_y}) with {electron} electrons\n"
        return desc
    
    def __str__(self):
        return super().__str__()
    
    def __repr__(self):
        return super().__repr__() + "\n" + self.describe()
    
    def order(self):
        return sum(self.electrons) / len(self.atoms)
    
    def dativity(self):
        return (self.electrons[0] - self.electrons[1]) / 2
    
    def electron_count(self):
        return sum(self.electrons)
    
    def atom_electrons(self, one_atom):
        if not one_atom in self.atoms:
            return None
        return self.electrons[self.atoms.index(one_atom)]
    
    def other_atoms(self, one_atom):
        if not one_atom in self.atoms:
            return False
        return list(atom for atom in self.atoms if atom != one_atom)
    
    def coords(self):
        return (self.atoms[0].coord_x, self.atoms[0].coord_y,
                self.atoms[1].coord_x, self.atoms[1].coord_y)
    
    def length(self):
        if len(self.atoms) < 2:
            return 0
        xdiff = self.atoms[0].coord_x - self.atoms[1].coord_x
        ydiff = self.atoms[0].coord_y - self.atoms[1].coord_y
        return sqrt(xdiff ** 2 + ydiff ** 2)

class Atom:
    def __init__(self, symbol, valence_el, coords, charge=0, fullshell=8,
                 hypervalent=False):
        self.symbol = symbol
        self.valence = valence_el
        self.charge = 0
        self.fullshell = fullshell
        self.hypervalent = hypervalent
        self.bonds = []
        self.coord_x = coords[0]
        self.coord_y = coords[1]
    
    def describe(self, short=False):
        desc  = f"{self.symbol} atom at ({self.coord_x}, {self.coord_y})"
        if short:
            return desc
        desc +=  "\n"
        desc += f"  Base valence electron number : {self.valence}\n"
        desc += f"  Base charge                  : {self.charge}\n"
        desc += f"  Electrons to fill shell      : {self.fullshell}\n"
        desc += f"  Can be hypervalent           : {self.hypervalent}\n"
        desc += f"  Number of covalent bonds     : {len(self.bonds)}\n"
        desc += f"  Total valence electron number: {self.electrons()}\n"
        for bond in self.bonds:
            desc += f"    {bond.order()} order bond with "
            for atom in bond.other_atoms(self):
                desc += f"{atom.symbol} at ({atom.coord_x}, {atom.coord_y})"
            desc += "\n"
        desc += f"  Non-bonding electrons      : {self.nonbonding_el()}\n"
        desc += f"  Unfilled valence orbitals  : {self.empty_valence()}\n"
        desc += f"  Radicals                   : {self.radicals()}\n"
        desc += f"  Non-bonding pairs          : {self.nonbonding_pairs()}\n"
        return desc
    
    def __str__(self):
        return super().__str__()
    
    def __repr__(self):
        return super().__repr__() + "\n" + self.describe(True)
    
    def bond(self, other_atom, order=1, dative=0):
        if not isinstance(other_atom, Atom):
            print("Only can bond to an Atom instance")
            return False
        if self.is_bonded(other_atom):
            assert other_atom.is_bonded(self), "Discrepancy in bond registration to atoms"
            print("Trying to bond with an atom already bonded with")
            return False
        assert not other_atom.is_bonded(self), "Discrepancy in bond registration to atoms"
        new_bond = Cov_bond([self, other_atom], [order-dative, order+dative])
        self.bonds.append(new_bond)
        other_atom.register_bond(new_bond)
        return new_bond
    
    def register_bond(self, new_bond):
        if self in new_bond.atoms:
            self.bonds.append(new_bond)
        else:
            assert False, "Trying to register a bond that does not work"
    
    def is_bonded(self, other_atom):
        bond_list_lists = list(bond.other_atoms(self) for bond in self.bonds)
        return other_atom in [item for sublist in bond_list_lists for
                              item in sublist]
    
    def bonded_atoms(self):
        bond_list_lists = list(bond.other_atoms(self) for bond in self.bonds)
        bond_list = [item for sublist in bond_list_lists for item in sublist]
        atomlist = []
        for atom in bond_list:
            if not atom in atomlist:
                atomlist.append(atom)
        return atomlist
    
    def bond_instance(self, other_atom):
        for bond in self.bonds:
            if other_atom in bond.other_atoms(self):
                return bond
        else:
            return False
    
    def electrons(self):
        from_bonds = 0
        for bond in self.bonds:
            from_bonds += bond.electron_count()
            from_bonds -= bond.atom_electrons(self)
        return from_bonds + self.valence
    
    def nonbonding_el(self):
        nonbonding = self.valence
        for bond in self.bonds:
            nonbonding -= bond.atom_electrons(self)
        return nonbonding
    
    def empty_valence(self):
        return self.fullshell - self.electrons()
    
    def radicals(self):
        radicals = min(self.empty_valence(), self.nonbonding_el())
        radicals = radicals % 2
        return radicals
    
    def nonbonding_pairs(self):
        nonbonding = self.nonbonding_el() - self.radicals()
        return nonbonding // 2

def add_atom_by_symbol(symbol, coords):
    if symbol == "H":
        return Atom("H", 1, coords, fullshell=2)
    if symbol == "B":
        return Atom("B", 3, coords)
    if symbol == "C":
        return Atom("C", 4, coords)
    if symbol == "N":
        return Atom("N", 5, coords)
    if symbol == "O":
        return Atom("O", 6, coords)
    if symbol == "F":
        return Atom("F", 7, coords)
    if symbol == "S":
        return Atom("S", 6, coords, hypervalent=True)
    return None

def find_molecule(one_atom):
    atomlist = [one_atom]
    bondlist = []
    current_num = 0
    while current_num < len(atomlist):
        current_atom = atomlist[current_num]
        for bonded_atom in current_atom.bonded_atoms():
            if not bonded_atom in atomlist:
                atomlist.append(bonded_atom)
            bond_instance = current_atom.bond_instance(bonded_atom)
            if not bond_instance in bondlist:
                bondlist.append(bond_instance)
        current_num += 1
    return atomlist, bondlist
