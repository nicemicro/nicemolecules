# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 15:52:19 2022

@author: nicemicro
"""

from math import sqrt
from typing import Optional, Sequence
from enum import IntEnum, auto

class BondingError(IntEnum):
    OK = auto()
    BOND_SELF = auto()
    BOND_EXISTS = auto()
    INSUFF_EL_FIRST = auto()
    INSUFF_EL_OTHER = auto()
    HYPERVALENT_FIRST = auto()
    HYPERVALENT_OTHER = auto()

class Cov_bond:
    def __init__(self, atoms, electrons: list[int]):
        self.atoms: list[Atom] = atoms
        self.electrons: list[int] = electrons
    
    def describe(self) -> str:
        desc  = f"Bond order: {self.order()}\n"
        for atom, electron in zip(self.atoms, self.electrons):
            desc += f"  Atom symbol: {atom.symbol}"
            desc += f" at ({atom.coord_x}, {atom.coord_y}) with {electron} electrons\n"
        return desc
    
    def __str__(self) -> str:
        return super().__str__()
    
    def __repr__(self) -> str:
        return super().__repr__() + "\n" + self.describe()
    
    def order(self) -> float:
        return sum(self.electrons) / len(self.atoms)
    
    def dativity(self) -> float:
        return (self.electrons[0] - self.electrons[1]) / 2
    
    def electron_count(self) -> int:
        return sum(self.electrons)
    
    def atom_electrons(self, one_atom):
        assert isinstance(one_atom, Atom), "Only atoms are involved in a bond"
        if not one_atom in self.atoms:
            return None
        return self.electrons[self.atoms.index(one_atom)]
    
    def other_atoms(self, one_atom):
        assert isinstance(one_atom, Atom), "Only atoms are involved in a bond"
        if not one_atom in self.atoms:
            return None
        return list(atom for atom in self.atoms if atom != one_atom)
    
    def coords(self) -> tuple[float, float, float, float]:
        return (self.atoms[0].coord_x, self.atoms[0].coord_y,
                self.atoms[1].coord_x, self.atoms[1].coord_y)
    
    def length(self) -> float:
        if len(self.atoms) < 2:
            return 0
        xdiff: float = self.atoms[0].coord_x - self.atoms[1].coord_x
        ydiff: float = self.atoms[0].coord_y - self.atoms[1].coord_y
        return sqrt(xdiff ** 2 + ydiff ** 2)

class Atom:
    def __init__(self, symbol: str, valence_el: int, coords: Sequence[float],
                 charge: int=0, fullshell: int=8, hypervalent:bool=False):
        self.symbol = symbol
        self.valence = valence_el
        self.charge = 0
        self.fullshell = fullshell
        self.hypervalent = hypervalent
        self.bonds: list[Cov_bond] = []
        self.coord_x: float = coords[0]
        self.coord_y:float = coords[1]
    
    def describe(self, short: bool=False) -> str:
        desc: str  = f"{self.symbol} atom at ({self.coord_x}, {self.coord_y})"
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
        desc += f"  Lone pairs                 : {self.lone_pairs()}\n"
        return desc
    
    def __str__(self) -> str:
        return super().__str__()
    
    def __repr__(self) -> str:
        return super().__repr__() + "\n" + self.describe(True)
    
    def bond(self, other_atom, order: int=1, dative: int=0) -> Cov_bond:
        assert isinstance(other_atom, Atom), "Only can bond to an Atom instance"
        bond_err = self.can_bond(other_atom, order, dative)
        if bond_err != BondingError.OK:
            raise RuntimeError(f"Bonding to atom failed. Error: {bond_err}")
        new_bond: Cov_bond = Cov_bond([self, other_atom], [order-dative, order+dative])
        self.bonds.append(new_bond)
        other_atom.register_bond(new_bond)
        return new_bond

    def can_bond(self, other_atom, order: int=1, dative: int=0,
            check_other: bool=True) -> BondingError:
        assert isinstance(other_atom, Atom), "Only can bond to an Atom instance"
        # We want to check whether the other atom in the list can accommodate the
        # new bond unless check_other is False.
        other_err = BondingError.OK
        if check_other:
            other_err = other_atom.can_bond(self, order, -dative, False)
        # If there is any error from other atoms, we report that error directly.
        if other_err != BondingError.OK:
            return other_err
        if other_atom == self:
            return BondingError.BOND_SELF
        if self.is_bonded(other_atom):
            return BondingError.BOND_EXISTS
        # If there are less nonbonding electrons are left than what we need
        if self.nonbonding_el() < order-dative:
            if check_other:
                return BondingError.INSUFF_EL_FIRST
            else:
                return BondingError.INSUFF_EL_OTHER # This atom is NOT the "first"
        if not self.hypervalent and self.empty_valence() < order+dative:
            if check_other:
                return BondingError.HYPERVALENT_FIRST
            else:
                return BondingError.HYPERVALENT_OTHER # This atom is NOT the "first"
        return BondingError.OK
    
    def register_bond(self, new_bond: Cov_bond):
        if self in new_bond.atoms:
            self.bonds.append(new_bond)
        else:
            assert False, "Trying to register a bond that does not list current atom"
    
    def is_bonded(self, other_atom) -> bool:
        assert isinstance(other_atom, Atom), "Only Atom instance can be bonded to an Atom"
        bond_list_lists = list(bond.other_atoms(self) for bond in self.bonds)
        return other_atom in [item for sublist in bond_list_lists for
                              item in sublist]
    
    def bonded_atoms(self):
        bond_list_lists = list(bond.other_atoms(self) for bond in self.bonds)
        bond_list = [item for sublist in bond_list_lists for item in sublist]
        atomlist: list[Atom] = []
        for atom in bond_list:
            if not atom in atomlist:
                atomlist.append(atom)
        return atomlist
    
    def bond_instance(self, other_atom) -> Optional[Cov_bond]:
        assert isinstance(other_atom, Atom), "Only Atom instance can be bonded to an Atom"
        for bond in self.bonds:
            if other_atom in bond.other_atoms(self):
                return bond
        else:
            return None
    
    def electrons(self) -> int:
        from_bonds: int = 0
        for bond in self.bonds:
            from_bonds += bond.electron_count()
            from_bonds -= bond.atom_electrons(self)
        return from_bonds + self.valence
    
    def nonbonding_el(self) -> int:
        nonbonding: int = self.valence
        for bond in self.bonds:
            nonbonding -= bond.atom_electrons(self)
        return nonbonding
    
    def empty_valence(self) -> int:
        return self.fullshell - self.electrons()
    
    def radicals(self) -> int:
        radicals: int = min(self.empty_valence(), self.nonbonding_el())
        radicals = radicals % 2
        return radicals
    
    def lone_pairs(self) -> int:
        nonbonding: int = self.nonbonding_el() - self.radicals()
        return nonbonding // 2

def add_atom_by_symbol(symbol: str, coords: Sequence[float]) -> Optional[Atom] :
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

def find_molecule(one_atom: Atom) -> tuple:
    atomlist: list[Atom] = [one_atom]
    bondlist: list[Cov_bond] = []
    current_num: int = 0
    while current_num < len(atomlist):
        current_atom: Atom = atomlist[current_num]
        for bonded_atom in current_atom.bonded_atoms():
            if not bonded_atom in atomlist:
                atomlist.append(bonded_atom)
            bond_instance: Optional[Cov_bond]
            bond_instance = current_atom.bond_instance(bonded_atom)
            if not (bond_instance is None) and (not bond_instance in bondlist):
                bondlist.append(bond_instance)
        current_num += 1
    return atomlist, bondlist
