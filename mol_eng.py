# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 15:52:19 2022

@author: nicemicro
"""

from math import sqrt
from typing import Optional, Sequence
from enum import IntEnum, auto
import elements as el


class BondingError(IntEnum):
    """Reasons a bond formation might be impossible."""
    OK = auto()
    BOND_SELF = auto()
    BOND_EXISTS = auto()
    INSUFF_EL_FIRST = auto()
    INSUFF_EL_OTHER = auto()
    HYPERVALENT_FIRST = auto()
    HYPERVALENT_OTHER = auto()


class CovBond:
    """
    Representing a covalent bond between atoms.
    """

    def __init__(self, atoms, electrons: list[int]):
        self.atoms: list[Atom] = atoms
        self.electrons: list[int] = electrons

    def describe(self) -> str:
        desc = f"Bond order: {self.order()}\n"
        for atom, electron in zip(self.atoms, self.electrons):
            desc += f"  Atom symbol: {atom.symbol}"
            desc += f" at ({atom.coord_x}, {atom.coord_y}) with {electron} electrons\n"
        return desc

    def __str__(self) -> str:
        return super().__str__()

    def __repr__(self) -> str:
        return super().__repr__() + "\n" + self.describe()

    def order(self) -> float:
        """Returns the bond order of this bond."""
        return sum(self.electrons) / len(self.atoms)

    def dativity(self) -> float:
        """Returns the dativity of this bond."""
        return (self.electrons[0] - self.electrons[1]) / 2

    def electron_count(self) -> int:
        """Returns the number of electrons participating in the formation of this bond."""
        return sum(self.electrons)

    def atom_electrons(self, one_atom):
        """Returns the number of electrons from a certain atom that was donated to this bond."""
        assert isinstance(one_atom, Atom), "Only atoms are involved in a bond."
        if not one_atom in self.atoms:
            return None
        return self.electrons[self.atoms.index(one_atom)]

    def other_atoms(self, one_atom):
        """Returns a list of other Atom instances, except the one received as parameter."""
        assert isinstance(one_atom, Atom), "Only atoms are involved in a bond."
        if not one_atom in self.atoms:
            return None
        return list(atom for atom in self.atoms if atom != one_atom)

    def coords(self) -> tuple[float, float, float, float]:
        """returns the coordinates of the bond."""
        return (
            self.atoms[0].coord_x,
            self.atoms[0].coord_y,
            self.atoms[1].coord_x,
            self.atoms[1].coord_y,
        )

    def length(self) -> float:
        if len(self.atoms) < 2:
            return 0
        xdiff: float = self.atoms[0].coord_x - self.atoms[1].coord_x
        ydiff: float = self.atoms[0].coord_y - self.atoms[1].coord_y
        return sqrt(xdiff**2 + ydiff**2)


class Atom:
    """
    Represents an atom symbol in 2D-space in a structural formula.
    """

    def get_symbol(self) -> str:
        return self._element.symbol

    def get_valence(self) -> int:
        return self._element.valence_el - self._charge

    def get_fullshell(self) -> int:
        return self._element.fullshell

    def get_hypervalent(self) -> bool:
        return self._element.hypervalent

    def get_charge(self) -> int:
        return self._charge

    def get_bonds(self) -> list[CovBond]:
        return self._bonds

    def get_bonded_atoms(self):
        bond_list_lists = list(bond.other_atoms(self) for bond in self.bonds)
        bond_list = [item for sublist in bond_list_lists for item in sublist]
        atomlist: list[Atom] = []
        for atom in bond_list:
            if not atom in atomlist:
                atomlist.append(atom)
        return atomlist

    def get_electrons(self) -> int:
        from_bonds: int = 0
        for bond in self.bonds:
            from_bonds += bond.electron_count()
            from_bonds -= bond.atom_electrons(self)
        return from_bonds + self.valence

    def get_nonbonding_el(self) -> int:
        nonbonding: int = self.valence
        for bond in self.bonds:
            nonbonding -= bond.atom_electrons(self)
        return nonbonding

    def get_empty_valence(self) -> int:
        return self.fullshell - self.get_electrons()

    def get_radicals(self) -> int:
        radicals: int = min(self.get_empty_valence(), self.get_nonbonding_el())
        radicals = radicals % 2
        return radicals

    def get_lone_pairs(self) -> int:
        nonbonding: int = self.get_nonbonding_el() - self.get_radicals()
        return nonbonding // 2

    def toss(self, new_value):
        raise AttributeError("Can't directly modify attribute")

    _element: el.Element
    _charge: int
    _bonds: list[CovBond]
    coord_x: float
    coord_y: float
    symbol = property(get_symbol, toss)
    valence = property(get_valence, toss)
    fullshell = property(get_fullshell, toss)
    hypervalent = property(get_hypervalent, toss)
    charge = property(get_charge, toss)
    bonds = property(get_bonds, toss)
    bonded_atoms = property(get_bonded_atoms, toss)
    electrons = property(get_electrons, toss)
    nonbonding_el = property(get_nonbonding_el, toss)
    empty_valence = property(get_empty_valence, toss)
    radicals = property(get_radicals, toss)
    lone_pairs = property(get_lone_pairs, toss)

    def __init__(self,
                 element: el.Element,
                 coords: Sequence[float],
                 charge: int = 0):
        self._element = element
        self._charge = charge
        self._bonds = []
        self.coord_x = coords[0]
        self.coord_y = coords[1]

    def describe(self, short: bool = False) -> str:
        desc: str = f"{self.symbol} atom at ({self.coord_x}, {self.coord_y})"
        if short:
            return desc
        desc += "\n"
        desc += f"  Base valence electron number : {self.valence}\n"
        desc += f"  Base charge                  : {self.charge}\n"
        desc += f"  Electrons to fill shell      : {self.fullshell}\n"
        desc += f"  Can be hypervalent           : {self.hypervalent}\n"
        desc += f"  Total valence electron number: {self.electrons}\n"
        desc += f"  Number of covalent bonds     : {len(self.bonds)}\n"
        for bond in self.bonds:
            desc += f"    {bond.order()} order bond with "
            for atom in bond.other_atoms(self):
                desc += f"{atom.symbol} at ({atom.coord_x}, {atom.coord_y})"
            desc += "\n"
        desc += f"  Non-bonding electrons        : {self.nonbonding_el}\n"
        desc += f"  Unfilled valence orbitals    : {self.empty_valence}\n"
        desc += f"  Radicals                     : {self.radicals}\n"
        desc += f"  Lone pairs                   : {self.lone_pairs}\n"
        return desc

    def __str__(self) -> str:
        return super().__str__()

    def __repr__(self) -> str:
        return super().__repr__() + "\n" + self.describe(True)

    def bond(self, other_atom, order: int = 1, dative: int = 0) -> CovBond:
        """
        Bonds the current atom to an other atom with the specified bond order and
        dativity (dative/coordinative bond).
        """
        assert isinstance(other_atom,
                          Atom), "Only can bond to an Atom instance"
        bond_err = self.can_bond(other_atom, order, dative)
        if bond_err != BondingError.OK:
            raise RuntimeError(f"Bonding to atom failed. Error: {bond_err}")
        new_bond: CovBond = CovBond([self, other_atom],
                                    [order - dative, order + dative])
        self.bonds.append(new_bond)
        other_atom.register_bond(new_bond)
        return new_bond

    def can_bond(self,
                 other_atom,
                 order: int = 1,
                 dative: int = 0,
                 check_other: bool = True) -> BondingError:
        """
        A check whether a bond can be created between the current atom and an other atom.
        """
        assert isinstance(other_atom,
                          Atom), "Only can bond to an Atom instance"
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
        if self.nonbonding_el < order - dative:
            if check_other:
                return BondingError.INSUFF_EL_FIRST
            else:
                return BondingError.INSUFF_EL_OTHER  # This atom is NOT the "first"
        if not self.hypervalent and self.empty_valence < order + dative:
            if check_other:
                return BondingError.HYPERVALENT_FIRST
            else:
                return BondingError.HYPERVALENT_OTHER  # This atom is NOT the "first"
        return BondingError.OK

    def register_bond(self, new_bond: CovBond):
        """
        Register a CovBond instance to this atom, created by the can_bond function
        of an other Atom instance.
        """
        if self in new_bond.atoms:
            self.bonds.append(new_bond)
        else:
            assert False, "Trying to register a bond that does not list current atom"

    def is_bonded(self, other_atom) -> bool:
        """
        Checks whether a bond between this, and an other Atom instance exists.
        """
        assert isinstance(other_atom,
                          Atom), "Only Atom instance can be bonded to an Atom"
        bond_list_lists = list(bond.other_atoms(self) for bond in self.bonds)
        return other_atom in [
            item for sublist in bond_list_lists for item in sublist
        ]

    def bond_instance(self, other_atom) -> Optional[CovBond]:
        assert isinstance(other_atom,
                          Atom), "Only Atom instance can be bonded to an Atom"
        for bond in self.bonds:
            if other_atom in bond.other_atoms(self):
                return bond
        else:
            return None


def add_atom_by_symbol(symbol: str, coords: Sequence[float]) -> Optional[Atom]:
    element_dict = {element.symbol: element for element in el.element_table}
    if not symbol in element_dict:
        return None
    sel_el = element_dict[symbol]
    return Atom(sel_el, coords)


def find_molecule(one_atom: Atom) -> tuple:
    atomlist: list[Atom] = [one_atom]
    bondlist: list[CovBond] = []
    current_num: int = 0
    while current_num < len(atomlist):
        current_atom: Atom = atomlist[current_num]
        for bonded_atom in current_atom.bonded_atoms():
            if not bonded_atom in atomlist:
                atomlist.append(bonded_atom)
            bond_instance: Optional[CovBond]
            bond_instance = current_atom.bond_instance(bonded_atom)
            if not (bond_instance is None) and (not bond_instance in bondlist):
                bondlist.append(bond_instance)
        current_num += 1
    return atomlist, bondlist
