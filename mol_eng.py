# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 15:52:19 2022

@author: nicemicro
"""

from math import sqrt, atan2, pi
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
    MISC_ERROR = auto()


class CovBond:
    """
    Representing a covalent bond between atoms.
    """

    def get_order(self) -> float:
        """Returns the bond order of this bond."""
        if not self.atoms:
            return 0
        return sum(self.electrons) / len(self.atoms)

    def set_order(self, new_order: int) -> None:
        """Changes the bond order."""
        new_electrons = self.electrons_calc(new_order, self.dativity)
        assert new_electrons is not None
        self.set_electrons(new_electrons)

    def get_dativity(self) -> float:
        """Returns the dativity of this bond."""
        return (self.electrons[0] - self.electrons[1]) / 2

    def set_dativity(self, new_dativity: int) -> None:
        """Changes the bond order."""
        new_electrons = self.electrons_calc(self.order, new_dativity)
        assert new_electrons is not None
        self.set_electrons(new_electrons)

    def get_electron_count(self) -> int:
        """Returns the number of electrons participating in the formation of
        this bond."""
        return sum(self.electrons)

    def get_coords(self) -> tuple[float, float, float, float]:
        """Returns the coordinates of the bond."""
        return (
            self.atoms[0].coord_x,
            self.atoms[0].coord_y,
            self.atoms[1].coord_x,
            self.atoms[1].coord_y,
        )

    def get_length(self) -> float:
        """Returns the length of the bond."""
        if len(self.atoms) < 2:
            return 0
        xdiff: float = self.atoms[0].coord_x - self.atoms[1].coord_x
        ydiff: float = self.atoms[0].coord_y - self.atoms[1].coord_y
        return sqrt(xdiff**2 + ydiff**2)

    def get_atoms(self):
        """Returns the atoms constituting this bond."""
        return self._atoms.copy()

    def get_electrons(self) -> list[int]:
        """Returns the number of electrons constituting this bond."""
        return self._electrons.copy()

    def set_electrons(self, new_electrons: list[int]) -> None:
        """Sets the electron configuration of the bond."""
        err = self.can_change_electrons(new_electrons)
        if err != BondingError.OK:
            raise RuntimeError(f"Changing bond configuration failed. Error: {err}")
        self._electrons = new_electrons.copy()

    def toss(self, new_value) -> None:
        """Prohibits changing attributes."""
        raise AttributeError("Can't directly modify attribute")

    order = property(get_order, set_order)
    dativity = property(get_dativity, set_dativity)
    electron_count = property(get_electron_count, toss)
    coords = property(get_coords, toss)
    length = property(get_length, toss)
    atoms = property(get_atoms, toss)
    electrons = property(get_electrons, set_electrons)

    def __init__(
        self, atoms, electrons: Optional[list[int]] = None, order: int = 1, dative=0
    ):
        self._atoms: list[Atom] = atoms
        if electrons is None:
            electrons = self.electrons_calc(order, dative)
        assert electrons is not None
        self._electrons = electrons

    def describe(self) -> str:
        """Returns a string description of the bond."""
        desc = f"Bond order: {self.order}\n"
        for atom, electron in zip(self.atoms, self.electrons):
            desc += f"  Atom symbol: {atom.symbol}"
            desc += f" at ({atom.coord_x}, {atom.coord_y}) with {electron} electrons\n"
        return desc

    def __repr__(self) -> str:
        return super().__repr__() + "\n" + self.describe()

    def electrons_calc(
        self, order: float = 1, dative: float = 0
    ) -> Optional[list[int]]:
        """Converts order and dativity values to electron configuration"""
        if order < 1:
            return None
        return [int(order + dative), int(order - dative)]

    def can_change_dativity(self, new_dativity: int) -> BondingError:
        """Checks whether the suggested dativity can be changed to for this
        bond."""
        new_electrons = self.electrons_calc(self.order, new_dativity)
        # print(new_electrons)
        if new_electrons is None:
            return BondingError.MISC_ERROR
        return self.can_change_electrons(new_electrons)

    def can_change_order(self, new_order: int) -> BondingError:
        """Checks whether the suggested order can be changed to for this
        bond."""
        if new_order < 1 or len(self.atoms) < 2:
            return BondingError.MISC_ERROR
        new_electrons = self.electrons_calc(new_order, self.dativity)
        if new_electrons is None:
            return BondingError.MISC_ERROR
        return self.can_change_electrons(new_electrons)

    def can_change_electrons(self, new_electrons: list[int]) -> BondingError:
        """Checks whether the suggested electron configuration for this
        bond is possible or not."""
        if len(new_electrons) != len(self.electrons):
            return BondingError.MISC_ERROR
        additional_all = sum(new_electrons) - self.electron_count
        # print(f"{self.electrons} -> {new_electrons}")
        # print(f"change: {additional_all}")
        for old_el, new_el, atom in zip(self.electrons, new_electrons, self.atoms):
            additional = new_el - old_el
            # print(
            #    f"  additional {additional} electrons to contribute, has {atom.nonbonding_el} electrons"
            # )
            if additional > atom.nonbonding_el:
                if atom == self.atoms[0]:
                    return BondingError.INSUFF_EL_FIRST
                return BondingError.INSUFF_EL_OTHER
            # print(
            #    f"  {additional_all-additional} electrons to receive, {atom.empty_valence} empty places"
            # )
            if (
                not atom.hypervalent
                and additional_all - additional > atom.empty_valence
            ):
                if atom == self.atoms[0]:
                    return BondingError.HYPERVALENT_FIRST
                return BondingError.HYPERVALENT_OTHER
        return BondingError.OK

    def atom_electrons(self, one_atom):
        """Returns the number of electrons from a certain atom that was donated
        to this bond."""
        assert isinstance(one_atom, Atom), "Only atoms are involved in a bond."
        if one_atom not in self.atoms:
            return None
        return self.electrons[self.atoms.index(one_atom)]

    def other_atoms(self, one_atom):
        """Returns a list of other Atom instances, except the one received as
        parameter."""
        assert isinstance(one_atom, Atom), "Only atoms are involved in a bond."
        if one_atom not in self.atoms:
            return None
        return list(atom for atom in self.atoms if atom != one_atom)

    def delete(self):
        """Removes itself from the atoms it was bonded to"""
        for connected_atom in self.atoms:
            connected_atom.remove_bond(self)
        self._atoms = []


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

    def set_charge(self, new_charge: int) -> None:
        err: BondingError = self.can_ionize(new_charge)
        if err != BondingError.OK:
            raise RuntimeError(f"Bonding to atom failed. Error: {err}")
        self._charge = new_charge

    def get_bonds(self) -> list[CovBond]:
        return self._bonds.copy()

    def get_bonded_atoms(self):
        bond_list_lists = list(bond.other_atoms(self) for bond in self._bonds)
        bond_list = [item for sublist in bond_list_lists for item in sublist]
        atomlist: list[Atom] = []
        for atom in bond_list:
            if not atom in atomlist:
                atomlist.append(atom)
        return atomlist

    def get_electrons(self) -> int:
        from_bonds: int = 0
        for bond in self._bonds:
            from_bonds += bond.electron_count
            from_bonds -= bond.atom_electrons(self)
        return from_bonds + self.valence

    def get_nonbonding_el(self) -> int:
        nonbonding: int = self.valence
        for bond in self._bonds:
            nonbonding -= bond.atom_electrons(self)
        return nonbonding

    def get_empty_valence(self) -> int:
        return int(self.fullshell - self.get_electrons())

    def get_radicals(self) -> int:
        radicals: int = min(self.get_empty_valence(), self.get_nonbonding_el())
        if radicals < 0:
            radicals = radicals % 2
        return radicals

    def get_lone_pairs(self) -> int:
        nonbonding: int = self.get_nonbonding_el() - self.get_radicals()
        return int(nonbonding // 2)

    def get_bond_angles(self) -> list[float]:
        """Returns the list of the angles of bonds relative to the x axis."""
        angle_list: list[float] = []
        for bond in self._bonds:
            xdiff: float = bond.other_atoms(self)[0].coord_x - self.coord_x
            ydiff: float = bond.other_atoms(self)[0].coord_y - self.coord_y
            angle: float = atan2(ydiff, xdiff)
            if angle < 0:
                angle = 2 * pi + angle
            angle_list.append(angle)
        return angle_list

    def get_electron_angles(self) -> list[float]:
        """Returns the list of the suggested angles of radicals and lone pairs."""
        angles: list[float] = self.get_bond_angles()
        angles.sort()
        el_num: int = self.get_lone_pairs() + self.get_radicals()
        el_angles: list[float] = []
        rel_ang: list[float] = relative_angles(angles)
        angle: float
        diff: float
        start: float
        repeat: int
        while el_num > 0:
            if len(angles) <= 1:
                diff = 2 * pi / (el_num + len(angles))
                start = 0
                repeat = el_num
                if len(angles) == 1:
                    start = angles[0] + diff
            else:
                largest: float = max(rel_ang)
                largest_loc: int = rel_ang.index(largest)
                start = angles[largest_loc]
                if largest_loc == len(rel_ang) - 1:
                    diff = angles[0] + 2 * pi - start
                else:
                    diff = angles[largest_loc + 1] - start
                rel_ang2 = rel_ang.copy()
                rel_ang2.pop(largest_loc)
                section: int = int(round(largest / max(rel_ang2) + 0.2))
                repeat = min(section, el_num)
                diff = diff / (repeat + 1)
                start += diff
            for index in range(repeat):
                angle = start + index * diff
                while angle > 2 * pi:
                    angle -= 2 * pi
                el_angles.append(angle)
                angles.append(angle)
                el_num -= 1
            el_angles.sort()
            angles.sort()
            rel_ang = relative_angles(angles)
        return el_angles

    def toss(self, new_value):
        """Throws an error if unmodifiable things are trying to get
        modified."""
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
    charge = property(get_charge, set_charge)
    bonds = property(get_bonds, toss)
    bonded_atoms = property(get_bonded_atoms, toss)
    electrons = property(get_electrons, toss)
    nonbonding_el = property(get_nonbonding_el, toss)
    empty_valence = property(get_empty_valence, toss)
    radicals = property(get_radicals, toss)
    lone_pairs = property(get_lone_pairs, toss)
    bond_angles = property(get_bond_angles, toss)
    electron_angles = property(get_electron_angles, toss)

    def __init__(self, element: el.Element, coords: Sequence[float], charge: int = 0):
        self._element = element
        self.charge = charge
        self._bonds = []
        self.coord_x = coords[0]
        self.coord_y = coords[1]

    def describe(self, short: bool = False) -> str:
        """Describes the properties of the atom in human-readable form"""
        charge_str: str = ""
        if self.charge > 0:
            charge_str = self.charge * "+"
        elif self.charge < 0:
            charge_str = (-self.charge) * "-"
        desc: str = f"{self.symbol}{charge_str} atom at "
        desc += f"({self.coord_x}, {self.coord_y})"
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
            desc += f"    {bond.order} order bond with "
            for atom in bond.other_atoms(self):
                desc += f"{atom.symbol} at ({atom.coord_x}, {atom.coord_y})"
            desc += "\n"
        desc += f"  Non-bonding electrons        : {self.nonbonding_el}\n"
        desc += f"  Unfilled valence orbitals    : {self.empty_valence}\n"
        desc += f"  Radicals                     : {self.radicals}\n"
        desc += f"  Lone pairs                   : {self.lone_pairs}\n"
        return desc

    def __repr__(self) -> str:
        return super().__repr__() + "\n" + self.describe(True)

    def can_change_element(self, element: el.Element) -> BondingError:
        """Checks whether the atom can be switched to be a different
        element."""
        in_bonds: int = 0
        donate_bond: int = 0
        for bond in self._bonds:
            in_bonds += bond.electron_count
            donate_bond += bond.atom_electrons(self)
        total_prosp_el = in_bonds - donate_bond + element.valence_el - self.charge
        if not element.hypervalent and total_prosp_el > element.fullshell:
            return BondingError.HYPERVALENT_FIRST
        if donate_bond + self.charge > element.valence_el:
            return BondingError.INSUFF_EL_FIRST
        return BondingError.OK

    def change_element(self, new_element: el.Element) -> None:
        """Changes the element of the atom to a new one."""
        err: BondingError = self.can_change_element(new_element)
        if err != BondingError.OK:
            raise RuntimeError(f"Bonding to atom failed. Error: {err}")
        self._element = new_element

    def can_ionize(self, charge: int = 0) -> BondingError:
        """Checks whether the charge of the atom can be changed to the
        specified value."""
        if self._element.valence_el - charge < 0:
            return BondingError.INSUFF_EL_FIRST
        if self._element.valence_el - charge > self.fullshell:
            return BondingError.HYPERVALENT_FIRST
        return BondingError.OK

    def bond(self, other_atom, order: int = 1, dative: int = 0) -> CovBond:
        """
        Bonds the current atom to an other atom with the specified bond order
        and dativity (dative/coordinative bond).
        """
        assert isinstance(other_atom, Atom), "Only can bond to an Atom instance"
        bond_err = self.can_bond(other_atom, order, dative)
        if bond_err != BondingError.OK:
            raise RuntimeError(f"Bonding to atom failed. Error: {bond_err}")
        new_bond: CovBond = CovBond([self, other_atom], order=order, dative=dative)
        self._bonds.append(new_bond)
        other_atom.register_bond(new_bond)
        return new_bond

    def can_bond(
        self, other_atom, order: int = 1, dative: int = 0, check_other: bool = True
    ) -> BondingError:
        """
        A check whether a bond can be created between the current atom and an
        other atom.
        """
        assert isinstance(other_atom, Atom), "Only can bond to an Atom instance"
        # We want to check whether the other atom in the list can accommodate
        # the new bond unless check_other is False.
        other_err = BondingError.OK
        if check_other:
            other_err = other_atom.can_bond(self, order, -dative, False)
        # If there is any error from other atoms, we report that error directly
        if other_err != BondingError.OK:
            return other_err
        if other_atom == self:
            return BondingError.BOND_SELF
        if self.is_bonded(other_atom):
            return BondingError.BOND_EXISTS
        # If there are less nonbonding electrons are left than what we need
        if self.nonbonding_el < order + dative:
            if check_other:
                return BondingError.INSUFF_EL_FIRST
            return BondingError.INSUFF_EL_OTHER  # This atom is NOT the "first"
        if not self.hypervalent and self.empty_valence < order - dative:
            if check_other:
                return BondingError.HYPERVALENT_FIRST
            return BondingError.HYPERVALENT_OTHER  # This atom is NOT the "first"
        return BondingError.OK

    def remove_bond(self, bond_instance: CovBond) -> None:
        """Removes a bond instance registered to this atom as a part of deleting
        a bond."""
        assert (
            bond_instance in self._bonds
        ), "Bond trying to be deleted is not registered to this atom"
        self._bonds.remove(bond_instance)

    def register_bond(self, new_bond: CovBond) -> None:
        """
        Register a CovBond instance to this atom, created by the can_bond
        function of an other Atom instance.
        """
        if self in new_bond.atoms:
            self._bonds.append(new_bond)
        else:
            assert False, "Trying to register a bond that does not list current atom"

    def is_bonded(self, other_atom) -> bool:
        """
        Checks whether a bond between this, and an other Atom instance exists.
        """
        assert isinstance(
            other_atom, Atom
        ), "Only Atom instance can be bonded to an Atom"
        bond_list_lists = list(bond.other_atoms(self) for bond in self._bonds)
        return other_atom in [item for sublist in bond_list_lists for item in sublist]

    def bond_instance(self, other_atom) -> Optional[CovBond]:
        """Returns the bond instance connecting the current atom to the other_atom.
        If there is no bond, returns None."""
        assert isinstance(
            other_atom, Atom
        ), "Only Atom instance can be bonded to an Atom"
        for bond in self._bonds:
            if other_atom in bond.other_atoms(self):
                return bond
        return None

    def unbond_all(self) -> list[CovBond]:
        """Removes all bonds from the atom. Returns the instances of the now
        unused bonds."""
        bonds_loop: list[CovBond] = self._bonds.copy()
        for bond in bonds_loop:
            bond.delete()
        return bonds_loop


def element_by_symbol(symbol: str) -> Optional[el.Element]:
    """Returns the corresponding element instance to the element symbol."""
    element_dict = {element.symbol: element for element in el.element_table}
    if symbol not in element_dict:
        return None
    return element_dict[symbol]


def add_atom_by_symbol(symbol: str, coords: Sequence[float]) -> Optional[Atom]:
    """Returns a new atom instance specified by its symbol and coordinates."""
    sel_el: Optional[el.Element] = element_by_symbol(symbol)
    if sel_el is None:
        return None
    assert isinstance(sel_el, el.Element)
    return Atom(sel_el, coords)


def relative_angles(angle_list: list[float]) -> list[float]:
    relative_angl: list[float] = []
    if len(angle_list) > 0:
        for a, b in zip(angle_list, angle_list[1:] + [angle_list[0] + 2 * pi]):
            relative_angl.append(b - a)
    return relative_angl


def find_molecule(one_atom: Atom) -> tuple:
    """Finds all CovBond and Atom instances that are somehow connected to the
    specified Atom instance, and hence form a molecule."""
    atomlist: list[Atom] = [one_atom]
    bondlist: list[CovBond] = []
    current_num: int = 0
    while current_num < len(atomlist):
        current_atom: Atom = atomlist[current_num]
        for bonded_atom in current_atom.bonded_atoms():
            if bonded_atom not in atomlist:
                atomlist.append(bonded_atom)
            bond_instance: Optional[CovBond]
            bond_instance = current_atom.bond_instance(bonded_atom)
            if not (bond_instance is None) and (bond_instance not in bondlist):
                bondlist.append(bond_instance)
        current_num += 1
    return atomlist, bondlist
