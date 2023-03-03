# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 15:52:19 2022

@author: nicemicro
"""

from math import sqrt, atan2, cos, sin, pi
from typing import Optional, Any, Union
from enum import IntEnum, auto
from xml.etree import ElementTree as ET
import elements as el


VERSION: str = "0.0.0"


class BondingError(IntEnum):
    """Reasons a bond formation might be impossible."""

    OK = auto()
    BOND_SELF = auto()
    BOND_EXISTS = auto()
    IMPROPER_BOND = auto()
    INSUFF_EL_FIRST = auto()
    INSUFF_EL_OTHER = auto()
    HYPERVALENT_FIRST = auto()
    HYPERVALENT_OTHER = auto()
    MISC_ERROR = auto()


class Bond:
    """
    Representing any type of bond between any number of atoms.
    """

    def get_atoms(self):
        """Returns the atoms constituting this bond."""
        return self._atoms.copy()

    def get_electrons(self) -> list[int]:
        """Returns the number of electrons constituting this bond per atom."""
        return self._electrons.copy()

    def get_electron_count(self) -> int:
        """Returns the number of electrons participating in the formation of
        this bond."""
        return sum(self.electrons)

    def toss(self, new_value) -> None:
        """Prohibits changing attributes."""
        raise AttributeError("Can't directly modify attribute")

    atoms = property(get_atoms, toss)
    electrons = property(get_electrons, toss)
    electron_count = property(get_electron_count, toss)

    def __init__(
        self, atoms, electrons: Optional[list[int]] = None
    ):
        assert isinstance(atoms, list)
        self._atoms: list[Atom] = atoms
        self._electrons: list[int]
        if electrons is None:
            self._electrons = [1] * len(atoms)
        else:
            if len(electrons) != len(atoms):
                raise RuntimeError("Atom and electon list has to be the same length")
            self._electrons = electrons

    def describe(self) -> str:
        """Returns a string description of the bond."""
        desc: str = ""
        for atom, electron in zip(self.atoms, self.electrons):
            desc += f"  Atom symbol: {atom.symbol}"
            desc += f" at ({atom.coord_x}, {atom.coord_y}) with {electron} electrons\n"
        return desc

    def __repr__(self) -> str:
        return super().__repr__() + "\n" + self.describe()

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


class CovBond(Bond):
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
        return atomdist(self.atoms[0], self.atoms[1])

    def set_electrons(self, new_electrons: list[int]) -> None:
        """Sets the electron configuration of the bond."""
        err = self.can_change_electrons(new_electrons)
        if err != BondingError.OK:
            raise RuntimeError(
                f"Changing bond configuration failed. Error: {err}"
            )
        self._electrons = new_electrons.copy()

    def get_electrons(self) -> list[int]:
        """Returns the number of electrons constituting this bond."""
        return super().get_electrons()

    def toss(self, new_value) -> None:
        """Prohibits changing attributes."""
        raise AttributeError("Can't directly modify attribute")

    order = property(get_order, set_order)
    dativity = property(get_dativity, set_dativity)
    coords = property(get_coords, toss)
    length = property(get_length, toss)
    electrons = property(get_electrons, set_electrons)

    def __init__(
        self, atoms, electrons: Optional[list[int]] = None, order: int = 1, dative: int = 0
    ):
        if electrons is None:
            electrons = self.electrons_calc(order, dative)
        assert electrons is not None
        super().__init__(atoms, electrons)

    def describe(self) -> str:
        desc: str = f"Bond order: {self.order}\n"
        desc += f"Bond dativity: {self.dativity}\n"
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


class Atom:
    """
    Represents an atom symbol in 2D-space in a structural formula.
    """

    def get_coords(self) -> tuple[float, float]:
        return (self.coord_x, self.coord_y)

    def set_coords(self, new_coords: tuple[float, float]) -> None:
        (self.coord_x, self.coord_y) = new_coords

    def get_symbol(self) -> str:
        return self._element.symbol

    def get_element(self) -> el.Element:
        return self._element

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

    def get_bonds(self) -> list[Bond]:
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
        return self.get_els_from_bonds() + self.valence

    def get_nonbonding_el(self) -> int:
        nonbonding: int = self.valence
        for bond in self._bonds:
            nonbonding -= bond.atom_electrons(self)
        return nonbonding

    def get_empty_valence(self) -> int:
        return int(self.fullshell - self.get_electrons())

    def get_els_from_bonds(self) -> int:
        from_bonds: int = 0
        for bond in self._bonds:
            from_bonds += bond.electron_count
            from_bonds -= bond.atom_electrons(self)
        return from_bonds

    def get_els_to_bonds(self) -> int:
        into_bonds: int = 0
        for bond in self._bonds:
            into_bonds += bond.atom_electrons(self)
        return into_bonds

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
            angle_list.append(atom_angle(self, bond.other_atoms(self)[0]))
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
    _bonds: list[Bond]
    coord_x: float
    coord_y: float
    coord_z: int
    coords = property(get_coords, set_coords)
    symbol = property(get_symbol, toss)
    element = property(get_element, toss)
    valence = property(get_valence, toss)
    fullshell = property(get_fullshell, toss)
    hypervalent = property(get_hypervalent, toss)
    charge = property(get_charge, set_charge)
    bonds = property(get_bonds, toss)
    bonded_atoms = property(get_bonded_atoms, toss)
    electrons = property(get_electrons, toss)
    nonbonding_el = property(get_nonbonding_el, toss)
    empty_valence = property(get_empty_valence, toss)
    els_from_bonds = property(get_els_from_bonds, toss)
    els_to_bonds = property(get_els_to_bonds, toss)
    radicals = property(get_radicals, toss)
    lone_pairs = property(get_lone_pairs, toss)
    bond_angles = property(get_bond_angles, toss)
    electron_angles = property(get_electron_angles, toss)

    def __init__(
        self,
        element: el.Element,
        coords: tuple[float, float],
        coord_z: int = 0,
        charge: int = 0
    ):
        """Initialize the atom with the element, coordinates and charge."""
        self._element = element
        self.charge = charge
        self._bonds = []
        (self.coord_x, self.coord_y) = coords
        self.coord_z = coord_z

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

    def to_dict(self) -> dict[str, str]:
        desc: dict[str, str] = {}
        desc["coord_x"] = str(self.coord_x)
        desc["coord_y"] = str(self.coord_y)
        desc["coord_z"] = str(self.coord_z)
        desc["charge"] = str(self._charge)
        return desc

    def can_change_element(self, element: el.Element) -> BondingError:
        """Checks whether the atom can be switched to be a different
        element."""
        in_bonds: int = 0
        donate_bond: int = 0
        for bond in self._bonds:
            in_bonds += bond.electron_count
            donate_bond += bond.atom_electrons(self)
        total_prosp_el: int = in_bonds - donate_bond + element.valence_el - self.charge
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

    def can_ionize(self, charge: int) -> BondingError:
        """Checks whether the charge of the atom can be changed to the
        specified value."""
        if self._element.valence_el - charge < 0:
            return BondingError.INSUFF_EL_FIRST
        if self._element.valence_el - charge > self.fullshell:
            return BondingError.HYPERVALENT_FIRST
        return BondingError.OK

    def bond(
        self,
        other_atom,
        order: int = 1,
        dative: int = 0,
        electrons: Optional[tuple[int, int]] = None
    ) -> CovBond:
        """
        Bonds the current atom to an other atom with the specified bond order
        and dativity (dative/coordinative bond).
        """
        assert isinstance(other_atom, Atom), "Only can bond to an Atom instance"
        bond_err = self.can_bond(
            other_atom, order=order, dative=dative, electrons=electrons
        )
        if bond_err != BondingError.OK:
            raise RuntimeError(f"Bonding to atom failed. Error: {bond_err}")
        new_bond: CovBond
        if electrons is None:
            new_bond = CovBond([self, other_atom], order=order, dative=dative)
        else:
            new_bond = CovBond([self, other_atom], electrons=list(electrons))
        self._bonds.append(new_bond)
        other_atom.register_bond(new_bond)
        return new_bond

    def can_bond(
        self,
        other_atom,
        order: int = 1,
        dative: int = 0,
        electrons: Optional[tuple[int, int]] = None,
        check_other: bool = True
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
            other_err = other_atom.can_bond(
                self, order=order, dative=-dative, electrons=electrons, check_other=False
            )
        # If there is any error from other atoms, we report that error directly
        if other_err != BondingError.OK:
            return other_err
        if other_atom == self:
            return BondingError.BOND_SELF
        if self.is_bonded(other_atom):
            return BondingError.BOND_EXISTS
        # If there are less nonbonding electrons are left than what we need
        if electrons is not None:
            if check_other:
                if self.nonbonding_el < electrons[0]:
                    return BondingError.INSUFF_EL_FIRST
                if not self.hypervalent and self.empty_valence < electrons[1]:
                    return BondingError.HYPERVALENT_FIRST
            else:
                if self.nonbonding_el < electrons[1]:
                    return BondingError.INSUFF_EL_OTHER  # This atom is NOT the "first"
                if not self.hypervalent and self.empty_valence < electrons[0]:
                    return BondingError.HYPERVALENT_OTHER
            return BondingError.OK
        if self.nonbonding_el < order + dative:
            if check_other:
                return BondingError.INSUFF_EL_FIRST
            return BondingError.INSUFF_EL_OTHER  # This atom is NOT the "first"
        if not self.hypervalent and self.empty_valence < order - dative:
            if check_other:
                return BondingError.HYPERVALENT_FIRST
            return BondingError.HYPERVALENT_OTHER  # This atom is NOT the "first"
        return BondingError.OK

    def remove_bond(self, bond_instance: Bond) -> None:
        """Removes a bond instance registered to this atom as a part of deleting
        a bond."""
        assert (
            bond_instance in self._bonds
        ), "Bond trying to be deleted is not registered to this atom"
        self._bonds.remove(bond_instance)

    def register_bond(self, new_bond: Bond) -> None:
        """
        Register a Bond instance to this atom, created by the bond()
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
        """Returns the covalent bond instance connecting the current atom to the other_atom.
        If there is no bond, returns None."""
        assert isinstance(
            other_atom, Atom
        ), "Only Atom instance can be bonded to an Atom"
        for bond in self._bonds:
            if isinstance(bond, CovBond) and other_atom in bond.other_atoms(self):
                return bond
        return None

    def unbond_all(self) -> list[Bond]:
        """Removes all bonds from the atom. Returns the instances of the now
        affected bonds."""
        bonds_loop: list[Bond] = self._bonds.copy()
        for bond in bonds_loop:
            bond.delete()
        return bonds_loop


class UnrestrictedAtom(Atom):
    """An atom that has no restrictions on bonding or ionization."""
    def get_nonbonding_el(self) -> int:
        return 0

    def get_empty_valence(self) -> int:
        return 0

    def get_radicals(self) -> int:
        return 0

    def get_lone_pairs(self) -> int:
        return 0

    def __init__(
        self,
        symbol: str,
        coords: tuple[float, float],
        coord_z: int = 0,
        charge: int = 0
    ) -> None:
        element = el.CustomElement(symbol)
        super().__init__(element, coords, coord_z=coord_z, charge=charge)

    def can_ionize(self, charge: int) -> BondingError:
        return BondingError.OK

    def can_bond(
        self,
        other_atom,
        order: int = 1,
        dative: int = 0,
        electrons: Optional[tuple[int, int]] = None,
        check_other: bool = True
    ) -> BondingError:
        assert isinstance(other_atom, Atom), "Only can bond to an Atom instance"
        other_err = BondingError.OK
        if check_other:
            other_err = other_atom.can_bond(self, order, -dative, electrons, False)
        # If there is any error from other atoms, we report that error directly
        if other_err != BondingError.OK:
            return other_err
        if other_atom == self:
            return BondingError.BOND_SELF
        if self.is_bonded(other_atom):
            return BondingError.BOND_EXISTS
        return BondingError.OK


def atom_from_dict(element: el.Element, params: dict[str, Any]) -> Atom:
    coord_x: float = 0.0
    coord_y: float = 0.0
    coord_z: int = 0
    charge: int = 0
    for key, val in params.items():
        value = str(val)
        if key == "coord_x":
            coord_x = float(value)
        if key == "coord_y":
            coord_y = float(value)
        if key == "coord_z":
            coord_z = int(value)
        if key == "charge":
            charge = int(value)
    return Atom(element, (coord_x, coord_y), coord_z, charge)


def undestr_atom_from_dict(symbol: str, params: dict[str, Any]) -> UnrestrictedAtom:
    coord_x: float = 0.0
    coord_y: float = 0.0
    coord_z: int = 0
    charge: int = 0
    for key, val in params.items():
        value = str(val)
        if key == "coord_x":
            coord_x = float(value)
        if key == "coord_y":
            coord_y = float(value)
        if key == "coord_z":
            coord_z = int(value)
        if key == "charge":
            charge = int(value)
    return UnrestrictedAtom(symbol, (coord_x, coord_y), coord_z, charge)


def element_by_symbol(symbol: str) -> Optional[el.Element]:
    """Returns the corresponding element instance to the element symbol."""
    element_dict = {element.symbol: element for element in el.element_table}
    if symbol not in element_dict:
        return None
    return element_dict[symbol]


def add_atom_by_symbol(symbol: str, coords: tuple[float, float]) -> Optional[Atom]:
    """Returns a new atom instance specified by its symbol and coordinates."""
    sel_el: Optional[el.Element] = element_by_symbol(symbol)
    if sel_el is None:
        return None
    assert isinstance(sel_el, el.Element)
    return Atom(sel_el, coords)


def relative_angles(angle_list: list[float]) -> list[float]:
    """Calculates the difference between angles that are consecutive
    in the angle list. Angles are always clcockwise, positive angles."""
    relative_angl: list[float] = []
    if len(angle_list) > 0:
        for a, b in zip(angle_list, angle_list[1:] + [angle_list[0] + 2 * pi]):
            relative_angl.append(b - a)
    return relative_angl


def angles_rel_angle(angle_list: list[float], fix_angle: float) -> list[float]:
    """Calculates the angles relative to a fixed angle from a
    list. Results are between -pi and +pi."""
    relative_angl: list[float] = []
    angle_calc: float = 0
    for angle in angle_list:
        angle_calc = angle - fix_angle
        while angle_calc > pi:
            angle_calc -= 2 * pi
        while angle_calc < -pi:
            angle_calc += 2 * pi
        relative_angl.append(angle_calc)
    return relative_angl


def angle_side_calc(angles: list[float], threshold: float = 0.0001) -> tuple[int, int]:
    """Returns the number of angles to the left and to the right with
    higher deviance from the center than the threshold. Returns:
    (angles_on_left, angles_on_right)"""
    left: int = 0
    right: int = 0
    for angle in angles:
        while angle > pi:
            angle -= 2 * pi
        while angle < -pi:
            angle += 2 * pi
        if threshold < angle < pi - threshold:
            right += 1
        if -pi + threshold < angle < -threshold:
            left += 1
    return left, right


def find_molecule(one_atom: Atom) -> tuple[list[Atom], list[CovBond]]:
    """Finds all CovBond and Atom instances that are somehow connected to the
    specified Atom instance, and hence form a molecule."""
    atomlist: list[Atom] = [one_atom]
    bondlist: list[CovBond] = []
    current_num: int = 0
    while current_num < len(atomlist):
        current_atom: Atom = atomlist[current_num]
        for bonded_atom in current_atom.bonded_atoms:
            if bonded_atom not in atomlist:
                atomlist.append(bonded_atom)
            bond_instance: Optional[CovBond]
            bond_instance = current_atom.bond_instance(bonded_atom)
            if not (bond_instance is None) and (bond_instance not in bondlist):
                bondlist.append(bond_instance)
        current_num += 1
    return atomlist, bondlist


def merge_molecules_atom(
    connect_to: Atom,
    to_replace: Atom,
    second_atom: Optional[tuple[int, int]] = None
) -> tuple[list[Atom], list[CovBond]]:
    added_atoms, added_bonds = find_molecule(to_replace)
    origin_x = connect_to.coord_x
    origin_y = connect_to.coord_y
    shift_x = -to_replace.coord_x
    shift_y = -to_replace.coord_y
    #print(f"shift: {shift_x},{shift_y}, origin:{origin_x},{origin_y}, second_atom:{second_atom}")
    if connect_to in added_atoms:
        raise RuntimeError("Can't merge a molecule with itself")
    added_atoms.remove(to_replace)
    if second_atom is not None:
        original_x = to_replace.bonded_atoms[0].coord_x + shift_x
        original_y = to_replace.bonded_atoms[0].coord_y + shift_y
        new_x = second_atom[0] - origin_x
        new_y = second_atom[1] - origin_y
        #print(f"{original_x},{original_y} -> {new_x},{new_y}")
        delta_angle = (
            atan2(new_y, new_x) -
            atan2(original_y, original_x)
        )
        multiplier = (
            sqrt(new_x ** 2 + new_y ** 2) /
            sqrt(original_x ** 2 + original_y ** 2)
        )
        #print(f"angle: {delta_angle}, multiplier: {multiplier}")
        for atom in added_atoms:
            atom.coord_x = (atom.coord_x + shift_x) * multiplier
            atom.coord_y = (atom.coord_y + shift_y) * multiplier
            atom.coords = (
                atom.coord_x * cos(delta_angle) - atom.coord_y * sin(delta_angle),
                atom.coord_x * sin(delta_angle) + atom.coord_y * cos(delta_angle)
            )
    else:
        for atom in added_atoms:
            atom.coord_x += shift_x
            atom.coord_y += shift_y
    for atom in added_atoms:
        atom.coord_x += origin_x
        atom.coord_y += origin_y
    bond_electrons: dict[Atom, tuple[int, int]] = {}
    other_atom: Atom
    for bond in to_replace.bonds:
        other_atom = bond.other_atoms(to_replace)[0]
        bond_electrons[other_atom] = (
            bond.atom_electrons(to_replace),
            bond.atom_electrons(other_atom)
        )
    removed_bonds = to_replace.unbond_all()
    for removed_bond in removed_bonds:
        assert isinstance(removed_bond, CovBond)
        added_bonds.remove(removed_bond)
    if connect_to is not None:
        for bond_to, electrons in bond_electrons.items():
            new_bond = connect_to.bond(
                bond_to, electrons=electrons
            )
            added_bonds.append(new_bond)
    return added_atoms, added_bonds


def merge_molecules_bonds(
    orig_bond: CovBond,
    template_bond: CovBond,
    direction: Optional[tuple[int, int]] = None
) -> tuple[list[Atom], list[CovBond]]:
    added_atoms, added_bonds = find_molecule(template_bond.atoms[0])
    if orig_bond.atoms[0] in added_atoms:
        raise RuntimeError("Can't merge a molecule with itself")
    if (
        abs(orig_bond.dativity) != abs(template_bond.dativity) or
        orig_bond.order != template_bond.order
    ):
        raise RuntimeError("Can't merge bonds of different dativity or order")
    origin_atom = orig_bond.atoms[0]
    directional_atom = orig_bond.atoms[1]
    if orig_bond.dativity == -template_bond.dativity:
        origin_atom, directional_atom = directional_atom, origin_atom
    shift_x = -template_bond.atoms[0].coord_x
    shift_y = -template_bond.atoms[0].coord_y
    rotate_to_x = directional_atom.coord_x - origin_atom.coord_x
    rotate_to_y = directional_atom.coord_y - origin_atom.coord_y
    template_x = template_bond.atoms[1].coord_x - template_bond.atoms[0].coord_x
    template_y = template_bond.atoms[1].coord_y - template_bond.atoms[0].coord_y
    delta_angle = (
        atan2(rotate_to_y, rotate_to_x) -
        atan2(template_y, template_x)
    )
    multiplier = (
        sqrt(rotate_to_x ** 2 + rotate_to_y ** 2) /
        sqrt(template_x ** 2 + template_y ** 2)
    )
    for atom in added_atoms:
        atom.coord_x = (atom.coord_x + shift_x) * multiplier
        atom.coord_y = (atom.coord_y + shift_y) * multiplier
        atom.coords = (
            atom.coord_x * cos(delta_angle) - atom.coord_y * sin(delta_angle),
            atom.coord_x * sin(delta_angle) + atom.coord_y * cos(delta_angle)
        )
    for atom in added_atoms:
        atom.coord_x += origin_atom.coord_x
        atom.coord_y += origin_atom.coord_y
    bond_electrons: list[dict[Atom, tuple[int, int]]] = []
    bonds_dict: dict[Atom, tuple[int, int]]
    for atom in template_bond.atoms:
        bonds_dict = {}
        for bond in atom.bonds:
            if bond == template_bond:
                continue
            other_atom = bond.other_atoms(atom)[0]
            bonds_dict[other_atom] = (
                bond.atom_electrons(atom),
                bond.atom_electrons(other_atom)
            )
        bond_electrons.append(bonds_dict)
    for atom in template_bond.atoms:
        removed_bonds = atom.unbond_all()
        for removed_bond in removed_bonds:
            assert isinstance(removed_bond, CovBond)
            added_bonds.remove(removed_bond)
        added_atoms.remove(atom)
    for atom, bonds_dict in zip([origin_atom, directional_atom], bond_electrons):
        for bond_to, electrons in bonds_dict.items():
            new_bond = atom.bond(
                bond_to, electrons=electrons
            )
            added_bonds.append(new_bond)
    return added_atoms, added_bonds


def merge_molecules(
    connect: Union[tuple[Atom, Atom], tuple[CovBond, CovBond]],
    second_atom: Optional[tuple[int, int]] = None
) -> tuple[list[Atom], list[CovBond]]:
    if isinstance(connect[0], Atom):
        return merge_molecules_atom(*connect, second_atom)
    if isinstance(connect[0], CovBond):
        return merge_molecules_bonds(*connect, second_atom)
    assert False, "Unreachable"
    return ([], [])


def can_merge_molecules(
    connect: Union[tuple[Atom, Atom], tuple[CovBond, CovBond]],
) -> BondingError:
    if isinstance(connect[0], Atom):
        connect_to: Atom = connect[0]
        to_replace: Atom = connect[1]
        if (
            connect_to.empty_valence < to_replace.els_from_bonds and
            not connect_to.hypervalent
        ):
            return BondingError.HYPERVALENT_FIRST
        if (
            connect_to.nonbonding_el < to_replace.els_to_bonds
        ):
            return BondingError.INSUFF_EL_FIRST
        return BondingError.OK
    if isinstance(connect[0], CovBond):
        orig_bond: CovBond = connect[0]
        template_bond: CovBond = connect[1]
        if (
            abs(orig_bond.dativity) != abs(template_bond.dativity) or
            orig_bond.order != template_bond.order
        ):
            return BondingError.IMPROPER_BOND
        orig_first: Atom = orig_bond.atoms[0]
        orig_second: Atom = orig_bond.atoms[1]
        template_first: Atom = template_bond.atoms[0]
        template_second: Atom = template_bond.atoms[1]
        if orig_bond.dativity == -template_bond.dativity:
            orig_first, orig_second = orig_second, orig_first
        if (
            orig_first.empty_valence < (
                template_first.els_from_bonds -
                template_bond.atom_electrons(template_second)
            ) and
            not orig_first.hypervalent
        ):
            return BondingError.HYPERVALENT_FIRST
        if (
            orig_first.nonbonding_el < (
                template_first.els_to_bonds -
                template_bond.atom_electrons(template_first)
            )
        ):
            return BondingError.INSUFF_EL_FIRST
        if (
            orig_second.empty_valence < (
                template_second.els_from_bonds -
                template_bond.atom_electrons(template_first)
            ) and
            not orig_second.hypervalent
        ):
            return BondingError.HYPERVALENT_OTHER
        if (
            orig_second.nonbonding_el < (
                template_second.els_to_bonds -
                template_bond.atom_electrons(template_second)
            )
        ):
            return BondingError.INSUFF_EL_OTHER
        return BondingError.OK
    return BondingError.MISC_ERROR


def atomdist(one_atom: Atom, other_atom: Atom) -> float:
    """Calculates the distance between two atoms."""
    return sqrt(
        (one_atom.coord_x - other_atom.coord_x) ** 2
        + (one_atom.coord_y - other_atom.coord_y) ** 2
    )


def atom_angle(one_atom: Atom, other_atom: Atom) -> float:
    """Returns the angle of the vector pointung from one_atom to
    other_atom, with the angle being between 0 and 2*pi."""
    xdiff: float = other_atom.coord_x - one_atom.coord_x
    ydiff: float = other_atom.coord_y - one_atom.coord_y
    angle: float = atan2(ydiff, xdiff)
    if angle < 0:
        angle = 2 * pi + angle
    return angle


def to_xml(atom_bond_list: list[Union[Atom, Bond]], filename: str) -> None:
    """
    Saves the list of atoms and bonds to an xml file.
    """
    root: ET.Element = ET.Element(
        "NiceMolecules_Data",
        attrib={"version": VERSION}
    )
    used_elements: list[el.Element] = []
    desc: dict[str, str]
    atom_list: list[Atom] = [
        thing for thing in atom_bond_list if isinstance(thing, Atom)
    ]
    bond_list: list[Bond] = [
        thing for thing in atom_bond_list if isinstance(thing, Bond)
    ]
    for atom in atom_list:
        if isinstance(atom, UnrestrictedAtom):
            continue
        new_element = atom.element
        if new_element not in used_elements:
            used_elements.append(new_element)
    for element in used_elements:
        desc = {"symbol": element.symbol}
        if element in el.element_table:
            desc["builtin"] = "y"
        data_unit = ET.SubElement(root, "Element", attrib=desc)
        for key, value in element.to_dict().items():
            ET.SubElement(data_unit, key).text = value
    for atom in atom_list:
        desc = {}
        if isinstance(atom, UnrestrictedAtom):
            desc["unrestricted"] = "y"
            desc["symbol"] = atom.symbol
        else:
            desc["element_index"] = str(used_elements.index(atom.element))
        data_unit = ET.SubElement(root, "Atom", attrib=desc)
        for key, value in atom.to_dict().items():
            ET.SubElement(data_unit, key).text = value
    for bond in bond_list:
        desc = {}
        if isinstance(bond, CovBond):
            desc["covalent"] = "y"
        data_unit = ET.SubElement(root, "Bond", attrib=desc)
        for atom, electrons in zip(bond.atoms, bond.electrons):
            ET.SubElement(
                data_unit,
                "Atom",
                attrib={"electrons": str(electrons)}
            ).text = str(atom_list.index(atom))
    tree = ET.ElementTree(root)
    tree.write(filename)


def read_xml(filename: str) -> list[Union[Atom, Bond]]:
    loaded: list[Union[Atom, Bond]] = []
    tree = ET.parse(filename)
    root = tree.getroot()
    assert root.tag == "NiceMolecules_Data", "The XML file has some issues."
    builtin_symbols: list[str] = [
        element.symbol for element in el.element_table
    ]
    element_list: list[el.Element] = []
    symbol: str
    if "version" in root.attrib:
        print("Loading file created with version ", root.attrib["version"])
    else:
        print("File version unknown.")
    for element in [thing for thing in root if thing.tag == "Element"]:
        symbol = element.attrib["symbol"]
        if "builtin" in element.attrib and element.attrib["builtin"] == "y":
            if symbol in builtin_symbols:
                element_list.append(
                    el.element_table[builtin_symbols.index(symbol)]
                )
            continue
        element_list.append(
            el.element_from_dict(
                symbol,
                {detail.tag: detail.text for detail in element}
            )
        )
    for atom in [thing for thing in root if thing.tag == "Atom"]:
        if "unrestricted" in atom.attrib and atom.attrib["unrestricted"] == "y":
            if "symbol" in atom.attrib:
                symbol = atom.attrib["symbol"]
            else:
                symbol = "R"
            loaded.append(
                undestr_atom_from_dict(
                    symbol,
                    {detail.tag: detail.text for detail in atom}
                )
            )
            continue
        loaded.append(
            atom_from_dict(
                element_list[int(atom.attrib["element_index"])],
                {detail.tag: detail.text for detail in atom}
            )
        )
    bond_instance: Bond
    for bond in [thing for thing in root if thing.tag == "Bond"]:
        if "covalent" in bond.attrib and bond.attrib["covalent"] == "y":
            bond_instance = CovBond(
                    [loaded[int(str(atom.text))] for atom in bond],
                    [int(atom.attrib["electrons"]) for atom in bond]
                )
            loaded.append(bond_instance)
            for atom in bond:
                atom_instance = loaded[int(str(atom.text))]
                assert isinstance(atom_instance, Atom)
                atom_instance.register_bond(bond_instance)
    return loaded
