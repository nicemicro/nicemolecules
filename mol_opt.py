# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 15:52:19 2022

@author: nicemicro
"""
from math import pi, sin, cos
from typing import Optional
import mol_eng as eng
from random import uniform

IDEAL_REL_ANGLES: list[list[float]] = [
    [0.0],
    [pi] * 2,
    [pi * 2 / 3] * 3,
    [pi / 2, pi / 3, pi / 2, pi * 2 / 3],
    [pi * 2 / 3, pi / 2, pi / 3, pi / 2],
    [pi * 2 / 3, pi / 3, pi / 3, pi * 2 / 3],
    [pi / 2, pi / 3, pi / 3, pi / 3, pi / 2],
    [pi / 3] * 6,
    [pi / 3, pi / 3, pi / 6, pi / 3, pi / 6, pi / 3, pi / 3],
    [pi / 3, pi / 6, pi / 3, pi / 6, pi / 6, pi / 3, pi / 6, pi / 3],
]


def sum_square(numbers: list[float]) -> float:
    summary: float = 0
    for number in numbers:
        summary += number**2
    return summary


def angle_diff(angles1: list[float], angles2: list[float]) -> list[float]:
    diff_list: list[float] = []
    diff: float
    for a, b in zip(angles1, angles2):
        diff = b - a
        while diff > pi:
            diff -= 2 * pi
        while b < -pi:
            diff += 2 * pi
        diff_list.append(diff)
    return diff_list


def angle_deltas(bond_angle_list: list[float], electron_num) -> list[float]:
    rel_angles: list[float] = eng.relative_angles(bond_angle_list)
    deltas: list[float] = []
    angles_comps: list[list[float]] = [
            angle_list
            for angle_list in IDEAL_REL_ANGLES
            if len(angle_list) == len(bond_angle_list) + electron_num
            ]
    angles_comps_all: list[list[float]]
    if electron_num == 0:
        angles_comps_all = angles_comps
    else:
        angles_comps_all = []
        for angles_comp in angles_comps:
            for start, _ in enumerate(angles_comp):
                if start == 0:
                    new_angles_comp = angles_comp
                else:
                    new_angles_comp = angles_comp[start:] + angles_comp[0:start]
                angles_comps_all.append(
                    new_angles_comp[:-electron_num - 1] +
                    [sum(new_angles_comp[-electron_num - 1:])]
                )
    best_diff: float = -1
    #print(f"angles: {[int(ang * 180 / pi) for ang in rel_angles]}")
    for angles_comp in angles_comps_all:
        for start, _ in enumerate(rel_angles):
            if start == 0:
                calc_angles = rel_angles
            else:
                calc_angles = rel_angles[start:] + rel_angles[0:start]
            #print(f"  shifted: {[int(ang * 180 / pi) for ang in calc_angles]}")
            #print(f"  comp to: {[int(ang * 180 / pi) for ang in angles_comp]}")
            new_diff = sum_square(angle_diff(calc_angles, angles_comp))
            #print(f"    diff: {new_diff}")
            if new_diff < best_diff or best_diff == -1:
                #print("      BEST DIFF")
                best_diff = new_diff
                deltas = angle_diff(calc_angles, angles_comp)
                if start != 0:
                    deltas = deltas[-start:] + deltas[:-start]
                #print(f"      deltas: {[int(ang * 180 / pi) for ang in deltas]}")
    return deltas


def bond_angle_instances(atom: eng.Atom) -> tuple[list[float], list[eng.CovBond]]:
    bond_angles: list[float] = atom.bond_angles
    bonds: list[eng.CovBond] = atom.bonds
    bond_angles_orig = bond_angles.copy()
    bond_angles.sort()
    bonds_ordered: list[eng.CovBond] = []
    for angle in bond_angles:
        index = bond_angles_orig.index(angle)
        bonds_ordered.append(bonds[index])
        bonds.pop(index)
        bond_angles_orig.pop(index)
    return bond_angles, bonds_ordered


def optimize_dist(
    atomlist: list[eng.Atom],
    target_len: float = 30.,
    alpha: float = 0.1,
) -> None:
    """Goes through all the atoms and moves them a bit so they are closer
    to the target bond length defined by the target_len parameter."""
    deltas_xs: list[float] = [0] * len(atomlist)
    deltas_ys: list[float] = [0] * len(atomlist)
    for atom_index, atom in enumerate(atomlist):
        for bond, angle in zip(atom.bonds, atom.bond_angles):
            delta_len = (bond.length - target_len) * alpha
            deltas_xs[atom_index] += delta_len * cos(angle)
            deltas_ys[atom_index] += delta_len * sin(angle)
    for atom, delta_x, delta_y in zip(atomlist, deltas_xs, deltas_ys):
        atom.coord_x += delta_x
        atom.coord_y += delta_y


def push_away_close(
    atomlist: list[eng.Atom],
    min_dist: float = 20.,
    alpha: float = 0.1,
) -> None:
    """Goes through all atos and pushes atoms that are too close, apart by
    a bit."""
    deltas_xs: list[float] = [0] * len(atomlist)
    deltas_ys: list[float] = [0] * len(atomlist)
    for atom_index, atom in enumerate(atomlist):
        for other_atom in atomlist:
            if other_atom == atom or atom.is_bonded(other_atom):
                continue
            if eng.atomdist(atom, other_atom) < min_dist:
                delta_len = min_dist - eng.atomdist(atom, other_atom)
                deltas_xs[atom_index] += (delta_len**2 * alpha) * cos(
                    eng.atom_angle(atom, other_atom)
                )
                deltas_ys[atom_index] += (delta_len**2 * alpha) * sin(
                    eng.atom_angle(atom, other_atom)
                )
                deltas_xs[atom_index] += uniform(
                    -min_dist * alpha, min_dist * alpha
                )
                deltas_ys[atom_index] += uniform(
                    -min_dist * alpha, min_dist * alpha
                )
    for atom, delta_x, delta_y in zip(atomlist, deltas_xs, deltas_ys):
        atom.coord_x += delta_x
        atom.coord_y += delta_y


def move_connected_atom(
    moved_atom: eng.Atom,
    orig_atom: eng.Atom,
    delta_x: float,
    delta_y: float,
) -> None:
    """Moving atoms and all its connections (except for those that are connected
    to the original atom) by some delta x and y."""
    atoms_to_move: list[eng.Atom] = [moved_atom]
    current_num: int = 0
    current_atom: eng.Atom
    while current_num < len(atoms_to_move):
        current_atom = atoms_to_move[current_num]
        for bonded_atom in current_atom.bonded_atoms:
            if (
                bonded_atom not in atoms_to_move and
                bonded_atom != orig_atom and
                bonded_atom not in orig_atom.bonded_atoms
            ):
                atoms_to_move.append(bonded_atom)
        current_num += 1
    for atom in atoms_to_move:
        atom.coord_x += delta_x
        atom.coord_y += delta_y


def optimize_relative_angles(
    atomlist: list[eng.Atom],
    alpha: float = 0.1,
) -> None:
    """Goes through the atoms and moves them towards a more favorable angle."""
    for atom in atomlist:
        #print(f"  Atom optimized {atom.symbol} ({atom.coord_x}, {atom.coord_y})")
        bond_angles: list[float]
        bonds: list[eng.CovBond]
        bond_angles, bonds = bond_angle_instances(atom)
        #print(f"    bond angles: {bond_angles}")
        delta_rel_angles: list[float] = angle_deltas(
            bond_angles, atom.radicals + atom.lone_pairs
        )
        delta_angles: list[float] = [0.] * len(delta_rel_angles)
        for index, delta_angle in enumerate(delta_rel_angles):
            delta_angles[index] -= delta_angle / 2
            if index + 1 == len(delta_rel_angles):
                delta_angles[0] += delta_angle / 2
            else:
                delta_angles[index + 1] += delta_angle / 2
        #print(f"    delta angles: {delta_angles}")
        #print(f"    Iterating through bonds")
        for bond, delta_angle in zip(bonds, delta_angles):
            delta_angle *= (alpha / 5)
            other_index = atomlist.index(bond.other_atoms(atom)[0])
            #print(f"      other atom index: {other_index}")
            rel_x = atomlist[other_index].coord_x - atom.coord_x
            rel_y = atomlist[other_index].coord_y - atom.coord_y
            #print(f"      relative position: ({rel_x}, {rel_y})")
            #print(f"         delta_x: {rel_x * cos(delta_angle) - rel_y * sin(delta_angle) - rel_x}")
            #print(f"         delta_y: {rel_x * sin(delta_angle) + rel_y * cos(delta_angle) - rel_y}")
            move_connected_atom(
                bond.other_atoms(atom)[0],
                atom,
                rel_x * cos(delta_angle) - rel_y * sin(delta_angle) - rel_x,
                rel_x * sin(delta_angle) + rel_y * cos(delta_angle) - rel_y,
            )


def optimize_tilt(atomlist: list[eng.Atom], alpha: float = 0.1) -> None:
    """Goes thrugh every atom, and moves them in the direction of an
    absolute angle on the screen that is a multiple of 30 degrees."""
    deltas_xs: list[float] = [0] * len(atomlist)
    deltas_ys: list[float] = [0] * len(atomlist)
    for atom_index, atom in enumerate(atomlist):
        for bond, angle in zip(atom.bonds, atom.bond_angles):
            delta_angle = sin(angle * 6 - pi) / 10 * alpha
            if cos(pi - angle * 6) > 0:
                delta_angle *= -1
            rel_x = atom.coord_x - bond.other_atoms(atom)[0].coord_x
            rel_y = atom.coord_y - bond.other_atoms(atom)[0].coord_y
            move_connected_atom(
                bond.other_atoms(atom)[0],
                atom,
                rel_x * cos(delta_angle) - rel_y * sin(delta_angle) - rel_x,
                rel_x * sin(delta_angle) + rel_y * cos(delta_angle) - rel_y,
            )


def optimize_2D(
    one_atom: eng.Atom,
    iterator: int = 3,
    segments: Optional[list[int]] = None,
    target_len: float = 30.,
    min_dist: float = 20.,
    alpha: float = 0.1,
) -> None:
    """Finds the molecule that contains the atom, then optimzies the 2D
    formula."""
    #print("--- Optimization Starts ---")
    if segments is None:
        segments = [1, 3, 10, 3, 1]
        #segments = [1]
    assert isinstance(segments, list)
    iter_pattern: list[int] = segments * iterator
    optim_pattern: list[str] = []
    for multiplier in iter_pattern:
        optim_pattern += ["len"] * multiplier
        optim_pattern += ["angles"] * multiplier
        optim_pattern += ["tilt"] * multiplier
    atomlist: list[eng.Atom]
    atomlist, _ = eng.find_molecule(one_atom)
    original_x: float = sum([atom.coord_x for atom in atomlist])
    original_y: float = sum([atom.coord_y for atom in atomlist])
    for optimizer in optim_pattern:
        if optimizer == "len":
            optimize_dist(atomlist, target_len, alpha)
        elif optimizer == "angles":
            optimize_relative_angles(atomlist, alpha)
        elif optimizer == "tilt":
            optimize_tilt(atomlist, alpha)
        push_away_close(atomlist, min_dist, alpha)
    new_x: float = sum([atom.coord_x for atom in atomlist])
    new_y: float = sum([atom.coord_y for atom in atomlist])
    delta_x = (new_x - original_x) / len(atomlist)
    delta_y = (new_y - original_y) / len(atomlist)
    for atom in atomlist:
        atom.coord_x -= delta_x
        atom.coord_x = round(atom.coord_x)
        atom.coord_y -= delta_y
        atom.coord_y = round(atom.coord_y)
    #for atom_index, atom in enumerate(atomlist):
    #    print(f"atom {atom_index}")
    #    for bond, angle in zip(atom.bonds, atom.bond_angles):
    #        print(f"  bond length {bond.length}, angle {angle / pi *180} deg")
