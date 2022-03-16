#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 07:43:08 2022

@author: nicemicro
"""

import mol_eng as eng

def carbon_monoxide():
    a = eng.Atom("C", 4, [20, 20], 0)
    b = eng.Atom("O", 6, [40, 20], 0)
    print(a)
    print(b)
    print("bond a to b: ", a.bond(b, 3, 1))
    print("bond b to a: ", b.bond(a, 3, -1))
    return [a, b], a.bonds

def XH(x = "C"):
    a = eng.add_atom_by_symbol(x, [20, 20])
    b = eng.add_atom_by_symbol("H", [40, 20])
    print(a)
    print(b)
    print("bond a to b: ", a.bond(b))
    return [a, b], a.bonds

def SF6(n: int=6, central: str="S", terminal: str="F") -> tuple[list[eng.Atom], list[eng.Cov_bond]]:
    loc_x = [60, 100, 40, 120, 60, 100]
    loc_y = [50, 50, 80, 80, 110, 110]
    molecules: list[eng.Atom] = []
    molecules.append(eng.add_atom_by_symbol(central, [80, 80]))
    for x, y in zip(loc_x[0:n], loc_y[0:n]):
        molecules.append(eng.add_atom_by_symbol(terminal, [x, y]))
        molecules[-1].bond(molecules[0])
    return molecules, molecules[0].bonds

def methane():
    a = eng.Atom("C", 4, [40, 40])
    b = eng.Atom("H", 1, [20, 40], fullshell=2)
    c = eng.Atom("H", 1, [40, 20], fullshell=2)
    d = eng.Atom("H", 1, [20, 60], fullshell=2)
    e = eng.Atom("H", 1, [60, 20], fullshell=2)
    a.bond(b)
    a.bond(c)
    a.bond(d)
    a.bond(e)
    return [a, b, c, d, e], a.bonds

def compl():
    molecules = []
    bonds = []
    molecules.append(eng.Atom("C", 4, [40, 40]))
    for index in range(4):
        molecules.append(eng.Atom("H", 1, [20 + (index % 2)*40, 20+(index // 2)*40], fullshell=2))
        new_bond = molecules[-1].bond(molecules[0])
        bonds.append(new_bond)
    molecules.append(eng.Atom("C", 4, [40, 150]))
    molecules.append(eng.Atom("O", 6, [70, 140]))
    new_bond = molecules[-2].bond(molecules[-1], 3, 1)
    bonds.append(new_bond)
    carbons = []
    hydrogens = []
    for index in range(3):
        carbons.append(eng.Atom("C", 4, [120+index*30, 60]))
        if len(carbons) > 1:
            new_bond = carbons[-1].bond(carbons[-2])
            bonds.append(new_bond)
        for index2 in range(2):
            hydrogens.append(eng.Atom("H", 1, [120+index*30, 30+index2*60], fullshell=2))
            new_bond = hydrogens[-1].bond(carbons[-1])
            bonds.append(new_bond)
    hydrogens.append(eng.Atom("H", 1, [210, 60], fullshell=2))
    new_bond = hydrogens[-1].bond(carbons[-1])
    bonds.append(new_bond)
    hydrogens.append(eng.Atom("H", 1, [90, 60], fullshell=2))
    new_bond = hydrogens[-1].bond(carbons[0])
    bonds.append(new_bond)
    return molecules + carbons + hydrogens, bonds

def description(atoms, to_show = None):
    if to_show is None:
        to_show = [0]
    for show_this in to_show:
        print(atoms[show_this].describe())
        print(atoms[show_this].bonds)
        print(atoms[show_this].bonds[0].dativity())

if __name__ == '__main__':
    #molecules, bonds = carbon_monoxide()
    #molecules, bonds = XH("B")
    molecules, bonds = SF6(3)
    #molecules, bonds = methane()
    #molecules, bonds = test3()
    description(molecules, [0])
