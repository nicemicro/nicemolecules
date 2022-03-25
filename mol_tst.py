#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 07:43:08 2022

@author: nicemicro
"""

import mol_eng as eng
import elements as el

def carbon_monoxide():
    a = eng.Atom(el.Carbon(), [20, 20])
    b = eng.Atom(el.Oxygen(), [40, 20])
    print(a)
    print(b)
    c = a.bond(b, 3, -1)
    print(f"bond a to b: {c}")
    print(c.electrons)
    print("dativity: ", c.dativity)
    print("set order to 1")
    c.order = 1
    print(c.electrons)
    print("set dativity to 0")
    c.dativity = 0
    print(c.electrons)
    print("set order to 2")
    c.order = 2
    print(c.electrons)
    print("set dativity to -1")
    c.dativity = -1
    print(c.electrons)
    return [a, b], a.bonds

def XH(x = "C"):
    a = eng.add_atom_by_symbol(x, [20, 20])
    b = eng.add_atom_by_symbol("H", [40, 20])
    assert isinstance(a, eng.Atom)
    assert isinstance(b, eng.Atom)
    print(a)
    print(b)
    c = a.bond(b)
    print(f"bond a to b: {c}")
    return [a, b], a.bonds

def SF6(n: int=6, central: str="S", terminal: str="F") -> tuple[list[eng.Atom], list[eng.CovBond]]:
    loc_x = [60, 100, 40, 120, 60, 100]
    loc_y = [50, 50, 80, 80, 110, 110]
    molecules: list[eng.Atom] = []
    new_atom = eng.add_atom_by_symbol(central, [80, 80])
    assert isinstance(new_atom, eng.Atom)
    molecules.append(new_atom)
    for x, y in zip(loc_x[0:n], loc_y[0:n]):
        new_atom = eng.add_atom_by_symbol(terminal, [x, y])
        assert isinstance(new_atom, eng.Atom)
        molecules.append(new_atom)
        molecules[-1].bond(molecules[0])
    return molecules, molecules[0].bonds

def methane():
    a = eng.Atom(el.Carbon(), [40, 40])
    b = eng.Atom(el.Hydrogen(), [20, 40])
    c = eng.Atom(el.Hydrogen(), [40, 40])
    d = eng.Atom(el.Hydrogen(), [20, 40])
    e = eng.Atom(el.Hydrogen(), [20, 40])
    a.bond(b)
    a.bond(c)
    a.bond(d)
    a.bond(e)
    return [a, b, c, d, e], a.bonds

def description(atoms, to_show = None):
    if to_show is None:
        to_show = [0]
    for show_this in to_show:
        print(atoms[show_this].describe())
        print(atoms[show_this].bonds)
        print(atoms[show_this].bonds[0].dativity)

if __name__ == '__main__':
    molecules, bonds = carbon_monoxide()
    #molecules, bonds = XH("B")
    #molecules, bonds = SF6(4, "C", "O")
    #molecules, bonds = methane()
    #molecules, bonds = test3()
    #description(molecules, [0])
