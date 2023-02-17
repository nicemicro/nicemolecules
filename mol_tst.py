#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 07:43:08 2022

@author: nicemicro
"""

import mol_eng as eng
import mol_opt as opt
import elements as el

def carbon_monoxide():
    a = eng.Atom(el.Carbon(), (20, 20))
    b = eng.Atom(el.Oxygen(), (40, 20))
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
    a = eng.add_atom_by_symbol(x, (20, 20))
    b = eng.add_atom_by_symbol("H", (40, 20))
    assert isinstance(a, eng.Atom)
    assert isinstance(b, eng.Atom)
    print(a)
    print(b)
    c = a.bond(b)
    print(f"bond a to b: {c}")
    return [a, b], a.bonds

def SF6(n: int=6, central: str="S", terminal: str="F") -> tuple[list[eng.Atom], list[eng.CovBond]]:
    loc_x = [60, 100, 120, 60, 100, 40]
    loc_y = [50, 50, 80, 110, 110, 80]
    atoms: list[eng.Atom] = []
    new_atom = eng.add_atom_by_symbol(central, (80, 80))
    assert isinstance(new_atom, eng.Atom)
    atoms.append(new_atom)
    for x, y in zip(loc_x[0:n], loc_y[0:n]):
        new_atom = eng.add_atom_by_symbol(terminal, (x, y))
        assert isinstance(new_atom, eng.Atom)
        atoms.append(new_atom)
        print(atoms[-1].can_bond(atoms[0], electrons=(1, 1)))
        print(atoms[0].can_bond(atoms[-1], electrons=(1, 1)))
        atoms[-1].bond(atoms[0], electrons=(1, 1))
    print(f"Bond angles: {atoms[0].bond_angles}")
    print(f"Electron angles: {atoms[0].electron_angles}")
    return atoms, atoms[0].bonds

def methane(x = 40, y = 40):
    a = eng.Atom(el.Carbon(), (x, y))
    b = eng.Atom(el.Hydrogen(), (x-20, y))
    c = eng.Atom(el.Hydrogen(), (x+20, y))
    d = eng.Atom(el.Hydrogen(), (x, y-20))
    e = eng.Atom(el.Hydrogen(), (x, y+20))
    a.bond(b)
    a.bond(c)
    a.bond(d)
    a.bond(e)
    print(f"Bond angles: {a.bond_angles}")
    print(f"Electron angles: {a.electron_angles}")
    return [a, b, c, d, e], a.bonds

def custom_element():
    custom_el = el.CustomElement("Cstm", 3, 4, False)
    return eng.Atom(custom_el, (150, 150))

def unrestricted_atom():
    return eng.UnrestrictedAtom("Unr", (200, 100))

def description(atoms, to_show = None):
    if to_show is None:
        to_show = [0]
    for show_this in to_show:
        print(atoms[show_this].describe())
        print(atoms[show_this].bonds)
        print(atoms[show_this].bonds[0].dativity)

def optimize_tst():
    atoms, bonds = SF6(4, "Xe", "F")
    print(atoms)
    print([bond.length for bond in bonds])
    opt.optimize_2D(atoms[0])
    print("== After optimization ==")
    print(atoms)
    print([bond.length for bond in bonds])

def join_molecules():
    atoms, bonds = methane()
    bonds.remove(atoms[4].bonds[0])
    atoms[0].bond_instance(atoms[4]).delete()
    atoms.pop(4)
    atoms2, _ = methane(40, 90)
    for atom in atoms + atoms2:
        print(f"{atom.symbol} at ({atom.coords})")
    atoms3, bonds3 = eng.merge_molecules(atoms[0], atoms2[3])
    print("----")
    #for atom in atoms3:
    #    print(f"{atom.symbol} at ({atom.coords})")
    print(atoms[0].describe())
    return atoms + atoms3, bonds + bonds3

if __name__ == '__main__':
    #atoms, bonds = carbon_monoxide()
    #atoms, bonds = XH("B")
    #atoms, bonds = SF6(6, "C", "H")
    #atoms, bonds = methane()
    atoms, bonds = join_molecules()
    #atoms += [custom_element()]
    #atoms += [unrestricted_atom()]
    #description(atoms, [0])
    #optimize_tst()
    #eng.to_xml(atoms + bonds + [custom_element()], "test.nm.xml")
