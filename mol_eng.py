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
            desc += f"  Atom letter: {atom.letter} with {electron} electrons\n"
        return desc
    
    def __str__(self):
        return super().__str__()
    
    def __repr__(self):
        return super().__repr__() + "\n" + self.describe()
    
    def order(self):
        return sum(self.electrons) / len(self.atoms)
    
    def electron_count(self):
        return sum(self.electrons)
    
    def atom_electrons(self, one_atom):
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
    def __init__(self, letter, valence_el, coords, charge=0, fullshell=8,
                 hypervalent=False):
        self.letter = letter
        self.valence = valence_el
        self.charge = 0
        self.fullshell = fullshell
        self.hypervalent = hypervalent
        self.bonds = []
        self.coord_x = coords[0]
        self.coord_y = coords[1]
    
    def describe(self, short=False):
        desc  = f"{self.letter} atom at ({self.coord_x}, {self.coord_y})"
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
                desc += f"{atom.letter} at ({atom.coord_x}, {atom.coord_y})"
            desc += "\n"
        return desc
    
    def __str__(self):
        return super().__str__()
    
    def __repr__(self):
        return super().__repr__() + "\n" + self.describe(True)
    
    def bond(self, other_atom, order=1, dative=0):
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
        return other_atom in [item for sublist in bond_list_lists for item in sublist]
    
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

def test1():
    print()
    print("--== Test 1 ==-")
    print()
    a = Atom("C", 4, [20, 20], 0)
    b = Atom("O", 6, [40, 20], 0)
    print(a)
    print(b)
    print("bond a to b: ", a.bond(b, 3, 1))
    print("bond b to a: ", b.bond(a, 3, -1))
    print(a)
    print(b)
    print(a.bonds)

def test2():
    print()
    print("--== Test 2 ==-")
    print()
    a = Atom("C", 4, [40, 40])
    b = Atom("H", 1, [20, 40], fullshell=2)
    c = Atom("H", 1, [40, 20], fullshell=2)
    d = Atom("H", 1, [20, 60], fullshell=2)
    e = Atom("H", 1, [60, 20], fullshell=2)
    a.bond(b)
    a.bond(c)
    a.bond(d)
    a.bond(e)
    print(a)
    print(a.bonds)

def test3():
    molecules = []
    bonds = []
    molecules.append(Atom("C", 4, [40, 40]))
    for index in range(4):
        molecules.append(Atom("H", 1, [20 + (index % 2)*40, 20+(index // 2)*40], fullshell=2))
        new_bond = molecules[-1].bond(molecules[0])
        bonds.append(new_bond)
    molecules.append(Atom("C", 4, [40, 150]))
    molecules.append(Atom("O", 6, [70, 120]))
    new_bond = molecules[-2].bond(molecules[-1], 3, 1)
    bonds.append(new_bond)
    carbons = []
    hydrogens = []
    for index in range(3):
        carbons.append(Atom("C", 4, [120+index*30, 60]))
        if len(carbons) > 1:
            new_bond = carbons[-1].bond(carbons[-2])
            bonds.append(new_bond)
        for index2 in range(2):
            hydrogens.append(Atom("H", 1, [120+index*30, 30+index2*60], fullshell=2))
            new_bond = hydrogens[-1].bond(carbons[-1])
            bonds.append(new_bond)
    hydrogens.append(Atom("H", 1, [210, 60], fullshell=2))
    new_bond = hydrogens[-1].bond(carbons[-1])
    bonds.append(new_bond)
    hydrogens.append(Atom("H", 1, [90, 60], fullshell=2))
    new_bond = hydrogens[-1].bond(carbons[0])
    bonds.append(new_bond)
    return molecules + carbons + hydrogens, bonds

if __name__ == '__main__':
    #test1()
    #test2()
    molecules, bonds = test3()
