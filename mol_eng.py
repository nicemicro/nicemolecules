# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 15:52:19 2022

@author: nicemicro
"""

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
    
    def describe(self):
        desc  = f"Location of atom at ({self.coord_x}, {self.coord_y})\n"
        desc += f"  Atom letter representation   : {self.letter}\n"
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
        return super().__repr__() + "\n" + self.describe()
    
    def bond(self, other_atom, order=1, dative=0):
        if self.is_bonded(other_atom):
            return False
        if other_atom.is_bonded(self):
            self.bonds.append(other_atom.bond_instance(self))
        else:
            new_bond = Cov_bond([self, other_atom], [order-dative, order+dative])
            self.bonds.append(new_bond)
            other_atom.bond(self, order, -dative)
        return True
    
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
    current_num = 0
    while current_num < len(atomlist):
        current_atom = atomlist[current_num]
        for bonded_atom in current_atom.bonded_atoms():
            if not bonded_atom in atomlist:
                atomlist.append(bonded_atom )
        current_num += 1
    return atomlist

def test1():
    a = Atom("C", 4, [0, 0], 0)
    b = Atom("O", 6, [1, 0], 0)
    print(a)
    print(b)
    print("bond a to b: ", a.bond(b, 3, 1))
    print("bond b to a: ", b.bond(a, 3, -1))
    print(a)
    print(b)
    print(a.bonds)

def test2():
    a = Atom("C", 4, [1, 1])
    b = Atom("H", 1, [0, 1], fullshell=2)
    c = Atom("H", 1, [1, 0], fullshell=2)
    d = Atom("H", 1, [1, 2], fullshell=2)
    e = Atom("H", 1, [2, 1], fullshell=2)
    a.bond(b)
    a.bond(c)
    a.bond(d)
    a.bond(e)
    print(a)
    print(a.bonds)

def test3():
    molecules = []
    molecules.append(Atom("C", 4, [1, 1]))
    for index in range(4):
        molecules.append(Atom("H", 1, [index % 2, index // 2], fullshell=2))
        molecules[-1].bond(molecules[0])
    carbons = []
    hydrogens = []
    for index in range(3):
        carbons.append(Atom("C", 4, [6 + index*2, 2]))
        if len(carbons) > 1:
            carbons[-1].bond(carbons[-2])
        for index2 in range(2):
            hydrogens.append(Atom("H", 1, [6 + index*2, index2 * 4], fullshell=2))
            hydrogens[-1].bond(carbons[-1])
    hydrogens.append(Atom("H", 1, [12, 2], fullshell=2))
    hydrogens[-1].bond(carbons[-1])
    hydrogens.append(Atom("H", 1, [4, 2], fullshell=2))
    hydrogens[-1].bond(carbons[0])
    for atom_member in find_molecule(hydrogens[3]):
        print(f"{atom_member.letter} ({atom_member.coord_x}, {atom_member.coord_y})")
    print()
    for atom_member in find_molecule(molecules[1]):
        print(f"{atom_member.letter} ({atom_member.coord_x}, {atom_member.coord_y})")
    print()
    for atom_member in find_molecule(carbons[0]):
        print(f"{atom_member.letter} ({atom_member.coord_x}, {atom_member.coord_y})")

test3()
