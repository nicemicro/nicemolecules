# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 15:52:19 2022

@author: nicemicro
"""

class atom:
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
    
    def bond(self, other_atom, order=1, dative=0):
        self.bonds.append({"atom": other_atom,
                           "order": order,
                           "dative": dative})
        if not other_atom.is_bonded(self):
            other_atom.bond(self, order, -dative)
        return True
    
    def is_bonded(self, other_atom):
        return other_atom in (bond["atom"] for bond in self.bonds)
    
    def electrons(self):
        from_bonds = 0
        for bond in self.bonds:
            from_bonds += bond["order"]
            from_bonds += bond["dative"]
        return from_bonds + self.valence