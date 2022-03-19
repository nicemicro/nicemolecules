#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 07:27:35 2022

@author: nicemicro
"""


class Element:
    """
    Abstract class that stores information about an element.
    Members:
        _symbol: str, the chemical symbol of an atom
        _valence_el: int, number of valence electrons
        _fullshell: int, the number of electrons in a filled valence shell
        _hypervalent: bool, whether the atom can accomodate more electrons
            then the full valence shell
    """
    _symbol: str = ""
    _valence_el: int = 0
    _fullshell: int = 8
    _hypervalent: bool = False

    def get_symbol(self):
        return self._symbol

    def get_valence_el(self):
        return self._valence_el

    def get_fullshell(self):
        return self._fullshell

    def get_hypervalent(self):
        return self._hypervalent

    def set_attr (self, new):
        raise AttributeError("Can't modify element properties")

    def del_attr(self):
        raise AttributeError("Can't delete element properties")

    symbol = property(get_symbol, set_attr, del_attr)
    valence_el = property(get_valence_el, set_attr, del_attr)
    fullshell = property(get_fullshell, set_attr, del_attr)
    hypervalent = property(get_hypervalent, set_attr, del_attr)

class Hydrogen(Element):
    _symbol = "H"
    _valence_el = 1
    _fullshell = 2


class Boron(Element):
    _symbol = "B"
    _valence_el = 3


class Carbon(Element):
    _symbol = "C"
    _valence_el = 4


class Nitrogen(Element):
    _symbol = "N"
    _valence_el = 5


class Oxygen(Element):
    _symbol = "O"
    _valence_el = 6


class Fluorine(Element):
    _symbol = "F"
    _valence_el = 7


class Silicon(Element):
    _symbol = "Si"
    _valence_el = 4
    _hypervalent = True


class Phosphorus(Element):
    _symbol = "P"
    _valence_el = 5
    _hypervalent = True


class Sulphur(Element):
    _symbol = "S"
    _valence_el = 6
    _hypervalent = True


element_table: list[Element] = []
element_table.append(Hydrogen())
element_table.append(Boron())
element_table.append(Carbon())
element_table.append(Nitrogen())
element_table.append(Oxygen())
element_table.append(Fluorine())
element_table.append(Silicon())
element_table.append(Phosphorus())
element_table.append(Sulphur())
