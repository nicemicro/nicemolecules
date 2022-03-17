#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 07:27:35 2022

@author: nicemicro
"""

class Element():
    symbol: str = ""
    valence_el: int = 0
    fullshell: int = 8
    hypervalent: bool = False

class Hydrogen(Element):
    symbol = "H"
    valence_el = 1
    fullshell = 2

class Boron(Element):
    symbol = "B"
    valence_el = 3

class Carbon(Element):
    symbol = "C"
    valence_el = 4

class Nitrogen(Element):
    symbol = "N"
    valence_el = 5

class Oxygen(Element):
    symbol = "O"
    valence_el = 6

class Fluorine(Element):
    symbol = "F"
    valence_el = 7

class Silicon(Element):
    symbol = "Si"
    valence_el = 4
    hypervalent = True

class Phosphorus(Element):
    symbol = "P"
    valence_el = 5
    hypervalent = True

class Sulphur(Element):
    symbol = "S"
    valence_el = 6
    hypervalent = True

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

