#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 07:27:35 2022

@author: nicemicro
"""
import mol_eng as eng
import tkinter as tk
from tkinter import ttk
from enum import IntEnum, auto

class TopToolbar(ttk.Frame):
    class Modes(IntEnum):
        SELECT = auto()
        ADD_ATOM = auto()
    
    def __init__(self, parent, controller):
        ttk.Frame.__init__(self, parent)
        self.controller = controller
        self.mode = self.Modes.SELECT
        self.atom = None
        ttk.Button(self, text="Select", command=self.set_select) \
            .grid(row=0, column=0, columnspan=1)
        ttk.Button(self, text="C", command=lambda: self.set_add("C")) \
            .grid(row=0, column=1, columnspan=1)
        ttk.Button(self, text="O", command=lambda: self.set_add("O")) \
            .grid(row=0, column=2, columnspan=1)
    
    def is_select(self):
        return self.mode == self.Modes.SELECT
    
    def what_add(self):
        if self.mode == self.Modes.ADD_ATOM:
            return self.atom
        return False
    
    def set_select(self):
        self.mode = self.Modes.SELECT
        self.atom = None
        
    def set_add(self, atom):
        self.mode = self.Modes.ADD_ATOM
        self.atom = atom

class AppContainer(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        self.atoms = []
        self.bonds = []
        self.graphics = {}
        self.selected = None
        self.event_listened = False
        self.title("Nice Molecules")
        
        container = ttk.Frame(self)
        container.grid(row=0, column=0, sticky="nsew")
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        
        self.toolbar = TopToolbar(container, self)
        self.toolbar.grid(row=0, column=0, sticky="nsew")
        
        self.mol_canvas = tk.Canvas(container, width=700, height=600)
        self.mol_canvas.grid(row=1, column=0, sticky="nsew")
        self.mol_canvas.bind("<Button-1>", self.leftclick_canvas)
    
    def leftclick_canvas(self, event):
        if self.event_listened:
            self.event_listened = False
            return
        if self.toolbar.is_select():
            return
        if self.toolbar.what_add():
            if not self.selected is None:
                self.mol_canvas.itemconfigure(self.selected, fill="black")
                self.selected = None
            closestitems = self.mol_canvas.find_closest(event.x, event.y)
            if len(closestitems) == 0:
                item = None
            else:
                item = closestitems[0]
            print(item)
            close = item in self.graphics and \
                isinstance(self.graphics[item], eng.Atom) and \
                -20 < self.graphics[item].coord_x - event.x < 20 and \
                -20 < self.graphics[item].coord_y - event.y < 20
            if close:
                return
            coords = [event.x, event.y]
            self.add_atom(self.toolbar.what_add(), coords)
    
    def leftclick_atom(self, atom_s, event):
        if self.toolbar.is_select():
            if not self.selected is None:
                self.mol_canvas.itemconfigure(self.selected, fill="black")
            self.selected = atom_s
            self.mol_canvas.itemconfigure(atom_s, fill="green")
            self.event_listened = True
            return
        if self.toolbar.what_add():
            self.event_listened = True
    
    def add_atom(self, atom_symbol, coords):
        new_atom = eng.add_atom_by_symbol(atom_symbol, coords)
        assert not new_atom is None, "Unknown symbol probably"
        self.atoms.append(new_atom)
        self.redraw_all_molecules(self.atoms, self.bonds)
    
    def redraw_all_molecules(self, atomlist, bondlist):
        self.mol_canvas.delete("all")
        self.graphics = {}
        obscure = 7
        shift = 2
        for bond in bondlist:
            x1, y1, x2, y2 = bond.coords()
            bondlen = bond.length()
            if bondlen == 0:
                continue
            bond_order = int(bond.order())
            dativity = bond.dativity()
            bondshifts = list(range(1-bond_order, bond_order+1, 2))
            for bondpart in bondshifts:
                xd1 = x1 + (obscure / bondlen) * (x2-x1)
                xd2 = x2 - (obscure / bondlen) * (x2-x1)
                yd1 = y1 + (obscure / bondlen) * (y2-y1)
                yd2 = y2 - (obscure / bondlen) * (y2-y1)
                xs = (y1 - y2) / bondlen * shift * bondpart
                ys = (x2 - x1) / bondlen * shift * bondpart
                if dativity == 0:
                    arrowtype = None
                elif dativity < 0:
                    arrowtype = "first"
                    dativity += 1
                elif dativity > 0:
                    arrowtype = "last"
                    dativity -= 1
                bondid = self.mol_canvas.create_line(xd1+xs, yd1+ys,
                                                     xd2+xs, yd2+ys,
                                                     arrow=arrowtype)
                self.graphics[bondid] = bond
        for atom in atomlist:
            atom_s = self.mol_canvas.create_text(atom.coord_x, atom.coord_y,
                                                 text=atom.symbol,
                                                 justify="center")
            self.graphics[atom_s] = atom
            self.mol_canvas.tag_bind(atom_s, "<Button-1>",
                                     lambda event, atom_s=atom_s:
                                         self.leftclick_atom(atom_s, event))
        
def main():
    app = AppContainer()
    app.mainloop()

if __name__ == '__main__':
    main()
