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

MINDIST = 20
DEFLEN = 30
cos30 = DEFLEN * 3 ** (1/2) / 2
XSHIFT = [0, 0, DEFLEN, -DEFLEN,
          DEFLEN / 2, DEFLEN / 2, -DEFLEN / 2, -DEFLEN / 2,
          cos30, -cos30, cos30, -cos30]
YSHIFT = [DEFLEN, -DEFLEN, 0, 0,
          cos30, -cos30, cos30, -cos30,
          DEFLEN / 2, DEFLEN / 2, -DEFLEN / 2, -DEFLEN / 2]

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
    class Modes(IntEnum):
        NORMAL = auto()
        ADD_LINKED_ATOM = auto()
    
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        self.atoms = []
        self.bonds = []
        self.graphics = {}
        self.selected = None
        self.event_listened = False
        self.mode = self.Modes.NORMAL
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
        
        #self.atoms, self.bonds = eng.test3()
        #self.redraw_all_molecules(self.atoms, self.bonds)
    
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
            close = item in self.graphics and \
                isinstance(self.graphics[item], eng.Atom) and \
                -MINDIST < self.graphics[item].coord_x - event.x < MINDIST and \
                -MINDIST < self.graphics[item].coord_y - event.y < MINDIST
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
            self.possible_atoms(self.graphics[atom_s])
            self.mode =self.Modes.ADD_LINKED_ATOM
        #TODO: if the mouse is up right away, the helper lines don't disappear
    
    def mouseup_atom(self, atom_s, event):
        if self.mode == self.Modes.ADD_LINKED_ATOM:
            closestitems = self.mol_canvas.find_closest(event.x, event.y)
            self.mol_canvas.itemconfigure("ui_help", fill="grey")
            if len(closestitems) == 0:
                return
            item = closestitems[0]
            tags = self.mol_canvas.gettags(item)
            if len(tags) == 0 or tags[0] != "ui_help":
                return
            new_x = int(self.mol_canvas.coords(item)[0])
            new_y = int(self.mol_canvas.coords(item)[1])
            atom_connect = self.graphics[atom_s]
            new_atom = self.add_atom(self.toolbar.what_add(), [new_x, new_y])
            new_bond = new_atom.bond(atom_connect)
            self.bonds.append(new_bond)
            self.redraw_all_molecules(self.atoms, self.bonds)
            #TODO: even if the cursor is on line, we should add the new atom at the end of the line!
        self.mol_canvas.delete("ui_help")
        self.redraw_all_molecules(self.atoms, self.bonds)
        self.mode = self.Modes.NORMAL
    
    def drag_atom(self, atom_s, event):
        if self.mode == self.Modes.ADD_LINKED_ATOM:
            closestitems = self.mol_canvas.find_closest(event.x, event.y)
            self.mol_canvas.itemconfigure("ui_help", fill="grey")
            if len(closestitems) == 0:
                return
            item = closestitems[0]
            tags = self.mol_canvas.gettags(item)
            if len(tags) == 0 or tags[0] != "ui_help":
                return
            self.mol_canvas.itemconfigure(tags[1], fill="blue")
    
    def possible_atoms(self, atom):
        num = 1
        for (deltax, deltay) in zip(XSHIFT, YSHIFT):
            self.draw_bond(atom.coord_x, atom.coord_y, atom.coord_x+deltax,
                           atom.coord_y+deltay, 30, color="grey",
                           tags=("ui_help", f"ui_help_{num}"))
            self.mol_canvas.create_text(atom.coord_x+deltax, atom.coord_y+deltay,
                                        text=self.toolbar.what_add(),
                                        justify="center", fill="grey",
                                        tags=("ui_help", f"ui_help_{num}"))
            num += 1
    
    def add_atom(self, atom_symbol, coords):
        new_atom = eng.add_atom_by_symbol(atom_symbol, coords)
        assert not new_atom is None, "Unknown symbol probably"
        self.atoms.append(new_atom)
        self.redraw_all_molecules(self.atoms, self.bonds)
        return new_atom

    def draw_bond(self, x1, y1, x2, y2, bondlen, order=1, dativity=0,
                  color="black", tags=None):
        if tags is None:
            tags = ()
        obscure = 7
        shift = 2
        bondshifts = list(range(1-order, order+1, 2))
        bond_objects = []
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
            bondid = self.mol_canvas.create_line(xd1+xs, yd1+ys, xd2+xs, yd2+ys,
                                                 arrow=arrowtype, fill=color,
                                                 tags=tags)
            bond_objects.append(bondid)
        return bond_objects
    
    def redraw_all_molecules(self, atomlist, bondlist):
        self.mol_canvas.delete("all")
        self.graphics = {}
        for bond in bondlist:
            x1, y1, x2, y2 = bond.coords()
            bondlen = bond.length()
            if bondlen == 0:
                continue
            bond_order = int(bond.order())
            dativity = bond.dativity()
            bond_drawings = self.draw_bond(x1, y1, x2, y2, bondlen, 
                                           bond_order, dativity)
            for bondid in bond_drawings:
                self.graphics[bondid] = bond
        for atom in atomlist:
            atom_s = self.mol_canvas.create_text(atom.coord_x, atom.coord_y,
                                                 text=atom.symbol,
                                                 justify="center")
            self.graphics[atom_s] = atom
            self.mol_canvas.tag_bind(atom_s, "<Button-1>",
                                     lambda event, atom_s=atom_s:
                                         self.leftclick_atom(atom_s, event))
            self.mol_canvas.tag_bind(atom_s, "<ButtonRelease-1>",
                                     lambda event, atom_s=atom_s:
                                         self.mouseup_atom(atom_s, event))
            self.mol_canvas.tag_bind(atom_s, "<B1-Motion>",
                                     lambda event, atom_s=atom_s:
                                         self.drag_atom(atom_s, event))
        
def main():
    app = AppContainer()
    app.mainloop()

if __name__ == '__main__':
    main()
