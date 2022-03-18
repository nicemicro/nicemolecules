#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 07:27:35 2022

@author: nicemicro
"""
import mol_eng as eng
import mol_tst as tst
import tkinter as tk
from tkinter import ttk
from math import sqrt
from enum import IntEnum, auto
from typing import Optional, Sequence, Union

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
        ADD_BOND = auto()
    
    def __init__(self, parent: ttk.Frame, controller: tk.Tk) -> None:
        ttk.Frame.__init__(self, parent)
        self.controller: tk.Tk = controller
        self.mode: int = self.Modes.ADD_ATOM
        self.symbol: str = "C"
        self.bond: list[int] = [1, 0]
        mode_row = ttk.Frame(self, padding="0 0 0 5")
        mode_row.grid(row=0, column=0, columnspan=2, sticky="nsew")
        symbol_sel = ttk.Frame(self, padding="0 0 0 5")
        symbol_sel.grid(row=1, column=0, sticky="nsw")
        bond_sel = ttk.Frame(self, padding="10 0 0 5")
        bond_sel.grid(row=1, column=1, sticky="nse")
        ttk.Button(mode_row, text="Add atom",
                command=lambda: self.set_mode(self.Modes.ADD_ATOM)) \
                .grid(row=0, column=0, columnspan=1)
        ttk.Button(mode_row, text="Connect atoms",
                command=lambda: self.set_mode(self.Modes.ADD_BOND)) \
                .grid(row=0, column=1, columnspan=1)
        ttk.Button(mode_row, text="Select",
                command=lambda: self.set_mode(self.Modes.SELECT)) \
                .grid(row=0, column=2, columnspan=1)
        ttk.Button(symbol_sel, text="C", command=lambda: self.set_symbol("C")) \
            .grid(row=0, column=0, columnspan=1)
        ttk.Button(symbol_sel, text="O", command=lambda: self.set_symbol("O")) \
            .grid(row=0, column=1, columnspan=1)
        ttk.Button(bond_sel, text="--", command=lambda: self.set_bond([1, None])) \
            .grid(row=0, column=0, columnspan=1)
        ttk.Button(bond_sel, text="==", command=lambda: self.set_bond([2, None])) \
            .grid(row=0, column=1, columnspan=1)
    
    def is_select(self) -> bool:
        return self.mode == self.Modes.SELECT
    
    def is_add(self) -> bool:
        return self.mode == self.Modes.ADD_ATOM
    
    def is_connect(self) -> bool:
        return self.mode == self.Modes.ADD_BOND
    
    def atom_symbol(self) -> str:
        return self.symbol
    
    def bond_order(self) -> int:
        return self.bond[0]
    
    def bond_dativity(self) -> int:
        return self.bond[1]
    
    def set_mode(self, new_mode: int = -1) -> None:
        if new_mode == -1:
            new_mode = self.Modes.SELECT
        self.mode = new_mode
        
    def set_symbol(self, new_symbol: str) -> None:
        self.symbol = new_symbol
    
    def set_bond(self, bondtype: list[Optional[int]]) -> None:
        for index, item in enumerate(bondtype):
            if not item is None:
                self.bond[index] = item

class AppContainer(tk.Tk):
    class Modes(IntEnum):
        NORMAL = auto()
        ADD_LINKED_ATOM = auto()
        LINK_ATOMS = auto()
    
    def __init__(self, *args, **kwargs) -> None:
        tk.Tk.__init__(self, *args, **kwargs)
        self.atoms: list[eng.Atom ]= []
        self.bonds: list[eng.Cov_bond] = []
        self.graphics: dict[int, Union[eng.Atom, eng.Cov_bond]] = {}
        self.selected: Optional[int] = None
        self.event_listened: bool = False
        self.mode = self.Modes.NORMAL
        self.title("Nice Molecules")
        
        container = ttk.Frame(self)
        container.grid(row=0, column=0, sticky="nsew")
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        
        self.toolbar = TopToolbar(container, self)
        self.toolbar.grid(row=0, column=0, sticky="nsew")
        
        self.mol_canvas = tk.Canvas(container, width=700, height=600,
                                    bg="white")
        self.mol_canvas.grid(row=1, column=0, sticky="nsew")
        self.mol_canvas.bind("<Button-1>", self.leftclick_canvas)
        
        #self.atoms, self.bonds = tst.SF6()
        #self.redraw_all_molecules(self.atoms, self.bonds)
    
    def leftclick_canvas(self, event: tk.Event) -> None:
        if self.event_listened:
            self.event_listened = False
            return
        if self.toolbar.is_select():
            return
        if self.toolbar.is_add():
            if not self.selected is None:
                self.mol_canvas.itemconfigure(self.selected, fill="black")
                self.selected = None
            closestitems: tuple[int, ...] = self.mol_canvas.find_closest(event.x, event.y)
            if len(closestitems) == 0:
                item = None
            else:
                item = closestitems[0]
            if item in self.graphics and isinstance(self.graphics[item], eng.Atom):
                item_obj = self.graphics[item]
                assert isinstance(item_obj, eng.Atom)
                if -MINDIST < item_obj.coord_x - event.x < MINDIST \
                        and -MINDIST < item_obj.coord_y - event.y < MINDIST:
                    return
            coords: list[int] = [event.x, event.y]
            self.add_atom(self.toolbar.atom_symbol(), coords)
    
    def leftdown_atom(self, sel_atom_num: int, event: tk.Event) -> None:
        sel_atom = self.graphics[sel_atom_num]
        assert isinstance(sel_atom, eng.Atom), "The selection should have been an atom"
        if self.toolbar.is_select():
            if not self.selected is None:
                self.mol_canvas.itemconfigure(self.selected, fill="black")
            self.selected = sel_atom_num
            self.mol_canvas.itemconfigure(sel_atom_num, fill="green")
            self.event_listened = True
            return
        if self.toolbar.is_add():
            self.event_listened = True
            test_atom = eng.add_atom_by_symbol(self.toolbar.atom_symbol(), [0, 0])
            if sel_atom.can_bond(test_atom, self.toolbar.bond_order(),
                    self.toolbar.bond_dativity()) != eng.BondingError.OK:
                return
            self.possible_atoms(sel_atom)
            self.mode = self.Modes.ADD_LINKED_ATOM
            return
        if self.toolbar.is_connect():
            self.event_listened = True
            self.possible_links(sel_atom)
            self.mode = self.Modes.LINK_ATOMS
    
    def close_selection(self, closestitems: tuple[int, ...]) \
            -> tuple[Optional[int], Optional[Sequence[str]]]:
        if len(closestitems) == 0:
            self.set_normal_mode()
            return None, None
        item: int = closestitems[0]
        tags: Sequence[str] = self.mol_canvas.gettags(item)
        if ((len(tags) < 2 or tags[0] != "ui_help") and 
                (len(tags) < 1 or tags[0] != "atom")):
            return None, None
        return item, tags
    
    def mouseup_atom(self, atom_s: int, event: tk.Event) -> None:
        if self.mode == self.Modes.ADD_LINKED_ATOM:
            closestitems = self.mol_canvas.find_closest(event.x, event.y)
            self.mol_canvas.itemconfigure("ui_help", fill="grey")
            item, tags = self.close_selection(closestitems)
            if item is None or tags is None or tags[0] != "ui_help":
                self.set_normal_mode()
                return
            atomplace = "atom_here-" + tags[1].split("-")[-1]
            new_x = int(self.mol_canvas.coords(atomplace)[0])
            new_y = int(self.mol_canvas.coords(atomplace)[1])
            atom_connect = self.graphics[atom_s]
            assert isinstance(atom_connect, eng.Atom)
            new_atom = self.add_atom(self.toolbar.atom_symbol(), [new_x, new_y])
            new_bond = atom_connect.bond(new_atom, order=self.toolbar.bond_order())
            self.bonds.append(new_bond)
            self.redraw_all_molecules(self.atoms, self.bonds)
            self.set_normal_mode()
        if self.mode == self.Modes.LINK_ATOMS:
            closestitems = self.mol_canvas.find_closest(event.x, event.y)
            self.mol_canvas.itemconfigure("ui_help", fill="#aaaaaa")
            self.mol_canvas.itemconfigure("atom", fill="black")
            item, tags = self.close_selection(closestitems)
            if item is None or tags is None:
                self.set_normal_mode()
                return
            selected_atom = -1
            if len(tags) >= 2 and tags[0] == "ui_help":
                selected_atom = int(tags[1].split("-")[-1])
            if len(tags) >= 1 and tags[0] == "atom":
                selected_atom = item
            if selected_atom == -1:
                self.set_normal_mode()
                return
            atom_from = self.graphics[atom_s]
            atom_to = self.graphics[selected_atom]
            assert isinstance(atom_from, eng.Atom), "Trying to bobd to a non-atom"
            assert isinstance(atom_to, eng.Atom), "Trying to bobd to a non-atom"
            if atom_from.can_bond(atom_to, self.toolbar.bond_order(),
                    self.toolbar.bond_dativity()) != eng.BondingError.OK:
                self.set_normal_mode()
                return
            new_bond = atom_from.bond(atom_to, order=self.toolbar.bond_order(),
                    dative=self.toolbar.bond_dativity())
            assert not new_bond is None, "Creating a bond failed"
            self.bonds.append(new_bond)
            self.redraw_all_molecules(self.atoms, self.bonds)
            self.set_normal_mode()
    
    def drag_atom(self, atom_s: int, event: tk.Event) -> None:
        if (self.mode == self.Modes.ADD_LINKED_ATOM or
                self.mode == self.Modes.LINK_ATOMS):
            closestitems = self.mol_canvas.find_closest(event.x, event.y)
            self.mol_canvas.itemconfigure("ui_help", fill="#aaaaaa")
            self.mol_canvas.itemconfigure("atom", fill="black")
            item, tags = self.close_selection(closestitems)
            if item is None or tags is None:
                return
            if self.mode == self.Modes.ADD_LINKED_ATOM:
                if len(tags) < 2 or tags[0] != "ui_help":
                    return
                self.mol_canvas.itemconfigure(tags[1], fill="blue")
            elif self.mode == self.Modes.LINK_ATOMS:
                selected_atom = -1
                if len(tags) >= 2 and tags[0] == "ui_help":
                    selected_atom = int(tags[1].split("-")[-1])
                if len(tags) >= 1 and tags[0] == "atom":
                    selected_atom = item
                if selected_atom == -1:
                    return
                atom_from = self.graphics[selected_atom]
                atom_to = self.graphics[atom_s]
                assert isinstance(atom_from, eng.Atom)
                assert isinstance(atom_to, eng.Atom)
                if atom_from.can_bond(atom_to, self.toolbar.bond_order(),
                        self.toolbar.bond_dativity()) != eng.BondingError.OK:
                    return
                self.mol_canvas.itemconfigure(f"ui_help-{selected_atom}",
                                              fill="blue")
                self.mol_canvas.itemconfigure(selected_atom, fill ="blue")
        
    def set_normal_mode(self) -> None:
        self.mol_canvas.delete("ui_help")
        self.mode = self.Modes.NORMAL
    
    def possible_links(self, atom_from: eng.Atom) -> None:
        for atom_link_num in self.graphics:
            atom_to = self.graphics[atom_link_num]
            if not isinstance(atom_to, eng.Atom):
                continue
            if atom_from.can_bond(atom_to, self.toolbar.bond_order(),
                    self.toolbar.bond_dativity()) != eng.BondingError.OK:
                continue
            self.draw_bond(atom_from.coord_x, atom_from.coord_y,
                    atom_to.coord_x, atom_to.coord_y, -1, color="#aaaaaa",
                    tags=("ui_help", f"ui_help-{atom_link_num}"))
    
    def possible_atoms(self, atom: eng.Atom) -> None:
        num = 1
        for (deltax, deltay) in zip(XSHIFT, YSHIFT):
            (x1, y1, x2, y2) = (atom.coord_x + deltax - MINDIST/2,
                                atom.coord_y + deltay - MINDIST/2,
                                atom.coord_x + deltax + MINDIST/2,
                                atom.coord_y + deltay + MINDIST/2)
            over = self.mol_canvas.find_overlapping(x1, y1, x2, y2)
            over_atom = False
            for obj in over:
                tags = self.mol_canvas.gettags(obj)
                if len(tags) > 0 and ("atom" in tags or "bond" in tags):
                    over_atom = True
                    break
            if over_atom:
                num += 1
                continue
            self.draw_bond(atom.coord_x, atom.coord_y, atom.coord_x+deltax,
                           atom.coord_y+deltay, 30, color="#aaaaaa",
                           tags=("ui_help", f"ui_help-{num}"))
            self.mol_canvas.create_text(atom.coord_x+deltax, atom.coord_y+deltay,
                                        text=self.toolbar.atom_symbol(),
                                        justify="center", fill="#aaaaaa",
                                        tags=("ui_help", f"ui_help-{num}",
                                              f"atom_here-{num}"))
            num += 1
    
    def add_atom(self, atom_symbol: str, coords: Sequence[float]) -> eng.Atom:
        new_atom = eng.add_atom_by_symbol(atom_symbol, coords)
        assert not new_atom is None, "Unknown symbol probably"
        self.atoms.append(new_atom)
        self.redraw_all_molecules(self.atoms, self.bonds)
        return new_atom

    def draw_bond(self, x1: float, y1: float, x2: float, y2: float,
                  bondlen: float=-1, order: float=-1, dativity: float=0,
                  color: str="black", tags: Optional[tuple[str, ...]]=None
                  ) -> list[int]:
        if tags is None:
            tags = ()
        if bondlen == -1:
            bondlen = sqrt((x2 - x1)**2 + (y2 - y1)**2)
        if order == -1:
            order = self.toolbar.bond_order()
        obscure = 7
        shift = 2
        bondshifts = list(range(1-int(order), int(order)+1, 2))
        bond_objects: list[int] = []
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
    
    def redraw_all_molecules(self, atomlist: list[eng.Atom],
                             bondlist: list[eng.Cov_bond]) -> None:
        self.mol_canvas.delete("all")
        self.graphics = {}
        for bond in bondlist:
            x1, y1, x2, y2 = bond.coords()
            bondlen = bond.length()
            if bondlen == 0:
                continue
            bond_order = bond.order()
            dativity = bond.dativity()
            bond_drawings = self.draw_bond(x1, y1, x2, y2, bondlen, 
                                           bond_order, dativity, tags=("bond", ))
            for bondid in bond_drawings:
                self.graphics[bondid] = bond
        for atom in atomlist:
            atom_s: int = self.mol_canvas.create_text(atom.coord_x, atom.coord_y,
                                                     text=atom.symbol,
                                                     justify="center",
                                                     tags=("atom"))
            self.graphics[atom_s] = atom
            #self.mol_canvas.addtag_above(f"atom-{atom_s}", atom_s)
            self.mol_canvas.tag_bind(atom_s, "<ButtonPress-1>",
                                     lambda event, atom_s=atom_s:
                                         self.leftdown_atom(atom_s, event))
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
