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
COS30 = DEFLEN * 3**(1 / 2) / 2
XSHIFT = [
    0, 0, DEFLEN, -DEFLEN,
    DEFLEN / 2, DEFLEN / 2, -DEFLEN / 2, -DEFLEN / 2,
    COS30, -COS30, COS30, -COS30
]
YSHIFT = [
    DEFLEN, -DEFLEN, 0, 0,
    COS30, -COS30, COS30, -COS30, DEFLEN / 2,
    DEFLEN / 2, -DEFLEN / 2, -DEFLEN / 2
]


class TopToolbar(ttk.Frame):
    """The top toolbar of the window: contains controls for the
    program, and stores editor information"""

    class Modes(IntEnum):
        """The selectable operations."""
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
        """Changes the operation mode of the application."""
        if new_mode == -1:
            new_mode = self.Modes.SELECT
        self.mode = new_mode

    def set_symbol(self, new_symbol: str) -> None:
        """Sets the default atom symbol for the application"""
        self.event_generate("<<AtomButtonPress>>", when="tail")
        self.symbol = new_symbol

    def set_bond(self, bondtype: list[Optional[int]]) -> None:
        """Sets the bond type: [order, dativity]"""
        for index, item in enumerate(bondtype):
            if item is not None:
                self.bond[index] = item


class AppContainer(tk.Tk):
    """The main application window."""

    class Modes(IntEnum):
        """UI input mdes of the canvas."""
        NORMAL = auto()
        ADD_LINKED_ATOM = auto()
        LINK_ATOMS = auto()

    def __init__(self, *args, **kwargs) -> None:
        tk.Tk.__init__(self, *args, **kwargs)
        self.atoms: list[eng.Atom] = []
        self.bonds: list[eng.CovBond] = []
        self.graphics: dict[int, Union[eng.Atom, eng.CovBond]] = {}
        self.selected: list[int] = []
        self.event_listened: bool = False
        self.mode = self.Modes.NORMAL
        self.title("Nice Molecules")

        container = ttk.Frame(self)
        container.grid(row=0, column=0, sticky="nsew")
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        self.toolbar = TopToolbar(container, self)
        self.toolbar.grid(row=0, column=0, sticky="nsew")

        self.mol_canvas = tk.Canvas(container,
                                    width=700,
                                    height=600,
                                    bg="white")
        self.mol_canvas.grid(row=1, column=0, sticky="nsew")
        self.mol_canvas.bind("<Button-1>", self.leftclick_canvas)
        self.bind("<Delete>", self.delkey_pressed)
        self.bind("<<AtomButtonPress>>", self.atom_button_pressed)

    def atom_button_pressed(self, event: tk.Event) -> None:
        """Handles the pressing of the atom symbol buttons in the toolbar:
            - in normal mode, if atoms are selected, those are changed to
              the symbol what has been selected on the toolbar."""
        if self.mode == self.Modes.NORMAL:
            changes: bool = False
            change_to: Optional[eng.el.Element]
            change_to = eng.element_by_symbol(self.toolbar.atom_symbol())
            if change_to is None:
                return
            for sel_num in self.selected:
                if sel_num not in self.graphics:
                    continue
                sel_item = self.graphics[sel_num]
                if isinstance(sel_item, eng.Atom):
                    err: eng.BondingError = sel_item.can_change_element(
                        change_to)
                    if err == eng.BondingError.OK:
                        sel_item.change_element(change_to)
                        changes = True
            if changes:
                self.selected = []
                self.redraw_all()

    def delkey_pressed(self, event: tk.Event) -> None:
        """Handles the event of the delete key being pressed:
            - in normal mode, if there are elements selected, the selected
              bonds removed, or atoms and bonds of atoms removed."""
        if self.mode == self.Modes.NORMAL and self.selected:
            for sel_num in self.selected:
                if sel_num not in self.graphics:
                    continue
                sel_item = self.graphics[sel_num]
                if isinstance(sel_item, eng.Atom):
                    removed_bonds = sel_item.unbond_all()
                    self.graphics = {
                        key: obj
                        for key, obj in self.graphics.items()
                        if obj != sel_item
                    }
                    self.atoms.remove(sel_item)
                    for removed_bond in removed_bonds:
                        self.graphics = {
                            key: obj
                            for key, obj in self.graphics.items()
                            if obj != removed_bond
                        }
                        self.bonds.remove(removed_bond)
                elif isinstance(sel_item, eng.CovBond):
                    sel_item.delete()
                    self.graphics = {
                        key: obj
                        for key, obj in self.graphics.items()
                        if obj != sel_item
                    }
                    self.bonds.remove(sel_item)
            self.selected = []
            self.redraw_all()

    def leftclick_canvas(self, event: tk.Event) -> None:
        """Handes left click on the canvas:
            - if the operation is adding an atom, adds a new atom at the
              position of mouse click"""
        if self.event_listened:
            self.event_listened = False
            return
        if self.toolbar.is_select():
            return
        if self.toolbar.is_add():
            if not len(self.selected) == 0:
                for sel_item in self.selected:
                    self.mol_canvas.itemconfigure(sel_item, fill="black")
                self.selected = []
            closestitems: tuple[int, ...] = self.mol_canvas.find_closest(
                event.x, event.y)
            if len(closestitems) == 0:
                item = None
            else:
                item = closestitems[0]
            if item in self.graphics and isinstance(self.graphics[item],
                                                    eng.Atom):
                item_obj = self.graphics[item]
                assert isinstance(item_obj, eng.Atom)
                if -MINDIST < item_obj.coord_x - event.x < MINDIST \
                        and -MINDIST < item_obj.coord_y - event.y < MINDIST:
                    return
            coords: list[int] = [event.x, event.y]
            self.add_atom(self.toolbar.atom_symbol(), coords)

    def leftdown_atom(self, sel_atom_num: int, _event: tk.Event) -> None:
        """Handles the event of mouse button press on atom, based on the
           operation:
           - if selection, current atom is selected,
           - if add atom, enters mode to add linked atom,
           - if connect atoms, enters mode to connect existing atoms."""
        sel_atom = self.graphics[sel_atom_num]
        assert isinstance(sel_atom,
                          eng.Atom), "The selection should have been an atom"
        if self.toolbar.is_select():
            if not len(self.selected) == 0:
                for sel_item in self.selected:
                    self.mol_canvas.itemconfigure(sel_item, fill="black")
            self.selected = [sel_atom_num]
            self.mol_canvas.itemconfigure(sel_atom_num, fill="green")
            self.event_listened = True
            return
        if self.toolbar.is_add():
            self.event_listened = True
            test_atom = eng.add_atom_by_symbol(self.toolbar.atom_symbol(),
                                               [0, 0])
            if sel_atom.can_bond(
                    test_atom, self.toolbar.bond_order(),
                    self.toolbar.bond_dativity()) != eng.BondingError.OK:
                return
            self.possible_atoms(sel_atom)
            self.mode = self.Modes.ADD_LINKED_ATOM
            return
        if self.toolbar.is_connect():
            self.event_listened = True
            self.possible_links(sel_atom)
            self.mode = self.Modes.LINK_ATOMS

    def mouseup_atom(self, atom_s: int, event: tk.Event) -> None:
        """Handles the events when the mouse is lifted after draging an atom:
            - adding linked atom: adds new atom and new bond,
            - linking atoms: adds new bond."""
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
            new_atom = self.add_atom(self.toolbar.atom_symbol(),
                                     [new_x, new_y])
            new_bond = atom_connect.bond(new_atom,
                                         order=self.toolbar.bond_order())
            self.bonds.append(new_bond)
            self.redraw_all()
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
            assert isinstance(atom_from,
                              eng.Atom), "Trying to bobd to a non-atom"
            assert isinstance(atom_to,
                              eng.Atom), "Trying to bobd to a non-atom"
            if atom_from.can_bond(
                    atom_to, self.toolbar.bond_order(),
                    self.toolbar.bond_dativity()) != eng.BondingError.OK:
                self.set_normal_mode()
                return
            new_bond = atom_from.bond(atom_to,
                                      order=self.toolbar.bond_order(),
                                      dative=self.toolbar.bond_dativity())
            assert new_bond is not None, "Creating a bond failed"
            self.bonds.append(new_bond)
            self.redraw_all()
            self.set_normal_mode()

    def drag_atom(self, atom_s: int, event: tk.Event) -> None:
        """Handling the dragging of an atom that was clicked on:
            - adding linked atom / linking atoms: highlights the
              placeholder where new atom / bond will be created."""
        if self.mode in (self.Modes.ADD_LINKED_ATOM,
                self.Modes.LINK_ATOMS):
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
                if atom_from.can_bond(
                        atom_to, self.toolbar.bond_order(),
                        self.toolbar.bond_dativity()) != eng.BondingError.OK:
                    return
                self.mol_canvas.itemconfigure(f"ui_help-{selected_atom}",
                                              fill="blue")
                self.mol_canvas.itemconfigure(selected_atom, fill="blue")

    def click_bond(self, bond_line_id: int, _event: tk.Event) -> None:
        """Handling the event of a bond being clicked on:
            - if operation is selection, selects the bond, and all lines
              that constitute the bond."""
        sel_bond = self.graphics[bond_line_id]
        assert isinstance(
            sel_bond, eng.CovBond), "The selection should have been an atom"
        if self.toolbar.is_select():
            tags: Sequence[str] = self.mol_canvas.gettags(bond_line_id)
            if len(tags) < 0 or tags[0] != "bond":
                return
            bond_id: int = int(tags[1].split("-")[1])
            if not len(self.selected) == 0:
                for sel_item in self.selected:
                    self.mol_canvas.itemconfigure(sel_item, fill="black")
            self.selected = list(
                self.mol_canvas.find_withtag(f"bond-{bond_id}"))
            self.mol_canvas.itemconfigure(f"bond-{bond_id}", fill="green")
            self.event_listened = True
            return

    def close_selection(self, closestitems: tuple[int, ...]) \
            -> tuple[Optional[int], Optional[Sequence[str]]]:
        """Finds the closest item on the canvas to the
        listed item on the canvas."""
        if len(closestitems) == 0:
            self.set_normal_mode()
            return None, None
        item: int = closestitems[0]
        tags: Sequence[str] = self.mol_canvas.gettags(item)
        if ((len(tags) < 2 or tags[0] != "ui_help")
                and (len(tags) < 1 or tags[0] != "atom")):
            return None, None
        return item, tags

    def set_normal_mode(self) -> None:
        """Returns the application to normal mode."""
        self.mol_canvas.delete("ui_help")
        self.mode = self.Modes.NORMAL

    def possible_links(self, atom_from: eng.Atom) -> None:
        """Draws possible links between atoms to create in
        atom connecting mode."""
        for atom_link_num in self.graphics:
            atom_to = self.graphics[atom_link_num]
            if not isinstance(atom_to, eng.Atom):
                continue
            if atom_from.can_bond(
                    atom_to, self.toolbar.bond_order(),
                    self.toolbar.bond_dativity()) != eng.BondingError.OK:
                continue
            self.draw_bond(atom_from.coord_x,
                           atom_from.coord_y,
                           atom_to.coord_x,
                           atom_to.coord_y,
                           -1,
                           color="#aaaaaa",
                           tags=("ui_help", f"ui_help-{atom_link_num}"))

    def possible_atoms(self, atom: eng.Atom) -> None:
        """Draws possible atoms that could be created and joined
        to the current atom in linked atom drawing mode."""
        num = 1
        for (deltax, deltay) in zip(XSHIFT, YSHIFT):
            (x1, y1, x2, y2) = (atom.coord_x + deltax - MINDIST / 2,
                                atom.coord_y + deltay - MINDIST / 2,
                                atom.coord_x + deltax + MINDIST / 2,
                                atom.coord_y + deltay + MINDIST / 2)
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
            self.draw_bond(atom.coord_x,
                           atom.coord_y,
                           atom.coord_x + deltax,
                           atom.coord_y + deltay,
                           30,
                           color="#aaaaaa",
                           tags=("ui_help", f"ui_help-{num}"))
            self.mol_canvas.create_text(atom.coord_x + deltax,
                                        atom.coord_y + deltay,
                                        text=self.toolbar.atom_symbol(),
                                        justify="center",
                                        fill="#aaaaaa",
                                        tags=("ui_help", f"ui_help-{num}",
                                              f"atom_here-{num}"))
            num += 1

    def add_atom(self, atom_symbol: str, coords: Sequence[float]) -> eng.Atom:
        """Adds an atom to the canvas with given symbol at
        given coordinates."""
        new_atom = eng.add_atom_by_symbol(atom_symbol, coords)
        assert new_atom is not None, "Unknown symbol probably"
        self.atoms.append(new_atom)
        self.redraw_all()
        return new_atom

    def draw_bond(self,
                  x1: float,
                  y1: float,
                  x2: float,
                  y2: float,
                  bondlen: float = -1,
                  order: float = -1,
                  dativity: float = 0,
                  color: str = "black",
                  tags: Optional[tuple[str, ...]] = None) -> list[int]:
        """Draws a bond between two atoms and returns the list of line
        IDs on the canvas that represent the chemical bond."""
        if tags is None:
            tags = ()
        if bondlen == -1:
            bondlen = sqrt((x2 - x1)**2 + (y2 - y1)**2)
        if order == -1:
            order = self.toolbar.bond_order()
        obscure = 7
        shift = 2
        bondshifts = list(range(1 - int(order), int(order) + 1, 2))
        bond_objects: list[int] = []
        for bondpart in bondshifts:
            xd1 = x1 + (obscure / bondlen) * (x2 - x1)
            xd2 = x2 - (obscure / bondlen) * (x2 - x1)
            yd1 = y1 + (obscure / bondlen) * (y2 - y1)
            yd2 = y2 - (obscure / bondlen) * (y2 - y1)
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
            bond_line_id = self.mol_canvas.create_line(xd1 + xs,
                                                       yd1 + ys,
                                                       xd2 + xs,
                                                       yd2 + ys,
                                                       arrow=arrowtype,
                                                       fill=color,
                                                       tags=tags)
            self.mol_canvas.tag_bind(bond_line_id,
                                     "<Button-1>",
                                     lambda event, bond_id=bond_line_id: self.
                                     click_bond(bond_line_id, event))
            bond_objects.append(bond_line_id)
        return bond_objects

    def redraw_all(self) -> None:
        """Redraws all atoms and bonds on a freshly cleared canvas."""
        self.mol_canvas.delete("all")
        self.graphics = {}
        bond_id: int = 0
        for bond in self.bonds:
            x1, y1, x2, y2 = bond.coords
            bondlen = bond.length
            if bondlen == 0:
                continue
            bond_order = bond.order
            dativity = bond.dativity
            bond_drawings = self.draw_bond(x1,
                                           y1,
                                           x2,
                                           y2,
                                           bondlen,
                                           bond_order,
                                           dativity,
                                           tags=("bond", f"bond-{bond_id}"))
            for bondid in bond_drawings:
                self.graphics[bondid] = bond
            bond_id += 1
        for atom in self.atoms:
            atom_s: int = self.mol_canvas.create_text(atom.coord_x,
                                                      atom.coord_y,
                                                      text=atom.symbol,
                                                      justify="center",
                                                      tags=("atom"))
            self.graphics[atom_s] = atom
            self.mol_canvas.tag_bind(
                atom_s,
                "<ButtonPress-1>",
                lambda event, atom_s=atom_s: self.leftdown_atom(atom_s, event))
            self.mol_canvas.tag_bind(
                atom_s,
                "<ButtonRelease-1>",
                lambda event, atom_s=atom_s: self.mouseup_atom(atom_s, event))
            self.mol_canvas.tag_bind(
                atom_s,
                "<B1-Motion>",
                lambda event, atom_s=atom_s: self.drag_atom(atom_s, event))


def main():
    app = AppContainer()
    app.mainloop()


if __name__ == '__main__':
    main()
