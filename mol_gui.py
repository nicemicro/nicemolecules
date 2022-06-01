#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 07:27:35 2022

@author: nicemicro
"""
import mol_eng as eng
import mol_tst as tst
import mol_opt as opt
import tkinter as tk
from tkinter import ttk
from tkinter import font as tkfont
from math import sqrt, cos, sin, pi, atan2, tan
from enum import IntEnum, auto
from typing import Optional, Sequence, Union, Any, Literal

MINDIST = 20
DEFLEN = 30
COS30 = DEFLEN * 3 ** (1 / 2) / 2
XSHIFT = [
    0,
    0,
    DEFLEN,
    -DEFLEN,
    DEFLEN / 2,
    DEFLEN / 2,
    -DEFLEN / 2,
    -DEFLEN / 2,
    COS30,
    -COS30,
    COS30,
    -COS30,
]
YSHIFT = [
    DEFLEN,
    -DEFLEN,
    0,
    0,
    COS30,
    -COS30,
    COS30,
    -COS30,
    DEFLEN / 2,
    DEFLEN / 2,
    -DEFLEN / 2,
    -DEFLEN / 2,
]


class TopToolbar(ttk.Frame):
    """The top toolbar of the window: contains controls for the
    program, and stores editor information"""

    class Modes(IntEnum):
        """The selectable operations."""

        SELECT = auto()
        ADD_ATOM = auto()
        ADD_BOND = auto()

    class EmptyValence(IntEnum):
        """Types of drawing atoms with empty valence shells"""

        NOTHING = auto()
        HYDROGENS = auto()
        ELECTRONS = auto()

    empty_val_names: dict[int, str]

    def get_is_select(self) -> bool:
        return self._mode == self.Modes.SELECT

    def get_is_add(self) -> bool:
        return self._mode == self.Modes.ADD_ATOM

    def get_is_connect(self) -> bool:
        return self._mode == self.Modes.ADD_BOND

    def get_atom_symbol(self) -> str:
        return self._symbol

    def get_bond_order(self) -> int:
        return self._bond[0]

    def get_bond_dativity(self) -> int:
        return self._bond[1]

    def get_empty_val_style(self) -> int:
        index: int = list(self.empty_val_names.values()).index(self._emptyvalence.get())
        return list(self.empty_val_names.keys())[index]

    def get_hide_carbon(self) -> bool:
        return self._hide_carbon.get() == 1

    def toss(self, new_value: Any) -> None:
        """Prohibits changing attributes."""
        raise AttributeError("This value can not be modified directly.")

    is_select = property(get_is_select, toss)
    is_add = property(get_is_add, toss)
    is_connect = property(get_is_connect, toss)
    atom_symbol = property(get_atom_symbol, toss)
    bond_order = property(get_bond_order, toss)
    bond_dativity = property(get_bond_dativity, toss)
    empty_val_style = property(get_empty_val_style, toss)
    hide_carbon = property(get_hide_carbon, toss)

    def __init__(self, parent: ttk.Frame, controller: tk.Tk) -> None:
        ttk.Frame.__init__(self, parent)
        #s = ttk.Style()
        #s.configure("test.TFrame", background="red")
        self.controller: tk.Tk = controller
        self._mode: int = self.Modes.ADD_ATOM
        self._symbol: str = "C"
        self._bond: list[int] = [1, 0]
        self.status_text = tk.StringVar()
        self.custom_atom_symbol = tk.StringVar()
        atom_symbols: list[list[str]] = [
            ["H", "C", "N", "O", "F", ""],
            ["", "", "P", "S", "Cl"],
            ["", "", "", "", "I", "Xe"],
        ]
        self.empty_val_names = {
            self.EmptyValence.NOTHING: "Nothing",
            self.EmptyValence.HYDROGENS: "Hydrogens",
            self.EmptyValence.ELECTRONS: "Electrons",
        }

        mode_row = ttk.Frame(self, padding="0 0 0 5")
        mode_row.grid(row=0, column=0, columnspan=1, sticky="nsew")
        mol_style = ttk.Frame(self, padding="10 0 0 5")
        mol_style.grid(row=0, column=1, columnspan=2, sticky="nse")
        symbol_sel = ttk.Frame(self, padding="0 0 0 5")
        symbol_sel.grid(row=1, column=0, columnspan=1, sticky="nsw")
        bond_sel = ttk.Frame(self, padding="10 0 0 5")
        bond_sel.grid(row=1, column=1, rowspan=2, sticky="nsw")
        ttk.Label(self, textvariable=self.status_text).grid(
            row=1, column=2, rowspan=1, sticky="nesw"
        )
        custom_add = ttk.Frame(self, padding="0 0 0 5")
        custom_add.grid(row=2, column=0, columnspan=1, sticky="nsew")
        self.status_text_update()

        ttk.Button(
            mode_row,
            text="Add atom",
            command=lambda: self.set_mode(self.Modes.ADD_ATOM),
        ).grid(row=0, column=0, columnspan=1)
        ttk.Button(
            mode_row,
            text="Connect atoms",
            command=lambda: self.set_mode(self.Modes.ADD_BOND),
        ).grid(row=0, column=1, columnspan=1)
        ttk.Button(
            mode_row, text="Select", command=lambda: self.set_mode(self.Modes.SELECT)
        ).grid(row=0, column=2, columnspan=1)

        self._hide_carbon = tk.IntVar()
        self._hide_carbon.set(0)
        ttk.Label(mol_style, text="Hide C atom symbols").grid(row=0, column=0)
        ttk.Checkbutton(
            mol_style,
            variable=self._hide_carbon,
            offvalue=0,
            onvalue=1,
            command=self.draw_mode_change,
        ).grid(row=0, column=1)
        self._emptyvalence = tk.StringVar()
        self._emptyvalence.set(self.empty_val_names[self.EmptyValence.NOTHING])
        ttk.Label(mol_style, text="Atoms with not full valence shells:").grid(
            row=0, column=2
        )
        ttk.OptionMenu(
            mol_style,
            self._emptyvalence,
            self.empty_val_names[self.EmptyValence.NOTHING],
            *self.empty_val_names.values(),
            command=self.draw_mode_change,
        ).grid(row=0, column=3)

        for row, elements in enumerate(atom_symbols):
            for index, symbol in enumerate(elements):
                if symbol == "":
                    continue
                ttk.Button(
                    symbol_sel,
                    text=symbol,
                    width=3,
                    command=lambda symbol=symbol: self.set_symbol(symbol),
                ).grid(row=row, column=index, columnspan=1)
        ttk.Entry(custom_add, textvariable=self.custom_atom_symbol).grid(row=0, column=0, sticky="nsew")
        ttk.Button(
            custom_add,
            text="Set",
            command=self.set_custom_symbol
        ).grid(row=0, column=1, sticky="nsew")

        ttk.Button(bond_sel, text="--", width=4, command=lambda: self.set_bond([1, None])).grid(
            row=0, column=0, columnspan=1
        )
        ttk.Button(bond_sel, text="==", width=4, command=lambda: self.set_bond([2, None])).grid(
            row=0, column=1, columnspan=1
        )
        ttk.Button(bond_sel, text="≡≡", width=4, command=lambda: self.set_bond([3, None])).grid(
            row=0, column=2, columnspan=1
        )
        ttk.Button(bond_sel, text="↮", width=4, command=lambda: self.set_bond([None, 0])).grid(
            row=1, column=0, columnspan=1
        )
        ttk.Button(bond_sel, text="⟶", width=4, command=lambda: self.set_bond([None, 1])).grid(
            row=1, column=1, columnspan=1
        )
        ttk.Button(bond_sel, text="⟵", width=4, command=lambda: self.set_bond([None, -1])).grid(
            row=1, column=2, columnspan=1
        )
        ttk.Button(bond_sel, text="(-)", width=4, command=self.minus_charge).grid(
            row=2, column=0, columnspan=1
        )
        ttk.Button(bond_sel, text="(+)", width=4, command=self.plus_charge).grid(
            row=2, column=1, columnspan=1
        )
        ttk.Button(bond_sel, text="UP", width=4, command=self.move_above).grid(
            row=3, column=0, columnspan=1
        )
        ttk.Button(bond_sel, text="DOWN", width=4, command=self.move_below).grid(
            row=3, column=1, columnspan=1
        )

    def draw_mode_change(self, _nothing: Optional[Any] = None) -> None:
        self.event_generate("<<RedrawAll>>", when="tail")

    def status_text_update(self) -> None:
        text: str = "MODE: "
        if self._mode == self.Modes.ADD_ATOM:
            text += "add atom\n"
        elif self._mode == self.Modes.ADD_BOND:
            text += "bond atoms\n"
        elif self._mode == self.Modes.SELECT:
            text += "select\n"
        text += f"ATOM SYMBOL: {self._symbol}\n"
        text += f"BOND ORDER: {self._bond[0]} "
        if self._bond[1] == 0:
            text += "--"
        elif self._bond[1] == -1:
            text += "⟵"
        elif self._bond[1] == 1:
            text += "⟶"
        self.status_text.set(text)

    def set_mode(self, new_mode: int = -1) -> None:
        """Changes the operation mode of the application."""
        if new_mode == -1:
            new_mode = self.Modes.SELECT
        self._mode = new_mode
        self.status_text_update()
        self.event_generate("<<ModeButtonPress>>", when="tail")

    def set_symbol(self, new_symbol: str) -> None:
        """Sets the default atom symbol for the application"""
        self.event_generate("<<AtomButtonPress>>", when="tail")
        self._symbol = new_symbol
        self.status_text_update()

    def set_custom_symbol(self) -> None:
        """Sets the default atom symbol to a custom one for the application"""
        new_symbol: str = self.custom_atom_symbol.get()
        if new_symbol == "":
            return
        self.set_symbol(new_symbol)

    def set_bond(self, bondtype: list[Optional[int]]) -> None:
        """Sets the bond type: [order, dativity]"""
        for index, item in enumerate(bondtype):
            if item is not None:
                self._bond[index] = item
                self.event_generate(f"<<BondButton{index}Press>>", when="tail")
        self.status_text_update()
    
    def plus_charge(self) -> None:
        self.event_generate("<<PlusChargeButtonPress>>", when="tail")

    def minus_charge(self) -> None:
        self.event_generate("<<MinusChargeButtonPress>>", when="tail")

    def move_above(self) -> None:
        self.event_generate("<<MoveAboveButtonPress>>", when="tail")

    def move_below(self) -> None:
        self.event_generate("<<MoveBelowButtonPress>>", when="tail")


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

        self.symbol_font: tkfont.Font = tkfont.nametofont("TkDefaultFont")
        self.symbol_index_font = tkfont.Font(
            font=self.symbol_font.copy(), name="symbol_index_font"
        )
        self.symbol_index_font.configure(size=8)
        self.symbol_font_sel = tkfont.Font(
            font=self.symbol_font.copy(), name="symbol_font_sel"
        )
        self.symbol_font_sel.configure(weight="bold")
        self.symbol_index_font_sel = tkfont.Font(
            font=self.symbol_index_font.copy(), name="symbol_index_font_sel"
        )
        self.symbol_index_font_sel.configure(weight="bold")

        container = ttk.Frame(self)
        container.grid(row=0, column=0, sticky="nsew")
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        self.toolbar = TopToolbar(container, self)
        self.toolbar.grid(row=0, column=0, sticky="nsew")

        self.mol_canvas = tk.Canvas(container, width=800, height=600, bg="white")
        self.mol_canvas.grid(row=1, column=0, sticky="nsew")
        self.mol_canvas.bind("<Button-1>", self.leftclick_canvas)
        self.bind("<Delete>", self.delkey_pressed)
        self.bind("<F5>", self.optimize_2D)
        self.bind("<<RedrawAll>>", lambda x: self.redraw_all())
        self.bind("<<ModeButtonPress>>", self.button_pressed_mode)
        self.bind("<<AtomButtonPress>>", self.button_pressed_atom)
        self.bind("<<BondButton0Press>>", self.button_pressed_order)
        self.bind("<<BondButton1Press>>", self.button_pressed_dativity)
        self.bind("<<PlusChargeButtonPress>>", lambda x: self.button_pressed_charge(1))
        self.bind("<<MinusChargeButtonPress>>", lambda x: self.button_pressed_charge(-1))
        self.bind("<<MoveAboveButtonPress>>", lambda x: self.button_pressed_zaxis(1))
        self.bind("<<MoveBelowButtonPress>>", lambda x: self.button_pressed_zaxis(-1))

    def button_pressed_mode(self, _event: tk.Event) -> None:
        if not self.toolbar.is_select:
            self.select_nothing()

    def button_pressed_atom(self, _event: tk.Event) -> None:
        """Handles the pressing of the atom symbol buttons in the toolbar:
        - in normal mode, if atoms are selected, those are changed to
          the symbol what has been selected on the toolbar."""
        if self.mode == self.Modes.NORMAL:
            changes: bool = False
            change_to: Optional[eng.el.Element]
            change_to = eng.element_by_symbol(self.toolbar.atom_symbol)
            if change_to is None:
                return
            for sel_num in self.selected:
                if sel_num not in self.graphics:
                    continue
                sel_item = self.graphics[sel_num]
                if isinstance(sel_item, eng.Atom):
                    err: eng.BondingError = sel_item.can_change_element(change_to)
                    if err == eng.BondingError.OK:
                        sel_item.change_element(change_to)
                        changes = True
            if changes:
                self.selected = []
                self.redraw_all()

    def button_pressed_order(self, _event: tk.Event) -> None:
        """Handles the pressing of the bond symbol buttons in the toolbar:
        - in normal mode, if bonds are selected, the order of those bonds
          is changed with regard to the toolbar setting."""
        if self.mode == self.Modes.NORMAL:
            changes: bool = False
            new_order = self.toolbar.bond_order
            for sel_num in self.selected:
                if sel_num not in self.graphics:
                    continue
                sel_item = self.graphics[sel_num]
                if isinstance(sel_item, eng.CovBond):
                    err: eng.BondingError = sel_item.can_change_order(new_order)
                    if err == eng.BondingError.OK:
                        sel_item.order = new_order
                        changes = True
            if changes:
                self.selected = []
                self.redraw_all()

    def button_pressed_dativity(self, _event: tk.Event) -> None:
        """Handles the pressing of the bond symbol buttons in the toolbar:
        - in normal mode, if bonds are selected, the dativity of those
          bonds is changed with regard to the toolbar setting."""
        if self.mode == self.Modes.NORMAL:
            changes: bool = False
            new_dativity = self.toolbar.bond_dativity
            for sel_num in self.selected:
                if sel_num not in self.graphics:
                    continue
                sel_item = self.graphics[sel_num]
                if isinstance(sel_item, eng.CovBond):
                    err: eng.BondingError = sel_item.can_change_dativity(new_dativity)
                    if err == eng.BondingError.OK:
                        sel_item.dativity = new_dativity
                        changes = True
            if changes:
                self.selected = []
                self.redraw_all()

    def delkey_pressed(self, _event: tk.Event) -> None:
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

    def button_pressed_charge(self, charge_change: int) -> None:
        """The charge change button has been pressed."""
        if self.mode == self.Modes.NORMAL:
            changes: list[int] = []
            change_to: int = 0
            for sel_num in self.selected:
                if sel_num not in self.graphics:
                    continue
                sel_item = self.graphics[sel_num]
                if isinstance(sel_item, eng.Atom):
                    if self.atoms.index(sel_item) in changes:
                        continue
                    change_to = sel_item.charge + charge_change
                    err: eng.BondingError = sel_item.can_ionize(change_to)
                    if err == eng.BondingError.OK:
                        sel_item.charge = change_to
                        changes.append(self.atoms.index(sel_item))
        self.redraw_all()

    def button_pressed_zaxis(self, change_z: int) -> None:
        """The Z coordinate of the selected atom is moved."""
        if self.mode == self.Modes.NORMAL:
            changes: list[int] = []
            for sel_num in self.selected:
                if sel_num not in self.graphics:
                    continue
                sel_item = self.graphics[sel_num]
                if isinstance(sel_item, eng.Atom):
                    if self.atoms.index(sel_item) in changes:
                        continue
                    sel_item.coord_z += change_z
                    changes.append(self.atoms.index(sel_item))
        self.redraw_all()

    def optimize_2D(self, _event: tk.Event) -> None:
        """Optimize the structural formula based on simple rules."""
        if len(self.selected) == 0:
            return
        for sel_num in self.selected:
            if sel_num not in self.graphics:
                continue
            sel_item = self.graphics[sel_num]
            if isinstance(sel_item, eng.Atom):
                opt.optimize_2D(sel_item, target_len=DEFLEN)
        self.redraw_all()

    def leftclick_canvas(self, event: tk.Event) -> None:
        """Handes left click on the canvas:
        - if the operation is adding an atom, adds a new atom at the
          position of mouse click"""
        if self.event_listened:
            self.event_listened = False
            return
        if self.toolbar.is_select:
            self.select_nothing()
            return
        if self.toolbar.is_add:
            closestitems: tuple[int, ...] = self.mol_canvas.find_overlapping(
                event.x - MINDIST / 2,
                event.y - MINDIST / 2,
                event.x + MINDIST / 2,
                event.y + MINDIST / 2,
            )
            if len(closestitems) == 0:
                item = None
            else:
                item = closestitems[0]
            if item in self.graphics and isinstance(self.graphics[item], eng.Atom):
                item_obj = self.graphics[item]
                assert isinstance(item_obj, eng.Atom)
                if (
                    -MINDIST < item_obj.coord_x - event.x < MINDIST
                    and -MINDIST < item_obj.coord_y - event.y < MINDIST
                ):
                    return
            coords = (event.x, event.y)
            self.add_atom(self.toolbar.atom_symbol, coords)

    def leftdown_atom(self, atom_id: int, _event: tk.Event) -> None:
        """Handles the event of mouse button press on atom, based on the
        operation:
        - if selection, current atom is selected,
        - if add atom, enters mode to add linked atom,
        - if connect atoms, enters mode to connect existing atoms."""
        sel_atom = self.atoms[atom_id]
        if self.toolbar.is_select:
            self.select_nothing()
            self.selected = list(self.mol_canvas.find_withtag(f"atom-{atom_id}"))
            self.selected += list(self.mol_canvas.find_withtag(f"atom_text-{atom_id}"))
            self.mol_canvas.itemconfigure(f"atom-{atom_id}", fill="green")
            self.mol_canvas.itemconfigure(f"atom_text-{atom_id}", fill="green")
            self.mol_canvas.itemconfigure(
                f"atom_rad-{atom_id}", outline="green", width=3
            )
            self.mol_canvas.itemconfigure(f"atom_lpair-{atom_id}", width=2)
            for sel_element in self.selected:
                if "atom_text" in self.mol_canvas.gettags(sel_element):
                    self.change_font_weight(sel_element, "bold")
            self.event_listened = True
            return
        if self.toolbar.is_add:
            self.event_listened = True
            test_atom = eng.add_atom_by_symbol(self.toolbar.atom_symbol, (0, 0))
            if test_atom is None:
                test_atom = eng.UnrestrictedAtom(self.toolbar.atom_symbol, (0, 0))
            assert test_atom is not None, "Test atom could not be created."
            if (
                sel_atom.can_bond(
                    test_atom, self.toolbar.bond_order, self.toolbar.bond_dativity
                )
                != eng.BondingError.OK
            ):
                return
            self.possible_atoms(sel_atom)
            self.mode = self.Modes.ADD_LINKED_ATOM
            return
        if self.toolbar.is_connect:
            self.event_listened = True
            self.possible_links(sel_atom)
            self.mode = self.Modes.LINK_ATOMS

    def mouseup_atom(self, atom_id: int, event: tk.Event) -> None:
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
            atom_connect = self.atoms[atom_id]
            new_atom = self.add_atom(self.toolbar.atom_symbol, (new_x, new_y))
            new_bond = atom_connect.bond(
                new_atom,
                order=self.toolbar.bond_order,
                dative=self.toolbar.bond_dativity,
            )
            self.bonds.append(new_bond)
            self.redraw_all()
            self.set_normal_mode()
        if self.mode == self.Modes.LINK_ATOMS:
            closestitems = self.mol_canvas.find_closest(event.x, event.y)
            self.mol_canvas.itemconfigure("ui_help", fill="#aaaaaa")
            self.mol_canvas.itemconfigure("atom", fill="black")
            self.mol_canvas.itemconfigure("empty", fill="")
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
            atom_from = self.atoms[atom_id]
            atom_to = self.graphics[selected_atom]
            assert isinstance(atom_to, eng.Atom), "Trying to bobd to a non-atom"
            if (
                atom_from.can_bond(
                    atom_to, self.toolbar.bond_order, self.toolbar.bond_dativity
                )
                != eng.BondingError.OK
            ):
                self.set_normal_mode()
                return
            new_bond = atom_from.bond(
                atom_to,
                order=self.toolbar.bond_order,
                dative=self.toolbar.bond_dativity,
            )
            assert new_bond is not None, "Creating a bond failed"
            self.bonds.append(new_bond)
            self.redraw_all()
            self.set_normal_mode()

    def drag_atom(self, atom_id: int, event: tk.Event) -> None:
        """Handling the dragging of an atom that was clicked on:
        - adding linked atom / linking atoms: highlights the
          placeholder where new atom / bond will be created."""
        if self.mode in (self.Modes.ADD_LINKED_ATOM, self.Modes.LINK_ATOMS):
            closestitems = self.mol_canvas.find_closest(event.x, event.y)
            self.mol_canvas.itemconfigure("ui_help", fill="#aaaaaa")
            self.mol_canvas.itemconfigure("atom", fill="black")
            self.mol_canvas.itemconfigure("empty", fill="")
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
                atom_from = self.atoms[atom_id]
                atom_to = self.graphics[selected_atom]
                assert isinstance(atom_to, eng.Atom)
                if (
                    atom_from.can_bond(
                        atom_to, self.toolbar.bond_order, self.toolbar.bond_dativity
                    )
                    != eng.BondingError.OK
                ):
                    return
                self.mol_canvas.itemconfigure(f"ui_help-{selected_atom}", fill="blue")
                self.mol_canvas.itemconfigure(selected_atom, fill="blue")

    def click_bond(self, bond_id: int, _event: tk.Event) -> None:
        """Handling the event of a bond being clicked on:
        - if operation is selection, selects the bond, and all lines
          that constitute the bond."""
        if self.toolbar.is_select:
            self.select_nothing()
            self.selected = list(self.mol_canvas.find_withtag(f"bond-{bond_id}"))
            self.mol_canvas.itemconfigure(f"bond-{bond_id}", fill="green", width=2)
            self.event_listened = True
            return

    def change_font_weight(
        self, item_num: int, new_weight: Literal["bold", "normal"]
    ) -> None:
        assert "atom_text" in self.mol_canvas.gettags(
            item_num
        ), "Only text items have fonts."
        if "index" in self.mol_canvas.gettags(item_num):
            if new_weight == "bold":
                self.mol_canvas.itemconfigure(item_num, font=self.symbol_index_font_sel)
                return
            self.mol_canvas.itemconfigure(item_num, font=self.symbol_index_font)
            return
        if new_weight == "bold":
            self.mol_canvas.itemconfigure(item_num, font=self.symbol_font_sel)
            return
        self.mol_canvas.itemconfigure(item_num, font=self.symbol_font)
        return

    def select_nothing(self) -> None:
        """Clears the selection and changes the previously selected objects
        back to their default look."""
        for sel_num in self.selected:
            tags: Sequence[str] = self.mol_canvas.gettags(sel_num)
            if "empty" in tags:
                self.mol_canvas.itemconfigure(sel_num, fill="")
            else:
                self.mol_canvas.itemconfigure(sel_num, fill="black")
            if "atom_text" in tags:
                self.change_font_weight(sel_num, "normal")
            if "bond" in tags:
                sel_item = self.graphics[sel_num]
                assert isinstance(sel_item, eng.CovBond)
                atom1: eng.Atom = sel_item.atoms[0]
                atom2: eng.Atom = sel_item.atoms[1]
                if atom1.coord_z > 0 and atom2.coord_z == atom1.coord_z:
                    self.mol_canvas.itemconfigure(sel_num, width=3)
                else:
                    self.mol_canvas.itemconfigure(sel_num, width=1)
            if "atom_rad" in tags:
                self.mol_canvas.itemconfigure(sel_num, outline="black", width=1)
            if "atom_lpair" in tags:
                self.mol_canvas.itemconfigure(sel_num, width=1)
        self.selected = []

    def close_selection(
        self, closestitems: tuple[int, ...]
    ) -> tuple[Optional[int], Optional[Sequence[str]]]:
        """Finds the closest item on the canvas to the
        listed item on the canvas."""
        if len(closestitems) == 0:
            self.set_normal_mode()
            return None, None
        item: int = closestitems[0]
        tags: Sequence[str] = self.mol_canvas.gettags(item)
        if (len(tags) < 2 or tags[0] != "ui_help") and (
            len(tags) < 1 or tags[0] != "atom"
        ):
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
            if (
                atom_from.can_bond(
                    atom_to, self.toolbar.bond_order, self.toolbar.bond_dativity
                )
                != eng.BondingError.OK
            ):
                continue
            self.draw_bond_atoms(
                atom_from,
                atom_to,
                -1,
                self.toolbar.bond_order,
                self.toolbar.bond_dativity,
                color="#aaaaaa",
                tags=("ui_help", f"ui_help-{atom_link_num}"),
            )

    def possible_atoms(self, atom: eng.Atom) -> None:
        """Draws possible atoms that could be created and joined
        to the current atom in linked atom drawing mode."""
        num = 1
        for (deltax, deltay) in zip(XSHIFT, YSHIFT):
            (x1, y1, x2, y2) = (
                atom.coord_x + deltax - MINDIST / 2,
                atom.coord_y + deltay - MINDIST / 2,
                atom.coord_x + deltax + MINDIST / 2,
                atom.coord_y + deltay + MINDIST / 2,
            )
            over = self.mol_canvas.find_overlapping(x1, y1, x2, y2)
            over_atom = False
            for obj in over:
                tags: Sequence[str] = self.mol_canvas.gettags(obj)
                if len(tags) > 0 and (
                    "atom" in tags or "bond" in tags or "atom_text" in tags
                ):
                    over_atom = True
                    break
            if over_atom:
                num += 1
                continue
            self.mol_canvas.create_text(
                atom.coord_x + deltax,
                atom.coord_y + deltay,
                text=self.toolbar.atom_symbol,
                font=self.symbol_font,
                justify="center",
                fill="#aaaaaa",
                tags=("ui_help", f"ui_help-{num}", f"atom_here-{num}"),
            )
            obscure: tuple[str, str]
            if self.hide_atom(atom):
                obscure = ("", f"atom_here-{num}")
            else:
                obscure = (f"atom_text-{self.atoms.index(atom)}", f"atom_here-{num}")
            self.draw_bond_lines(
                (
                    atom.coord_x,
                    atom.coord_y,
                    atom.coord_x + deltax,
                    atom.coord_y + deltay,
                ),
                30,
                self.toolbar.bond_order,
                self.toolbar.bond_dativity,
                obscure_area=obscure,
                color="#aaaaaa",
                tags=("ui_help", f"ui_help-{num}"),
            )
            num += 1

    def add_atom(self, atom_symbol: str, coords: tuple[float, float]) -> eng.Atom:
        """Adds an atom to the canvas with given symbol at
        given coordinates."""
        new_atom = eng.add_atom_by_symbol(atom_symbol, coords)
        if new_atom is None:
            new_atom = eng.UnrestrictedAtom(atom_symbol, coords)
        assert new_atom is not None, "Unsuccessful at creating a new atom."
        self.atoms.append(new_atom)
        self.redraw_all()
        return new_atom

    def draw_bond_atoms(
        self,
        atom1: eng.Atom,
        atom2: eng.Atom,
        bondlen: float = -1,
        order: float = -1,
        dativity: float = 0,
        color: str = "black",
        tags: Optional[tuple[str, ...]] = None,
    ) -> list[int]:
        """Draws the bond between two atom instances."""
        x1, y1, x2, y2 = atom1.coord_x, atom1.coord_y, atom2.coord_x, atom2.coord_y
        line_thickness: int = 1
        line_dash: Optional[tuple[int, int]] = None
        if bondlen == -1:
            bondlen = sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        angle: float = atan2(
            atom2.coord_y - atom1.coord_y, atom2.coord_x - atom1.coord_x
        )
        reverse_second_bond: bool = False
        if order == 2:
            left1, right1 = eng.angle_side_calc(
                eng.angles_rel_angle(atom1.bond_angles, angle), 0.01
            )
            left2, right2 = eng.angle_side_calc(
                eng.angles_rel_angle(atom2.bond_angles, angle), 0.01
            )
            reverse_second_bond = left1 + left2 > right1 + right2
        atomtag1: str = ""
        atomtag2: str = ""
        if not self.hide_atom(atom1):
            atomtag1 = f"atom_text-{self.atoms.index(atom1)}"
        if not self.hide_atom(atom2):
            atomtag2 = f"atom_text-{self.atoms.index(atom2)}"
        if atom1.coord_z > 0 and atom1.coord_z == atom2.coord_z:
            line_thickness = 3
        if atom1.coord_z < 0 and atom1.coord_z == atom2.coord_z:
            line_dash = (2, 1)
        triangle: int = atom2.coord_z - atom1.coord_z
        if order > 0 and dativity != 0:
            triangle = 0
        bond_objects = self.draw_bond_lines(
            (x1, y1, x2, y2),
            bondlen,
            order,
            dativity,
            obscure_area=(atomtag1, atomtag2),
            reverse_second_bond=reverse_second_bond,
            color=color,
            triangle=triangle,
            width=line_thickness,
            dashed=line_dash,
            tags=tags,
        )
        return bond_objects

    def draw_bond(
        self,
        bond: eng.CovBond,
        color: str = "black",
        tags: Optional[tuple[str, ...]] = None,
    ) -> list[int]:
        """Calls the bond line drawing function with parameters extracted from the
        instance of CobBond it receives."""
        atom1: eng.Atom = bond.atoms[0]
        atom2: eng.Atom = bond.atoms[1]
        bond_objects = self.draw_bond_atoms(
            atom1, atom2, bond.length, bond.order, bond.dativity, color=color, tags=tags
        )
        return bond_objects

    def find_edge(
        self, coords: tuple[float, float], tag: str, angle: float
    ) -> tuple[float, float]:
        """Finds the coordinates at the edge of canvas objects with tag, starting from
        coords, moving in the direction of angle."""
        while angle > pi:
            angle -= 2 * pi
        while angle < -pi:
            angle += 2 * pi
        x_look, y_look = coords
        while True:
            item_ids = self.mol_canvas.find_overlapping(x_look, y_look, x_look, y_look)
            found_id: int = -1
            for item_id in item_ids:
                if tag in self.mol_canvas.gettags(item_id):
                    found_id = item_id
                    break
            if found_id == -1:
                return x_look, y_look
            (xb1, yb1, xb2, yb2) = self.mol_canvas.bbox(found_id)
            if angle == 0:
                x_look = xb2 + 1
                continue
            if angle == pi / 2:
                y_look = yb2 + 1
                continue
            if angle == pi or angle == -pi:
                x_look = xb1 - 1
                continue
            if angle == -pi / 2:
                y_look = yb1 - 1
                continue
            dx: float = 0
            dy: float = 0
            if -pi < angle < 0:
                dx = (yb1 - 1 - y_look) / tan(angle)
            if -pi / 2 < angle < pi / 2:
                dy = (xb2 + 1 - x_look) * tan(angle)
            if 0 < angle < pi:
                dx = (yb2 + 1 - y_look) / tan(angle)
            if angle > pi / 2 or angle < -pi / 2:
                dy = (xb1 - 1 - x_look) * tan(angle)
            x_look = max(min(x_look + dx, xb2 + 1), xb1 - 1)
            y_look = max(min(y_look + dy, yb2 + 1), yb1 - 1)

    def draw_bond_lines(
        self,
        coords: tuple[float, float, float, float],
        bondlen: float = -1,
        order: float = -1,
        dativity: float = 0,
        obscure_area: Optional[tuple[str, str]] = None,
        reverse_second_bond: bool = False,
        triangle: int = 0,
        color: str = "black",
        width: int = 1,
        dashed: Optional[tuple[int, int]] = None,
        tags: Optional[tuple[str, ...]] = None,
    ) -> list[int]:
        """Draws a bond between two atoms and returns the list of line
        IDs on the canvas that represent the chemical bond."""
        x1, y1, x2, y2 = coords
        angle: float = atan2(y2 - y1, x2 - x1)
        if tags is None:
            tags = ()
        if bondlen == -1:
            bondlen = sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        if order == -1:
            order = self.toolbar.bond_order
        if obscure_area is None:
            obscure_area = ("", "")
        assert obscure_area is not None
        shift = 2
        order_i: int = int(order)
        if (obscure_area[0] == "") + (obscure_area[1] == "") > 0:
            bondshifts = list(
                range(2 - order_i % 2 - order_i, order_i + 2 - order_i % 2, 2)
            )
        else:
            bondshifts = list(range(1 - order_i, order_i + 1, 2))
        if obscure_area[0] != "":
            xd1, yd1 = self.find_edge((x1, y1), obscure_area[0], angle)
        else:
            xd1, yd1 = x1, y1
        if obscure_area[1] != "":
            xd2, yd2 = self.find_edge((x2, y2), obscure_area[1], pi + angle)
        else:
            xd2, yd2 = x2, y2
        bond_objects: list[int] = []
        for bondpart in bondshifts:
            if reverse_second_bond:
                bondpart *= -1
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
            if triangle == 0:
                bond_line_id = self.mol_canvas.create_line(
                    (xd1 + xs, yd1 + ys, xd2 + xs, yd2 + ys),
                    arrow=arrowtype,
                    fill=color,
                    width=width,
                    dash=dashed,
                    tags=tags,
                )
            elif triangle > 0:
                angle: float = atan2(yd2 - yd1, xd2 - xd1)
                bond_line_id = self.mol_canvas.create_polygon(
                    (
                        xd1 + xs,
                        yd1 + ys,
                        xd2 + xs + 2 * cos(angle - pi / 2),
                        yd2 + ys + 2 * sin(angle - pi / 2),
                        xd2 + xs - 2 * cos(angle - pi / 2),
                        yd2 + ys - 2 * sin(angle - pi / 2),
                    ),
                    fill=color,
                    tags=tags
                )
            else:
                angle: float = atan2(yd2 - yd1, xd2 - xd1)
                bond_line_id = self.mol_canvas.create_polygon(
                    (
                        xd1 + xs + 2 * cos(angle - pi / 2),
                        yd1 + ys + 2 * sin(angle - pi / 2),
                        xd1 + xs - 2 * cos(angle - pi / 2),
                        yd1 + ys - 2 * sin(angle - pi / 2),
                        xd2 + xs,
                        yd2 + ys,
                    ),
                    fill=color,
                    tags=tags
                )
            bond_objects.append(bond_line_id)
        return bond_objects

    def hide_atom(self, atom: eng.Atom) -> bool:
        """Determines if an atom should be hidden (and symbol not shown on the formula)"""
        if atom.symbol == "C" and self.toolbar.hide_carbon and len(atom.bonds) > 0:
            if len(atom.bonds) > 2:
                return True
            if len(atom.bonds) == 2:
                angles = atom.bond_angles
                angles.sort()
                return not (pi - 0.1 < angles[1] - angles[0] < pi + 0.1)
            linked_atoms_sym = [bond.other_atoms(atom)[0].symbol for bond in atom.bonds]
            link_at_linkum = [
                len(bond.other_atoms(atom)[0].bonds) for bond in atom.bonds
            ]
            secondlinks = sum(link_at_linkum)
            if "C" in linked_atoms_sym and secondlinks > 1:
                return True
        return False

    def draw_atom_electrons(self, atom: eng.Atom, atom_id: int) -> list[int]:
        """Draws lone pairs and radicals around atoms."""
        atom_s: list[int] = []
        lone_pair: int = atom.lone_pairs
        radical: int = atom.radicals
        el_angles: list[float] = atom.electron_angles
        assert lone_pair + radical == len(
            el_angles
        ), "Error in electron number calculation"
        for angle in el_angles:
            xcent: float = atom.coord_x + MINDIST / 2 * cos(angle)
            ycent: float = atom.coord_y + MINDIST / 2 * sin(angle)
            if lone_pair > 0:
                line_coord: tuple[float, float, float, float] = (
                    xcent + MINDIST / 4 * cos(angle - pi / 2),
                    ycent + MINDIST / 4 * sin(angle - pi / 2),
                    xcent - MINDIST / 4 * cos(angle - pi / 2),
                    ycent - MINDIST / 4 * sin(angle - pi / 2),
                )
                atom_s.append(
                    self.mol_canvas.create_line(
                        line_coord,
                        tags=(
                            f"atom-{atom_id}",
                            "atom_lpair",
                            f"atom_lpair-{atom_id}",
                        ),
                    )
                )
                lone_pair -= 1
            else:
                atom_s.append(
                    self.mol_canvas.create_oval(
                        (xcent, ycent, xcent + 1, ycent + 1),
                        tags=(
                            f"atom-{atom_id}",
                            "atom_rad",
                            f"atom_rad-{atom_id}",
                        ),
                    )
                )
                radical -= 1
        return atom_s

    def draw_atom_with_H(self, atom: eng.Atom, atom_id: int) -> list[int]:
        """Draws the atoms with H next to it if needed"""
        atom_s: list[int] = []
        left_n, right_n = eng.angle_side_calc(
            eng.angles_rel_angle(atom.bond_angles, -pi / 2), pi / 6 + 0.1
        )
        right: bool = right_n > 0
        left: bool = left_n > 0
        if left == right:
            atom_s.append(
                self.mol_canvas.create_text(
                    atom.coord_x,
                    atom.coord_y,
                    text=f"{atom.symbol}H",
                    justify="center",
                    font=self.symbol_font,
                    tags=(
                        "atom",
                        f"atom-{atom_id}",
                        "atom_text",
                        f"atom_text-{atom_id}",
                        f"atom_core-{atom_id}",
                    ),
                )
            )
            if atom.radicals == 1:
                return atom_s
            (_, bby1, bbx2, bby2) = self.mol_canvas.bbox(f"atom-{atom_id}")
            atom_s.append(
                self.mol_canvas.create_text(
                    bbx2 + 2,
                    (bby1 + bby2 * 3) / 4,
                    text=f"{atom.radicals}",
                    justify="left",
                    font=self.symbol_index_font,
                    tags=(
                        f"atom-{atom_id}",
                        "atom_text",
                        f"atom_text-{atom_id}",
                        "index",
                    ),
                )
            )
            return atom_s
        atom_s.append(
            self.mol_canvas.create_text(
                atom.coord_x,
                atom.coord_y,
                text=atom.symbol,
                justify="center",
                font=self.symbol_font,
                tags=(
                    "atom",
                    f"atom-{atom_id}",
                    "atom_text",
                    f"atom_text-{atom_id}",
                    f"atom_core-{atom_id}",
                ),
            )
        )
        (bbx1, bby1, bbx2, bby2) = self.mol_canvas.bbox(f"atom-{atom_id}")
        if left:
            atom_s.append(
                self.mol_canvas.create_text(
                    bbx2 + 3,
                    atom.coord_y,
                    text="H",
                    justify="left",
                    font=self.symbol_font,
                    tags=(
                        f"atom-{atom_id}",
                        "atom_text",
                        f"atom_text-{atom_id}",
                        f"atom_core-{atom_id}",
                    ),
                )
            )
            if atom.radicals == 1:
                return atom_s
            (bbx1, bby1, bbx2, bby2) = self.mol_canvas.bbox(f"atom-{atom_id}")
            atom_s.append(
                self.mol_canvas.create_text(
                    bbx2 + 2,
                    (bby1 + bby2 * 3) / 4,
                    text=f"{atom.radicals}",
                    justify="left",
                    font=self.symbol_index_font,
                    tags=(
                        f"atom-{atom_id}",
                        "atom_text",
                        f"atom_text-{atom_id}",
                        "index",
                    ),
                )
            )
            return atom_s
        if not right:
            assert False, "Unreachable"
        if atom.radicals > 1:
            atom_s.append(
                self.mol_canvas.create_text(
                    bbx1 - 2,
                    (bby1 + bby2 * 3) / 4,
                    text=f"{atom.radicals}",
                    justify="right",
                    font=self.symbol_index_font,
                    tags=(
                        f"atom-{atom_id}",
                        "atom_text",
                        f"atom_text-{atom_id}",
                        "index",
                    ),
                )
            )
            (bbx1, bby1, bbx2, bby2) = self.mol_canvas.bbox(f"atom-{atom_id}")
        atom_s.append(
            self.mol_canvas.create_text(
                bbx1 - 3,
                atom.coord_y,
                text="H",
                justify="right",
                font=self.symbol_font,
                tags=(
                    f"atom-{atom_id}",
                    "atom_text",
                    f"atom_text-{atom_id}",
                    f"atom_core-{atom_id}",
                ),
            )
        )
        return atom_s

    def draw_atom_charge(
        self,
        coord_x: float,
        coord_y: float,
        charge: int,
        atom_id: int
    ) -> int:
        """Draws the charge symbol next to an atom."""
        charge_text: str
        if charge == 0:
            charge_text = "0"
        elif charge == -1:
            charge_text = "−"
        elif charge == 1:
            charge_text = "+"
        elif charge > 1:
            charge_text = f"{charge}+"
        elif charge < 1:
            charge_text = f"{-charge}−"
        return self.mol_canvas.create_text(
            coord_x,
            coord_y,
            text=charge_text,
            font=self.symbol_index_font,
            justify="left",
            tags=(
                f"atom-{atom_id}",
                "atom_text",
                f"atom_text-{atom_id}",
                "index",
            ),
        )

    def draw_atom(self, atom: eng.Atom, atom_id: int) -> list[int]:
        """Draws the representation of the atom on the screen."""
        atom_s: list[int] = []
        if self.hide_atom(atom):
            atom_s.append(
                self.mol_canvas.create_rectangle(
                    atom.coord_x - MINDIST / 3,
                    atom.coord_y - MINDIST / 3,
                    atom.coord_x + MINDIST / 3,
                    atom.coord_y + MINDIST / 3,
                    width=0,
                    tags=("atom", f"atom-{atom_id}", "empty"),
                )
            )
            if atom.charge != 0:
                atom_s.append(
                    self.draw_atom_charge(
                        atom.coord_x + 10,
                        atom.coord_y - 10,
                        atom.charge,
                        atom_id
                    )
                )
            return atom_s
        if (
            self.toolbar.empty_val_style == self.toolbar.EmptyValence.HYDROGENS
            and atom.radicals > 0
        ):
            atom_s += self.draw_atom_with_H(atom, atom_id)
        else:
            atom_s.append(
                self.mol_canvas.create_text(
                    atom.coord_x,
                    atom.coord_y,
                    text=atom.symbol,
                    justify="center",
                    font=self.symbol_font,
                    tags=(
                        "atom",
                        f"atom-{atom_id}",
                        "atom_text",
                        f"atom_text-{atom_id}",
                        f"atom_core-{atom_id}",
                    ),
                )
            )
            if self.toolbar.empty_val_style == self.toolbar.EmptyValence.ELECTRONS:
                atom_s += self.draw_atom_electrons(atom, atom_id)
        if atom.charge != 0:
            (_, bby1, bbx2, bby2) = self.mol_canvas.bbox(f"atom_core-{atom_id}")
            atom_s.append(
                self.draw_atom_charge(
                    bbx2 + 2,
                    (5 * bby1 - bby2) / 4,
                    atom.charge,
                    atom_id
                )
            )
        return atom_s

    def redraw_all(self) -> None:
        """Redraws all atoms and bonds on a freshly cleared canvas."""
        self.mol_canvas.delete("all")
        selected_objects: list[Union[eng.Atom, eng.CovBond]]
        selected_objects = [self.graphics[item] for item in self.selected]
        self.selected = []
        self.graphics = {}
        for atom_id, atom in enumerate(self.atoms):
            atom_s: list[int] = self.draw_atom(atom, atom_id)
            for atom_graph_ind in atom_s:
                self.graphics[atom_graph_ind] = atom
            if atom in selected_objects:
                self.selected += atom_s
            self.mol_canvas.tag_bind(
                f"atom-{atom_id}",
                "<ButtonPress-1>",
                lambda event, atom_id=atom_id: self.leftdown_atom(atom_id, event),
            )
            self.mol_canvas.tag_bind(
                f"atom-{atom_id}",
                "<ButtonRelease-1>",
                lambda event, atom_id=atom_id: self.mouseup_atom(atom_id, event),
            )
            self.mol_canvas.tag_bind(
                f"atom-{atom_id}",
                "<B1-Motion>",
                lambda event, atom_id=atom_id: self.drag_atom(atom_id, event),
            )
        for bond_id, bond in enumerate(self.bonds):
            bondlen = bond.length
            if bondlen == 0:
                continue
            bond_drawings: list[int] = self.draw_bond(
                bond,
                tags=("bond", f"bond-{bond_id}"),
            )
            for bondid in bond_drawings:
                self.graphics[bondid] = bond
            if bond in selected_objects:
                self.selected += bond_drawings
            self.mol_canvas.tag_bind(
                f"bond-{bond_id}",
                "<Button-1>",
                lambda event, bond_id=bond_id: self.click_bond(bond_id, event),
            )
        for item_num in self.selected:
            tags: Sequence[str] = self.mol_canvas.gettags(item_num)
            self.mol_canvas.itemconfigure(item_num, fill="green")
            gui_element = self.graphics[item_num]
            if isinstance(gui_element, eng.CovBond):
                self.mol_canvas.itemconfigure(item_num, width=2)
            if isinstance(gui_element, eng.Atom):
                if "atom_text" in tags:
                    self.change_font_weight(item_num, "bold")


def main() -> None:
    app = AppContainer()
    app.mainloop()


if __name__ == "__main__":
    main()
