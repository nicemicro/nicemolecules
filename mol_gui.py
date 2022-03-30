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
from tkinter import font as tkfont
from math import sqrt
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
        index: int = list(
                self.empty_val_names.values()).index(
                    self._emptyvalence.get())
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
        self.controller: tk.Tk = controller
        self._mode: int = self.Modes.ADD_ATOM
        self._symbol: str = "C"
        self._bond: list[int] = [1, 0]
        self.status_text = tk.StringVar()
        atom_symbols: list[list[str]] = [
                ["H", "C", "N", "O", "F"],
                ["", "", "P", "S"]]
        self.empty_val_names = {
                self.EmptyValence.NOTHING: "Nothing",
                self.EmptyValence.HYDROGENS: "Hydrogens",
                self.EmptyValence.ELECTRONS: "Electrons"
        }

        mode_row = ttk.Frame(self, padding="0 0 0 5")
        mode_row.grid(row=0, column=0, sticky="nsew")
        mol_style = ttk.Frame(self, padding="10 0 0 5")
        mol_style.grid(row=0, column=1, columnspan=2, sticky="nse")
        symbol_sel = ttk.Frame(self, padding="0 0 0 5")
        symbol_sel.grid(row=1, column=0, columnspan=2, sticky="nsw")
        bond_sel = ttk.Frame(self, padding="10 0 0 5")
        bond_sel.grid(row=1, column=2, sticky="nse")
        ttk.Label(self, textvariable=self.status_text).grid(
            row=0, column=3, rowspan=2, sticky="nesw"
        )
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
                    command=lambda symbol=symbol: self.set_symbol(symbol),
                ).grid(row=row, column=index, columnspan=1)

        ttk.Button(bond_sel, text="--", command=lambda: self.set_bond([1, None])).grid(
            row=0, column=0, columnspan=1
        )
        ttk.Button(bond_sel, text="==", command=lambda: self.set_bond([2, None])).grid(
            row=0, column=1, columnspan=1
        )
        ttk.Button(bond_sel, text="≡≡", command=lambda: self.set_bond([3, None])).grid(
            row=0, column=2, columnspan=1
        )
        ttk.Button(bond_sel, text="↮", command=lambda: self.set_bond([None, 0])).grid(
            row=1, column=0, columnspan=1
        )
        ttk.Button(bond_sel, text="⟶", command=lambda: self.set_bond([None, 1])).grid(
            row=1, column=1, columnspan=1
        )
        ttk.Button(bond_sel, text="⟵", command=lambda: self.set_bond([None, -1])).grid(
            row=1, column=2, columnspan=1
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

    def set_bond(self, bondtype: list[Optional[int]]) -> None:
        """Sets the bond type: [order, dativity]"""
        for index, item in enumerate(bondtype):
            if item is not None:
                self._bond[index] = item
                self.event_generate(f"<<BondButton{index}Press>>", when="tail")
        self.status_text_update()


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
            font=self.symbol_font.copy(),
            name="symbol_index_font"
        )
        self.symbol_index_font.configure(size=8)
        self.symbol_font_sel = tkfont.Font(
            font=self.symbol_font.copy(),
            name="symbol_font_sel"
        )
        self.symbol_font_sel.configure(weight="bold")
        self.symbol_index_font_sel = tkfont.Font(
            font=self.symbol_index_font.copy(),
            name="symbol_index_font_sel"
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
        self.bind("<<RedrawAll>>", lambda x: self.redraw_all())
        self.bind("<<ModeButtonPress>>", self.button_pressed_mode)
        self.bind("<<AtomButtonPress>>", self.button_pressed_atom)
        self.bind("<<BondButton0Press>>", self.button_pressed_order)
        self.bind("<<BondButton1Press>>", self.button_pressed_dativity)
    
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

    def button_pressed_order(self, event: tk.Event) -> None:
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

    def button_pressed_dativity(self, event: tk.Event) -> None:
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
        if self.toolbar.is_select:
            self.select_nothing()
            return
        if self.toolbar.is_add:
            closestitems: tuple[int, ...] = self.mol_canvas.find_closest(
                event.x, event.y
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
            coords: list[int] = [event.x, event.y]
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
            for sel_element in self.selected:
                if "atom_text" in self.mol_canvas.gettags(sel_element):
                    self.change_font_weight(sel_element, "bold")
            self.event_listened = True
            return
        if self.toolbar.is_add:
            self.event_listened = True
            test_atom = eng.add_atom_by_symbol(self.toolbar.atom_symbol, [0, 0])
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
            new_atom = self.add_atom(self.toolbar.atom_symbol, [new_x, new_y])
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

    def change_font_weight(self, item_num: int, new_weight: Literal["bold", "normal"]) -> None:
        assert "atom_text" in self.mol_canvas.gettags(item_num), \
            "Only text items have fonts."
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
        for sel_item in self.selected:
            tags: Sequence[str] = self.mol_canvas.gettags(sel_item)
            if "empty" in tags:
                self.mol_canvas.itemconfigure(sel_item, fill="")
            else:
                self.mol_canvas.itemconfigure(sel_item, fill="black")
            if "atom_text" in tags:
                self.change_font_weight(sel_item, "normal")
            if "bond" in tags:
                self.mol_canvas.itemconfigure(sel_item, width=1)
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
            cut_bond: tuple[bool, bool]
            cut_bond = (not self.hide_atom(atom_from), not self.hide_atom(atom_to))
            self.draw_bond(
                (
                    atom_from.coord_x,
                    atom_from.coord_y,
                    atom_to.coord_x,
                    atom_to.coord_y,
                ),
                -1,
                self.toolbar.bond_order,
                self.toolbar.bond_dativity,
                cut=cut_bond,
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
                if len(tags) > 0 and ("atom" in tags or "bond" in tags):
                    over_atom = True
                    break
            if over_atom:
                num += 1
                continue
            cut_bond: tuple[bool, bool]
            if self.hide_atom(atom):
                cut_bond = (False, True)
            else:
                cut_bond = (True, True)
            self.draw_bond(
                (
                    atom.coord_x,
                    atom.coord_y,
                    atom.coord_x + deltax,
                    atom.coord_y + deltay,
                ),
                30,
                self.toolbar.bond_order,
                self.toolbar.bond_dativity,
                cut=cut_bond,
                color="#aaaaaa",
                tags=("ui_help", f"ui_help-{num}"),
            )
            self.mol_canvas.create_text(
                atom.coord_x + deltax,
                atom.coord_y + deltay,
                text=self.toolbar.atom_symbol,
                justify="center",
                fill="#aaaaaa",
                tags=("ui_help", f"ui_help-{num}", f"atom_here-{num}"),
            )
            num += 1

    def add_atom(self, atom_symbol: str, coords: Sequence[float]) -> eng.Atom:
        """Adds an atom to the canvas with given symbol at
        given coordinates."""
        new_atom = eng.add_atom_by_symbol(atom_symbol, coords)
        assert new_atom is not None, "Unknown symbol probably"
        self.atoms.append(new_atom)
        self.redraw_all()
        return new_atom

    def draw_bond(
        self,
        coords: tuple[float, float, float, float],
        bondlen: float = -1,
        order: float = -1,
        dativity: float = 0,
        cut: Optional[tuple[bool, bool]] = None,
        color: str = "black",
        tags: Optional[tuple[str, ...]] = None,
    ) -> list[int]:
        """Draws a bond between two atoms and returns the list of line
        IDs on the canvas that represent the chemical bond."""
        x1, y1, x2, y2 = coords
        if tags is None:
            tags = ()
        if bondlen == -1:
            bondlen = sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        if order == -1:
            order = self.toolbar.bond_order
        if cut is None:
            cut = (True, True)
        assert cut is not None
        obscure = 7
        shift = 2
        bondshifts = list(range(1 - int(order), int(order) + 1, 2))
        bond_objects: list[int] = []
        for bondpart in bondshifts:
            if cut[0]:
                xd1 = x1 + (obscure / bondlen) * (x2 - x1)
                yd1 = y1 + (obscure / bondlen) * (y2 - y1)
            else:
                xd1, yd1 = x1, y1
            if cut[1]:
                xd2 = x2 - (obscure / bondlen) * (x2 - x1)
                yd2 = y2 - (obscure / bondlen) * (y2 - y1)
            else:
                xd2, yd2 = x2, y2
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
            bond_line_id = self.mol_canvas.create_line(
                xd1 + xs,
                yd1 + ys,
                xd2 + xs,
                yd2 + ys,
                arrow=arrowtype,
                fill=color,
                tags=tags,
            )
            bond_objects.append(bond_line_id)
        return bond_objects

    def hide_atom(self, atom: eng.Atom) -> bool:
        """Determines if an atom should be hidden (and symbol not shown on the formula)"""
        if atom.symbol == "C" and self.toolbar.hide_carbon and len(atom.bonds) > 0:
            linked_atoms_sym = [bond.other_atoms(atom)[0].symbol for bond in atom.bonds]
            if "C" in linked_atoms_sym:
                return True
        return False

    def redraw_all(self) -> None:
        """Redraws all atoms and bonds on a freshly cleared canvas."""
        self.mol_canvas.delete("all")
        selected_objects: list[Union[eng.Atom, eng.CovBond]]
        selected_objects = [self.graphics[item] for item in self.selected]
        self.selected = []
        self.graphics = {}
        for bond_id, bond in enumerate(self.bonds):
            x1, y1, x2, y2 = bond.coords
            bondlen = bond.length
            if bondlen == 0:
                continue
            bond_order = bond.order
            dativity = bond.dativity
            cut_bond = (
                not self.hide_atom(bond.atoms[0]),
                not self.hide_atom(bond.atoms[1]),
            )
            bond_drawings: list[int] = self.draw_bond(
                (x1, y1, x2, y2),
                bondlen,
                bond_order,
                dativity,
                cut=cut_bond,
                tags=("bond", f"bond-{bond_id}"),
            )
            for bondid in bond_drawings:
                self.graphics[bondid] = bond
            if bond in selected_objects:
                self.selected += bond_drawings
            self.mol_canvas.tag_bind(
                f"bond-{bond_id}",
                "<Button-1>",
                lambda event, bond_id=bond_id: self.click_bond(
                    bond_id, event
                ),
            )
        for atom_id, atom in enumerate(self.atoms):
            atom_s: list[int] = []
            if self.hide_atom(atom):
                atom_s.append(self.mol_canvas.create_rectangle(
                    atom.coord_x - MINDIST / 3,
                    atom.coord_y - MINDIST / 3,
                    atom.coord_x + MINDIST / 3,
                    atom.coord_y + MINDIST / 3,
                    width=0,
                    tags=("atom", f"atom-{atom_id}", "empty"),
                ))
            else:
                textrep: str = atom.symbol
                if (
                    self.toolbar.empty_val_style ==
                    self.toolbar.EmptyValence.HYDROGENS
                    and atom.empty_valence > 0
                ):
                    textrep += "H"
                atom_s.append(self.mol_canvas.create_text(
                    atom.coord_x,
                    atom.coord_y,
                    text=textrep,
                    justify="center",
                    font=self.symbol_font,
                    tags=("atom", f"atom-{atom_id}", "atom_text", f"atom_text-{atom_id}"),
                    ))
                if (
                    self.toolbar.empty_val_style ==
                    self.toolbar.EmptyValence.HYDROGENS
                    and atom.empty_valence > 1
                ):
                    (_, bby1, bbx2, bby2) = self.mol_canvas.bbox(atom_s[0])
                    atom_s.append(self.mol_canvas.create_text(
                        bbx2 + 2,
                        (bby1 + bby2 * 3) / 4,
                        text=f"{atom.empty_valence}",
                        justify="left",
                        font=self.symbol_index_font,
                        tags=("atom_text", f"atom_text-{atom_id}", "index")
                        ))
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
