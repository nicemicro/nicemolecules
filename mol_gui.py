#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 07:27:35 2022

@author: nicemicro
"""
import mol_eng as eng
import tkinter as tk
from tkinter import ttk

class AppContainer(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        atoms, bonds = eng.test3()
        self.title("Nice Molecules")
        
        container = ttk.Frame(self)
        container.grid(row=0, column=0, sticky="nsew")
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.mol_canvas = tk.Canvas(container, width=700, height=600)
        self.mol_canvas.grid(row=0, column=0, sticky="nsew")
        self.draw_all_molecules(atoms, bonds)
    
    def draw_all_molecules(self, atomlist, bondlist):
        for bond in bondlist:
            x1, y1, x2, y2 = bond.coords()
            bondlen = bond.length()
            bond_order = bond.order()
            obscure = 7
            xd1 = x1 + (obscure / bondlen) * (x2-x1)
            xd2 = x2 - (obscure / bondlen) * (x2-x1)
            yd1 = y1 + (obscure / bondlen) * (y2-y1)
            yd2 = y2 - (obscure / bondlen) * (y2-y1)
            self.mol_canvas.create_line(xd1, yd1, xd2, yd2)
        for atom in atomlist:
            self.mol_canvas.create_text(atom.coord_x, atom.coord_y,
                                        text=atom.letter, justify="center")
        
def main():
    app = AppContainer()
    app.mainloop()

if __name__ == '__main__':
    main()
