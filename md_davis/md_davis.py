import os
from tkinter import Tk
from tkinter import ttk

dirname = os.path.dirname(__file__)

root = Tk()
root.title('MD DaVis')
root.iconbitmap(os.path.join(dirname, 'MD_DaVis.ico'))

frm = ttk.Frame(root, padding=10)
frm.grid()

ttk.Label(frm, text="Hello World!").grid(column=0, row=0)
ttk.Button(frm, text="Quit", command=root.destroy).grid(column=1, row=0)

root.mainloop()

