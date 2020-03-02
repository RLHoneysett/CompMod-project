"""
Tool to convert from metric units to the reduced units used in
Particle3D and particle_manybody
"""

import tkinter as tk


# Declare global variables
mol_vol_metric = None
mol_vol_reduced = None

temp_metric = None
temp_reduced = None

length_metric = None
length_reduced = None

# Function called when button is pressed
def convert_mol():
    
    global mol_vol_metric
    global mol_vol_reduced
    
    #Convert units
    try:
        val = mol_vol_metric.get()
        mol_vol_reduced.set(val * 0.04206275737)
    except:
        pass

def convert_temp():
    
    global temp_metric
    global temp_reduced
    
    try:
        val = temp_metric.get()
        temp_reduced.set(val * 8.347245409e-3)
    except:
        pass

def convert_length():
    
    global length_metric
    global length_reduced

    try:
        val = length_metric.get()
        length_reduced.set(val * 0.2936857562)
    except:
        pass

# Create the main window
root = tk.Tk()
root.title("LJ Reduced Unit Converter")



# Create the main container
frame = tk.Frame(root, bg ="#33353F")

# Lay out main container, specify we want it to grow with window size
frame.pack(fill=tk.BOTH, expand = True)

# Allow middle cell to grow with window size
frame.columnconfigure(1, weight=0)
frame.rowconfigure(4, weight=0)

# Variables for holding data
mol_vol_metric = tk.DoubleVar()
mol_vol_reduced = tk.DoubleVar()

temp_metric = tk.DoubleVar()
temp_reduced = tk.DoubleVar()

length_metric = tk.DoubleVar()
length_reduced = tk.DoubleVar()

# Create widgets (formatting is messy)
entry_mol_metric = tk.Entry(frame, width=7, textvariable = mol_vol_metric, relief=tk.FLAT, bg="#DCDCDC")
label_mol = tk.Label(frame, text = "Molar Volume", bg = "#33353F", fg="#FFFFFF", font="Lucidasans 12")
label_mol_unit_metric = tk.Label(frame, text = "cm³/mol", bg = "#33353F", fg="#FFFFFF", font="Lucidasans 10 italic")
label_mol_equal = tk.Label(frame, text = "converts to", bg = "#33353F", fg="#FFFFFF", font="Lucidasans 12")
label_mol_reduced = tk.Label(frame, textvariable = mol_vol_reduced, bg = "#33353F", fg="#FFFFFF", font="Lucidasans 12")
label_mol_unit_reduced = tk.Label(frame, text = "in Reduced Units", bg = "#33353F", fg="#FFFFFF", font="Lucidasans 12")
button_mol_convert = tk.Button(frame, text="Convert", command=convert_mol, bg="#F48C22",activebackground="#CE751A", activeforeground="#000000",fg="#000000", relief=tk.FLAT, font="Lucidasans", highlightthickness=0)

entry_temp_metric = tk.Entry(frame, width=7, textvariable = temp_metric, relief=tk.FLAT, bg="#DCDCDC")
label_temp = tk.Label(frame, text = "Temperature", bg = "#33353F", fg="#FFFFFF", font="Lucidasans 12")
label_temp_unit_metric = tk.Label(frame, text = "K", bg = "#33353F", fg="#FFFFFF", font="Lucidasans 10 italic")
label_temp_equal = tk.Label(frame, text = "converts to", bg = "#33353F", fg="#FFFFFF", font="Lucidasans 12")
label_temp_reduced = tk.Label(frame, textvariable = temp_reduced, bg = "#33353F", fg="#FFFFFF", font="Lucidasans 12")
label_temp_unit_reduced = tk.Label(frame, text = "in Reduced Units", bg = "#33353F", fg="#FFFFFF", font="Lucidasans 12")
button_temp_convert = tk.Button(frame, text = "Convert", command = convert_temp, bg="#F48C22",activebackground="#CE751A", activeforeground="#000000",fg="#000000", relief=tk.FLAT, font="Lucidasans", highlightthickness=0)

entry_length_metric = tk.Entry(frame, width=7, textvariable = length_metric, relief=tk.FLAT, bg="#DCDCDC")
label_length = tk.Label(frame, text = "Length", bg = "#33353F", fg="#FFFFFF", font="Lucidasans 12")
label_length_unit_metric = tk.Label(frame, text = "Å", bg = "#33353F", fg="#FFFFFF", font="Lucidasans 10 italic")
label_length_equal = tk.Label(frame, text = "converts to", bg = "#33353F", fg="#FFFFFF", font="Lucidasans 12")
label_length_reduced = tk.Label(frame, textvariable = length_reduced, bg = "#33353F", fg="#FFFFFF", font="Lucidasans 12")
label_length_unit_reduced = tk.Label(frame, text = "in Reduced Units", bg = "#33353F", fg="#FFFFFF", font="Lucidasans 12")
button_length_convert = tk.Button(frame, text = "Convert", command = convert_length, bg="#F48C22",activebackground="#CE751A", activeforeground="#000000",fg="#000000", relief=tk.FLAT, font="Lucidasans", highlightthickness=0)

# Lay out widgets
entry_mol_metric.grid(row=0, column=1, padx=5, pady=5)
label_mol.grid(row=0, column=0, padx=5, pady=5, sticky=tk.E)
label_mol_unit_metric.grid(row=0, column=2, padx=5, pady=5, sticky=tk.W)
label_mol_equal.grid(row=1, column=0, padx=5, pady=5, sticky=tk.E)
label_mol_reduced.grid(row=1, column=1, padx=5, pady=5)
label_mol_unit_reduced.grid(row=1, column=2, padx=5, pady=5, sticky=tk.W)
button_mol_convert.grid(row=2, column=1, columnspan=2, padx=5, pady=5, sticky=tk.E)

entry_temp_metric.grid(row=3, column=1, padx=5, pady=5)
label_temp.grid(row=3, column=0, padx=5, pady=5, sticky=tk.E)
label_temp_unit_metric.grid(row=3, column=2, padx=5, pady=5, sticky=tk.W)
label_temp_equal.grid(row=4, column=0, padx=5, pady=5, sticky=tk.E)
label_temp_reduced.grid(row=4, column=1, padx=5, pady=5)
label_temp_unit_reduced.grid(row=4, column=2, padx=5, pady=5, sticky=tk.W)
button_temp_convert.grid(row=5, column=1, columnspan=2, padx=5, pady=5, sticky=tk.E)

entry_length_metric.grid(row=6, column=1, padx=5, pady=5)
label_length.grid(row=6, column=0, padx=5, pady=5, sticky=tk.E)
label_length_unit_metric.grid(row=6, column=2, padx=5, pady=5, sticky=tk.W)
label_length_equal.grid(row=7, column=0, padx=5, pady=5, sticky=tk.E)
label_length_reduced.grid(row=7, column=1, padx=5, pady=5)
label_length_unit_reduced.grid(row=7, column=2, padx=5, pady=5, sticky=tk.W)
button_length_convert.grid(row=8, column=1, columnspan=2, padx=5, pady=5, sticky=tk.E)


# Place cursor in entry box by default
entry_mol_metric.focus()


# Run forever
root.mainloop()

    
