import tkinter as tk
from tkinter import ttk
import sys

import analysis
from unit_converter import *

head_offset = 5
N_entries = 0
abs_column = 0

entries = []
text_entries = []
entry_unit_vars = []

def is_run_from_idle():
    return bool("idlelib" in sys.modules)

if is_run_from_idle():
    print("Please do not run the program on IDLE. Use the terminal/cmd instead.")
    quit()

class entry:
    def __init__(self, entry_field, entry_label, unit_field, unit_type):
        self.entry_field = entry_field
        self.entry_label = entry_label
        self.unit_field = unit_field
        self.unit_type = unit_type

    def get_value(self):
        if not self.entry_field.get("1.0", "end-1c") == "":
            if self.unit_type == "float":
                return float(self.entry_field.get("1.0", "end-1c"))
            elif self.unit_type == "int":
                return int(float(self.entry_field.get("1.0", "end-1c")))
            elif self.unit_type == "string":
                return self.entry_field.get("1.0", "end-1c")
        else:
            return 0

def focus_next_widget(event):
    event.widget.tk_focusNext().focus()
    return("break")

def show_about():
    about_win = tk.Toplevel()
    about_win.title("About LETALIS")
    about_text = "Liquid propellant rocket Engine Thermal AnaLysIS. (It sounded nice.)\n"
    about_text += "Thermal analysis and cooling system design tool for liquid propellant rocket engines.\n"
    about_text += "Written by H. A. Güler (arda-guler @ Github). All rights reserved.\n\n"
    about_text += "The software is provided as-is, without warranty of any kind, express or implied.\n"
    about_text += "Analysis results are preliminary, and variations from experimental results should be expected."
    about_label = tk.Label(about_win, text=about_text)
    about_label.pack()

def create_entry(label, units, unit_type):
    global N_entries, entries, head_offset, text_entries, entry_unit_vars

    abs_offset = head_offset + N_entries

    entry_unit = tk.StringVar()
    entry_unit.set(units[0])
    entry_unit_vars.append(entry_unit)
    
    new_label = tk.Label(mw, text=label, anchor='e')
    new_label.grid(row = abs_offset, column = 0 + abs_column*3, sticky="e")
    
    new_textfield = tk.Text(mw, height=1, width=10)
    new_textfield.grid(row = abs_offset, column = 1 + abs_column*3)
    new_textfield.bind("<Tab>", focus_next_widget)
    text_entries.append(new_textfield)
    
    new_unitpicker = tk.OptionMenu(mw, entry_unit, *units)
    new_unitpicker.config(width=10)
    new_unitpicker.grid(row = abs_offset, column=2 + abs_column*3)

    new_entry = entry(new_textfield, label, entry_unit, unit_type)
    entries.append(new_entry)

    N_entries += 1

def create_label(label):
    global N_entries, head_offset

    abs_offset = head_offset + N_entries

    new_label = tk.Label(mw, text=label)
    new_label.grid(row = abs_offset, column = 0 + abs_column*3, columnspan=3)

    N_entries += 1

def export_file():
    global entries

    save_filename = filename_field.get("1.0", "end-1c")
    if not save_filename.endswith(".lpre"):
        save_filename = save_filename + ".lpre"

    with open(save_filename, "w") as cf:
        for entry in entries:
            cf.write(str(entry.entry_label) + " " + str(entry.get_value()) + " "  + str(entry.unit_field.get()) + "\n")

def import_file():
    global entries, entry_unit_vars, text_entries
    global material_unit, no_unit, length_units, angle_units, mass_flow_units, pressure_units, temperature_units,\
           velocity_units, conductivity_units, molecular_mass_units, viscosity_units, time_units

    import_filename = filename_field.get("1.0", "end-1c")
    if not import_filename.endswith(".lpre"):
        import_filename = import_filename + ".lpre"

    import_file = open(import_filename, "r")
    import_lines = import_file.readlines()
    import_file.close()

    for n_line in range(len(entries)):
        line = import_lines[n_line]
        line = line[:-1]
        line = line.split(" ")

        cval = ""
        for element in line:
            if not element == line[-1]:
                try:
                    cval = float(element)
                except:
                    pass

                if cval == "":
                    if element == "SS" or element == "CCZ" or element == "Jet_A1" or element == "bell" or element == "conic":
                        cval = element
                
                text_entries[n_line].delete('1.0',"end")
                text_entries[n_line].insert('1.0',str(cval))
                
            else:
                entry_unit_vars[n_line].set(element)

def start_analysis():
    global material_unit, no_unit, length_units, angle_units, mass_flow_units,\
           pressure_units, temperature_units, velocity_units, conductivity_units,\
           molecular_mass_units, viscosity_units, time_units
    
    params = []
    for entry in entries:
        cvalue = entry.get_value()
        cunit = str(entry.unit_field.get())

        if type(cvalue) == float:
            if cunit in length_units:
                converted_value = convert_unit(cvalue, cunit, "m")
            elif cunit in angle_units:
                converted_value = convert_unit(cvalue, cunit, "deg")
            elif cunit in mass_flow_units:
                converted_value = convert_unit(cvalue, cunit, "kg/s")
            elif cunit in pressure_units:
                converted_value = convert_unit(cvalue, cunit, "Pa")
            elif cunit in temperature_units:
                converted_value = convert_unit(cvalue, cunit, "K")
            elif cunit in velocity_units:
                converted_value = convert_unit(cvalue, cunit, "m/s")
            elif cunit in conductivity_units:
                converted_value = convert_unit(cvalue, cunit, "W/(m*K)")
            elif cunit in molecular_mass_units:
                converted_value = convert_unit(cvalue, cunit, "g/mol")
            elif cunit in viscosity_units:
                converted_value = convert_unit(cvalue, cunit, "millipoise")
            elif cunit in time_units:
                converted_value = convert_unit(cvalue, cunit, "s")
            else:
                converted_value = cvalue
        else:
            converted_value = cvalue

        params.append(converted_value)
        
    analysis.perform(params)

print("Please don't close this window while working with LETALIS.")

mw = tk.Tk()
mw.title("LETALIS")
mw.iconbitmap('icon.ico')

material_unit = ["[Material]"]
nozzle_type_unit = ["[NozzleType]"]
no_unit = ["# (unitless)"]
percentage_unit = ["%"]
length_units = ["mm", "cm", "m", "km"]
angle_units = ["deg", "rad"]
mass_flow_units = ["kg/s"]
pressure_units = ["MPa", "kPa", "Pa", "bar", "atm"]
temperature_units = ["K", "C"]
velocity_units = ["m/s", "km/s"]
conductivity_units = ["W/(m*K)"]
molecular_mass_units = ["g/mol", "kg/kmol"]
viscosity_units = ["millipoise", "kg/(m*s)"]
time_units = ["s", "min", "hr"]

import_button = tk.Button(mw, text="Import Design", width=15, command=import_file)
import_button.grid(row=0, column=0, sticky="w")

export_button = tk.Button(mw, text="Export Design", width=15, command=export_file)
export_button.grid(row=1, column=0, sticky="w")

filename_field_label = tk.Label(mw, text="Filename/path:")
filename_field_label.grid(row=0, column=1, columnspan=2, sticky="w")
filename_field = tk.Text(mw, height=1, width=25)
filename_field.grid(row=1, column=1, columnspan=2)

analyze_button = tk.Button(mw, text="PERFORM ANALYSIS", width=25, height=1, font=("Arial",17), command=start_analysis, bg="red", fg="white")
analyze_button.grid(row=0, column=6, columnspan=3, rowspan=3)

about_button = tk.Button(mw, text="About", command=show_about)
about_button.grid(row=0, column=3, columnspan=3, rowspan=2)

hsep1 = ttk.Separator(mw, orient='horizontal')
hsep1.place(x=0, y=60, relwidth=1, relheight=0.2)

inputs_label = tk.Label(mw, text="DESIGN PARAMETERS", font=("Arial", 13))
inputs_label.grid(row=3, column=0, columnspan=9)

create_label("Thrust Chamber Geometry")
create_entry("Engine Length", length_units, "float")
create_entry("Combustion Chamber Diameter", length_units, "float")
create_entry("Throat Diameter", length_units, "float")
create_entry("Exit Diameter", length_units, "float")
create_entry("Chamber Contraction Angle", angle_units, "float")
create_entry("Chamber Radius of Curvature", length_units, "float")
create_entry("Nozzle Type: 'bell' or 'conic'?", nozzle_type_unit, "string")
create_label("Conic Nozzle Params")
create_entry("Nozzle Expansion (Half-)Angle", angle_units, "float")
create_entry("Throat Upstream Radius of Curvature", length_units, "float")
create_entry("Throat Downstream Radius of Curvature", length_units, "float")
create_label("Bell Nozzle Params")
create_entry("%Length to Conic (optional)", percentage_unit, "float")
create_entry("Initial Divergence (optional)", angle_units, "float")
create_entry("Exit Angle (optional)", angle_units, "float")
abs_column += 1
N_entries = 0

create_label("Regenerative Cooling System")
create_entry("Number of Coolant Channels", no_unit, "int")
create_entry("Inner Wall Thickness (ribs excluded)", length_units, "float")
create_entry("Coolant Channel Tangential Width", length_units, "float")
create_entry("Coolant Channel Radial Depth", length_units, "float")

create_label("Film Cooling System (WIP!)")
create_entry("Film Cooling Injection Point 1", length_units, "float")
create_entry("Film Cooling Mass Flow 1", mass_flow_units, "float")
create_entry("Film Cooling Injection Point 2", length_units, "float")
create_entry("Film Cooling Mass Flow 2", mass_flow_units, "float")

create_label("Thrust Chamber Thermodynamics")
create_entry("Combustion Chamber Mass Flow", mass_flow_units, "float")
create_entry("Chamber Pressure (p_c)", pressure_units, "float")
create_entry("Chamber Stagnation Temperature (T_c)", temperature_units, "float")
create_entry("Characteristic Velocity (C*)", velocity_units, "float")
create_entry("Gas Conductivity of Combustion Mixture", conductivity_units, "float")
create_entry("Avg. Molecular Mass of Combustion Mix.", molecular_mass_units, "float")
create_entry("Initial Wall Temperature", temperature_units, "float")
create_entry("Viscosity at Combustion Chamber", viscosity_units, "float")
create_entry("Gamma at Combustion Chamber", no_unit, "float")
create_entry("Viscosity at Throat", viscosity_units, "float")
create_entry("Gamma at Throat", no_unit, "float")
abs_column += 1
N_entries = 0

create_label("Materials and Material Properties")
create_entry("Inner Wall Material", material_unit, "string")
create_entry("Outer Shell Material", material_unit, "string")
create_entry("Coolant Fluid", material_unit, "string")
create_entry("Coolant Mass Flow Rate", mass_flow_units, "float")
create_entry("Coolant Manifold Temperature", temperature_units, "float")
create_entry("Coolant Manifold Pressure", pressure_units, "float")

create_label("Analysis Setup")
create_entry("Analysis Vertical Fineness", no_unit, "int")
create_entry("Analysis End Time", time_units, "float")
create_entry("Time Steps", time_units, "float")

mw.mainloop()
