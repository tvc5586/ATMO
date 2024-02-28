import PySimpleGUI as sg
import os
import datetime
import re
import pandas as pd
import numpy as np
from netCDF4 import Dataset
from os import system, name
import time

import traceback 

import parallel
import plot
import convert
import boxmod as bm

# TODO: Customize conversion folder name
# TODO: Allow users to check scenario settings
# TODO: Allow users to customize height setting
# TODO: Allow users to save scenario settings to harddrive

# ====================================================================================================================================================

sg.theme("default1") # Set GUI Theme

# Heights for plotting
z = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
     1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 
     5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 
     10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 
     14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 
     18.0, 18.5, 19.0, 19.5, 20, 21, 22, 23, 24, 25, 
     26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 
     38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 
     50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 
     62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 
     74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 
     86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 
     98, 99, 100, 105, 110, 115, 120, 125, 130, 135, 
     140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 
     190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 
     240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 
     290, 295, 300, 305, 310, 315, 320, 325, 330, 335, 
     340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 
     390, 395, 400, 420, 440, 460, 480, 500, 520, 540, 
     560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 
     760, 780, 800, 820, 840, 860, 880, 900, 920, 940, 
     960, 980, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 
     1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 
     2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 
     3500, 3600, 3700, 3800, 3900, 4000]

# ====================================================================================================================================================



# ====================================================================================================================================================
   
def main():
    
    # =============================== Computation Mode Elements ===============================
                             
    env_page = [[sg.Button("Select Equations", size = (25, 1), key = "select_eqn")],
                [sg.Text("Selected equation files:", key = "eqn_file_show")],
                [sg.Button("Generate Specie Config", disabled = True, size = (25, 1), key = "specie_config_gen")],
                [sg.Button("Select Sppecie Config File", size = (25, 1), key = "specie_config_sel")],
                [sg.Text("Selected specie config file:", key = "specie_file_show")],]
    
    binary_switches = [[sg.Text("Box Model Mode", size = (35, 1)), sg.Button('Off', size = (5, 1), button_color = "white on red", key = "box_model")],
                       [sg.Text("Altitude Calculation", size = (35, 1)), 
                        sg.Button('On', size = (5, 1), button_color = "white on green", key = "altitude")],
                       [sg.Text("Vertical Mix", size = (35, 1)), 
                        sg.Button('On', size = (5, 1), button_color = "white on green", key = "vertical_mix")],
                       [sg.Text("Surface Decomposition", size = (35, 1)), 
                        sg.Button('On', size = (5, 1), button_color = "white on green", key = "decomposition")],
                       [sg.Text("Temperature Calculation", size = (35, 1)), 
                        sg.Button('Off', size = (5, 1), button_color = "white on red", key = "temp"), 
                        sg.Button('?', size = (1, 1), key = "temp_info")],
                       ]

    numerical_values = [[sg.Text("Output Folder", size = (25, 1)),
                         sg.In(size = (13, 1), enable_events = True, key = "-FOLDER-"),
                         sg.FolderBrowse(),],
                        [sg.Text("Name of NetCDF Output File", size = (25, 1)), sg.In("Output", size = (13, 1), key = "name")],
                        [sg.Text("Generate Other Formats", size = (25, 1)), 
                         sg.Checkbox('CSV (.csv)', default=False, key = "CSV"), 
                         sg.Checkbox('Text (.txt)', default=False, key = "TXT")],
                        [sg.Text("Running Time (Day)", size = (25, 1)), sg.In("1", size = (5, 1), key = "running_time")],
                        [sg.Text("Time Stamp (Second)", size = (25, 1)), sg.In("300", size = (5, 1), key = "deltim")],]
    
    date_location = [[sg.Text("Year", size = (35, 1)), sg.In("2012", size = (8, 1), key = "DL1")],
                     [sg.Text("Month", size = (35, 1)), sg.In("3", size = (8, 1), key = "DL2")],
                     [sg.Text("Day", size = (35, 1)), sg.In("24", size = (8, 1), key = "DL3")],
                     [sg.Text("Latitude", size = (35, 1)), sg.In("71.2906", size = (8, 1), key = "DL4")],
                     [sg.Text("Longitude", size = (35, 1)), sg.In("156.7886", size = (8, 1), key = "DL5")],
                     [sg.Text("STD Longitude", size = (35, 1)), 
                      sg.In("156.7886", size = (8, 1), key = "DL6"), 
                      sg.Button('?', size = (1, 1), key = "STDLong_info")]]
    
    # =============================== Conversion Mode Elements ===============================
    
    convert_func = [[sg.Text("Input Location"), sg.In(size = (25, 1), enable_events = True, key = "-INNETCDF-"), sg.FolderBrowse(),],
                    [sg.Listbox(values=[], enable_events = True, size = (50, 18), key = "-NETCDF LIST-")],
                    [sg.Text("Output Location"), sg.In(size = (25, 1), key = "-OUTLOCATION-"), sg.FolderBrowse(),],
                    [sg.Text("Select Formats"), 
                     sg.Checkbox('CSV (.csv)', default=False, key = "TOCSV"), 
                     sg.Checkbox('Text (.txt)', default=False, key = "TOTXT")]]
    
    # =============================== Plot Mode Elements ===============================
    
    file_places = [[sg.Text("Input Location"), sg.In(size = (25, 1), enable_events = True, key = "-INFILE-"), sg.FolderBrowse(),],
                   [sg.Listbox(values=[], enable_events = True, size = (50, 18), key = "-FILE LIST-")]]
    
    line_plot = [[sg.Text("Specie", size = (18, 1)), sg.Combo([], size = (15, 1), key = "line_specie")],
                 [sg.Text("Height: 0.1m", size = (18, 1), key = "height_display"), 
                  sg.Slider(range=(1, 248), 
                  resolution = 1, 
                  orientation ='horizontal', 
                  disable_number_display = True, 
                  enable_events = True, 
                  key='line_height')]]
    
    contourf_plot = [[sg.Text("Specie", size = (18, 1)), 
                      sg.Combo([], size = (15, 1),
                      key = "contourf_specie")],
                      [sg.Text("Smoothness: 5", size = (18, 1), key = "smooth"), 
                       sg.Slider(range=(5, 100), 
                       resolution = 5, 
                       orientation ='horizontal', 
                       disable_number_display = True, 
                       enable_events = True, 
                       key='contourf_level'),
                       sg.Button('?', size = (1, 1), key = "contourf_info")]]
    
    # =============================== Sub Menu Layouts ===============================
    
    no_task = [[sg.VPush()],
               [sg.Text("No active task", size = (23, 1), key = "Check")],
               [sg.VPush()]]
        
    csv_prog = [[sg.VPush()],
                [sg.Text("CSV conversion progress:", size = (23, 1)), sg.Text("0%", key = "csvper"),
                 sg.ProgressBar(50, orientation = 'horizontal', size = (44, 15), border_width = 2, key = "csvbar")],
                [sg.VPush()]]
                
    txt_prog = [[sg.VPush()],
                [sg.Text("Text conversion progress:", size = (23, 1)), sg.Text("0%", key = "txtper"),
                 sg.ProgressBar(50, orientation = 'horizontal', size = (44, 15), border_width = 2, key = "txtbar")],
                [sg.VPush()]]
    
    plot_button = [[sg.VPush()],
                   [sg.Button("Select File", size = (15, 1), key = "file_in")],
                   [sg.Button("Line Plot", size = (15, 1), disabled = True, key = "Line_Plot")],
                   [sg.Button("Contourf Plot", size = (15, 1), disabled = True, key = "Contourf_Plot")],
                   [sg.Button("Plot", size = (15, 1), disabled = True, key = "show")],
                   [sg.VPush()]]
    
    convert_button = [[sg.VPush()],
                      [sg.Button("Convert", size = (15, 1), disabled = True, key = "convert")],
                      [sg.VPush()]]
    
    compute_button = [[sg.VPush()],
                      [sg.Button("Env Config", size = (15, 1), key = "env_config")],
                      [sg.Button("Binary Switch", size = (15, 1), key = "binary_switch")],
                      [sg.Button("Numerical Input", size = (15, 1), key = "numerical_input")],
                      [sg.Button("Date & Location", size = (15, 1), key = "D&L")],
                      [sg.Button("Save Scenario", size = (15, 1), key = "save_scene")],
                      [sg.Button("Execute", disabled = True, size = (15, 1), key = "run")],
                      [sg.VPush()]]

    # =============================== Main Layouts ===============================
    
    plot_mode_layout = [[sg.VPush()],
                        [sg.Text("Plot Mode Menu", font = ("Arial", 14))],
                        [sg.HorizontalSeparator(color = "black")],
                        [sg.Column(plot_button, size = (150, 435)),
                         sg.VSeparator(color = "black"),
                         sg.Column(file_places, size = (470, 435), pad = ((0, 0),(10, 0)), visible = True, key = "-FML-"),
                         sg.Column(line_plot, size = (470, 435), pad = ((0, 0),(10, 0)), visible = False, key = "-LS-"),
                         sg.Column(contourf_plot, size = (470, 435), pad = ((0, 0),(10, 0)), visible = False, key = "-CS-")],
                        [sg.VPush()]]
    
    conversion_mode_layout = [[sg.VPush()],
                              [sg.Text("Conversion Mode Menu", font = ("Arial", 14))],
                              [sg.HorizontalSeparator(color = "black")],
                              [sg.Column(convert_button, size = (150, 435)),
                               sg.VSeparator(color = "black"),
                               sg.Column(convert_func, size = (470, 435), pad = ((0, 0),(10, 0)))],
                              [sg.VPush()]]
                                 
    computation_mode_layout = [[sg.VPush()],
                               [sg.Text("Computation Mode Menu", font = ("Arial", 14))],
                               [sg.HorizontalSeparator(color = "black")],
                               [sg.Column(compute_button, size = (150, 435)),
                                sg.VSeparator(color = "black"),
                                sg.Column(env_page, size = (470, 435), pad = ((0, 0),(10, 0)), visible = True, key = "env_face"),
                                sg.Column(binary_switches, size = (470, 435), pad = ((0, 0),(10, 0)), visible = True, key = "binary_face"),
                                sg.Column(numerical_values, size = (470, 435), pad = ((0, 0),(10, 0)), visible = False, key = "numerical_face"),
                                sg.Column(date_location, size = (470, 435), pad = ((0, 0),(10, 0)), visible = False, key = "date_loc_face")],
                               [sg.VPush()]]
    
    main_layout = [[sg.Text("Main Menu", font = ("Arial", 15))],
                   [sg.Button("Computation Mode", size = (15, 2), key = "Compute")],
                   [sg.Button("Conversion Mode", size = (15, 2), key = "Convert")],
                   [sg.Button("Plot Mode", size = (15, 2), key = "Plot")],
                   [sg.Exit(size = (15, 2))]]
                   
    layout = [[sg.VPush()],
              [sg.Text("ATmospheric Model - One-dimensional", font = ("Arial", 20))],
              [sg.Text("Version 2"),],
              [sg.HorizontalSeparator(color = "black")],
              [sg.Column(main_layout, size = (150, 490), element_justification = "center", key = "-ML-"),
               sg.VSeparator(color = "black"),
               sg.Column(computation_mode_layout, size = (620, 490), visible = True, element_justification = "left", key = "-CML-"),
               sg.Column(plot_mode_layout, size = (620, 490), visible = False, element_justification = "left", key = "-PML-"),
               sg.Column(conversion_mode_layout, size = (620, 490), visible = False, element_justification = "left", key = "-OML-"),],
              [sg.HorizontalSeparator(color = "black")],
              [sg.Column(no_task, size = (770, 40), visible = True, element_justification = "center", key = "no_task"),
               sg.Column(csv_prog, size = (770, 40), visible = False, element_justification = "center", key = "-CSV-"), 
               sg.Column(txt_prog, size = (770, 40), visible = False, element_justification = "center", key = "-TEXT-")],
              [sg.VPush()],]
              
    window = sg.Window("ATMO", layout, 
                       size = (770, 665), 
                       resizable = False, 
                       enable_close_attempted_event = True, 
                       element_justification = "center").finalize()
    
    # =============================== Local Variables ===============================
    
    # Values for buttons in Computation Mode 
    box_model = True
    altitude  = True
    vmix      = True
    decomp    = True
    temp      = True
        
    # Values for buttons to change plot type
    line     = False
    contourf = False
    
    # Create plot object
    plt = plot.Plot()
    
    # Name holders
    filename   = ""
    eqn_file_list = []
    all_spc = []
    specie_file_name = ""
    spec_csv_only_check = False
    input_array = []
    scene_num = 0
    netcdfname = ""
    file_selected = False
    specie_list = []
        
    # =============================== Run the Event Loop ===============================
    
    while True:
    
        event, values = window.read()
    
        # =============================== Main Menu Functions ===============================
    
        if ((event == "Exit" or event == sg.WINDOW_CLOSE_ATTEMPTED_EVENT) and 
             sg.PopupOKCancel('Are you sure to exit the program?', title = "Exit")) == "OK":
            # Quit program
            break
            
        elif event == "Compute":
            # Switch Layout
            window['-OML-'].update(visible = False)
            window['-FML-'].update(visible = False)
            window['-PML-'].update(visible = False)
            window['-LS-'].update(visible = False)
            window['-CS-'].update(visible = False)
            window['-CML-'].update(visible = True)
        
        elif event == "Convert":
            # Switch Layout
            window['-CML-'].update(visible = False)
            window['-FML-'].update(visible = False)
            window['-PML-'].update(visible = False)
            window['-LS-'].update(visible = False)
            window['-CS-'].update(visible = False)
            window['-OML-'].update(visible = True)
        
        elif event == "Plot":
            # Switch Layout
            window['-OML-'].update(visible = False)
            window['-CML-'].update(visible = False)
            window['-LS-'].update(visible = False)
            window['-CS-'].update(visible = False)
            window['-FML-'].update(visible = True)
            window['-PML-'].update(visible = True)
            
        # =============================== Computation Mode Functions ===============================

        # ------------------------------- Switch sub-layouts -------------------------------

        elif event == "env_config":
            window["env_face"].update(visible = True)
            window["binary_face"].update(visible = False)
            window["numerical_face"].update(visible = False)
            window["date_loc_face"].update(visible = False)

        elif event == "binary_switch":
            window["env_face"].update(visible = False)
            window["binary_face"].update(visible = True)
            window["numerical_face"].update(visible = False)
            window["date_loc_face"].update(visible = False)
        
        elif event == "numerical_input":
            window["env_face"].update(visible = False)
            window["binary_face"].update(visible = False)
            window["numerical_face"].update(visible = True)
            window["date_loc_face"].update(visible = False)
                    
        elif event == "D&L":
            window["env_face"].update(visible = False)
            window["binary_face"].update(visible = False)
            window["numerical_face"].update(visible = False)
            window["date_loc_face"].update(visible = True)
        
        # ------------------------------- Select Equation Files -------------------------------
        
        elif event == "select_eqn":

            file_list = ""
            
            file_path_string = sg.popup_get_file('Select Equation Files', multiple_files = True)
            
            if file_path_string != None:
              eqn_file_list = list(file_path_string.split(";"))
             
            try:
              eqn_func_check = bm.function_not_exist(eqn_file_list) # Check if all functions are implemented
              
              species = bm.load.read_csvs(eqn_file_list, True)
              mech = bm.EqnSet(species, name="Specie_Gen", long_name="Specie_Generation")
              all_spc = mech.spcs # List of all species
                
              for i in eqn_file_list:              
                # for windows
                if name == 'nt':
                  file_list += i.split('\\')[-1] + "\n" 
                # for mac and linux(here, os.name is 'posix')
                else:
                  file_list += i.split('/')[-1] + "\n"
            
              window['eqn_file_show'].update("Select equation files:\n----------------------\n" + file_list)
              window['specie_config_gen'].update(disabled = False)  

            except Exception as e:
              traceback.print_exc()
              sg.popup("Could not read the equation files selected.\nPlease select the correct equation files", title = "Error")
                
        # ------------------------------- Generate Specie File -------------------------------

        elif event == "specie_config_gen":
            try:
              specie_file_name = sg.popup_get_text("Enter a name for the specie configuration file.") + ".csv"
            
              if specie_file_name == "": # If user enters nothing
                specie_file_name = "Specie Config"
            
              temp = np.full((len(all_spc), 5), "", dtype = object)
              temp[:, 0] = all_spc
            
              pd.DataFrame(temp, columns = ["Specie", "Initial Value", "Unit", "henry", "f0"]).to_csv(os.path.join(os.getcwd(), specie_file_name), index = False, header = True) 
            
              sg.popup("A file has been generated for configuring initial values.\nThe unit for each specie must be either 'ppbv' or 'molec_cm3'.\nPlease fill the file before continuing.", title = "Notice")

              window['specie_file_show'].update("Selected specie config file:\n------------------------------\n" + specie_file_name)
            
            except: # If user cancels
              pass
        
        # ------------------------------- Select Specie File -------------------------------
        
        elif event == "specie_config_sel":
            try:
                specie_file_name = sg.popup_get_file('Select Specie Config Files', multiple_files = False)
                spec_csv_only_check = True if specie_file_name[-4:] == ".csv" else False
                
                temp_specie_file_name = ""
                
                # for windows
                if name == 'nt':
                  temp_specie_file_name = specie_file_name.split('\\')[-1] 
                # for mac and linux(here, os.name is 'posix')
                else:
                  temp_specie_file_name = specie_file_name.split('/')[-1]
                
                window['specie_file_show'].update("Selected specie config file:\n------------------------------\n" + temp_specie_file_name)

            except: # If user cancels
                pass
            
        # ------------------------------- Save Scenario -------------------------------
            
        elif event == "save_scene":
            if values['-FOLDER-'] == "" or values['name'] == "":
                sg.popup("Output location cannot be empty!", title = "Error")
                        
            elif specie_file_name == "":
                sg.popup("A specie file must be presented!", title = "Error")
                
            elif eqn_file_list == []:
                sg.popup("An equation file must be presented!", title = "Error")
                
            elif not spec_csv_only_check:
                sg.popup("Specie file can only be of CSV type!", title = "Error")
                
            else:    
                try:
                    # Check if inputs are numbers
                    Running_days = float(values['running_time'])
                    Deltim = float(values['deltim'])
                    
                    # Check if date and location are in correct format   
                    Date = re.split('-', re.split('\s+', str(datetime.datetime(int(values['DL1']), int(values['DL2']), int(values['DL3']))))[0])
                    values['DL1'] = Date[0]
                    values['DL2'] = Date[1]
                    values['DL3'] = Date[2]
                    
                    # Check if inputs are negative
                    if Running_days < 0:
                        sg.popup("Running_days cannot be negative!", title = "Error")
                    
                    elif Deltim <= 0:
                        sg.popup("Time stamp cannot be non-positive!", title = "Error")
                    
                    # Check if latitude is valid
                    elif abs(float(values['DL4'])) > 90:
                        sg.popup("Invalid Latitude!", title = "Error")
                
                    # Check if longitude is valid
                    elif abs(float(values['DL5'])) > 180:
                        sg.popup("Invalid Longitutde!", title = "Error")
             
                    # Check if STD longitude is valid
                    elif abs(float(values['DL6'])) > 180:
                        sg.popup("Invalid STD Longitutde!", title = "Error")
                            
                    else: 
                        DL_value  = [0] * 6
                    
                        for i in range(6):
                            if DL_value[i] != float(values[f'DL{i + 1}']):
                                DL_value[i] = float(values[f'DL{i + 1}']) # Assign configured initial values
                        
                        # Save input arguments to an array
                        input_array.append([values['-FOLDER-'], values['name'], 
                                            window.Element("box_model").get_text(),
                                            window.Element("altitude").get_text(),
                                            window.Element("vertical_mix").get_text(),
                                            window.Element("decomposition").get_text(),
                                            window.Element("temp").get_text(),
                                            Running_days, Deltim, DL_value, all_spc,
                                            eqn_file_list, os.path.join(os.getcwd(), specie_file_name)])
                        
                        # Unlock execution button
                        window['run'].update(disabled = False)
                        
                        # Update number of scenarios stored
                        scene_num += 1
                        window["Check"].update(str(scene_num) + " scenarios waiting")
                
                except Exception as e: 
                        if "day" in str(e):
                            sg.popup("Day is out of range for month!", title = "Error")
                    
                        elif "month" in str(e):
                            sg.popup("Month is out of range!", title = "Error")
                        
                        elif "year" in str(e):
                            sg.popup("Year is out of range!", title = "Error")
                            
                        else:
                            print(e)
                            sg.popup("Value should be numbers!", title = "Error") 

        # ------------------------------- Run Computation -------------------------------
        
        elif event == "run":
            #TODO: Confirmation window
                                                                          
              try:
                max_val = 0 # Reset max value
                                
                window["Check"].update("Processing " + str(scene_num) + " scenarios")
                   
                time.sleep(5) # Ensure GUI update taskbar
                   
                # Run CLI Program
                for i in (parallel.multiprocess(input_array)):
                
                  event, values = window.read(10)
                
                  # Update number of processing scenarios
                  scene_num -= 1
                  window["Check"].update("Processing " + str(scene_num) + " scenarios")
                                    
                # Conversion
                if values['CSV']:
                # Show CSV Conversion Progression Bar
                  window['-CSV-'].update(visible = True)
                  CSV_max = True
                  
                  filename = os.path.join(values['-FOLDER-'], values['name'] + ".nc")
                  specie_list = list(Dataset(filename, mode = "r").variables['Specie'][:])
                                    
                  for i in (convert.Convert.To_CSV(filename), values['-FOLDER-'], specie_list):
                    event, values = window.read(10)
                    
                    # Get maximum value
                    if CSV_max:
                      CSV_max = False
                      max_val = i
                      window["csvbar"].update(max = max_val, current_count = 0)
                        
                    # Update Progression Bar
                    else:
                      window["csvper"].update(str(int(i/max_val*100)) + "%")
                      window["csvbar"].update(max = max_val, current_count = i)
                                    
                  # Hide CSV Conversion Progression Bar
                  window['-CSV-'].update(visible = False)
                                    
                if values['TXT']:
                  # Show Text Conversion Progression Bar
                  window['-TEXT-'].update(visible = True)
                  Text_max = True
                  
                  filename = os.path.join(values['-FOLDER-'], values['name'] + ".nc")
                  specie_list = list(Dataset(filename, mode = "r").variables['Specie'][:])
                           
                  for i in (convert.Convert.To_Text(filename), values['-FOLDER-'], specie_list):
                    event, values = window.read(10)
                
                    # Get maximum value
                    if Text_max:
                      Text_max = False
                      max_val = i
                      window["txtbar"].update(max = max_val, current_count = 0)
                        
                    # Update Progression Bar
                    else:
                      window["txtper"].update(str(int(i/max_val*100)) + "%")
                      window["txtbar"].update(max = max_val, current_count = i)
                                    
                # Hide Text Conversion Progression Bar
                window['-TEXT-'].update(visible = False)
            
                window["Check"].update(str("No active task"))
                                
              except Exception as e: 
                traceback.print_exc() 
                sg.popup("Unknown error. Please contact the developers.", title = "Error") 
                      
        # ------------------------------- Toggle Binary Switches -------------------------------
          
        elif event == "box_model":
            box_model = not box_model
            window.Element("box_model").Update(("On", "Off")[box_model], button_color = (("white", ("green", "red")[box_model])))
            
        elif event == "altitude":
            altitude = not altitude
            window.Element("altitude").Update(("Off", "On")[altitude], button_color = (("white", ("red", "green")[altitude])))
            
        elif event == "vertical_mix":
            vmix = not vmix
            window.Element("vertical_mix").Update(("Off", "On")[vmix], button_color = (("white", ("red", "green")[vmix])))
            
        elif event == "decomposition":
            decomp = not decomp
            window.Element("decomposition").Update(("Off", "On")[decomp], button_color = (("white", ("red", "green")[decomp])))
            
        elif event == "temp":
            temp = not temp
            window.Element("temp").Update(("On", "Off")[temp], button_color = (("white", ("green", "red")[temp])))
                    
        # Information for some switches
        elif event == "temp_info":
            sg.popup(("Determine if temperature calculation is done by using a function. "
                      "Default (Off) means temperature is the same across all heights."), title = "Temperature Calculation")
                                  
        elif event == "STDLong_info":
            sg.popup("STD Longtitude is 75 plus 15 for daylight saving time", title = "STD Longtitude")
                        
        # =============================== Conversion Mode Functions ===============================
                    
        elif event == "-INNETCDF-":
            # Select folder where input file is located
            folder = values["-INNETCDF-"]
            
            try:
                # Get list of files in folder
                file_list = os.listdir(folder)
            except FileNotFoundError:
                file_list = []
            
            fnames = [f for f in file_list if os.path.isfile(os.path.join(folder, f)) and f.lower().endswith((".nc"))] # Put pathes of all netCDF files in a list
            fnames.sort(key=str.lower) # Sort list alphabetically
            window["-NETCDF LIST-"].update(fnames) # Show avaiable input files
        
        elif event == "-NETCDF LIST-":
            # Select input file
            try:
                netcdfname = os.path.join(values["-INNETCDF-"], values["-NETCDF LIST-"][0]) # Update input file path
                window['convert'].update(disabled = False) # Enable convert button
                            
            except:
                pass 
                
        elif event == "convert":
            # Check if output location is selected
            if values["-OUTLOCATION-"] == "":
                sg.popup("Output location cannot be empty!", title = "Error")
            
            # No formats are selected
            elif not values['TOCSV'] and not values['TOTXT']:
                sg.popup("No formats are selected!", title = "Error")
            
            else:
                max_val = 0 # Reset max value
                specie_list = list(Dataset(netcdfname, mode = "r").variables['Specie'][:]) # Get specie names
                
                window["txtper"].update(str(max_val) + "%") # Reset progress percentage
                window["txtbar"].update(max = max_val, current_count = 0) # Reset progress bar
                window["csvper"].update(str(max_val) + "%") # Reset progress percentage
                window["csvbar"].update(max = max_val, current_count = 0) # Reset progress bar
                
                window["no_task"].update(visible = False)
                
                # Switch Layout
                if values['TOTXT'] == True and values['TOCSV'] == False:
                    window['-CSV-'].update(visible = False)
                    window['-TEXT-'].update(visible = True)
            
                else:
                    window['-CSV-'].update(visible = True)
                    window['-TEXT-'].update(visible = False)
             
                # Run Conversion
                if values['TOCSV'] == values['TOTXT']:
                    try:    
                        # Convert to csv
                        for i in (convert.Convert.To_CSV(netcdfname, values["-OUTLOCATION-"], specie_list)):
                            event, values = window.read(10)
                    
                            # Get maximum value
                            if max_val == 0:
                                max_val = i
                                window["csvbar"].update(max = max_val, current_count = 0)
                        
                            # Update Progression Bar
                            else:
                                window["csvper"].update(str(int(i/max_val*100)) + "%")
                                window["csvbar"].update(max = max_val, current_count = i)
                
                        # Switch layout
                        window['-CSV-'].update(visible = False)
                        window['-TEXT-'].update(visible = True)
                     
                        max_val = 0 # Reset max value
                
                        # Convert to txt
                        for i in (convert.Convert.To_Text(netcdfname, values["-OUTLOCATION-"], specie_list)):
                            event, values = window.read(10)
                
                            # Get maximum value
                            if max_val == 0:
                                max_val = i
                                window["txtbar"].update(max = max_val, current_count = 0)
                        
                            # Update Progression Bar
                            else:
                                window["txtper"].update(str(int(i/max_val*100)) + "%")
                                window["txtbar"].update(max = max_val, current_count = i)
                        
                    except Exception as e: 
                        print(e)
                    
                    # Hide CSV and Text progress bars
                    window['-CSV-'].update(visible = False)
                    window['-TEXT-'].update(visible = False)
        
                else:
                    # Run Conversion
                    try:
                        if values['TOCSV'] == True:
                            for i in (convert.Convert.To_CSV(netcdfname, values["-OUTLOCATION-"], specie_list)):
                                event, values = window.read(10)
                    
                                # Get maximum value
                                if max_val == 0:
                                    max_val = i
                                    window["csvbar"].update(max = max_val, current_count = 0)
                        
                                # Update Progression Bar
                                else:
                                    window["csvper"].update(str(int(i/max_val*100)) + "%")
                                    window["csvbar"].update(max = max_val, current_count = i)
                
                        if values['TOTXT'] == True:
                            for i in (convert.Convert.To_Text(netcdfname, values["-OUTLOCATION-"], specie_list)):
                                event, values = window.read(10)
                
                                # Get maximum value
                                if max_val == 0:
                                    max_val = i
                                    window["txtbar"].update(max = max_val, current_count = 0)
                        
                                # Update Progression Bar
                                else:
                                    window["txtper"].update(str(int(i/max_val*100)) + "%")
                                    window["txtbar"].update(max = max_val, current_count = i)
                            
                    except Exception as e: 
                        print(e)
                    
                    # Hide CSV and Text progress bars
                    window['-CSV-'].update(visible = False)
                    window['-TEXT-'].update(visible = False)
            
                # Hide CSV and Text progress bars
                window["no_task"].update(visible = True)
        
        # =============================== Plot Mode Functions ===============================
        
        elif event == "file_in":
            # Switch Layout
            window['-OML-'].update(visible = False)
            window['-CML-'].update(visible = False)
            window['-LS-'].update(visible = False)
            window['-CS-'].update(visible = False)
            window['-PML-'].update(visible = True)
            window['-FML-'].update(visible = True)
                    
        elif event == "-INFILE-":
            # Select folder where input file is located
            folder = values["-INFILE-"]
            
            try:
                # Get list of files in folder
                file_list = os.listdir(folder)
            except FileNotFoundError:
                file_list = []
            
            fnames = [f for f in file_list if os.path.isfile(os.path.join(folder, f)) and f.lower().endswith((".nc", ".txt", ".csv"))] # Put pathes of all netCDF, text and CSV files in a list
            fnames.sort(key=str.lower) # Sort list alphabetically
            window["-FILE LIST-"].update(fnames) # Show avaiable input files
        
        elif event == "-FILE LIST-":
            # Select input file
            try:
                filename = os.path.join(values["-INFILE-"], values["-FILE LIST-"][0]) # Update input file path
                file_selected = True
                
                # Check filetype and update content
                if os.path.splitext(filename)[1] != ".nc":
                                    
                    window['line_specie'].update(disabled = True)
                    window['contourf_specie'].update(disabled = True)
                
                else:
                    specie_list = list(Dataset(filename, mode = "r").variables['Specie'][:])
                    window['line_specie'].update(values = specie_list, value=values['line_specie'])
                    window['contourf_specie'].update(values = specie_list, value=values['contourf_specie'])
                    
                    window['line_specie'].update(disabled = False)
                    window['contourf_specie'].update(disabled = False)
                            
            except:
                pass
                
            if file_selected:
                # Unlock functions
                window['Line_Plot'].update(disabled = False)
                window['Contourf_Plot'].update(disabled = False)
                window["show"].update(disabled = False)
        
        elif event == 'Line_Plot':
            # Switch Layout
            window['-OML-'].update(visible = False)
            window['-CML-'].update(visible = False)
            window['-FML-'].update(visible = False)
            window['-CS-'].update(visible = False)
            window['-PML-'].update(visible = True)
            window['-LS-'].update(visible = True)
            
            line = True
            contourf = False
        
        elif event == 'Contourf_Plot':
            # Switch Layout
            window['-OML-'].update(visible = False)
            window['-CML-'].update(visible = False)
            window['-FML-'].update(visible = False)
            window['-LS-'].update(visible = False)
            window['-PML-'].update(visible = True)
            window['-CS-'].update(visible = True)
            
            line = False
            contourf = True
        
        elif event == "line_height":
            # Display height
            window['height_display'].update("Height: " + str(z[int(values["line_height"] - 1)]) + "m")
        
        elif event == "contourf_level":
            # Display contourf level
            window["smooth"].update("Smoothness: " + str(int(values["contourf_level"])))
       
        elif event == "show":
            # Display Plots
            
            # Show line plot
            if line == True:
                
                # If input file is an netCDF file
                if os.path.splitext(filename)[1] == ".nc":
                    plt.Line_NC(z, filename, specie_list.index(values['line_specie']), values['line_specie'], values['line_height'] - 1)
                    
                # If input file is not an netCDF file
                else:
                    plt.Line(z, filename, values['line_height'] - 1)
            
            # Show contour plot                    
            elif contourf == True:
                
                # If input file is an netCDF file
                if os.path.splitext(filename)[1] == ".nc":
                    plt.Contourf_NC(filename, specie_list.index(values['contourf_specie']), values['contourf_specie'], int(values["contourf_level"]))
                    
                # If input file is not an netCDF file
                else:
                    plt.Contourf(z, filename, int(values["contourf_level"]))
        # Information for some switches
        elif event == "contourf_info":
            sg.popup("Determine how smooth the plot will be.", title = "Smoothness")
        
        
        # =============================== Mode Functions End ===============================
        
    # Close Window
    window.close()

# ====================================================================================================================================================

if __name__ == "__main__":
    main()

