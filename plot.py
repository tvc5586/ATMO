from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os.path

class Plot():

  def Line_NC(self, z, input_file, specie, specie_name, height):
  
    # Get Data
    fh = Dataset(input_file, mode = "r")
    gas = fh.variables['Concentration']
    Kh = fh.variables["Kh"]
    
    # Calculate Time
    X, Y, Z = gas.shape
    day = (X-1)*2/24/(60/5)*60/60
    x_range = np.arange(0, day*24, (day*24/(X-1)))
    x_range = np.append(x_range, [day*24])
    
    # Plot
    if specie_name == "Kh":
      plt.plot(x_range, Kh[:, height], label = specie_name)
      
    else:
      plt.plot(x_range, gas[:, height, specie], label = specie_name)
    
    # Put attributes
    if specie_name == "Kh":
      plt.title(specie_name + ' at height ' + str(z[int(height)]) + "m")
      plt.ylabel('m\N{SUPERSCRIPT TWO}s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}')
    
    else:
      plt.title('Concentration of '+ specie_name + ' at height ' + str(z[int(height)]) + "m")
      plt.ylabel('ppbv')
    
    plt.xlabel('Time(Hour)')
    plt.legend(loc = 'best')
    plt.show()
  
  
  def Line(self, z, input_file, height):
  
    # Preprocess
    file_name = os.path.splitext(input_file) # Get extension
    specie_name = file_name[0].split("/")[-1]
  
    # Check data
    if file_name[1] == '.txt':
      df = pd.read_csv(input_file, sep = " ")
    
    else:
      df = pd.read_csv(input_file)
      
    # Calculate Time
    X, Y = df.shape
    day = (X-1)*2/24/(60/5)*60/60
    x_range = np.arange(0, day*24, (day*24/(X-1)))
    x_range = np.append(x_range, [day*24])
    
    # Plot
    plt.plot(x_range, df.iloc[:, int(height)], label = specie_name)
    
    # Put attributes
    if specie_name == "Kh":
      plt.title(specie_name + ' at height ' + str(z[int(height)]) + "m")
      plt.ylabel('m\N{SUPERSCRIPT TWO}s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}')
    
    else:
      plt.title('Concentration of '+ specie_name + ' at height ' + str(z[int(height)]) + "m")
      plt.ylabel('ppbv')
    
    plt.xlabel('Time(Hour)')
    plt.legend(loc = 'best')
    plt.show()
  
  def Contourf_NC(self, input_file, specie, specie_name, pre):
  
    # Get Data
    fh = Dataset(input_file, mode = "r")
    gas = fh.variables['Concentration']
    Kh = fh.variables["Kh"]
    
    # Calculate Time
    X, Y, Z = gas.shape
    day = (X-1)*2/24/(60/5)*60/60
    x_range = np.arange(0, day*24, (day*24/(X-1)))
    x_range = np.append(x_range, [day*24])
    
    # Get Time, Height and Precision
    y, x = np.meshgrid(fh.variables['z_coordinate'], x_range)
    max_val = 0
    contourf_levels = 0
    
    if specie_name == "Kh": 
      max_val = Kh[:, :].max() if Kh[:, :].max() != 0 else 1
      contourf_levels = np.arange(Kh[:, :].min(), max_val + (max_val / 1e10), (max_val - Kh[:, :].min()) / 10)
    
    else:
      max_val = gas[:, :, specie].max() if gas[:, :, specie].max() != 0 else 1
      contourf_levels = np.arange(gas[:, :, specie].min(), max_val + (max_val / 1e10), (max_val - gas[:, :, specie].min()) / pre)
    
    # Plot
    if specie_name == "Kh":
      fig = plt.contourf(x, y, Kh[:, :], contourf_levels)
      
    else:
      fig = plt.contourf(x, y, gas[:, :, specie], contourf_levels)
    
    # Put attributes
    if specie_name == "Kh":
      plt.colorbar(fig, label = "m\N{SUPERSCRIPT TWO}s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}", orientation = "vertical")
    
    else:
      plt.colorbar(fig, label = "ppbv", orientation = "vertical")
    
    plt.title(specie_name)
    plt.xlabel('Time (Hour)')
    plt.ylabel('Height (Meter)')
    plt.show()
  
  
  def Contourf(self, z, input_file, pre):
    
    # Preprocess
    file_name = os.path.splitext(input_file) # Get extension
    specie_name = file_name[0].split("/")[-1]
  
    # Check data
    if file_name[1] == '.txt':
      df = pd.read_csv(input_file, sep = " ")
    
    else:
      df = pd.read_csv(input_file)
      
    # Calculate Time
    X, Y = df.shape
    day = (X-1)*2/24/(60/5)*60/60
    x_range = np.arange(0, day*24, (day*24/(X-1)))
    x_range = np.append(x_range, [day*24])
    
    # Get Time, Height and Precision
    y, x = np.meshgrid(z, x_range)
    max_val = df.to_numpy().max() if df.to_numpy().max() != 0 else 1
    contourf_levels = np.arange(df.to_numpy().min(), max_val + (max_val / 1e10), (max_val - df.to_numpy().min()) / pre)
    
    # Plot
    fig = plt.contourf(x, y, df.iloc[:, :], contourf_levels)
    
    # Put attributes
    if specie_name == "Kh":
      plt.colorbar(fig, label = "m\N{SUPERSCRIPT TWO}s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}", orientation = "vertical")
    
    else:
      plt.colorbar(fig, label = "ppbv", orientation = "vertical")
    
    plt.title(specie_name)
    plt.xlabel('Time (Hour)')
    plt.ylabel('Height (Meter)')
    plt.show()
