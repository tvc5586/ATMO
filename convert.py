from netCDF4 import Dataset
import numpy as np
import pandas as pd
import xarray as xr
import os

class Convert():

  def To_CSV(in_location, out_location, specie_list):
    
    # Create Folder
    try:
      os.mkdir(out_location + "/CSV Files")
    except:
      pass
  
    # Get Data
    fh = Dataset(in_location, mode = "r")
    gas = fh.variables['Concentration']
    Kh = fh.variables["Kh"]
    
    ds = xr.DataArray(np.asarray(gas)) # Convert netCDF variable to xarray DataArray
    specie_dim = list(list(dict(ds.sizes).items())[-1])[-1] # # Convert dimensions of DataArray to dictionary, then access value of last dimension
  
    yield (specie_dim + 1)
  
    for i in range(specie_dim):
      df = pd.DataFrame(np.asarray(ds[:, :, i])) # Convert a slice of DataArray to numpy array, then to DataFrame
      df.to_csv(out_location + "/CSV Files/" + specie_list[i] + ".csv", index = False) # Save DataFrame as csv
      yield (i + 1)
  
    df = pd.DataFrame(np.asarray(Kh)) # Convert Kh to numpy array, then to DataFrame
    df.to_csv(out_location + "/CSV Files/Kh.csv", index = False) # Save DataFrame as csv
    yield (specie_dim + 1)
  
  def To_Text(in_location, out_location, specie_list):
  
    # Create Folder
    try:
      os.mkdir(out_location + "/Text Files")
    except:
      pass
      
    # Get Data
    fh = Dataset(in_location, mode = "r")
    gas = fh.variables['Concentration']
    Kh = fh.variables["Kh"]
    
    ds = xr.DataArray(np.asarray(gas)) # Convert netCDF variable to xarray DataArray
    specie_dim = list(list(dict(ds.sizes).items())[-1])[-1] # # Convert dimensions of DataArray to dictionary, then access value of last dimension
    
    yield (specie_dim + 1)
    
    for i in range(specie_dim):
      np.savetxt(out_location + "/Text Files/" + specie_list[i] + ".txt", np.asarray(ds[:, :, i])) # Convert a slice of DataArray to numpy array, then save numpy array as txt
      yield (i + 1)
      
    np.savetxt(out_location + "/Text Files/Kh.txt", np.asarray(Kh)) # Convert Kh to numpy array, then save numpy array as txt
    yield (specie_dim + 1)
