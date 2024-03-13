'''
ATMO Driver File
'''
from os import remove, getpid, getcwd
from os.path import join
import pandas as pd
import numpy as np
import netCDF4 as nc

import boxmod as bm

def main(Input):
  """
  Run one scenario using the box model
  
  Parameters:
  -----------
  Input: List
    List of input variables gathered from either GUI or selected files that includes:
    
    location: String
      Path to the directory/folder for the output file
      
    name: String
      Name of the output file
      
    box: Boolean
      Binary switch for choosing to run the model as a box model or not
      
    altitude: Boolean
      Binary switch for running vertical mixing coefficient or not
      
    vimx: Boolean
      Binary switch for running vertical mixing or not
      
    decomp: Boolean
      Binary switch for running surface deposition or not
      
    temp: Boolean
      Binary switch for calculating temperature using a formula or not
      
    Running_days: float
      Used to calculate total running time and iterations
      
    deltim: float
      Time stamp (in seconds)
      
    DL_values: List
      Date & Location values including: 1) Year
                                        2) Month
                                        3) Day
                                        4) Latitude
                                        5) Longitude
                                        6) Standard Longitude
                                        
    all_spc: List
      List containing all specie names
      
    eqn_list: List
      List containing paths to all selected equation files
      
    specie_file: String
      Path to specie file
    
  Return:
  -----------
  1: int
    Used to trigger GUI to update taskbar
  """

  # %% Load mechanism and examine
  #eqns = bm.load.read_csvs(eqns_list, True) # True - specie_only
  
  # Create mechanism (EqnSet)
  #mech = bm.EqnSet(eqns, name="Scenario_Name", long_name="Scenario_Description")
  
  #---------------------------
  # Code for list all species
  #all_spc = mech.spcs # List
  #n_spc = mech.n_spc
  #assert n_spc == len(all_spc)
  
  #print(all_spc)
  #---------------------------

  (location, name,
   box, altitude, vmix, decomp, temp, 
   Running_days, deltim,  
   DL_values, all_spc,
   eqn_list, specie_file) = [Input[i] for i in range(len(Input))]

  I_mix = True if vmix == 'On' else False
  I_deposition = True if decomp == 'On' else False
  I_atk        = True if altitude == 'On' else False
  Box_Model = False if box == 'Off' else True # True to run the program as a box model, False to run the complete program
  Temp      = False if temp == 'Off' else True # True to compute temperature using a function, False to set temperature as a constant
  
  # Set level
  env_level = 248 # TODO: Let user define height levels themselves
  
  if(Box_Model):
    ops_level = 1
  else:
    ops_level = env_level 
  
  # Configure location & date
  normal = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
  abnormal = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    
  iyear = int(DL_values[0]) # Assign year
    
  # Calculate day
  jday = 0
    
  for i in range(int(DL_values[1]) - 1):
    if iyear % 4 == 0:
      jday += abnormal[i]
    else:
      jday += normal[i]
    
  jday += int(DL_values[2])
  
  ntimes = 60 * 60 / deltim * 24 * Running_days
  
  # Clear old output file
  path = join(location, name + '.nc')
  
  try:
    remove(path)
  except:
    pass
  
  # Prepare to write netCDF file  
  ncfile = nc.Dataset(path, mode='w', format='NETCDF4')
  
  # Define netCDF file
  output = bm.netCDF()
  
  # Create var
  timloc = 0
  y = 0  # For output
  updated_eqns_list = []
  
  df = pd.read_csv(specie_file)
  all_spc = list(df.iloc[:, 0])
  init_values = df.iloc[:, 1]
  unit_spc = df.iloc[:, 2]
  henry = df.iloc[:, 3]
  f0 = df.iloc[:, 4]
  
  akh = np.zeros(env_level)
  Vdep = np.zeros(len(all_spc))
  
  dens2con_a = 1e-3 * (1 / 28.97) * 6.022e23
  dens2con_w = 1e-3 * (1 / (2 * 1.0079 + 15.9994)) * 6.022e23
  
  # Initialize all variables
  (dlong, sinlat, coslat, 
   sindec, cosdec, eqtm, 
   t_prof, theta, p_prof, 
   z, dz, dzf, wind, relh, vcp) = bm.init(Temp, iyear, jday, 
                                          DL_values[5], 
                                          DL_values[4], 
                                          DL_values[3], 
                                          env_level, init_values,
                                          all_spc.index("H2O"))
  
  output.define(248, all_spc, z) # TODO: same as above
  
  # Time loop
  for i in range(1, int(ntimes + 2)):
    t0 = timloc * 3600
    t1 = t0 + deltim
  
    if((i % (600 / int(deltim))) == 1):
      parameter = [""] * 8
      
      parameter[0] = y
      parameter[1] = timloc + DL_values[4] / 15
      parameter[2] = Vdep
      parameter[3] = akh
      parameter[4] = t_prof
      parameter[5] = theta
      parameter[6] = p_prof
      parameter[7] = vcp
        
      # Pass value
      output.pass_param(ncfile, parameter)
      # Write to output file
      output.output_netcdf_file()
        
      y = y + 1
        
    # Sun angel
    coszen = bm.zenith(timloc, eqtm, dlong, sinlat, sindec, coslat, cosdec)
  
    source = np.zeros((env_level, len(all_spc)))
    
    # Vertical mixing coefficient
    if(I_atk):
      akh, theta = bm.atk(all_spc.index("H2O"), env_level, theta, t_prof, p_prof, vcp, wind, z)
    
    # Surface deposition
    if(I_deposition):
      Vdep, source = bm.sinks(akh, f0, henry, all_spc.index("O3"), env_level, len(all_spc), relh, source, vcp, dz, z)
        
    # Vertical mixing
    if(I_mix):
      vcp = bm.newc(env_level, akh, dzf, dz, deltim, len(all_spc), vcp, source, all_spc.index("H2O"), all_spc)
        
      source = np.zeros((env_level, len(all_spc)))
        
      vcp = bm.newc(env_level, akh, dzf, dz, deltim, len(all_spc), vcp, source, all_spc.index("H2O"), all_spc) # Need vertical mixing to smooth near surface since extreme high resolution, 0.1m
  
    # For each height level
    for j in range(ops_level):
        
      # Update rate expression
      rho_phy = (p_prof[j] / (t_prof[j] * 287.058))
      C_M = dens2con_a * (p_prof[j] / (t_prof[j] * 287.058))
      C_H2O = dens2con_w * ((vcp[j, all_spc.index("O3")] * 18 / 28.8)) * rho_phy
      updated_eqns_list = bm.kinetic_function_update(eqn_list, t_prof[j], C_M, C_H2O, getpid())
      
      # Load mechanism and examine
      eqns = bm.load.read_csvs(updated_eqns_list, False)
  
      # Create mechanism (EqnSet)
      mech = bm.EqnSet(eqns, name = "", long_name = "")
  
      # Assign specie concentrations to the correct dictionary based on defined unit
      conv = dens2con_a * rho_phy
      oconv = 1 / conv
      
      ppbv = {all_spc[i]: conv * max(vcp[j, i], 0) for i in range(len(all_spc)) if unit_spc[i] == "ppbv"}
      molec_cm3 = {all_spc[i]: conv * max(vcp[j, i], 0) for i in range(len(all_spc)) if unit_spc[i] == "molec_cm3"}
  
      # Run box model
      df = bm.run_exp(mech,
                      c0_ppbv = ppbv,
                      c0_molec_cm3 = molec_cm3,
                      p = p_prof[j], T = t_prof[j],
                      t_tot = deltim,
                      dt = deltim,
                     )
      
      for m in range(len(all_spc)):
        vcp[j, m] = oconv * max(df.iloc[1][m], 0)
  
    timloc = t1 / 3600
      
  ncfile.close()

  # Auto delete updated equation files after finishing computation
  for i in updated_eqns_list:
    remove(i)

  return 1

if __name__ == "__main__":
  
    generate_specie_file = True
    specie_file_name = ""
  
    # Select folder containing all equations
    eqn_file_list = [join(getcwd(), "Equations", "racm.csv")]
    species = bm.load.read_csvs(eqn_file_list, True)
    mech = bm.EqnSet(species, name="Specie_Gen", long_name="Specie_Generation")
    all_spc = mech.spcs # List of all species

    # Generate specie file
    if(generate_specie_file):
      temp = np.full((len(all_spc), 5), "", dtype = object)
      temp[:, 0] = all_spc
      pd.DataFrame(temp, columns = ["Specie", "Initial Value", "Unit", "henry", "f0"]).to_csv(os.path.join(os.getcwd(), specie_file_name), index = False, header = True) 

    # Run program
    else
    #           output folder, output file name
      _ = main([getcwd(),      "test_output",
    #           box model, altitude, vmix, decomp, temperature
                "Off",     "On",     "On", "On",   "Off",
    #           running days, timestamp (in seconds)
                1,            300,
    #           year,    month, day,  latitude, longitude, STD longitude
                [2012.0, 3.0,   24.0, 71.2906,  156.7886,  156.7886],
                all_spc, eqn_file_list, 
    #           specie file
                join(getcwd(), specie_file_name)])
