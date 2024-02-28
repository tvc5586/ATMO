"""
Generate netCDF output file
"""

import numpy as np

__all__ = ("netCDF")

class netCDF():

  specie = None
  time = None
  z_coordinate = None
  Concentration = None
  Kh = None
  Vdep = None
  temperature = None
  potential_temperature = None
  pressure = None

  def define(self, n_level, specie_list, z_coordinate):

    self.nlev = n_level
    self.spec_list = specie_list
    self.nspec = len(specie_list)
    self.z_coordinate = z_coordinate

  def pass_param(self,ncfile,parameter):
    
    self.ncfile = ncfile
    self.parameter = parameter

  def output_netcdf_file(self):

    ncfile = self.ncfile
    
    t = self.parameter[0]

    print(t)

    if t == 0:

      #! Define the dimensions. 

      ncfile.createDimension('level', self.nlev)
      ncfile.createDimension('species', self.nspec)
      ncfile.createDimension('time', None)

      #! Assign create Variables`

      specie = ncfile.createVariable("Specie", str, ('species',))
      time = ncfile.createVariable('time', np.float32, ('time',))
      z_coordinate = ncfile.createVariable('z_coordinate', np.float32, ('level',))
      Concentration = ncfile.createVariable('Concentration',np.float64,('time','level','species'))
      Kh = ncfile.createVariable('Kh',np.float64,('time','level'))
      Vdep = ncfile.createVariable('Vdep',np.float64,('time','species'))
      temperature = ncfile.createVariable('temperature',np.float64,('time','level'))
      potential_temperature = ncfile.createVariable('potential_temperature',np.float64,('time','level'))
      pressure = ncfile.createVariable('pressure',np.float64,('time','level'))

      #! Assign units attributes to the netCDF variables.

      Kh.units = 'm2/s'
      Vdep.description = 'deposition velocity'
      Vdep.units = 'm/s'
      temperature.units = 'K'
      potential_temperature.units = 'K'
      pressure.units = 'Pa'
      Concentration.units = 'molec cm^-3'
      Concentration.SpeciesList = self.spec_list
      time.units = 'UTC'
      z_coordinate.units = 'm'

      ### DATA CREATION ###

      netCDF.specie = ncfile.variables['Specie']
      netCDF.time = ncfile.variables['time']
      netCDF.z_coordinate = ncfile.variables['z_coordinate']
      netCDF.Concentration = ncfile.variables['Concentration']
      netCDF.Kh = ncfile.variables['Kh']
      netCDF.Vdep = ncfile.variables['Vdep']
      netCDF.temperature = ncfile.variables['temperature']
      netCDF.potential_temperature = ncfile.variables['potential_temperature']
      netCDF.pressure = ncfile.variables['pressure']
      
      netCDF.z_coordinate[:] = self.z_coordinate
      netCDF.specie[:] = np.array(self.spec_list, dtype='object')

    ### DATA WRITING ###

    netCDF.time[t] = self.parameter[1]
    netCDF.Vdep[t] = self.parameter[2]
    netCDF.Kh[t] = self.parameter[3]
    netCDF.temperature[t] = self.parameter[4]
    netCDF.potential_temperature[t] = self.parameter[5]
    netCDF.pressure[t] = self.parameter[6]
    netCDF.Concentration[t] = self.parameter[7]
