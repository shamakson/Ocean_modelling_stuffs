#=======================================================================#
#=======================================================================#
# Create MPI-OM inputs from OMIP input files prepared by                #
# Frank Roeske and the German OMIP Group.                               #
# My input variable are from HIPPO-simulations                          #
#=======================================================================#
# Author: Eric Samakinwa                                                #
#=======================================================================#

from netCDF4 import Dataset
import numpy as np
from numpy import dtype

#========================================================================#
# read netCDF file                                                       #
#========================================================================#
# open a netCDF file to read   OMIP       	 			 #
#========================================================================#

filename = "OMIP_input/2m_dewpoint_temp_arctic_corr.nc"

#========================================================================#
# Read variable from HIPPO-simulations                                   #
#========================================================================#
filename2 = '/scratch3/eric/HIPPO_r001/1780/2m_dewpoint_temp_arctic_corr.nc'

ncin = Dataset(filename, 'r', format='NETCDF4')
ncin2 = Dataset(filename2, 'r', format='NETCDF4')

# Read input data
vin = ncin2.variables['dpt_temp_2m']

#------------------
# write netCDF file
#------------------

# open a netCDF file to write
ncout = Dataset('/scratch3/eric/prep_mpidata/1780.2m_dewpoint_temp_arctic_corr.nc', 'w', format='NETCDF4')

#========================================================================#
# Copy longitude axis                                                    #
#========================================================================#

lonpt_in = ncin.variables['lonpt']
lon_edge_in = ncin.variables['lon_edge']

latpt_in = ncin.variables['latpt']
lat_edge_in = ncin.variables['lat_edge']

Time_in = ncin.variables['Time']

# get length of axis data
#ntime = len(tin)
nlat = len(latpt_in)
nlon = len(lonpt_in)

nlat2 = len(lat_edge_in)
nlon2 = len(lon_edge_in)


# define axis size
ncout.createDimension('lonpt', nlon)
ncout.createDimension('lon_edge', nlon2)

ncout.createDimension('latpt', nlat)
ncout.createDimension('lat_edge', nlat2)

ncout.createDimension('Time', None)  # unlimited

##########################################################################
# create longitude axis
#========================================================================#
lonpt = ncout.createVariable('lonpt', dtype('float32').char, ('lonpt'))
#lonpt.standard_name = 'longitude'
lonpt.long_name = 'longitude at Points'
lonpt.units = 'degrees'

#========================================================================#
# create longitude axis2                                                 #
#========================================================================#
lon_edge = ncout.createVariable('lon_edge', dtype('float32').char, ('lon_edge'))
#lon_edge.standard_name = 'longitude_edge'
lon_edge.long_name = 'Longitude at Edges of Cells'
lon_edge.units = 'degrees'

#=======================================================================#
# create latitude axis                                                  #
#=======================================================================#
latpt = ncout.createVariable('latpt', dtype('float32').char, ('latpt'))
#latpt.standard_name = 'latitude'
latpt.long_name = 'Latitude at Points'
latpt.units = 'degrees'

#=======================================================================#
# create latitude axis2                                                 #
#=======================================================================#
lat_edge = ncout.createVariable('lat_edge', dtype('float32').char, ('lat_edge'))
lat_edge.standard_name = 'latitude_edge'
lat_edge.long_name = 'Latitude at Edges of Cells'
lat_edge.units = 'degrees'

#=======================================================================#
# create time axis                                                      #
#=======================================================================#
Time = ncout.createVariable('Time', dtype('float32').char, ('Time',))
Time.long_name = 'Year Day'
Time.units = 'days'

#=======================================================================#
# create variable array                                                 #
#=======================================================================#
vout = ncout.createVariable('dpt_temp_2m', dtype('float32').char, ('Time', 'latpt', 'lonpt'))
vout.long_name = '2m Dewpoint Temperature'
vout.units = 'K'
vout.valid_range=np.array((150. , 350. ))

#=======================================================================#
# copy axis from original dataset                                       #
#=======================================================================#
lonpt[:] = lonpt_in[:]
lon_edge[:] = lon_edge_in[:]

latpt[:] = latpt_in[:]
lat_edge[:] = lat_edge_in[:]

Time[:] = Time_in[:]

vout[:] = vin[:]
#=======================================================================#
# close files                                                           #
#=======================================================================#
ncin.close()
ncout.close()

exit()
