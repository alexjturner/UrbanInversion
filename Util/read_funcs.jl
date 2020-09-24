### =======================================================================
### = read_funcs.jl
### = Alex Turner
### = 03/09/2020
### =----------------------------------------------------------------------
### = NOTES:
### =  ( 1) Functions for reading data
### =  ( 2) All functions take a single input to facilitate parallel reads.
### =----------------------------------------------------------------------
### = INPUTS:
### =  (  ) N/A
### =----------------------------------------------------------------------
### = OUTPUTS:
### =  (  ) N/A
### =======================================================================

### Libraries
using NetCDF
using Dates
using SparseArrays


### ========================================
### Functions to read the netCDF files
### ========================================

### Make functions to read in parallel with pmap
@everywhere function read_lat(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"obs_lat")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
@everywhere function read_lon(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"obs_lon")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
@everywhere function read_agl(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"obs_agl")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
@everywhere function read_co2(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"co2")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
@everywhere function read_co2_err(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"co2_err")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
@everywhere function read_obs_time(file_name)
   try
      ncid = NetCDF.open(file_name)
      out  = Dates.DateTime(
        NetCDF.readvar(ncid,"yr")[1],
        NetCDF.readvar(ncid,"mon")[1],
        NetCDF.readvar(ncid,"day")[1],
        NetCDF.readvar(ncid,"hr")[1])
      NetCDF.close(ncid)
      out
   catch
      Dates.DateTime(0,1,1)
   end
end
# NOAA Background
@everywhere function read_bkg_co2_NOAA(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"bkg_co2_NOAA")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
@everywhere function read_bkg_err_NOAA(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"bkg_err_NOAA")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
# NASA Background
@everywhere function read_bkg_co2_NASA(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"bkg_co2_NASA")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
@everywhere function read_bkg_err_NASA(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"bkg_err_NASA")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
# AmeriFlux Background
@everywhere function read_bkg_co2_AmeriFlux(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"ameriflux_co2")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
@everywhere function read_bkg_err_AmeriFlux(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"ameriflux_err")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
# Trajectory end location (AmeriFlux)
@everywhere function read_amf_lat(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"ameriflux_lat")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
@everywhere function read_amf_lon(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"ameriflux_lon")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
@everywhere function read_amf_agl(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"ameriflux_agl")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
@everywhere function read_amf_time(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = Dates.unix2datetime(NetCDF.readvar(ncid,"ameriflux_julian")[1])
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
# Trajectory end location
@everywhere function read_end_lat(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"end_lat")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
@everywhere function read_end_lon(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"end_lon")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
@everywhere function read_end_agl(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"end_agl")[1]
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end


### ========================================
### Functions to read the footprints
### ========================================

# Function to read the time points for the footprint
@everywhere function read_ftime(file_name)
   try
      ncid = NetCDF.open(file_name)
      out  = map(x -> Dates.unix2datetime(x),
        NetCDF.readvar(ncid,"time"))
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
# Define the function to read the footprints
@everywhere function read_footprint(footDat)
   (file_name, lon_ind, lat_ind) = footDat
   try
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"foot",start=[findfirst(lon_ind),findfirst(lat_ind),1],count=[sum(lon_ind),sum(lat_ind),-1])
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
# Define the function to read the footprints into sparse arrays
@everywhere function read_sparse_footprint(footDat)
   (file_name, lon_ind, lat_ind) = footDat
   try
      ncid  = NetCDF.open(file_name)
      foot  = NetCDF.readvar(ncid,"foot",start=[findfirst(lon_ind),findfirst(lat_ind),1],count=[sum(lon_ind),sum(lat_ind),-1])
      ftime = map(x -> Dates.unix2datetime(x),NetCDF.readvar(ncid,"time"))
      NetCDF.close(ncid)
      #sparse(squeeze(sum(foot,3),3))
      iS    = Array{Any,1}(undef,size(foot)[3])
      jS    = Array{Any,1}(undef,size(foot)[3])
      vS    = Array{Any,1}(undef,size(foot)[3])
      for i = 1:size(foot)[3]
         (iS[i],jS[i],vS[i]) = findnz(sparse(foot[:,:,i]))
      end
      [ftime,iS,jS,vS]
   catch
      NaN
   end
end


### ========================================
### Functions to read the emissions
### ========================================

# Function to read emissions
@everywhere function read_ems_time(file_name)
   try 
      ncid = NetCDF.open(file_name)
      out  = Dates.DateTime(
        NetCDF.readvar(ncid,"yr")[1],
        NetCDF.readvar(ncid,"mon")[1],
        NetCDF.readvar(ncid,"day")[1],
        NetCDF.readvar(ncid,"hr")[1])
      NetCDF.close(ncid)
      out 
   catch
      Dates.DateTime(0,1,1)
   end 
end
@everywhere function read_emissions(emsDat)
   (file_name, lon_ind, lat_ind) = emsDat
   try
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"flx_total",start=[findfirst(lon_ind),findfirst(lat_ind)],count=[sum(lon_ind),sum(lat_ind)])
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end
@everywhere function read_emissions_bio(emsDat)
   (file_name, lon_ind, lat_ind) = emsDat
   try
      ncid = NetCDF.open(file_name)
      out  = NetCDF.readvar(ncid,"flx_bio",start=[findfirst(lon_ind),findfirst(lat_ind)],count=[sum(lon_ind),sum(lat_ind)])
      NetCDF.close(ncid)
      out
   catch
      NaN
   end
end


### =======================================================================
### =                            E   N   D                                =
### =======================================================================
