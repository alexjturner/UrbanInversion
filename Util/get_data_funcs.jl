### =======================================================================
### = get_data_funcs.jl
### = Alex Turner
### = 03/09/2020
### =----------------------------------------------------------------------
### = NOTES:
### =  ( 1) Functions for reading data for the inversion
### =----------------------------------------------------------------------
### = SUBFUNCTIONS:
### =  ( 1) get_footprint_data :: Reads footprint data.
### =  ( 2) get_emission_data  :: Reads emission data.
### =  ( 3) save_emissions     :: Saves emissions to netcdf files.
### =  ( 3) save_observations  :: Saves the resulting concentrations.
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
using Printf
using DataFrames

### Load my own functions
include("./read_funcs.jl")
include("./reshape_funcs.jl")

### Globals
global origin
global IntType
global FltType


### ========================================
### Read the footprint data
### ========================================

### Read the data in
function get_footprint_data(lat_range,lon_range,time_range)
   # Get gridinfo from the first file
   searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
   path  = "./obs/"
   key   = "obs_"
   flist = searchdir(path,key)
   nObs  = length(flist)
   # Create a vector of file names
   obs_names = Array{String,1}(undef, nObs)
   for i in 1:nObs
      obs_names[i] = "$path$(flist[i])"
   end
   ncid  = NetCDF.open(obs_names[1])
   lat   = NetCDF.readvar(ncid,"lat")
   lon   = NetCDF.readvar(ncid,"lon")
   NetCDF.close(ncid)
   # Get info about the observations
   obsLat      = pmap(read_lat,obs_names)
   obsLon      = pmap(read_lon,obs_names)
   obsAGL      = pmap(read_agl,obs_names)
   obsCO2      = pmap(read_co2,obs_names)
   obsErr      = pmap(read_co2_err,obs_names)
   bkgCO2_NOAA = pmap(read_bkg_co2_NOAA,obs_names)
   bkgErr_NOAA = pmap(read_bkg_err_NOAA,obs_names)
   bkgCO2_NASA = pmap(read_bkg_co2_NASA,obs_names)
   bkgErr_NASA = pmap(read_bkg_err_NASA,obs_names)
   bkgCO2_AMFX = pmap(read_bkg_co2_AmeriFlux,obs_names)
   bkgErr_AMFX = pmap(read_bkg_err_AmeriFlux,obs_names)
   endLat      = pmap(read_end_lat,obs_names)
   endLon      = pmap(read_end_lon,obs_names)
   endAGL      = pmap(read_end_agl,obs_names)
   amfLat      = pmap(read_amf_lat,obs_names)
   amfLon      = pmap(read_amf_lon,obs_names)
   amfAGL      = pmap(read_amf_agl,obs_names)
   amfTime     = pmap(read_amf_time,obs_names)
   obsTime     = pmap(read_obs_time,obs_names)
   # Check for missing values in the background data
   bkgCO2_NOAA[bkgCO2_NOAA .<= -999] .= NaN
   bkgErr_NOAA[bkgErr_NOAA .<= -999] .= NaN
   bkgCO2_NASA[bkgCO2_NASA .<= -999] .= NaN
   bkgErr_NASA[bkgErr_NASA .<= -999] .= NaN
   bkgCO2_AMFX[bkgCO2_AMFX .<= -999] .= NaN
   bkgErr_AMFX[bkgErr_AMFX .<= -999] .= NaN
   # Trim the observations down if they're out of our spatio-temporal domain
   obsInd = map(Bool,ones(size(obsLat)))
   for i in 1:length(obsLat)
      obsInd[i] = lat_range[1]  <= obsLat[i]  && obsLat[i]  <= lat_range[2] &&
                  lon_range[1]  <= obsLon[i]  && obsLon[i]  <= lon_range[2] &&
                  time_range[1] <= obsTime[i] && obsTime[i] <  time_range[2]
   end
   obs_names   = obs_names[obsInd]
   obsLat      = obsLat[obsInd]
   obsLon      = obsLon[obsInd]
   obsAGL      = obsAGL[obsInd]
   obsCO2      = obsCO2[obsInd]
   obsErr      = obsErr[obsInd]
   bkgCO2_NOAA = bkgCO2_NOAA[obsInd]
   bkgErr_NOAA = bkgErr_NOAA[obsInd]
   bkgCO2_NASA = bkgCO2_NASA[obsInd]
   bkgErr_NASA = bkgErr_NASA[obsInd]
   bkgCO2_AMFX = bkgCO2_AMFX[obsInd]
   bkgErr_AMFX = bkgErr_AMFX[obsInd]
   endLon      = endLon[obsInd]
   endLat      = endLat[obsInd]
   endAGL      = endAGL[obsInd]
   amfLon      = amfLon[obsInd]
   amfLat      = amfLat[obsInd]
   amfAGL      = amfAGL[obsInd]
   amfTime     = amfTime[obsInd]
   obsTime     = obsTime[obsInd]
   trajDat     = (endLon,endLat,endAGL,amfLon,amfLat,amfAGL,amfTime)
   obsInfo     = (obsTime,obsLat,obsLon,obsAGL,obsCO2,obsErr,bkgCO2_NOAA,bkgErr_NOAA,bkgCO2_NASA,bkgErr_NASA,bkgCO2_AMFX,bkgErr_AMFX,trajDat)
   # Get the lat and lon indicies
   lat_ind = map(Bool,ones(size(lat)))
   lon_ind = map(Bool,ones(size(lon)))
   for i in 1:length(lat)
      lat_ind[i] = lat_range[1] <= lat[i] && lat[i] <= lat_range[2]
   end
   for i in 1:length(lon)
      lon_ind[i] = lon_range[1] <= lon[i] && lon[i] <= lon_range[2]
   end
   lat = lat[lat_ind]
   lon = lon[lon_ind]
   # Get the footprints
   footDat = [(n,lon_ind,lat_ind) for n in obs_names]
   footprint = pmap(read_sparse_footprint,footDat)
   # Return data
   return (lat,lon,obsInfo,footprint)
end


### ========================================
### Read the emission data
### ========================================

### Check if we've read it before
# Standard
function get_emission_data(lat_range,lon_range,time_range)
   # Get gridinfo from the first file
   searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
   path  = "./ems/"
   key   = "BEACON_"
   flist = searchdir(path,key)
   nEms  = length(flist)
   # Create a vector of file names
   ems_names = Array{String,1}(undef, nEms)
   for i in 1:nEms
      ems_names[i] = "$path$(flist[i])"
   end 
   # Get timestamps for the emissions 
   emsTime = pmap(read_ems_time,ems_names)
   # Trim the observations down if they're out of our spatio-temporal domain
   emsInd = map(Bool,ones(size(emsTime)))
   for i in 1:length(emsTime)
      emsInd[i] = time_range[1] <= emsTime[i] && emsTime[i] <= time_range[2]
   end
   ems_names = ems_names[emsInd]
   emsTime   = emsTime[emsInd]
   # Read the lats and lons
   ncid = NetCDF.open(ems_names[1])
   lat  = NetCDF.readvar(ncid,"lat")
   lon  = NetCDF.readvar(ncid,"lon")
   NetCDF.close(ncid)
   # Get the lat and lon indicies
   lat_ind = map(Bool,ones(size(lat)))
   lon_ind = map(Bool,ones(size(lon)))
   for i in 1:length(lat)
      lat_ind[i] = lat_range[1] <= lat[i] && lat[i] <= lat_range[2]
   end
   for i in 1:length(lon)
      lon_ind[i] = lon_range[1] <= lon[i] && lon[i] <= lon_range[2]
   end
   lat = lat[lat_ind]
   lon = lon[lon_ind]
   # Read the emissions
   emsDat = [(n,lon_ind,lat_ind) for n in ems_names]
   emsF   = pmap(read_emissions,emsDat)
   # Put all the data into a single array
   ems = Array{FltType,3}(undef,(size(emsF[1])[1],size(emsF[1])[2],length(emsF)))
   for i in 1:size(ems)[3]
      ems[:,:,i] = emsF[i]
   end
   # Return data
   return (lat,lon,ems,emsTime)
end


### ========================================
### Write the emissions
### ========================================

### Function to save the emissions
function save_emissions(outDir,ems,lon,lat,emsTimes)
   # Make the directory if it doesn't exist
   ispath(outDir) ? nothing : mkdir(outDir)
   # Define the variables we're saving
   dFormat = Dates.DateFormat("mm/dd/yyyy HH:MM:SS")
   londim  = NetCDF.NcDim("lon",  size(ems)[1])
   latdim  = NetCDF.NcDim("lat",  size(ems)[2])
   infdim  = NetCDF.NcDim("info",            1)
   v1      = NetCDF.NcVar("lon",               [londim],        compress=9)
   v2      = NetCDF.NcVar("lat",               [latdim],        compress=9)
   v3      = NetCDF.NcVar("yr",                [infdim],        compress=9)
   v4      = NetCDF.NcVar("mon",               [infdim],        compress=9)
   v5      = NetCDF.NcVar("day",               [infdim],        compress=9)
   v6      = NetCDF.NcVar("hr",                [infdim],        compress=9)
   v7      = NetCDF.NcVar("ems_posterior",     [londim,latdim], compress=9)
   v8      = NetCDF.NcVar("flx_posterior",     [londim,latdim], compress=9)
   v9      = NetCDF.NcVar("ems_prior",         [londim,latdim], compress=9)
   v10     = NetCDF.NcVar("flx_total_prior",   [londim,latdim], compress=9)
   v11     = NetCDF.NcVar("flx_traffic_prior", [londim,latdim], compress=9)
   v12     = NetCDF.NcVar("flx_point_prior",   [londim,latdim], compress=9)
   v13     = NetCDF.NcVar("flx_bio_prior",     [londim,latdim], compress=9)
   # Loop through the timesteps
   for i = 1:size(ems)[3]
      # Define the filename and info for this timestep
      filename = @sprintf( "%s/ems_%s.nc", outDir, Dates.format( emsTimes[i], Dates.DateFormat( "yyyyxmmxddxHH" ) ) )
      priname  = @sprintf( "./ems/BEACON_%s.ncdf", Dates.format( emsTimes[i], Dates.DateFormat( "yyyyxmmxddxHH" ) ) )
      fInfo    = Dict("emission_date" => @sprintf("%s UTC",Dates.format(emsTimes[i],dFormat)),
                      "creation_date" => Dates.format(Dates.now(),dFormat) )
      # Read the prior
      lon_prior   = NetCDF.ncread(priname, "lon")
      lat_prior   = NetCDF.ncread(priname, "lat")
      emsPrior    = NetCDF.ncread(priname, "ems_total")
      flx_total   = NetCDF.ncread(priname, "flx_total")
      flx_traffic = NetCDF.ncread(priname, "flx_traffic")
      flx_point   = NetCDF.ncread(priname, "flx_point")
      flx_bio     = NetCDF.ncread(priname, "flx_bio")
      # Trim the prior to our domain
      lon_ind     = (lon[1] .<= lon_prior) .& (lon_prior .<= lon[end])
      lat_ind     = (lat[1] .<= lat_prior) .& (lat_prior .<= lat[end])
      emsPrior    = emsPrior[lon_ind,lat_ind]
      flx_total   = flx_total[lon_ind,lat_ind]
      flx_traffic = flx_traffic[lon_ind,lat_ind]
      flx_point   = flx_point[lon_ind,lat_ind]
      flx_bio     = flx_bio[lon_ind,lat_ind]
      flxRatio    = emsPrior ./ flx_total
      for j = 1:size(flxRatio)[2]
         iNaN              = isnan.(flxRatio[:,j])
         flxRatio[iNaN,j] .= mean(flxRatio[.!iNaN,j])
      end
      # Save the data
      #@printf("  * Saving: %s\n",filename)
      ncid = NetCDF.create(filename,[v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13],mode=NC_NETCDF4)
      NetCDF.putvar(ncid, "lon",               convert.(Float32,lon)      )
      NetCDF.putvar(ncid, "lat",               convert.(Float32,lat)      )
      NetCDF.putvar(ncid, "yr",                [Dates.year(emsTimes[i])]  )
      NetCDF.putvar(ncid, "mon",               [Dates.month(emsTimes[i])] )
      NetCDF.putvar(ncid, "day",               [Dates.day(emsTimes[i])]   )
      NetCDF.putvar(ncid, "hr",                [Dates.hour(emsTimes[i])]  )
      NetCDF.putvar(ncid, "ems_posterior",     convert.(Float32,ems[:,:,i].*flxRatio) )
      NetCDF.putvar(ncid, "flx_posterior",     convert.(Float32,ems[:,:,i])           )
      NetCDF.putvar(ncid, "ems_prior",         convert.(Float32,emsPrior)             )
      NetCDF.putvar(ncid, "flx_total_prior",   convert.(Float32,flx_total)            )
      NetCDF.putvar(ncid, "flx_traffic_prior", convert.(Float32,flx_traffic)          )
      NetCDF.putvar(ncid, "flx_point_prior",   convert.(Float32,flx_point)            )
      NetCDF.putvar(ncid, "flx_bio_prior",     convert.(Float32,flx_bio)              )
      NetCDF.putatt(ncid, "flx_bio_prior",     Dict("Unit"=>"umol/m2/s") )
      NetCDF.putatt(ncid, "ems_posterior",     Dict("Unit"=>"tC/hr")     )
      NetCDF.putatt(ncid, "flx_posterior",     Dict("Unit"=>"umol/m2/s") )
      NetCDF.putatt(ncid, "ems_prior",         Dict("Unit"=>"tC/hr")     )
      NetCDF.putatt(ncid, "flx_total_prior",   Dict("Unit"=>"umol/m2/s") )
      NetCDF.putatt(ncid, "flx_traffic_prior", Dict("Unit"=>"umol/m2/s") )
      NetCDF.putatt(ncid, "flx_point_prior",   Dict("Unit"=>"umol/m2/s") )
      NetCDF.putatt(ncid, "attributes",        fInfo )
      NetCDF.close( ncid )
   end
end


### ========================================
### Write the observations
### ========================================

### Function to save the observations
function save_observations(obsCSVname,obs_compare)
   # Load data from the tuple
   (obs_time,obs_lon,obs_lat,obs_agl,obs_co2,prior_co2,posterior_co2,bkg_co2,obs_co2_err,model_err,bkg_co2_err,iTrain) = obs_compare
   # Define the variables we're saving
   df = DataFrame([:iTrain        => iTrain,                                 # Training or eval? 
                   :obs_time      => Dates.value.(obs_time .- origin)/1000,  # Seconds since origin
                   :obs_lon       => obs_lon,                                # Degrees
                   :obs_lat       => obs_lat,                                # Degrees
                   :obs_agl       => obs_agl,                                # mAGL
                   :obs_co2       => obs_co2,                                # ppm
                   :prior_co2     => prior_co2,                              # ppm
                   :posterior_co2 => posterior_co2,                          # ppm
                   :bkg_co2       => bkg_co2,                                # ppm
                   :obs_co2_err   => obs_co2_err,                            # ppm
                   :model_err     => model_err,                              # ppm
                   :bkg_co2_err   => bkg_co2_err])                           # ppm
   CSV.write(obsCSVname, df, writeheader=true)
end


### =======================================================================
### =                            E   N   D                                =
### =======================================================================
