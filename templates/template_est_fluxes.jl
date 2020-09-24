### =======================================================================
### = template_est_fluxes.jl
### = Alex Turner
### = 03/09/2020
### =----------------------------------------------------------------------
### = NOTES:
### =  ( 1) Runs the real inversion.
### =======================================================================

### Offline plotting
using PyCall
pyimport("matplotlib").use("Agg")
pyplot = pyimport("matplotlib.pyplot")

### Libraries
using Distributed
@everywhere using Printf
@everywhere using NetCDF
@everywhere using CSV
@everywhere using Dates
@everywhere using DataFrames
using PyPlot
plt = PyPlot

### Load my functions
include("./Util/plot_funcs.jl")
include("./Util/read_funcs.jl")
include("./Util/get_data_funcs.jl")
include("./Util/reshape_funcs.jl")
include("./Util/solve_inv_funcs.jl")


### ========================================
### Parameters
### ========================================

### Globals needed for other files
global origin     = Dates.DateTime(1970,1,1,0,0,0)
global IntType    = Int64   # Precision we are using for integers
global FltType    = Float32 # Precision we are using for floats
global lowBound   = 1e-5    # Smallest number to include in our correlation matrices
global diag_prior = false   # Are we using a diagonal prior?
global minUncert  = 1.0     # Lower bound on the prior uncertainty (umol/m2/s)
global minObsErr  = 1.0     # Lower bound on the measurement uncertainty (ppm)
global obsFreq    = 60.0    # Observations per hour

### Time period for the inversion
startTime  = Dates.DateTime(YearYear,MonthMonth,DayDay) # Start of inversion period
invWindow  = Dates.Day(WindowWindow)               # Time window for inversion
nHr        = BackHoursBackHours                         # Number of backhours
UTC_to_PST = Dates.Hour(-8)

### Different grids to use
# Full grid used for footprints
full_xLim = [ -125.0, -120.0 ]
full_yLim = [   36.0,   40.0 ]
big_xLim,big_yLim = full_xLim,full_yLim
# Medium sized grid
medium_xLim = [-123.60,-121.60]
medium_yLim = [  36.80,  38.60]
med_xLim,med_yLim = medium_xLim,medium_yLim
# Bay Area domain (smallest grid)
BayArea_xLim = [-123.10,-121.80]
BayArea_yLim = [  37.35,  38.40]
small_xLim,small_yLim = BayArea_xLim,BayArea_yLim
# Inversion grid to use
Inv_lonLim = small_xLim
Inv_latLim = small_yLim

### Specify the uncertainties
ems_uncert =  50/100	# %
# Model error at different hours of the day
mod_err = [[00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23]  # hour
           [ 3  3  3  3  3  3  4  5  8  6  4  2  1  1  1  1  2  4  6  8  5  4  3  3]] # error (ppm)

### Information for the k-fold cross validation
cross_valid = CrossValidCrossValid # Set to true for cross validation
k_fold      = kFoldkFold  # Percentage to leave out
k_ind       = kIndkInd  # Which index are we running?
k_fold_dat  = (cross_valid, k_fold, k_ind)
# Set k_ind to zero if we are not doing cross_validation
k_ind = cross_valid ? k_ind : 0

### Save data?
save_plot = false
save_ems  = true

### Store the inversion data
startTime  = startTime - UTC_to_PST # Convert from PST to UTC
endTime    = startTime + invWindow
Inv_timLim = [ startTime - Dates.Hour(nHr), endTime + Dates.Hour(nHr) ]

### Make output folders?
ispath("./stored_data")   ? nothing : mkdir("./stored_data")
ispath("./output")        ? nothing : mkdir("./output")
ispath("./output/ems")    ? nothing : mkdir("./output/ems")
ispath("./output/obs")    ? nothing : mkdir("./output/obs")
ispath("./output/images") ? nothing : mkdir("./output/images")


### ========================================
### Get the data 
### ========================================

### Read the observation information
println("* READING OBSERVATION AND FOOTPRINT DATA")
(latF,lonF,obsInfo,foot) = get_footprint_data(Inv_latLim,Inv_lonLim,Inv_timLim)

### Read the emission information
println("* READING EMISSIONS DATA")
(latE,lonE,ems,emsTime) = get_emission_data(Inv_latLim,Inv_lonLim,Inv_timLim)
#(latE,lonE,ems,emsTime) = get_emission_data_2x(Inv_latLim,Inv_lonLim,Inv_timLim)

### Diagnostic
println("* ALL DATA READ")

### Make sure the lats/lons match
if (sum(abs.(lonE - lonF)) > 1e-3) || (sum(abs.(latE - latF)) > 1e-3)
   error("LATS AND LONS DON'T MATCH!")
else
   lon, lat = lonE, latE
end


### ========================================
### Estimate the emissions
### ========================================

### Remove the extra processors
rmprocs([2:24])

### Construct the input tuples
inv_data   = (lon,lat,Inv_latLim,Inv_lonLim,startTime,endTime,ems_uncert,mod_err,k_fold_dat)
foot_data  = (obsInfo,foot,nHr)
ems_data   = (ems,emsTime)
input_data = (inv_data,foot_data,ems_data)

### Estimate emissions
println("* ESTIMATING EMISSIONS")
@time (lon_red,lat_red,ems_Times,x_pri,x_hat,obs_compare,dofs) = estimate_ems(input_data)


### ========================================
### Store the solution
### ========================================

### Save emissions?
if save_ems

   ### Estimated emissions
   println("* STORING EMISSIONS")
   outDir = cross_valid ? @sprintf("./output/ems_%04i",k_ind) : "./output/ems"
   save_emissions(outDir,x_hat,lon_red,lat_red,ems_Times)
   #save_emissions_2x(outDir,x_hat,lon_red,lat_red,ems_Times)

   ### Estimated emissions
   println("* STORING OBSERVATIONS")
   dateName   = Dates.format(startTime+UTC_to_PST,Dates.DateFormat("yyyymmdd"))
   if cross_valid
      obsCSVname = @sprintf("./output/obs/%s-%04iSuffixSuffix.csv",dateName,k_ind)
   else
      obsCSVname = @sprintf("./output/obs/%sSuffixSuffix.csv",dateName)
   end
   save_observations(obsCSVname,obs_compare)

end


### ========================================
### Plot the solution
### ========================================

### Save plots?
if save_plot

   ### Get the average emissions
   avg_pri = sum(x_pri,dims=3)[:,:,1]/size(x_pri)[3]
   avg_hat = sum(x_hat,dims=3)[:,:,1]/size(x_hat)[3]
   avg_dif = sum(x_hat - x_pri,dims=3)[:,:,1]/size(x_hat)[3]

   ### Plot the average emissions
   println("* PLOTTING")
   zLims  = [-10,30]
   dLims  = [-10,10]
   zLabel = L"CO$_2$ Flux $(\mu\mathrm{mol}\, \mathrm{s}^{-1}\, \mathrm{m}^{-2})$"
   dLabel = L"$\Delta$CO$_2$ Flux $(\mu\mathrm{mol}\, \mathrm{s}^{-1}\, \mathrm{m}^{-2})$"
   # Plot name info
   outDir   = "./output/images"
   dateName = Dates.format(startTime+UTC_to_PST,Dates.DateFormat("yyyymmdd"))
   dateName = cross_valid ? @sprintf("%s-%04i",dateName,k_ind) : dateName
   dateName = @sprintf("%sSuffixSuffix",dateName)
   # Prior
   plt.figure(figsize=(5,5))
   plotEMS_neg(lon_red,lat_red,transpose(avg_pri),BayArea_xLim,BayArea_yLim,zLims,zLabel)
   plt.title("Prior Fluxes")
   plt.savefig(@sprintf("%s/AvgPrior_%s.png",outDir,dateName))
   # Posterior
   plt.figure(figsize=(5,5))
   plotEMS_neg(lon_red,lat_red,transpose(avg_hat),BayArea_xLim,BayArea_yLim,zLims,zLabel)
   plt.title("Posterior Fluxes")
   plt.savefig(@sprintf("%s/AvgPosterior_%s.png",outDir,dateName))

   ### Plot the difference
   # Posterior minus Prior
   plt.figure(figsize=(5,5))
   plotEMS_diff(lon_red,lat_red,transpose(avg_dif),BayArea_xLim,BayArea_yLim,dLims,dLabel)
   plt.title("Difference")
   plt.savefig(@sprintf("%s/AvgDifference_%s.png",outDir,dateName))

   ### All three together
   plt.figure(figsize=(15,5))
   plt.subplot(1,3,1);plotEMS_neg(lon_red,lat_red,transpose(avg_pri),BayArea_xLim,BayArea_yLim,zLims,zLabel);plt.title("Prior Fluxes")
   plt.subplot(1,3,2);plotEMS_neg(lon_red,lat_red,transpose(avg_hat),BayArea_xLim,BayArea_yLim,zLims,zLabel);plt.title("Posterior Fluxes")
   plt.subplot(1,3,3);plotEMS_diff(lon_red,lat_red,transpose(avg_dif),BayArea_xLim,BayArea_yLim,dLims,dLabel);plt.title("Difference")
   plt.savefig(@sprintf("%s/SummaryFigure_%s.png",outDir,dateName))

end

println("=======================================================================")
println("=                            E   N   D                                =")
println("=======================================================================")


### =======================================================================
### =                            E   N   D                                =
### =======================================================================
