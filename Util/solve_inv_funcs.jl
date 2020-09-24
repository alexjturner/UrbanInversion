### =======================================================================
### = solve_inv_funcs.jl
### = Alex Turner
### = 03/09/2020
### =----------------------------------------------------------------------
### = NOTES:
### =  ( 1) Function for inverting for emissions with obs
### =----------------------------------------------------------------------
### = SUBFUNCTIONS:
### =  ( 1) estimate_ems         :: Computes emissions.
### =  ( 2) build_t              :: Builds the temporal correlation matrix.
### =  ( 3) build_xy             :: Builds the spatial correlation matrix.
### =  ( 4) invert_ems           :: Function to do the inversion.
### =  ( 5) debug_invert_ems     :: Same as above but with some debugging.
### =  ( 6) hq_parallel          :: Computes the HQ matrix product.
### =  ( 7) hq_parallel_function :: Utility function for computing HQ.
### =  ( 8) vsplit               :: Utility function for splitting arrays.
### =  ( 9) hq                   :: Computes the HQ matrix product.
### =  (10) hqht                 :: Computes the HQH^T matrix product.
### =  (11) hqht_fast            :: Computes the HQH^T matrix product.
### =----------------------------------------------------------------------
### = INPUTS:
### =  (  ) N/A
### =----------------------------------------------------------------------
### = OUTPUTS:
### =  (  ) N/A
### =======================================================================

### Libraries
using Distributed
@everywhere using NetCDF
@everywhere using CSV
@everywhere using Dates
@everywhere using LinearAlgebra
@everywhere using Random
@everywhere using StatsBase
@everywhere using Statistics
@everywhere using SparseArrays
@everywhere using MKLSparse
@everywhere using Optim
@everywhere using LBFGSB
@everywhere using GeometricalPredicates

### Load my own functions
include("./read_funcs.jl")
include("./reshape_funcs.jl")

### Globals
global origin
global IntType    # Precision we are using for integers
global FltType    # Precision we are using for floats
global lowBound   # Smallest number to include in our sparse arrays
global diag_prior # Are we using a diagonal prior?
global minUncert  # Lower bound on the prior uncertainty (umol/m2/s)
global minObsErr  # Lower bound on the measurement uncertainty (ppm)
global obsFreq    # Observations per hour


### ========================================
### Estimate the emissions
### ========================================

### Function to estimate emissions
@everywhere function estimate_ems(input_data)
   ### Get data we need
   # Unpack data from the tuples
   (inv_data,foot_data,ems_data)                                                                                                = input_data
   (lon,lat,latLim,lonLim,startTime,endTime,ems_uncert,modErr,k_fold_dat)                                                       = inv_data
   (obsInfo,foot,nHr)                                                                                                           = foot_data
   (ems,emsTimes)                                                                                                               = ems_data
   (obsTime,obsLat,obsLon,obsAGL,obsCO2,obsErr,bkgCO2_NOAA,bkgErr_NOAA,bkgCO2_NASA,bkgErr_NASA,bkgCO2_AMFX,bkgErr_AMFX,trajDat) = obsInfo
   # Clear the tuples
   input_data,foot_data,ems_data = nothing,nothing,nothing
   Base.GC.gc() # Force garbage collection
   # Get sizes
   gridSize = size(ems)[1:2]
   emsSize  = size(ems)
   nObs     = length(foot)
   nEms     = length(emsTimes)
   # Choose the background to use
   bkgCO2 = bkgCO2_NOAA
   bkgErr = [maximum([bkgErr_NOAA[i],bkgErr_NASA[i]]) for i in 1:length(bkgCO2)]
   for i = 1:length(bkgErr)
      if isnan(bkgErr[i])
         if isnan(bkgErr_NOAA[i]) && isnan(bkgErr_NASA[i])
            bkgErr[i] = NaN
         elseif isnan(bkgErr_NASA[i])
            bkgErr[i] = bkgErr_NOAA[i]
         elseif isnan(bkgErr_NOAA[i])
            bkgErr[i] = bkgErr_NASA[i]
         end
      end
   end
   # Define regions for using the AmeriFlux data
   ameriflux_lon = -121.8
   ameriflux_lat =   38.2
   maxDist       =   25.0 # maximum allowable distance from the AmeriFlux site
   # Find trajectories coming from the east
   (endLon,endLat,endAGL,amfLon,amfLat,amfAGL,amfTime) = trajDat
   amfDist = 110 .* sqrt.( (amfLon .- ameriflux_lon).^2 + (amfLat .- ameriflux_lat).^2)
   amfUSE  = amfDist .< maxDist
   for i in 1:length(bkgCO2)
      if amfUSE[i]
         bkgCO2[i] = bkgCO2_AMFX[i]
         bkgErr[i] = maximum([bkgErr[i],bkgErr_AMFX[i]])
      end
   end
   # Throw out nighttime data
   #for i in 1:length(bkgCO2)
   #   if 4 < Dates.value(Dates.Hour(obsTime[i])) && Dates.value(Dates.Hour(obsTime[i])) < 13 # Throw out data between 10pm and 5am
   #      obsLon[i] = -9999
   #      obsLat[i] = -9999
   #   end
   #end

   ### Reduce the matrix sizes
   # Index for the observations
   ind_OBS = map(Bool,ones(size(foot)))
   for i in 1:nObs
      ind_OBS[i] = ind_OBS[i] &&
                   startTime <= obsTime[i] && obsTime[i] <  (endTime + Dates.Hour(nHr)) &&
                   latLim[1] <= obsLat[i]  && obsLat[i]  <= latLim[2]                   &&
                   lonLim[1] <= obsLon[i]  && obsLon[i]  <= lonLim[2]
   end
   # Index for the emissions
   ind_EMS      = map(Bool,ones(length(emsTimes)))
   ind_EMS_calc = map(Bool,ones(length(emsTimes)))
   for i in 1:nEms
      ind_EMS[i] = ind_EMS[i] &&
               (startTime - Dates.Hour(nHr)) <= emsTimes[i] && emsTimes[i] < (endTime + Dates.Hour(nHr))
      ind_EMS_calc[i] = ind_EMS_calc[i] &&
               startTime <= emsTimes[i] && emsTimes[i] < endTime
   end
   ind_EMS_calc = ind_EMS_calc[ind_EMS]
   # Index for reduced comparison region
   latEval_ind = map(Bool,ones(size(lat)))
   lonEval_ind = map(Bool,ones(size(lon)))
   for i in 1:length(lat)
      latEval_ind[i] = latLim[1] <= lat[i] && lat[i] <= latLim[2]
   end 
   for i in 1:length(lon)
      lonEval_ind[i] = lonLim[1] <= lon[i] && lon[i] <= lonLim[2]
   end 
   # Reduce dimensions
   foot     = foot[ind_OBS]
   obsTime  = obsTime[ind_OBS]
   obsLat   = obsLat[ind_OBS]
   obsLon   = obsLon[ind_OBS]
   obsAGL   = obsAGL[ind_OBS]
   obsCO2   = obsCO2[ind_OBS]
   obsErr   = obsErr[ind_OBS]
   bkgCO2   = bkgCO2[ind_OBS]
   bkgErr   = bkgErr[ind_OBS]
   emsTimes = emsTimes[ind_EMS]
   ems      = ems[:,:,ind_EMS]
   # Get new dimensions
   nObs = length(foot)
   nEms = length(emsTimes)
   # Build the model error
   modErr = Array{FltType,1}(undef,nObs)
   for i = 1:nObs
      modErr[i] = mod_err[2,findmin(abs.(mod_err[1,:].-Dates.hour(obsTime[i])))[2]]
   end

   ### Construct "K" and "x"
   # Initialize the empty arrays
   nX = prod(size(ems))
   nG = length(lon)*length(lat)
   # Build an array with start times and indices
   ind_arr = Array{IntType,2}(undef,(nEms,2))
   icount  = 1
   for i = 1:nEms
      ind_arr[i,:] = [icount,Dates.value(emsTimes[i] - origin)]
      icount      += nG
   end
   # Construct the arrays
   x_pri = flatten_ems(ems,emsTimes,ind_arr,gridSize)
   K_mat = flatten_Jacobian(foot,nHr,gridSize,ind_arr,nX,nObs)
   # Make matrices for the inversion and clear up space
   foot  = nothing
   K_mat = convert(SparseMatrixCSC{FltType,IntType},K_mat)
   x_pri = convert(SparseMatrixCSC{FltType,IntType},sparse(x_pri))
   Base.GC.gc() # Force garbage collection

   ### Are we doing k-fold cross-validation?
   (cross_valid, k_fold, k_ind) = k_fold_dat
   iTrain                       = map(Bool,ones(nObs))
   if cross_valid
      obs_to_pick_from  = [i for i in 1:length(obsTime) if startTime <= obsTime[i] && obsTime[i] <  endTime]
      nDraw             = convert(IntType,floor(length(obs_to_pick_from)*(k_fold/100)))
      obs_eval          = StatsBase.sample(MersenneTwister(k_ind),obs_to_pick_from,nDraw,replace=false)
      iTrain[obs_eval] .= false
   end

   ### Divide the observational error by sqrt(n) and keep a lower bound of 1 ppm
   obsErr                       = obsErr ./ sqrt(obsFreq)
   obsErr[obsErr .< minObsErr] .= minObsErr;

   ### Create the observations
   y_obs = convert(SparseMatrixCSC{FltType,IntType},sparse(obsCO2))
   y_bkg = convert(SparseMatrixCSC{FltType,IntType},sparse(bkgCO2))
   y_err = convert(SparseMatrixCSC{FltType,IntType},sparse(sqrt.( (obsErr).^2 + (bkgErr).^2 + (modErr).^2 )))

   
   ### Define more arrays for the inversion
   y_est    = K_mat * x_pri + y_bkg
   mismatch = (y_obs - y_est)

   ### Create the covariance matrices
   # Observational
   if diag_prior # Diagonal prior
      So_d = (y_err[iTrain]).^2
   else # Full covariance
      # Parameters
      tau_time, tau_space = 1, 2 # hour, km
      # Build diagonals
      So_d = (y_err).^2
      So   = sparse(convert(FltType,1.0)I,size(So_d)[1],size(So_d)[1])
      for i = 1:size(So_d)[1]
         So[i,i]  = So_d[i]
      end 
      # Add the off-diagonals
      for i = 1:nObs
         for j = 1:i 
            if j < i
               time_val   = Dates.value(Dates.Hour(obsTime[i] - obsTime[j]))
               dist_val   = sqrt( ( (obsLon[i]-obsLon[j])*110 )^2 + ( (obsLat[i]-obsLat[j])*110 )^2 )
               time_decay = exp( -abs(time_val) / tau_time  )
               dist_decay = exp( -abs(dist_val) / tau_space )
               sig_val    = time_decay * dist_decay * sqrt(So_d[i]*So_d[j])
               if (time_decay * dist_decay) > lowBound
                  So[i,j] = sig_val
                  So[j,i] = sig_val
               end 
            end 
         end 
      end 
      # Remove the non-training data
      So = So[iTrain,iTrain]
   end
   # Prior
   if diag_prior # Diagonal prior
      Sa_d = (ems_uncert .* x_pri).^2
      Sa_d = [ maximum([x,minUncert^2]) for x in Sa_d ]
      Sa_d = convert(SparseMatrixCSC{FltType,IntType},sparse(Sa_d))
   else # Full covariance
      # Parameters
      tau_day, tau_hr, tau_len = 1, 5, 5 # days, hours, & km
      mVal = minUncert/ems_uncert # Minimum value for the covariance
      # Filenames for the covariance matrices
      fname_t  = @sprintf("./stored_data/corr_T_%i_%id_%ih.csv",nEms,tau_day,tau_hr)
      fname_xy = @sprintf("./stored_data/corr_XY_%i_%ikm.csv",nG,tau_len)
      # Have we read this temporal matrix before?
      if isfile(fname_t)
         # Load the data and store the indices
         f = open(fname_t);  dat_t  = readdlm(f,','); close(f)
         # Pull out the indices
         I_t , J_t,  V_t  = dat_t[:,1],  dat_t[:,2],  dat_t[:,3]
         # Clear
         dat_t = nothing
         Base.GC.gc()
      else # Get the data
         Sa_t  = build_t(tau_day,tau_hr,x_pri,ind_arr)
         # Get the non-zero elements and save them
         I_t ,J_t ,V_t  = findnz(Sa_t);
         # Save the non-zero indices
         df_t  = DataFrame([:I => I_t,  :J => J_t,  :V => V_t])
         CSV.write(fname_t,  df_t,  writeheader=false)
         # Clear
         Sa_t, df_t = nothing, nothing
         Base.GC.gc()
      end
      # Have we read this spatial matrix before?
      if isfile(fname_xy)
         # Load the data and store the indices
         f = open(fname_xy); dat_xy = readdlm(f,','); close(f)
         # Pull out the indices
         I_xy, J_xy, V_xy = dat_xy[:,1], dat_xy[:,2], dat_xy[:,3]
         # Clear
         dat_xy = nothing
         Base.GC.gc()
      else # Get the data
         # Define the grid
         lat_grid = zeros(FltType,nG)
         lon_grid = zeros(FltType,nG)
         for lonInd = 1:gridSize[1]
            for latInd = 1:gridSize[2]
               ij           = 1 .+ (latInd .- 1) .+ (lonInd .- 1).*gridSize[2]
               lon_grid[ij] = lon[lonInd].*110 # the 110 is to convert from deg to km
               lat_grid[ij] = lat[latInd].*110
            end
         end
         # Build it
         Sa_xy = build_xy(tau_len,x_pri,ind_arr,lon_grid,lat_grid)
         # Get the non-zero elements and save them
         I_xy,J_xy,V_xy = findnz(Sa_xy);
         # Save the non-zero indices
         df_xy = DataFrame([:I => I_xy, :J => J_xy, :V => V_xy])
         CSV.write(fname_xy, df_xy, writeheader=false)
         # Clear
         Sa_xy, df_xy = nothing, nothing
         Base.GC.gc()
      end
      # Define the grid
      lat_grid = zeros(IntType,nG)
      lon_grid = zeros(IntType,nG)
      tim_grid = zeros(IntType,nEms)
      for lonInd = 1:gridSize[1]
         for latInd = 1:gridSize[2]
            ij           = 1 .+ (latInd .- 1) .+ (lonInd .- 1).*gridSize[2]
            lon_grid[ij] = lonInd
            lat_grid[ij] = latInd
         end
      end
      # Build a vector with the mean emissions at each timestep
      tEms = [maximum([mean(ems[:,:,i]),mVal]) for i = 1:nEms]
      # Build a vector with the mean emissions at each spatial location
      sEms = zeros(nG)
      for lonInd = 1:gridSize[1]
         for latInd = 1:gridSize[2]
            ij       = 1 .+ (latInd .- 1) .+ (lonInd .- 1).*gridSize[2]
            sEms[ij] = maximum([mean(ems[lonInd,latInd,:]),mVal])
         end
      end
      # Get temporal covariance
      iInd = convert.(IntType,mod.(I_t .- 1,nEms) .+ 1)
      jInd = convert.(IntType,mod.(J_t .- 1,nEms) .+ 1)
      C_t  = sqrt.(sEms[iInd].*sEms[jInd])
      # Get spatial covariance
      iLon = lon_grid[convert.(IntType,mod.(I_xy .- 1,nG) .+ 1)]
      iLat = lat_grid[convert.(IntType,mod.(I_xy .- 1,nG) .+ 1)]
      jLon = lon_grid[convert.(IntType,mod.(J_xy .- 1,nG) .+ 1)]
      jLat = lat_grid[convert.(IntType,mod.(J_xy .- 1,nG) .+ 1)]
      iInd = 1 .+ (iLat .- 1) + (iLon .- 1).*gridSize[2]
      jInd = 1 .+ (jLat .- 1) + (jLon .- 1).*gridSize[2]
      C_xy = sqrt.(sEms[iInd].*sEms[jInd])
      # Apply the emission uncertainty and correlation
      C_t  = ems_uncert .* V_t  .* C_t
      C_xy = ems_uncert .* V_xy .* C_xy
      # Make the sparse matrices
      Sa_t  = sparse(I_t, J_t, C_t, nEms,nEms)
      Sa_xy = sparse(I_xy,J_xy,C_xy,nG,  nG)
      Base.GC.gc() # Garbage collect
      ## Save the covariances
      #df_t  = DataFrame([:I => I_t,  :J => J_t,  :V => C_t])
      #df_xy = DataFrame([:I => I_xy, :J => J_xy, :V => C_xy])
      #CSV.write(@sprintf("./stored_data/cov_T_%i_%id_%ih.csv",nEms,tau_day,tau_hr), df_t, writeheader=false)
      #CSV.write(@sprintf("./stored_data/cov_XY_%i_%ikm.csv",nG,tau_len), df_xy, writeheader=false)
   end

   ### Invert
   if diag_prior # Diagonal
      (x_hat,dofs) = invert_ems(K_mat[iTrain,:],mismatch[iTrain],Sa_d,So_d,x_pri)
   else # Full covariance
      #(x_hat,dofs) = invert_ems(K_mat[iTrain,:],mismatch[iTrain],Sa_t,Sa_xy,So,x_pri)
      (x_hat,dofs) = debug_invert_ems(K_mat[iTrain,:],mismatch[iTrain],Sa_t,Sa_xy,So,x_pri)
   end
   Base.GC.gc() # Garbage collect

   ### Store the comparison to the observations
   # Get the data
   obs_co2       = obsCO2
   obs_co2_err   = obsErr
   bkg_co2       = bkgCO2
   bkg_co2_err   = bkgErr
   model_err     = modErr
   prior_co2     = Array(K_mat * x_pri)[:] + bkg_co2
   posterior_co2 = Array(K_mat * x_hat)[:] + bkg_co2
   # Convert types
   obs_co2       = convert(Array{FltType,1},obs_co2)
   obs_co2_err   = convert(Array{FltType,1},obs_co2_err)
   bkg_co2       = convert(Array{FltType,1},bkg_co2)
   bkg_co2_err   = convert(Array{FltType,1},bkg_co2_err)
   model_err     = convert(Array{FltType,1},model_err)
   prior_co2     = convert(Array{FltType,1},prior_co2)
   posterior_co2 = convert(Array{FltType,1},posterior_co2)
   # Filter for observations that are out of our time window
   ind_OBS = map(Bool,ones(size(posterior_co2)))
   for i in 1:nObs
      ind_OBS[i] = ind_OBS[i] &&
                   startTime <= obsTime[i] && obsTime[i] <  endTime   &&
                   latLim[1] <= obsLat[i]  && obsLat[i]  <= latLim[2] &&
                   lonLim[1] <= obsLon[i]  && obsLon[i]  <= lonLim[2]
   end
   # Create the output structure
   obs_time      = obsTime[ind_OBS]
   obs_lat       = obsLat[ind_OBS]
   obs_lon       = obsLon[ind_OBS]
   obs_agl       = obsAGL[ind_OBS]
   obs_co2       = obs_co2[ind_OBS]
   prior_co2     = prior_co2[ind_OBS]
   posterior_co2 = posterior_co2[ind_OBS]
   bkg_co2       = bkg_co2[ind_OBS]
   obs_co2_err   = obs_co2_err[ind_OBS]
   model_err     = model_err[ind_OBS]
   bkg_co2_err   = bkg_co2_err[ind_OBS]
   iTrain        = iTrain[ind_OBS]
   obs_compare   = (obs_time,obs_lon,obs_lat,obs_agl,obs_co2,prior_co2,posterior_co2,bkg_co2,obs_co2_err,model_err,bkg_co2_err,iTrain)

   ### Evaluate the solution in our region of interest
   # Expand the posterior
   x_hat    = expand_ems(x_hat,emsTimes,ind_arr,gridSize)
   x_pri    = expand_ems(x_pri,emsTimes,ind_arr,gridSize)
   # Reduced domain for comparison
   x_pri    = x_pri[lonEval_ind,latEval_ind,ind_EMS_calc]
   x_hat    = x_hat[lonEval_ind,latEval_ind,ind_EMS_calc]
   lat_red  = lat[latEval_ind]
   lon_red  = lon[lonEval_ind]
   emsTimes = emsTimes[ind_EMS_calc]
   
   ### Return the results
   return (lon_red,lat_red,emsTimes,x_pri,x_hat,obs_compare,dofs)
end


### ========================================
### Build the covariance matrices
### ========================================

### Build correlation matrices
# Temporal
function build_t(tau_day,tau_hr,x_pri,ind_arr)
   # Get dimensions
   nEms = convert(IntType,size(ind_arr)[1])
   nG   = convert(IntType,ind_arr[2,1] - ind_arr[1,1])
   # Initialize arrays
   tmpA     = Array{FltType,1}(undef,nG)
   tmpB     = Array{FltType,1}(undef,nG)
   Sa_t     = spzeros(FltType,nEms,nEms)
   ind_time = convert(Array{IntType,1},ind_arr[:,1])
   # Build the matrix
   for i = 1:nEms
      tmpA[:] .= 0
      tmpA[:]  = vec(x_pri[ind_time[i]:ind_time[i]+nG-1])
      for j = 1:i
         tmpB[:]   .= 0
         tmpB[:]    = vec(x_pri[ind_time[j]:ind_time[j]+nG-1])
         hour_apart = (ind_arr[i,2] - ind_arr[j,2])/(60*60*1000)    # hours apart
         days_apart = (ind_arr[i,2] - ind_arr[j,2])/(60*60*1000)/24 # days apart
         hour_apart = 12   - abs(12   - mod(hour_apart,24)) # need to mod the hours
         temp_hours = exp( -abs(hour_apart) / tau_hr  )
         temp_days  = exp( -abs(days_apart) / tau_day )
         cor_val    = cor(tmpA[:],tmpB[:])
         if isnan(cor_val)
            cor_val = 1
         end
         sig_val = cor_val * temp_hours * temp_days
         if sig_val > lowBound
            Sa_t[i,j] = sig_val
            Sa_t[j,i] = sig_val
         end
      end
   end
   return (Sa_t)
end
# Spatial
function build_xy(tau_len,x_pri,ind_arr,lon_grid,lat_grid)
   # Function to check if values are all zero
   @inline function allzero(x)
      @inbounds for i=1:length(x)
         x[i] == 0 || return false
      end
      return true
   end
   # Parameters
   min_dist = 30
   # Get dimensions
   nEms = convert(IntType,size(ind_arr)[1])
   nG   = convert(IntType,ind_arr[2,1] - ind_arr[1,1])
   # Initialize arrays
   emsALL   = zeros(FltType,nEms,nG)
   dist_mat = zeros(FltType,nG,nG)
   ocean    = map(Bool,ones(nG,1))
   Sa_xy    = zeros(FltType,nG,nG)
   ind_time = convert(Array{IntType,1},ind_arr[:,1])
   # First pass to get data
   for i = 1:nG
      #@printf("Pass 1 #%i/%i (%5.2f):     \n",i,nG,(i/nG*100))
      emsALL[:,i]   = vec(x_pri[ind_time .- 1 .+ i])
      ocean[i]      = allzero(emsALL[:,i])
      dist_mat[i,:] = sqrt.( (lon_grid[i] .- lon_grid).^2 .+ (lat_grid[i] .- lat_grid).^2 )
   end
   # Second pass to fill
   for i = 1:nG
      #@printf("Pass 2 #%i/%i (%5.2f):     \n",i,nG,(i/nG*100))
      for j = 1:i
         if dist_mat[i,j] < min_dist
            # don't allow ocean-land interactions
            if (ocean[i] && ocean[j]) || (!ocean[i] && !ocean[j])
               if (ocean[i] && ocean[j]) # ocean
                  cor_val = 1.0
               else # land
                  cor_val = nEms > 1 ? cor(emsALL[:,i],emsALL[:,j]) : NaN 
                  if isnan(cor_val)
                     cor_val = 1.0
                  end
               end
               dist_decay = exp( -abs(dist_mat[i,j]) / tau_len )
               sig_val    = cor_val * dist_decay
               if sig_val > lowBound
                  Sa_xy[i,j] = sig_val
                  Sa_xy[j,i] = sig_val
               end
            end
         end
      end 
   end 
   # Third pass to ensure there are no ocean-land interactions
   for i = 1:nG
      #@printf("Pass 3 #%i/%i (%5.2f):     \n",i,nG,(i/nG*100))
      for j = 1:i
         if dist_mat[i,j] < min_dist
            # don't allow ocean-land interactions
            if (ocean[i] && !ocean[j]) || (!ocean[i] && ocean[j])
               Sa_xy[i,j] = 0
               Sa_xy[j,i] = 0
            end
         end
      end 
   end 
   Sa_xy = sparse(Sa_xy)
   return (Sa_xy)
end


### ========================================
### Inversion routines
### ========================================

### Invert
# Sparse with full prior covariance
function invert_ems(K_mat::SparseMatrixCSC,mismatch::AbstractArray,Sa_t::SparseMatrixCSC,Sa_xy::SparseMatrixCSC,So::SparseMatrixCSC,x_pri)
   # Initialize DOFs
   dofs = 0
   # Invert: xp = xa + (KB)'(KBK' + R)^-1(y - Kxa)
   KSa      = hq(K_mat,Sa_t,Sa_xy);            Base.GC.gc()
   G        = hqht_fast(KSa,K_mat,Sa_t,Sa_xy); Base.GC.gc()
   G        = G + So;                          Base.GC.gc()
   mismatch = G \ mismatch;                    Base.GC.gc()
   x_dif    = KSa' * sparse(mismatch);         Base.GC.gc()
   x_hat    = x_pri + x_dif;                   Base.GC.gc()
   # Return data
   return (x_hat,dofs)
end
# Sparse with full prior covariance
function debug_invert_ems(K_mat::SparseMatrixCSC,mismatch::AbstractArray,Sa_t::SparseMatrixCSC,Sa_xy::SparseMatrixCSC,So::SparseMatrixCSC,x_pri)
   # Initialize DOFs
   dofs = 0
   # Invert: xp = xa + (KB)'(KBK' + R)^-1(y - Kxa)
@printf("HQ:     ")
   @time KSa      = hq(K_mat,Sa_t,Sa_xy);            Base.GC.gc()
@printf("G:      ")
   @time G        = hqht_fast(KSa,K_mat,Sa_t,Sa_xy); Base.GC.gc()
@printf("G:      ")
   @time G        = G + So;                          Base.GC.gc()
@printf("inv(G): ")
   @time mismatch = G \ Array(mismatch[:]);             Base.GC.gc()
@printf("KSaTG:  ")
   @time x_dif    = KSa' * sparse(mismatch);         Base.GC.gc()
@printf("x_hat:  ")
   @time x_hat    = x_pri + x_dif;                   Base.GC.gc()
   # Return data
   return (x_hat,dofs)
end
# Sparse with diagonal covariance matrices
function invert_ems(K_mat::SparseMatrixCSC,mismatch,Sa_d,So_d,x_pri)
   # Make covariance matrices and initialize DOFs
   dofs = 0
   Sa   = sparse(convert(FltType,1.0)I,size(Sa_d)[1],size(Sa_d)[1])
   So   = sparse(convert(FltType,1.0)I,size(So_d)[1],size(So_d)[1])
   for i = 1:size(Sa_d)[1]
      Sa[i,i]  = Sa_d[i]
   end
   for i = 1:size(So_d)[1]
      So[i,i]  = So_d[i]
   end
   # Invert: xp = xa + BK'(KBK' + R)^-1(y - Kxa)
   KSa      = K_mat * Sa;                   Base.GC.gc()
   G        = KSa * transpose(K_mat) + So;  Base.GC.gc()
   mismatch = G \ mismatch;                 Base.GC.gc()
   x_dif    = transpose(KSa) * mismatch;    Base.GC.gc()
   x_hat    = x_pri + x_dif;                Base.GC.gc()
   # Return data
   return (x_hat,dofs)
end


### ========================================
### Yadav & Michalak matrix operations
### ========================================

### Compute HQ and HQHT
# HQ indirect method (parallel)
@everywhere function hq_parallel(H::SparseMatrixCSC,D::SparseMatrixCSC,E::SparseMatrixCSC)
   # Get important dimensions
   n    = size(H)[1]
   p, q = size(D)
   r, t = size(E)
   # Initialize the output matrix
   HQ_out = convert(typeof(H),spzeros(n,q*t))
   # Break up the work
   jvec = vsplit([1:p],nprocs())
   nnzA = nnz(H)
   for i = 1:q
      HQsum = mapreduce( fetch, +, Any[ @spawnat pp hq_parallel_function(H,D[:,i],jvec[pp],n,r,nnzA) for pp = procs() ] )
      HQ_out[:,i*t-t+1:i*t] = HQsum * E
   end
   return (HQ_out)
end
@everywhere function hq_parallel_function(H,D,jvec,n,r,nnzA)
   HQsum  = convert(typeof(D),spzeros(n,r))
   rowInd = Array(IntType,nnzA)
   colInd = Array(IntType,nnzA)
   ValPtr = Array(eltype(H),nnzA)
   for j = jvec
      if D[j] != 0
         rowInd[:] = 0
         colInd[:] = 0
         ValPtr[:] = 
         counter   = 0
         jstart    = (j-1)*r+1
         jend      = j*r
         for ii = 1:n
            for jj = jstart:jend
               if H[ii,jj] != 0 # Skip if H is zero
                  counter        += 1
                  rowInd[counter] = ii
                  colInd[counter] = jj-jstart+1
                  if D[j] == 1 # Skip multiplication if D is 1
                     ValPtr[counter] = HQsum[ii,colInd[counter]] + H[ii,jj]
                  else
                     ValPtr[counter] = HQsum[ii,colInd[counter]] + H[ii,jj] .* D[j]
                  end
               end
            end
         end
         # Save the data
         HQsum = sparse(rowInd[1:counter],colInd[1:counter],ValPtr[1:counter],n,r)
      end
   end
   return HQsum
end
@everywhere function vsplit(v,r)
   a = Array{Any,1}(undef,0)
   for i = 1:r
      push!(a,v[i:r:length(v)])
   end
   return a
end
# HQ indirect method
function hq(H::SparseMatrixCSC,D::SparseMatrixCSC,E::SparseMatrixCSC)
   # Get important dimensions
   n    = size(H)[1]
   p, q = size(D)
   r, t = size(E)
   # Initialize the output matrix
   HQ_out = convert(typeof(H),spzeros(n,q*t))
   for i = 1:q
      HQsum = convert(typeof(H),spzeros(n,r))
      for j = 1:p
         if D[j,i] != 0 && D[j,i] != 1 # Avoid multiplying by blocks that are zero or one
            HQsum = HQsum + H[:,(j-1)*r+1:j*r] .* D[j,i]
         elseif D[j,i] == 1 # Only add if it is one
            HQsum = HQsum + H[:,(j-1)*r+1:j*r]
         end
      end
      HQ_out[:,i*t-t+1:i*t] = HQsum * E
      #HQ_out[:,i*t-t+1:i*t] = spmm(HQsum, E)
   end
   return (HQ_out)
end
# HQHT indirect method
function hqht(H::SparseMatrixCSC,D::SparseMatrixCSC,E::SparseMatrixCSC)
   # Get important dimensions
   n    = size(H)[1]
   p, q = size(D)
   r, t = size(E)
   # Initialize the output matrix
   HQHT_out = convert(typeof(H),spzeros(n,n))
   for i = 1:q
      HQsum = convert(typeof(H),spzeros(n,r))
      for j = 1:p
         if D[j,i] != 0 && D[j,i] != 1 # Avoid multiplying by blocks that are zero or one
            HQsum = HQsum + H[:,(j-1)*r+1:j*r] .* D[j,i]
         elseif D[j,i] == 1 # Only add if it is one
            HQsum = HQsum + H[:,(j-1)*r+1:j*r]
         end
      end
      HQHT_out = HQHT_out + HQsum * E * H[:,(i-1)*r+1:i*r]'
   end
   return (HQHT_out)
end

# HQHT with a precomputed HQ
function hqht_fast(HQ::SparseMatrixCSC,H::SparseMatrixCSC,D::SparseMatrixCSC,E::SparseMatrixCSC)
   # Get important dimensions
   n    = size(H)[1]
   p, q = size(D)
   r, t = size(E)
   # Initialize the output matrix
   HQHT_out = convert(typeof(H),spzeros(n,n))
   for i = 1:q 
      HQHT_out = HQHT_out + HQ[:,i*t-t+1:i*t]*H[:,(i-1)*r+1:i*r]'
   end 
   return (HQHT_out)
end


### =======================================================================
### =                            E   N   D                                =
### =======================================================================
