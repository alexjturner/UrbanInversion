### =======================================================================
### = reshape_funcs.jl
### = Alex Turner
### = 03/09/2020
### =----------------------------------------------------------------------
### = NOTES:
### =  ( 1) Functions to reshape emissions and the jacobian
### =----------------------------------------------------------------------
### = SUBFUNCTIONS:
### =  ( 1) flatten_ems      :: Reshapes the 3D emissions to a 1D vector.
### =  ( 2) expand_ems       :: Reshapes the 1D emissions to a 3D array.
### =  ( 3) flatten_jacobian :: Reshapes the structure to a sparse matrix.
### =======================================================================

### Libraries
using Dates

### Globals
global origin
global lowBound
global IntType
global FltType


### ========================================
### Emissions
### ========================================

### Flatten the emissions
@everywhere function flatten_ems(ems,emsTimes,ind_arr,gridSize)
    nEms  = length(emsTimes)
    x_mat = Array{FltType,1}(undef,nEms*prod(gridSize))
    for k = 1:nEms
       timInd     = Dates.value((emsTimes[k]) - origin)
       minval,pos = findmin(abs.(timInd .- ind_arr[:,2]))
       for lonInd = 1:gridSize[1]
          for latInd = 1:gridSize[2]
             ij        = ind_arr[pos,1] .+ (latInd .- 1) .+ (lonInd .- 1).*gridSize[2]
             x_mat[ij] = ems[lonInd,latInd,k]
          end
       end
    end
    return(x_mat)
end

### Unflatten the emissions
@everywhere function expand_ems(x_mat,emsTimes,ind_arr,gridSize)
    nEms = length(emsTimes)
    ems  = Array{FltType,3}(undef,(gridSize[1],gridSize[2],nEms))
    for k = 1:nEms
       timInd     = Dates.value((emsTimes[k]) - origin)
       minval,pos = findmin(abs.(timInd .- ind_arr[:,2]))
       for lonInd = 1:gridSize[1]
          for latInd = 1:gridSize[2]
             ij                   = ind_arr[pos,1] .+ (latInd .- 1) .+ (lonInd .- 1).*gridSize[2]
             ems[lonInd,latInd,k] = x_mat[ij]
          end
       end
    end
    return(ems)
end


### ========================================
### Jacobian
### ========================================

### Flatten the jacobian
@everywhere function flatten_Jacobian(foot,nHr,gridSize,ind_arr,nX,nObs)
    # Get the total number of elements
    nG   = Int(ind_arr[2,1] - ind_arr[1,1])
    nnzA = 0
    for i = 1:nObs
       nnzA += sum(map(x -> length(x),(foot[i])[4]))
    end
    # Initialize the sparse arrays
    iStore = map(Bool,zeros(nnzA))
    rowInd = Array{IntType,1}(undef,nnzA)
    colInd = Array{IntType,1}(undef,nnzA)
    ValPtr = Array{FltType,1}(undef,nnzA)
    ValPtr = zeros(FltType,nnzA)
    # Populate the arrays
    ij = 0
    for i = 1:nObs
       ftime = (foot[i])[1] # Get the timesteps for this footprint
       for j = 1:length((foot[i])[1])
          hrDiff = Dates.value(Dates.Hour(maximum(ftime) - ftime[j]))
          if hrDiff <= nHr # Only include if we're within "nHr" of the observation
             # Get information from the structure
             II = ((foot[i])[2])[j]
             JJ = ((foot[i])[3])[j]
             VV = ((foot[i])[4])[j]
             # Number of non-zero elements
             nnzC       = length(VV)
             # Fill the row indices and actual values
             rowInd[(ij+1):(ij+nnzC)] .= i
             ValPtr[(ij+1):(ij+nnzC)]  = VV
             # Get the column index
             lonInd     = II
             latInd     = JJ
             timInd     = try Dates.value(ftime[j] - origin); catch; -9999; end
             minval,pos = findmin(abs.(timInd .- ind_arr[:,2]))
             if minval < (30*60*1000) # Is this within 30 minutes of our time?
                colInd[(ij+1):(ij+nnzC)]  = ind_arr[pos,1] .+ (latInd .- 1) .+ (lonInd .- 1).*gridSize[2]
                iStore[(ij+1):(ij+nnzC)] .= true
             else
                ValPtr[(ij+1):(ij+nnzC)] .= 0
             end
             ij += nnzC # Keep track of our location in the array
          end
       end
    end
    # Now construct our sparse matrix
    foot  = nothing; Base.GC.gc()
    inds  = iStore  # 0 .< ValPtr
    K_mat = sparse(rowInd[inds],colInd[inds],ValPtr[inds],nObs,nX)
    return(K_mat)
end


### =======================================================================
### =                            E   N   D                                =
### =======================================================================
