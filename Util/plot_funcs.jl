### =======================================================================
### = plot_funcs.jl
### = Alex Turner
### = 03/09/2020
### =----------------------------------------------------------------------
### = NOTES:
### =  ( 1) Functions for plotting spatial data in the Bay Area
### =----------------------------------------------------------------------
### = SUBFUNCTIONS:
### =  (  ) N/A
### =----------------------------------------------------------------------
### = INPUTS:
### =  (  ) N/A
### =----------------------------------------------------------------------
### = OUTPUTS:
### =  (  ) N/A
### =======================================================================

### Libraries
using PyPlot
using DataFrames
using DelimitedFiles


### ========================================
### Function to add the coastline
### ========================================

function addWorld()
   ### Get the coastline
   f = open("./Util/World_Coastline.csv")
   dat = readdlm(f,',')
   close(f)
   for iter in eachindex(dat)
      if dat[iter] <= -9999
         dat[iter] = NaN
      end
   end

   ### Add the dat
   #PyPlot.hold(true)
   PyPlot.plot(dat[:,1],dat[:,2],color=[.4,.4,.4],linewidth=0.4)
end
function addBayArea()
   ### Get the coastline
   #f = open("./Util/BayArea_Coastline.csv")
   f = open("./Util/NOS80k.csv")
   dat = readdlm(f,',')
   close(f)
   for iter in eachindex(dat)
      if dat[iter] <= -9999
         dat[iter] = NaN
      end
   end

   ### Add the dat
   #PyPlot.hold(true)
   PyPlot.plot(dat[:,1],dat[:,2],color=[.4,.4,.4],linewidth=0.4)
end


### ========================================
### Generic plotting routine
### ========================================

### Define my Colormap
# Negative
my_cmap_neg = PyPlot.ColorMap([PyPlot.RGB(0.6745,0.8275,0.6196),
                           PyPlot.RGB(0.9686,0.9686,0.9686),
                           PyPlot.RGB(0.8000,0.7007,0.8170),
                           PyPlot.RGB(0.6314,0.4327,0.6654),
                           PyPlot.RGB(0.4627,0.1647,0.5137)],
                           256,1.)
# Positive
my_cmap_pos = PyPlot.ColorMap([PyPlot.RGB(0.9686,0.9686,0.9686),
                           PyPlot.RGB(0.4627,0.1647,0.5137)],
                           256,1.)
# Blue
my_cmap_blue = PyPlot.ColorMap([PyPlot.RGB(0.9686,0.9686,0.9686),
                           PyPlot.RGB(0.0980,0.0980,0.4392)],
                           256,1.)
# Blue (reverse)
my_cmap_blue_r = PyPlot.ColorMap([PyPlot.RGB(0.0980,0.0980,0.4392),
                           PyPlot.RGB(0.9686,0.9686,0.9686)],
                           256,1.)
# My redblue
my_redblue = PyPlot.ColorMap([PyPlot.RGB(0.0196,0.1882,0.3804),
			      PyPlot.RGB(0.1294,0.4000,0.6745),
			      PyPlot.RGB(0.2627,0.5765,0.7647),
			      PyPlot.RGB(0.5725,0.7725,0.8706),
			      PyPlot.RGB(0.8196,0.8980,0.9412),
			      PyPlot.RGB(0.9686,0.9686,0.9686),
			      PyPlot.RGB(0.9922,0.8588,0.7804),
			      PyPlot.RGB(0.9569,0.6471,0.5098),
			      PyPlot.RGB(0.8392,0.3765,0.3020),
			      PyPlot.RGB(0.6980,0.0941,0.1686),
                              PyPlot.RGB(0.4039,0.0000,0.1216)],
                              256,1.)


### Standard plotting
function plotEMS(lon,lat,Z,lonLims,latLims,zlims,zLabel)
   h = PyPlot.pcolor(lon,lat,Z,cmap=my_cmap_pos,vmin=zlims[1],vmax=zlims[2])
   PyPlot.ylim(latLims)
   PyPlot.xlim(lonLims)
   cb = PyPlot.colorbar(h)
   cb.set_label(zLabel)
   PyPlot.ticklabel_format(useOffset=false)
   PyPlot.tick_params(axis="both",which="both",bottom="off",top="off",left="off",right="off",labelbottom="off",labelleft="off")
   addBayArea()
end
function plotEMS_neg(lon,lat,Z,lonLims,latLims,zlims,zLabel)
   h = PyPlot.pcolor(lon,lat,Z,cmap=my_cmap_neg,vmin=zlims[1],vmax=zlims[2])
   PyPlot.ylim(latLims)
   PyPlot.xlim(lonLims)
   cb = PyPlot.colorbar(h)
   cb.set_label(zLabel)
   PyPlot.ticklabel_format(useOffset=false)
   PyPlot.tick_params(axis="both",which="both",bottom="off",top="off",left="off",right="off",labelbottom="off",labelleft="off")
   addBayArea()
end
function plotEMS_diff(lon,lat,Z,lonLims,latLims,zlims,zLabel)
   zlims = [-maximum(abs.(zlims)),maximum(abs.(zlims))]
   h = PyPlot.pcolor(lon,lat,Z,cmap=my_redblue,vmin=zlims[1],vmax=zlims[2])
   PyPlot.ylim(latLims)
   PyPlot.xlim(lonLims)
   cb = PyPlot.colorbar(h)
   cb.set_label(zLabel)
   PyPlot.ticklabel_format(useOffset=false)
   PyPlot.tick_params(axis="both",which="both",bottom="off",top="off",left="off",right="off",labelbottom="off",labelleft="off")
   addBayArea()
end
function plotFoot(lon,lat,Z,lonLims,latLims,zlims,zLabel)
   h = PyPlot.pcolor(lon,lat,Z,cmap=my_cmap_blue,vmin=zlims[1],vmax=zlims[2])
   PyPlot.ylim(latLims)
   PyPlot.xlim(lonLims)
   cb = PyPlot.colorbar(h)
   cb.set_label(zLabel)
   PyPlot.ticklabel_format(useOffset=false)
   PyPlot.tick_params(axis="both",which="both",bottom="off",top="off",left="off",right="off",labelbottom="off",labelleft="off")
   addBayArea()
end
function plotFoot_rev(lon,lat,Z,lonLims,latLims,zlims,zLabel)
   h = PyPlot.pcolor(lon,lat,Z,cmap=my_cmap_blue_r,vmin=zlims[1],vmax=zlims[2])
   PyPlot.ylim(latLims)
   PyPlot.xlim(lonLims)
   cb = PyPlot.colorbar(h)
   cb.set_label(zLabel)
   PyPlot.ticklabel_format(useOffset=false)
   PyPlot.tick_params(axis="both",which="both",bottom="off",top="off",left="off",right="off",labelbottom="off",labelleft="off")
   addBayArea()
end


### =======================================================================
### =                            E   N   D                                =
### =======================================================================
