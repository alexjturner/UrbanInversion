#!/bin/csh -f

### Specify the base directory
set baseDir = "/global/home/users/aturner/BEACON_Inv"
set obsDir  = "/global/scratch/aturner/STILT/BEACON/out"
set emsDir  = "/global/scratch/aturner/BEACON_Emissions/ems"

### Define the link pattern
set yrPat  = "20??"
set monPat = "0[1-6]"
set dayPat = "??"
set hrPat  = "??"
# Obs
set timPat = "${yrPat}${monPat}${dayPat}${hrPat}"
set obsVar = "obs_${timPat}_*.nc"
# Emissions
set timPat = "${yrPat}x${monPat}x${dayPat}x${hrPat}"
set emsVar = "BEACON_${timPat}.ncdf"

### Remove the files
rm $baseDir/obs/obs_*
rm $baseDir/ems/BEACON_*
rm $baseDir/ems_2x/BEACON_*

### Link the obs
ln -s ${obsDir}/obs_level2/$obsVar ${baseDir}/obs/.
##ln -s ${obsDir}/obs_level1/$obsVar ${baseDir}/obs/.
echo 'There are '`ls ${baseDir}/obs/obs_* | wc -l`' obs files'

### Link the emissions
ln -s ${emsDir}/turner/$emsVar ${baseDir}/ems/.
echo 'There are '`ls ${baseDir}/ems/*.ncdf | wc -l`' ems files'

exit(0)
