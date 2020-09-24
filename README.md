# UrbanInversion
README file for the urban inversion using BEACO2N data from Turner et al. (2020)

Alex Turner

September 24, 2020


# Info
 * The published paper can be found here: "XX".
 * BEACO2N observations can be found on the project website here: "http://beacon.berkeley.edu".
 * Users will need to link the obs files into the obs directory and the emissions files into the ems directory: see `linkDat.csh`.
 * Users should then edit the time window and run parameters in `make_scripts.csh` before building the driver scripts.
 * Additional inversion parameters can be modified in `templates/template_est_fluxes.jl`.
 * The actual inversion code is located in `Util/solve_inv_funcs.jl`.
 * Input data (obs and emission files) and output will be uploaded to the ORNL DAAC upon paper acceptance.

