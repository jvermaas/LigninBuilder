# Remove overlapping/pathological structures
# Use source step3.tcl from Tk console to run
package require ligninbuilder
::ligninbuilder::minimizestructure . namd2 +p8   "parameters extraparameters.prm 
 parameters ../par_all36_cgenff.prm 
" 
exit
