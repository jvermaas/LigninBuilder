# To generate combined psf/pdb file..
# Use source step2.tcl from Tk console to run
package require psfgen
package require topotools
set name switchgrass_nch_10
topology top_lignin.top
resetpsf

for {set i 1} {$i <=  10  }  {incr i} {
readpsf switchgrass_chnum_$i.psf
coordpdb switchgrass_chnum_$i.pdb
}
writepdb $name.pdb
writepsf $name.psf

# Generate extraparameters.prm
exec python3 ../findmissingterms.py
# Generate GROMACS *.top file 
mol new $name.psf
mol addfile $name.pdb
topo writegmxtop switchgrass_nch_10.top [list  par_lignin.prm  extraparameters.prm ]
exit
