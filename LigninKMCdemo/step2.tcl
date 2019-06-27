package require psfgen
package require ligninbuilder
source psfgen.tcl
#Remember the arguments are in the order of where do I look for the .psf files, and where do I want to put the pdbs relative to that path.
::ligninbuilder::makelignincoordinates . .
