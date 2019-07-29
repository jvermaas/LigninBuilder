package require ligninbuilder
#Just to see what changed, we'll load the un-minimized structure, and compare with the result after minimization
mol load psf L.psf pdb L.pdb
::ligninbuilder::minimizestructures . namd2 "+p8"
mol addfile L.pdb
