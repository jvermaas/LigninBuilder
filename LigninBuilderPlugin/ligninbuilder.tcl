package provide ligninbuilder 1.1

namespace eval ::ligninbuilder:: {
package require psfgen
package require readcharmmpar
package require exectool
package require topotools

proc rotatebond {idxlist idx idx2 rot {mid top} {f now}} {
	#idxlist is the list of indices that will get rotated. idx is bonded to an element in idxlist, and that is the bond around which we will be rotating.
	set rotatesel [atomselect $mid "index $idxlist" frame $f]
	set basesel [atomselect $mid "index $idx"]
	set partner [atomselect $mid "index $idx2"]
	if {[$basesel num] != 1} {
		puts "Second argument results in a atomselection of the wrong size."
		puts "Atomselection size: [$basesel num]"
		return
	}
	if {[$partner num] != 1} {
		puts "Third argument results in a atomselection of the wrong size."
		puts "Atomselection size: [$parter num]"
		return
	}
	set M [trans bond  [lindex [$basesel get {x y z}] 0] [lindex [$partner get {x y z}] 0] $rot deg]
	$rotatesel move $M
	$partner delete
	$rotatesel delete
	$basesel delete
}
#------------------------------------------------------------------------------
#Split a bond, and see which side is bigger. We'll be moving the smaller of the two sides.
proc selectsmaller { idx1 idx2 {mid top} } {
	set sel1 [atomselect $mid "index $idx1"]
	set sel2 [atomselect $mid "index $idx2"]
	#Expand the tree.
	set newsel1 [atomselect $mid "index $idx1 [join [$sel1 getbonds]] and not index $idx2"]
	set nidx1 [$newsel1 get index]
	set newsel2 [atomselect $mid "index $idx2 [join [$sel2 getbonds]] and not index $nidx1"]
	if { [$sel1 num] == [$newsel1 num] || [$sel2 num] == [$newsel2 num]} {
		puts "Bond not rotateable"
		return
	}
	#Keep expanding the selections until we are out of stuff to expand to on one side.
	while { [$sel1 num] != [$newsel1 num] && [$sel2 num] != [$newsel2 num] } {
		$sel1 delete
		$sel2 delete
		set sel1 $newsel1
		set sel2 $newsel2
		set nidx2 [$newsel2 get index]
		set newsel1 [atomselect $mid "index $nidx1 [join [$sel1 getbonds]] and not index $nidx2"]
		set nidx1 [$newsel1 get index]
		set newsel2 [atomselect $mid "index $nidx2 [join [$sel2 getbonds]] and not index $nidx1"]
	}
	set nidx2 [$newsel2 get index]
	if { [$sel1 num] == [$newsel1 num] } {
		set retval [list $idx2 $idx1 $nidx1]
	} else {
		set retval [list $idx1 $idx2 $nidx2]
	}
	#Cleanup.
	$sel1 delete
	$sel2 delete
	$newsel1 delete
	$newsel2 delete
	return $retval
}
#------------------------------------------------------------------------------
#This is the guy that rotates a bond based on two indices.
proc autorotatebond {idx1 idx2 rot {mid top} {f now}} {
	#Designed to be the "lazy" way of rotating a bond. Just give it two bonded atoms and a rotation, and it'll figure out the rest!
	lassign [selectsmaller $idx1 $idx2 $mid] idx i ilist
	rotatebond $ilist $idx $i $rot $mid $f
}
#Generate coordinates for monomers. They will be on top of each other, but that is OK!
proc placemonomers {mid monomerlist} {
	global env
	animate dup $mid
	#This gets you starting points for each residue. They will be on top of each other, and not at all optimized.
	for { set i 1 } { $i <= [llength $monomerlist] } { incr i } {
		set sel [atomselect $mid "resid $i"]
		set resname [string tolower [lsort -unique [$sel get resname]]]
		if { $resname == "guas" } {
			set resname guai
		} elseif { $resname == "phps" } {
			set resname php
		} elseif { $resname == "syrs" } {
			set resname syr
		}
		set newmol [mol new [file join $env(LIGNINBUILDERDIR) startingstructure $resname.js] waitfor all]
		set newsel [atomselect top "name [$sel get name]" frame [expr {floor(rand()*[molinfo top get numframes])}]]
		$sel set {x y z} [$newsel get {x y z}]
		$sel delete
		$newsel delete
		mol delete $newmol
	}
}
proc applyfitcommands {mid fitcommands} {
	global env
	set counter 1
	foreach fitcommand $fitcommands {
		#animate dup $mid
		set sel1 [atomselect $mid "name \"C\[1-6\]\" and resid [lindex $fitcommand 1]"]
		set sel2 [atomselect $mid "name \"C\[1-6\]\" and resid [lindex $fitcommand 2]"]
		lassign [expandselections $sel1 $sel2 $mid] mtxt1 mtxt2
		set msel1 [atomselect $mid $mtxt1]
		set msel2 [atomselect $mid $mtxt2]
		set name [string tolower [string range [lindex [$sel1 get resname] 0] 0 0][lindex $fitcommand 0][string range [lindex [$sel2 get resname] 0] 0 0]]
		if { [file exists [file join $env(LIGNINBUILDERDIR) startingstructure $name.js]] == 0 && [string range [lindex $fitcommand 0] 0 1] == "B5" } {
			set name [string tolower P[lindex $fitcommand 0]]
		} elseif { [file exists [file join $env(LIGNINBUILDERDIR) startingstructure $name.js]] == 0 && ([lindex $fitcommand 0] in [list "BB" "AOG" "GOG"] || [lindex [$sel2 get resname] 0] == "TRCN") } {
			set name [string tolower P[lindex $fitcommand 0]P]
		} elseif { [lindex $fitcommand 0] == "SPIR" } {
			set name [string tolower PB1[string range [lindex [$sel2 get resname] 0] 0 0]]
		} elseif { [file exists [file join $env(LIGNINBUILDERDIR) startingstructure $name.js]] == 0 } {
			set name [string tolower P[lindex $fitcommand 0][string range [lindex [$sel2 get resname] 0] 0 0]]
		} 
		#puts $name
		if { [file exists [file join $env(LIGNINBUILDERDIR) startingstructure $name.js]] == 0 } {
			error "Can't find appropriate template for $fitcommand"
		}
		set newmol [mol new [file join $env(LIGNINBUILDERDIR) startingstructure $name.js] waitfor all]
		set frame [expr {floor(rand()*[molinfo $newmol get numframes])}]
		set ref1 [atomselect $newmol "name \"C\[1-6\]\" and resid 1" frame $frame]
		set ref2 [atomselect $newmol "name \"C\[1-6\]\" and resid 2" frame $frame]
		$msel1 move [measure fit $sel1 $ref1]
		$msel2 move [measure fit $sel2 $ref2]
		#puts "$counter of [llength $fitcommands]"
		alignmiddleselection $mid [lindex $fitcommand 1] [lindex $fitcommand 2] $newmol $frame
		$sel1 delete
		$sel2 delete
		$msel1 delete
		$msel2 delete
		$ref1 delete
		$ref2 delete
		mol delete $newmol
		# if { $segname == "L12" } {
		# 	animate write js [format "$outdir/initial-$segname-%03d.js" $commandnumber] waitfor all $mid
		# 	incr commandnumber
		# }
		incr counter
	}
}
proc applyterminalpatches {mid outdir segname} {
	#This applies a patch to the C1 end of lignin.
	set makedbbl [atomselect $mid "(same residue as name HO7) and (same residue as name HO8)"]
	set realsegname [lindex [$makedbbl get segname] 0]
	set reslist [lsort -unique -integer [$makedbbl get resid]]
	foreach res $reslist {
		#Make the double bond carbons have a dihedral of 180
		set tmpsel [atomselect $mid "resid $res and name C1 C7 C8 C9"]
		set dihedval [measure dihed [$tmpsel get index] molid $mid]
		autorotatebond [lindex [$tmpsel get index] 1] [lindex [$tmpsel get index] 2] [expr {180-$dihedval}] $mid
		$tmpsel delete
	}
	set allsel [atomselect $mid "all"]
	$allsel writepdb [file join $outdir $segname.pdb]
	resetpsf
	readpsf [file join $outdir $segname.psf]
	coordpdb [file join $outdir $segname.pdb]
	foreach res $reslist {
		patch DBBL $realsegname:$res
	}
	regenerate angles dihedrals
	guesscoord
	$allsel delete
	writepsf [file join $outdir $segname.psf]
	writepdb [file join $outdir $segname.pdb]
}
proc readpsfremarks {psf} {
	set fin [open $psf r]
	set contents [read $fin]
	close $fin
	set lines [split $contents "\n"]
	set fitcommands [list ]
	foreach l $lines {
		if { [string first "REMARKS patch" $l] != -1 } {
			set patchline [regexp -all -inline {\S+} $l]
			set patchname [lindex $patchline 2]
			#puts "[llength $patchline] $patchname"
			#These are the patches that have only two components. "Easy" to build.
			if {[llength $patchline] == 5} {
				lappend fitcommands [list $patchname [lindex [split [lindex $patchline 3] ":"] 1] [lindex [split [lindex $patchline 4] ":"] 1]]
			}
		}
	}
	return $fitcommands
}
#If there are multiple lignins in a single psf, make some attempt to space them out.
proc separatelignins {mid} {
	set asel [atomselect $mid "all"]
	set mw [vecsum [$asel get mass]]
	set radius [expr {6 * (($mw * 0.001) ** (1.0/3.0))}]
	set nfrags [llength [lsort -unique [$asel get fragment]]]
	set dimension [expr {int(ceil( $nfrags ** (1.0/3.0)))}]
	set offset [expr {2 * $radius / $dimension}]
	for { set f 0 } { $f < $nfrags } { incr f } {
		set sel [atomselect $mid "fragment $f"]
		$sel moveby [vecscale -1 [measure center $sel]]
		set multmatrix [list [expr { $f % $dimension }] [ expr { ($f / $dimension) % $dimension }] [ expr { ($f / $dimension) / $dimension }] ]
		$sel moveby [vecscale $offset $multmatrix]
		$sel delete
	}
	$asel delete
}
proc makelignincoordinates {inputdir {outdir .}} {
	global env
	topology [file join $env(LIGNINBUILDERDIR) top_lignin.top]
	set psflist [lsort [glob [file join $inputdir "*psf"]]]

	foreach psf $psflist {
		set mid [mol new $psf]
		set allsel [atomselect top "all"]
		set chargesum [vecsum [$allsel get charge]]
		puts $chargesum
		if { [expr {abs($chargesum - round($chargesum))} > 1e-3] } {
			puts "This isn't an integer..."
			puts "PSF charge sum not an integer for file $psf"
			exit
		}
		set c1sel [atomselect top "name C1"]; #This works because all lignin monomers have a C1 atom that remains unmodified.
		set monomerlist [$c1sel get resname]
		placemonomers $mid $monomerlist
		set fitcommands [readpsfremarks $psf]
		applyfitcommands $mid $fitcommands
		separatelignins $mid
		applyterminalpatches $mid [file join $inputdir $outdir] [file rootname [file tail $psf]]
	}
}
proc makelignin {monomerlist patchdescription {segname L} {outdir .}} {
	global env
	resetpsf
	set i 1
	segment $segname {
		foreach monomer $monomerlist {
			switch $monomer {
				1 {
					set mononame SYR
				}
				2 {
					set mononame PHP
				}
				3 {
					set mononame GUAI
				}
				default {
					error "Something has gone wrong here: $monomer"
				}
			}
			residue $i $mononame
			incr i
		}
	}
	set fitcommands [list ]
	foreach patch $patchdescription {
		lassign $patch m1 m2 c1 c2
		if {$c2 == 1 || $c2 == 8} { 
			if {$c1 == 1} {
				patch BO4 $segname:$m1 $segname:$m2 
				lappend fitcommands [list BO4 $m1 $m2]
			} elseif {$c1 == 2} {
				patch BO4 $segname:$m2 $segname:$m1
				lappend fitcommands [list BO4 $m2 $m1]
			} else {
				error "Something has gone wrong here: B-O-4 $patch"
			}
		} elseif {$c2 == 2} {
			if {$c1 == 1} {
				set monomer [lindex $monomerlist [expr {$m2-1}]]
				switch $monomer {
					1 {
						error "Something has gone wrong here: $patch. B5 cannot bond to SYR"
					}
					2 {
						set patchname B5P
					}
					3 {
						set patchname B5G
					}
					default {
						error "Something has gone wrong here: $monomer"
					}
				}
				patch $patchname $segname:$m1 $segname:$m2
				lappend fitcommands [list $patchname $m1 $m2]
			} else {
				error "Something has gone wrong here: B-5 $patch"
			}
		} elseif {$c2 == 3} {
			if {$c1 == 3} {
				patch 55 $segname:$m1 $segname:$m2
				lappend fitcommands [list 55 $m1 $m2]
			} else {
				error "Something has gone wrong here: 5-5 $patch"
			}
		} elseif {$c2 == 4} {
			if {$c1 == 2} {
				patch 4O5 $segname:$m1 $segname:$m2
				lappend fitcommands [list 4O5 $m1 $m2]
			} elseif {$c1 == 3} {
				patch 4O5 $segname:$m2 $segname:$m1
				lappend fitcommands [list 4O5 $m2 $m1]
			} else {
				error "Something has gone wrong here: 4-O-5 $patch"
			}
		} elseif {$c2 == 5} {
			if {$c1 == 1} {
				patch B1 $segname:$m1 $segname:$m2
				lappend fitcommands [list B1 $m1 $m2]
			} else {
				error "Something has gone wrong here: B-1 $patch"
			}
		} elseif {$c2 == 6} {
			if {$c1 == 4} {
				patch AO4 $segname:$m1 $segname:$m2
				lappend fitcommands [list AO4 $m1 $m2]
			} else {
				error "Something has gone wrong here: A-O-4 $patch"
			}
		} elseif {$c2 == 7} {
			if {$c1 == 1} {
				patch BB $segname:$m1 $segname:$m2
				lappend fitcommands [list BB $m1 $m2]
			} else {
				error "Something has gone wrong here: B-B $patch"
			}
		}
	}
	regenerate angles dihedrals
	writepsf [file join $outdir $segname.psf]
	set mid [mol new [file join $outdir $segname.psf]]
	set allsel [atomselect top "all"]
	set chargesum [vecsum [$allsel get charge]]
	puts $chargesum
	if { [expr {abs($chargesum - round($chargesum))} > 1e-3] } {
		puts "This isn't an integer..."
		exit
	}
	placemonomers $mid $monomerlist
	
	# if { $segname == "L12" } {
	# 	animate write js "$outdir/initial-$segname.js" waitfor all $mid
	# 	set commandnumber 0
	# }
	#Now we look back at the patches we issued, and try to get the local geometry to match!
	applyfitcommands $mid $fitcommands
	
	#Apply patches to C1
	applyterminalpatches $mid $outdir $segname
	#animate write dcd "$segname.dcd" waitfor all $mid
}
proc expandselections {sel1 sel2 {mid top}} {
	set frag [lsort -integer -unique [$sel1 get fragment]]
	set idx1 [$sel1 get index]
	set idx2 [$sel2 get index]
	set reslist2 [lsort -integer -unique [$sel2 get residue]]
	#Expand the tree.
	set newsel1 [atomselect $mid "fragment $frag and index $idx1 [join [$sel1 getbonds]] and not residue $reslist2"]
	set nidx1 [$newsel1 get index]
	set reslist1 [lsort -integer -unique [$newsel1 get residue]]
	set newsel2 [atomselect $mid "fragment $frag and index $idx2 [join [$sel2 getbonds]] and not residue $reslist1"]
	set iteration 0
	#Keep expanding the selections until we are out of stuff to expand to on both sides.
	while { [$sel1 num] != [$newsel1 num] || [$sel2 num] != [$newsel2 num] } {
		#These next two commands short-circuit alot of needless atom-selecting when one of the two selections is complete.
		if { [$sel1 num] == [$newsel1 num] } {
			$newsel2 delete
			set newsel2 [atomselect $mid "fragment $frag and not residue $reslist1"]
		}
		if { [$sel2 num] == [$newsel2 num] } {
			$newsel1 delete
			set newsel1 [atomselect $mid "fragment $frag and not residue $reslist2"]
		}
		if { $iteration } {
			$sel1 delete
			$sel2 delete
		}
		incr iteration
		set sel1 $newsel1
		set sel2 $newsel2
		set nidx2 [$newsel2 get index]
		set reslist2 [lsort -integer -unique [$sel2 get residue]]
		set newsel1 [atomselect $mid "index $nidx1 [join [$sel1 getbonds]] and not residue $reslist2"]
		set nidx1 [$newsel1 get index]
		set reslist1 [lsort -integer -unique [$newsel1 get residue]]
		set newsel2 [atomselect $mid "index $nidx2 [join [$sel2 getbonds]] and not residue $reslist1"]
	}
	set txt1 "index [$newsel1 get index]"
	set txt2 "index [$newsel2 get index]"
	#Cleanup.
	$sel1 delete
	$sel2 delete
	$newsel1 delete
	$newsel2 delete
	return [list $txt1 $txt2]
}
proc alignmiddleselection {mid r1 r2 newmol f} {
	set midsel [atomselect $newmol "(not name \"C\[1-6\]\") and ((resid 1 and withinbonds 1 of resid 2) or (resid 2 and withinbonds 1 of resid 1))"]
	if { [$midsel num] == 0 } {
		$midsel delete
		return
	}
	set idx [$midsel get index]
	set newmidsel [atomselect $newmol "(not name \"C\[1-6\]\") and index $idx [join [$midsel getbonds]]"]
	while { [$midsel num] != [$newmidsel num] } {
		$midsel delete
		set midsel $newmidsel
		set idx [$midsel get index]
		set newmidsel [atomselect $newmol "(not name \"C\[1-6\]\") and index $idx [join [$midsel getbonds]]"]
	}
	set idx [$newmidsel get index]
	$midsel delete
	$newmidsel delete
	set lowsel [atomselect $newmol "resid 1 and index $idx"]
	set highsel [atomselect $newmol "resid 2 and index $idx"]
	if { [$lowsel num] } {
		set lsel [atomselect $mid "resid $r1 and name [$lowsel get name]"]
		if { [$lsel num] } {
			$lowsel delete
			set lowsel [atomselect $newmol "resid 1 and name [$lsel get name]" frame $f]
			$lsel set {x y z} [$lowsel get {x y z}]
		}
		$lsel delete
	}
	$lowsel delete
	if { [$highsel num] } {
		set hsel [atomselect $mid "resid $r2 and name [$highsel get name]"]
		if { [$hsel num] } {
			$highsel delete
			set highsel [atomselect $newmol "resid 2 and name [$hsel get name]" frame $f]
			$hsel set {x y z} [$highsel get {x y z}]
		}
		$hsel delete
	}
	$highsel delete
}
proc buildfromlibrary {inputlibrary outputdirectory} {
	global env
	topology [file join $env(LIGNINBUILDERDIR) top_lignin.top]
	set name [lindex [split [file tail $inputlibrary]] 0]
	file mkdir $outputdirectory
	set fp [open $inputlibrary r]
	set file_data [read $fp]
	close $fp
	set data [split $file_data "\n"]
	set lcounter 0
	for { set lnum 0 } { $lnum < [llength $data] } { incr lnum } {
		set line [lindex $data $lnum]
		if { [string range $line 0 2] == "***" } {
			incr lnum 4
			set linedata [regexp -inline -all -- {\w+}  [lindex $data $lnum] ]
			set seqlength [lrange $linedata 0 2]
			set patches [lrange $linedata 3 end]
			puts $linedata
			set nummonomers [expr {int([vecsum $seqlength])}]
			set numpatches [expr {int([vecsum $patches])}]
			set monolist [list ]
			incr lnum 2
			for { set i 0 } { $i < $nummonomers } { incr i } {
				lappend monolist [lindex [regexp -inline -all -- {\w+}  [lindex $data $lnum] ] 1]
				incr lnum
			}
			incr lnum 3
			set patchlist [list ]
			for { set i 0 } { $i < $numpatches} { incr i } {
				lappend patchlist [regexp -inline -all -- {\w+}  [lindex $data $lnum] ]
				incr lnum
			}
			puts $monolist
			puts $patchlist
			puts $nummonomers
			puts $numpatches
			makelignin $monolist $patchlist L$lcounter $outputdirectory
			incr lcounter
		}
	}
}
# proc findmissing {directory} {
# 	global env
# 	set paramlist [list [file join $env(CHARMMPARDIR) par_all36_carb.prm] [file join $env(LIGNINBUILDERDIR) par_lignin.prm]]
# 	foreach t {bond angle dihed improper atom} {
# 		set ${t}params [list ]
# 	}
# 	foreach p $paramlist {
# 		set parameters [::Pararead::read_charmm_parameters $p]
# 		foreach t {bond angle dihed improper atom} plist [lrange $parameters 0 5] {
# 			lappend ${t}params $plist
# 		}
# 	}
# 	foreach t {bond angle dihed improper atom} {
# 		set ${t}params [join ${t}params]
# 		puts $t
# 	}
# }

proc minimizestructures {directory namdbin namdargs {namdextraconf ""}} {
	global env
	#Copy over parameters from their respective directories.
	set paramlist [list [file join $env(CHARMMPARDIR) par_all36_carb.prm] [file join $env(LIGNINBUILDERDIR) par_lignin.prm] [file join $env(LIGNINBUILDERDIR) extraterms-par_lignin.prm] [file join $env(LIGNINBUILDERDIR) minimize.namd] [file join $env(LIGNINBUILDERDIR) colvars.conf]]
	foreach par $paramlist {
		set f [open $par "r"]
		set out [open [file join $directory [file tail $par]] w ]
		set txt [read -nonewline ${f}]
		puts $out $txt
		close $f
		close $out
	}
	#Minimization command needs to be last.
	set psflist [lsort [glob [file join $directory "*psf"]]]
	foreach psf $psflist {
		set tail [file tail $psf]
		set name [file rootname $tail]
		set mid [mol new $psf]
		puts [molinfo list]
		mol addfile [file join $directory $name.pdb] waitfor all
		puts [molinfo list]
		set asel [atomselect $mid "all"]
		$asel set beta 0
		animate write psf [file join $directory L.psf]
		animate write pdb [file join $directory L.pdb]
		puts [molinfo list]
		mdffi sim $asel -o [file join $directory grid.dx] -res 10 -spacing 1 -mol $mid
		puts [molinfo list]
		set finished 0
		set counter 0
		set fout [open [file join $directory "addenda.namd"] w ]
		puts $fout $namdextraconf
		close $fout
		set othersel [atomselect $mid "within 4 of occupancy > 0 and not withinbonds 3 of occupancy > 0"]
		set badbeta [atomselect $mid "beta > 0"]
		puts "$namdbin $namdargs [file join $directory minimize.namd] $namdextraconf"
		while { ! $finished } {
			::ExecTool::exec $namdbin $namdargs [file join $directory minimize.namd] > [file join $directory $name.log]
			animate delete all $mid
			incr counter
			incr finished
			mol addfile [file join $directory out.coor] type namdbin waitfor all
			set unfinished [vecsum [$asel get beta]]
			$asel set beta 0
			$asel set occupancy 0
			foreach bond [topo getbondlist -molid $mid] { 
				if { [measure bond $bond molid $mid] > 1.65 } {
					set ssel [atomselect $mid "same residue as index $bond"]
					$ssel set beta 1
					$ssel set occupancy 1
					set finished 0
					$ssel delete
				}
			}
			$othersel update
			set fout [open [file join $directory "addenda.namd"] w ]
			puts $fout $namdextraconf
			if { [$othersel num] } {
				puts $fout "colvars on\ncolvarsconfig colvars.conf\n"
			}
			close $fout
			if { ! $finished && $counter > 100 } {
				error "I'm sorry, I don't seem to be able to minimize $namd, and I tried 100 times!"
			}
			if { ! $finished } {
				$badbeta update
				mdffi sim $badbeta -o [file join $directory grid.dx] -res 10 -spacing 1 -mol $mid
			}
			if { $unfinished } {
				set finished 0
			}
			$asel writepdb [file join $directory L.pdb]
		}
		$asel writepdb [file join $directory $name.pdb]
		$asel delete
		$othersel delete
		$badbeta delete
		mol delete $mid
	}
}

#This was the original test reader.
# foreach library [glob libraries/*txt] {
# 	puts $library
# 	puts [file tail $library]
# 	set name [lindex [split [file tail $library]] 0]
# 	file mkdir $name
# 	set fp [open $library r]
# 	set file_data [read $fp]
# 	close $fp
# 	set data [split $file_data "\n"]
# 	set lcounter 0
# 	for { set lnum 0 } { $lnum < [llength $data] } { incr lnum } {
# 		set line [lindex $data $lnum]
# 		#puts $line
# 		#exit
# 		if { [string range $line 0 2] == "***" } {
# 			incr lnum 4
# 			set linedata [regexp -inline -all -- {\w+}  [lindex $data $lnum] ]
# 			set seqlength [lrange $linedata 0 2]
# 			set patches [lrange $linedata 3 end]
# 			puts $linedata
# 			set nummonomers [expr {int([vecsum $seqlength])}]
# 			set numpatches [expr {int([vecsum $patches])}]
# 			set monolist [list ]
# 			incr lnum 2
# 			for { set i 0 } { $i < $nummonomers } { incr i } {
# 				lappend monolist [lindex [regexp -inline -all -- {\w+}  [lindex $data $lnum] ] 1]
# 				incr lnum
# 			}
# 			incr lnum 3
# 			set patchlist [list ]
# 			for { set i 0 } { $i < $numpatches} { incr i } {
# 				lappend patchlist [regexp -inline -all -- {\w+}  [lindex $data $lnum] ]
# 				incr lnum
# 			}
# 			puts $monolist
# 			puts $patchlist
# 			puts $nummonomers
# 			puts $numpatches
# 			makelignin $monolist $patchlist L$lcounter $name
# 			incr lcounter
# 			#exit
# 		}
# 	}
# }

}
