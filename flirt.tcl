# FLexible InteRaction Tool (FlIRT)
# Written by D.K. Weber, 2014

namespace eval ::flirt:: {
    namespace export flirt
}

proc flirt { sel1 sel2 { args } } {
    set timei [expr double([clock clicks -milliseconds])]
    set mode 0
    if {[lsearch $args "-contact"] != -1} { set mode contacts }
    if {[lsearch $args "-hbond"] != -1} { set mode hbonds }
    if {[lsearch $args "-cationpi"] != -1} { set mode cationpi }
    if {[lsearch $args "-nonpolar"] != -1} { set mode nonpolar }
    if { $mode == 0 } { return [puts "No mode selected!"] }

    # Fix values and switches
    set c_cut [ flirt::variable_assign $args "-contact" 1 null ]
    set h_cut [ flirt::variable_assign $args "-hbond" 1 null ]
    set h_ang [ flirt::variable_assign $args "-hbond" 2 null ]
    set p_cut [ flirt::variable_assign $args "-cationpi" 1 null ]
    set p_tol [ flirt::variable_assign $args "-cationpi" 2 null ]
    set np_cut [ flirt::variable_assign $args "-nonpolar" 1 null ]
    set mol [ flirt::variable_assign $args "-mol" 1 top ]
    set start [ flirt::variable_assign $args "-start" 1 0 ]
    set stop [ flirt::variable_assign $args "-stop" 1 [molinfo $mol get numframes] ]
    set skip [ flirt::variable_assign $args "-skip" 1 1 ]
    set attr1 [ flirt::variable_assign $args "-attr1" 1 resid ]
    set attr2 [ flirt::variable_assign $args "-attr2" 1 resid ]
    set prefix [ flirt::variable_assign $args "-prefix" 1 $mode ]
    set exc_sel [ flirt::variable_assign $args "-exclude" 1 null ]
    set exc_atr [ flirt::variable_assign $args "-exclude" 2 index ]
    set betas [ flirt::switch_on $args "-betas" off "" ]
    set xy [ flirt::switch_on $args "-xy" off "Retreiving mean/dev $mode by $attr1 (sel1)" ]
    set xyz [ flirt::switch_on $args "-xyz" off "Retreiving $mode every frame by $attr1 (sel1)" ]
    set pairxy [ flirt::switch_on $args "-pairxy" off "Retreiving mean/dev $mode by $attr1 (sel1) to $attr2 (sel2) pair" ]
    set pairxyz [ flirt::switch_on $args "-pairxyz" off "Retreiving $mode every frame by $attr1 (sel1) to $attr2 (sel2) pair" ]
    set total [ flirt::switch_on $args "-total" off "Retreiving total $mode each frame" ]
    set zeros [ flirt::switch_on $args "-zeros" off "Writing zeros" ] 
    set combine [ flirt::switch_on $args "-combine" off "Combining donor and acceptor interactions for hbonds and cationpi" ]
    set self [ flirt::switch_on $args "-self" off "Computing hydrogen bonds for same selection." ]
    
    # Create output files
    if { $xy == on } { set o_xy [open "$prefix.xy.dat" w] }
    if { $xyz == on } { set o_xyz [open "$prefix.xyz.dat" w] }
    if { $pairxy == on } { set o_pxy [open "$prefix.pairxy.dat" w] }
    if { $pairxyz == on } { set o_pxyz [open "$prefix.pairxyz.dat" w] }
    if { $total == on } { set o_total [open "$prefix.total.dat" w] }

    # create list pointers fir sel1
    set s1list1 [$sel1 get index]
    set s1list2 [$sel1 get name]
    set s1list3 [$sel1 get resid]
    set s1list4 [$sel1 get resname]
    #set s1list5 [$sel1 get "x y z"]
    set s1list6 [$sel1 get $attr1]
    set s1list7 [$sel1 get $exc_atr]
    set c 0
    foreach index $s1list1 { set sel1ref($index) $c; incr c }

    # create list pointers fir sel2
    set s2list1 [$sel2 get index]
    set s2list2 [$sel2 get name]
    set s2list3 [$sel2 get resid]
    set s2list4 [$sel2 get resname]
    #set s2list5 [$sel2 get "x y z"]
    set s2list6 [$sel2 get $attr2]
    set s2list7 [$sel2 get $exc_atr]
    set c 0
    foreach index $s2list1 { set sel2ref($index) $c; incr c }

    # Determine if attributes are string or numeric. Makes outputs cleaner by lsort.
    set num_attrs { "residue" "resid" "index" "serial" }
    set str_attrs { "resname" "type" "name" "structure" }
    if {[lsearch $num_attrs $attr1] != -1} { set s1u [lsort -unique -integer -increasing $s1list6] }
    if {[lsearch $str_attrs $attr1] != -1} { set s1u [lsort -unique -ascii $s1list6] } 

    #####################################
    ############ Main loop ##############
    #####################################
    # Get number of frames for progess metre
    set nf [expr (($stop + 1) - $start) / $skip]
    if { $stop <= $start } { set $nf 1; puts "Measing one frame"}
    if { $stop > $start } { puts "Looping frames $start to $stop every $skip frames"}
    set f 0
    for {set j $start} {$j <= $stop} {incr j $skip} {
	$sel1 frame $j
	$sel2 frame $j
	# Put in blank list to prevent errors if no contacts detected
	set sel1contr {}
	set sel2contr {}

	### Contacts Mode ###
	if { $mode == "contacts" } {
	    set result [measure contacts $c_cut $sel1 $sel2]
	    set sel1contr [lindex $result 0]
	    set sel2contr [lindex $result 1]
	}

	### Nonpolar Mode ###
	if { $mode == "nonpolar" } {
	    set result [measure_nonpolar $np_cut $sel1 $sel2]
	    set sel1contr [lindex $result 0]
	    set sel2contr [lindex $result 1]
	}

	### H-bond Mode ###
	if { $mode == "hbonds" } {
	    set result1 [measure hbonds $h_cut $h_ang $sel1 $sel2]
	    set sel1contr [lindex $result1 0]
	    set sel2contr [lindex $result1 1]
	    if { $combine == on } {
		set result2 [measure hbonds $h_cut $h_ang $sel2 $sel1]
		foreach i [lindex $result2 1] { lappend sel1contr $i }
		foreach i [lindex $result2 0] { lappend sel2contr $i }
	    }
	}

	### cationpi Mode ###
	if { $mode == "cationpi" } {
	    set result1 [measure_cationpi $p_cut $p_tol $sel1 $sel2]
	    set sel1contr [lindex $result1 0]
	    set sel2contr [lindex $result1 1]
	    if { $combine == on } {
		set result2 [measure_cationpi $p_cut $p_tol $sel2 $sel1]
		foreach i [lindex $result2 1] { lappend sel1contr $i }
		foreach i [lindex $result2 0] { lappend sel2contr $i }
	    }
	}

	# Contacts sorting section
	# A lot of doubling up of code. Ideally have each
	# subection in own precudure, but arrays are
	# difficult to return.
	# Accumulate index hits into desired attribute
	set c 0

	# Update v1.1 - Initialise tc to prevent fail if frame has no contacts
	set tc 0
 
	set pairlist {}
	foreach i $sel1contr { 
	    set i2 [lindex $sel2contr $c]
	    set at1 [lindex $s1list6 $sel1ref($i)]
	    set at2 [lindex $s2list6 $sel2ref($i2)]
	    set check1 [lindex $s1list7 $sel1ref($i)]
	    set check2 [lindex $s2list7 $sel2ref($i2)]

	    # If exclusion option is set:
	    if { $exc_sel != "null" } {
		# Exclude intra-atrribute contacts
		if { $exc_sel == "intra" } { 
		    if { $check1 != $check2 } {
			incr tc
			incr r($at1)
			if {[info exist r($at1.$at2)] == 0} { lappend pairlist "$at1.$at2" }
			if {[info exist r2($at1.$at2)] == 0} { 
			    lappend pairlistx "$at1" 
			    lappend pairlisty "$at2"
			}
			incr r2($at1.$at2)
			incr r($at1.$at2) 
			incr c
		    }
		}
		# Exclude inter-residue contacts
		if { $exc_sel == "inter" } { 
		    if { $check1 == $check2 } {
			incr tc
			incr r($at1)
			if {[info exist r($at1.$at2)] == 0} { lappend pairlist "$at1.$at2" }
			if {[info exist r2($at1.$at2)] == 0} { 
			    lappend pairlistx "$at1" 
			    lappend pairlisty "$at2"
			}
			incr r2($at1.$at2)
			incr r($at1.$at2) 
			incr c
		    }
		}	
	    }
	    # End exclude option

	    # If the exclude option is not set
	    if { $exc_sel == "null" } {
		incr tc
		incr r($at1)
		if {[info exist r($at1.$at2)] == 0} { lappend pairlist "$at1.$at2" }
		if {[info exist r2($at1.$at2)] == 0} { 
		    lappend pairlistx "$at1" 
		    lappend pairlisty "$at2"
		}
		incr r2($at1.$at2)
		incr r($at1.$at2) 
		incr c
	    }
	}
	# End of Sorting section
	
	# Maintain running statistics
	# Selection 1 statistics
	foreach i $s1u {
	    if { $zeros == on } {
		if {[info exist r($i)] == 0} { set r($i) 0 }
	    }
	    if {[info exist r($i)] == 1} {
		incr xy_cs($i) $r($i)
		incr xy_cs2($i) [expr int(pow($r($i),2))]
		# Output time-evolution information so
		# it can be removed from memory
		if { $xyz == on } { puts $o_xyz "$j\t\t$i\t\t$r($i)" }
	    }
	}

	# Pairwise statistics
	foreach i $pairlist {
	    if { $pairxy == on } { 
		incr pxy_cs($i) $r($i)
		incr pxy_cs2($i) [expr int(pow($r($i),2))]
	    }
	    if { $pairxyz == on } { puts $o_pxyz "$j\t\t$i\t\t$r($i)" }
	}

	# Housekeeping at end of frame
	# Put space after frame block for pm3d plotting
	if { $xyz == on } { puts $o_xyz "" }
	# Output total number of contacts for the frame
	if { $total == on } { puts $o_total "$j\t\t$tc" }
	# Remove frame-specific pairlist and contact 
	# information from memory
	unset -nocomplain pairlist tc r
	incr f
	if { $nf >= 10 } { flirt::progress $f $nf }
    }
    puts "Done!"
    ############ End Main loop ##############

    # Output overall information
    if { $xy == on } {
	foreach i $s1u {
	    if {[info exist xy_cs($i)] == 1} {
		set istat [flirt::stats $xy_cs($i) $xy_cs2($i) $f]
		puts $o_xy "$i\t\t[lindex $istat 0]\t\t[lindex $istat 1]"
	    }}}

    if { $pairxy == on } {
	if {[info exist pairlistx] == 1} {
	    set c 0
	    foreach i $pairlistx {
		set yval [lindex $pairlisty $c]
		set istat [flirt::stats $pxy_cs($i.$yval) $pxy_cs2($i.$yval) $f]
		puts $o_pxy "$i\t\t$yval\t\t[lindex $istat 0]\t\t[lindex $istat 1]"
		incr c
	    }}}
    
    # Assign beta values
    set c 0
    set betalist {}
    if { $betas == on } { 
	foreach i $s1list1 {
	    set at1 [lindex $s1list6 $c]
	    if {[info exist xy_cs($at1)] == 1} {
		lappend betalist [lindex [flirt::stats $xy_cs($at1) $xy_cs2($at1) $f] 0]
	    } else { lappend betalist 0 }
	    incr c
	}
	$sel1 set beta $betalist
    }

    # Final housekeeping
    if { $xy == on } { close $o_xy }
    if { $xyz == on } { close $o_xyz }
    if { $pairxy == on } { close $o_pxy }
    if { $pairxyz == on } { close $o_pxyz }
    if { $total == on } { close $o_total }

    # Output time taken
    set timen [expr double([clock clicks -milliseconds])]
    return [puts "[format %.3f [expr ($timen - $timei) / 60000 ]] minutes"]
}


################################################################################
################################ Libraries #####################################
################################################################################

# Cation-pi library
proc flirt::check_cationpi { resname name } {
    # Library of pi systems
    set pisystems(TYR) { CG CD2 CE2 CZ CE1 CD1 }
    set pisystems(PHE) { CG CD2 CE2 CZ CE1 CD1 }
    set pisystems(TRP) { CE2 CD2 CE3 CZ3 CH2 CZ2 }

    # Library of cations
    set cations(DPC) { N4 }
    set cations(ASM) { N1 }
    set cations(ARG) { CZ }
    set cations(LYS) { NZ }
    set cations("Na+") { "Na+" }
    set cations("K+") { "K+" }

    # Check values
    if { [info exist pisystems($resname)] == 1 } {
	if { [lsearch $pisystems($resname) $name] != -1} {
	    return "pi"
	}}

    if {[info exist cations($resname)] == 1} {
	if {[lsearch $cations($resname) $name] != -1} {
      return "cation"
	}}
    return "null"
}

# Library of non-polar atom names
proc flirt::check_hydrophobic { resname name } {
    set hphobe(ALA) { CA CB HA HB1 HB2 HB3 }
    set hphobe(ARG) { CA CB CG HA HB2 HB3 HG2 HG3 CD HD2 HD3 }
    set hphobe(ASN) { CA CB HA HB2 HB3 }
    set hphobe(ASP) { CA CB HA HB2 HB3 }
    set hphobe(CYS) { CA CB HA HB2 HB3 HG }
    set hphobe(GLN) { CA CB CG HA HB2 HB3 HG2 HG3 }
    set hphobe(GLU) { CA CB CG HA HB2 HB3 HG2 HG3 }
    set hphobe(GLY) { CA HA2 HA3 }

    # Include AMBER protonation states HID (HD2 proton)
    # HIE (HE1 protein) and HIP (both HD2 and HE1 protons)
    set hphobe(HID) { CA CB CG CD2 HA HB2 HB3 HD2 CE1 HE1 }
    set hphobe(HIE) { CA CB CG CD2 HA HB2 HB3 HD2 CE1 HE1 }
    set hphobe(HIP) { CA CB CG CD2 HA HB2 HB3 HD2 CE1 HE1 }
    set hphobe(HIS) { CA CB CG CD2 HA HB2 HB3 HD2 CE1 HE1 }
    set hphobe(ILE) { CA CB CG1 CG2 CD1 HA HB HG21 HG22 HG23 HG12 HG13 HD11 HD12 HD13 }
    set hphobe(LEU) { CA CB CG CD1 CD2 HA HB2 HB3 HG HD21 HD22 HD23 HD11 HD12 HD13 }
    set hphobe(LYS) { CA CB CG CD HA HB2 HB3 HG2 HG3 HD2 HD3 CE HE2 HE3 }
    set hphobe(MET) { CA CB CG HA HB2 HB3 HG2 HG3 CE HE1 HE2 HE3 }
    set hphobe(PHE) { CA CB CG CD1 CD2 CE1 CE2 CZ HA HB2 HB3 HE1 HE2 HD1 HD2 HZ }
    set hphobe(PRO) { CA CB CG CD HA HB2 HB3 HG2 HG3 HD2 HD3 }
    set hphobe(SER) { CA CB HA HB2 HB3 }
    set hphobe(THR) { CA CB CG2 HA HB HG21 HG22 HG23 }
    set hphobe(TRP) { CA CB CG CD1 CD2 CE3 CZ2 CZ3 HA HB2 HB3 HD1 HE3 HZ2 HZ3 CE2 CH2 HH2 }	
    set hphobe(TYR) { CA CB CG CD1 CD2 CE1 CE2 CZ HA HB2 HB3 HE1 HE2 HD1 HD2  }
    set hphobe(VAL) { CA CB CG1 CG2 HA HB HG11 HG12 HG13 HG21 HG22 HG23 }

    # User defined
    set hphobe(DPC) { 
	C1 H1A H1B H1C C2 H2A H2B H2C C3 H3A H3B H3C C5 H5A H5B C6 H6A H6B C12 
	H12A H12B C13 H13A H13B C14 H14A H14B C15 H15A H15B C16 H16A H16B C17 H17A H17B C18 
	H18A H18B C19 H19A H19B C20 H20A H20B C21 H21A H21B C22 H22A H22B C23 H23A H23B H23C 
    }

    set hphobe(ASM) { 
	C2 H2A H2B H2C C3 H3A H3B H3C C4 H4A H4B H4C C5 H5A H5B C6 H6A H6B C12 
	H12A H12B C13 H13A C14 H14A C16 H16A C17 H17A C18 H18A H18B C19 H19A H19B C20 H20A H20B 
	C21 H21A H21B C22 H22A H22B C23 H23A H23B C24 H24A H24B C25 H25A H25B C26 H26A H26B C27 
	H27A H27B C28 H28A H28B C29 H29A H29B C30 H30A H30B H30C C34 H34A H34B H34C 
    }

    if {[info exist hphobe($resname)] == 1} {
	if {[lsearch $hphobe($resname) $name] != -1} {
	    return "hphobe"
	}}
    return "null"
}

################################################################################
############# New Procedure Introduced into Global Namespace ###################
################################################################################

# Selection 1 contains pi-system and sel2 contains cations
proc measure_cationpi { p_cut p_tol sel1 sel2 } {
    set sel1index [$sel1 get index]
    set sel1resid [$sel1 get resid]
    set sel1rname [$sel1 get resname]
    set sel1aname [$sel1 get name]
    set coords1 [$sel1 get "x y z"]
    set sel2index [$sel2 get index]
    set sel2resid [$sel2 get resid]
    set sel2rname [$sel2 get resname]
    set sel2aname [$sel2 get name]
    set coords2 [$sel2 get "x y z"]
    
    # Make reference arrays for selection1
    set c 0
    foreach i $sel1index {
	set aname [lindex $sel1aname $c]
	set resid [lindex $sel1resid $c]
	set rname [lindex $sel1rname $c]
	if { [flirt::check_cationpi $rname $aname] == "pi" } {
	    # Make array for indices making up pi-system
	    set sel1ref($i) $c
	    lappend resid_pi($resid) $i
	}
	incr c
    }
    # Same for cation in selection 2
    set c 0
    foreach i $sel2index {
	set aname [lindex $sel2aname $c]
	set resid [lindex $sel2resid $c]
	set rname [lindex $sel2rname $c]
	if { [flirt::check_cationpi $rname $aname] == "cation" } {
	    # Make array for indices making up pi-system
	    set sel2ref($i) $c
	    lappend resid_cation($resid) $i
	}
	incr c
    }

    set sel1contr {}
    set sel2contr {}
    # All cotacts that show up within p_cut - filtering step to
    # limit pairwise distance measurements to make.
    set result [measure contacts $p_cut $sel1 $sel2]
    set s1contr [lindex $result 0]
    set s2contr [lindex $result 1]
    set c 0
    ### Find unique pairs of indices
    foreach i $s1contr {
	set i2 [lindex $s2contr $c]
	if {[info exist sel1ref($i)] != 0} { 
	    if {[info exist sel2ref($i2)] != 0} { 
		set rid1 [lindex $sel1resid $sel1ref($i)]
		incr cppairs($rid1.$i2)
    	
		# If the interaction is up to 6, then
		# check if distances are within tolarence
		if { $cppairs($rid1.$i2) == 6 } {
		    set coord2 [lindex $coords2 $sel2ref($i2)]
		    foreach atom $resid_pi($rid1) {
			set coord1 [lindex $coords1 $sel1ref($atom)]
			lappend tempresult [flirt::dist $coord1 $coord2]
		    }
		    set diff [expr [flirt::max $tempresult] - [flirt::min $tempresult]]
		    set mindist [flirt::min $tempresult]
		    if { $diff <= $p_tol } {
			foreach atom $resid_pi($rid1) {
			    lappend sel1contr $atom
			    lappend sel2contr $i2
			}
		    }
		}
	    }}
	incr c
    }
    return [list $sel1contr $sel2contr]
}


proc measure_nonpolar { cut sel1 sel2 } {
    set sel1index [$sel1 get index]
    set sel1rname [$sel1 get resname]
    set sel1aname [$sel1 get name]
    set sel2index [$sel2 get index]
    set sel2rname [$sel2 get resname]
    set sel2aname [$sel2 get name]
  
    # Make reference arrays for selection1
    set c 0
    foreach i $sel1index {
	set aname [lindex $sel1aname $c]
	set rname [lindex $sel1rname $c]
	if { [flirt::check_hydrophobic $rname $aname] == "hphobe" } {
	    set sel1ref($i) $c
	}
	incr c
    }

    set c 0
    foreach i $sel2index {
	set aname [lindex $sel2aname $c]
	set rname [lindex $sel2rname $c]
	if { [flirt::check_hydrophobic $rname $aname] == "hphobe" } {
      set sel2ref($i) $c
	}
	incr c
    }

    set sel1contr {}
    set sel2contr {}
    
    set result [measure contacts $cut $sel1 $sel2]
    set s1contr [lindex $result 0]
    set s2contr [lindex $result 1]
    set c 0
    ### Find unique pairs of indices
    foreach i $s1contr {
	set i2 [lindex $s2contr $c]
	if {[info exist sel1ref($i)] != 0} { 
	    if {[info exist sel2ref($i2)] != 0} { 
		lappend sel1contr $i
		lappend sel2contr $i2
	    }}
	incr c
    }
    return [list $sel1contr $sel2contr]
}


################################################################################
################################ Dependencies ##################################
################################################################################

proc flirt::switch_on { args flag default message } {
    set value $default
    if { [ lsearch $args $flag ] != -1 } { 
	set value on
	puts $message
    }
    return $value
}

proc flirt::variable_assign { args flag input_index default } {
    set value $default
    if { [ lsearch $args $flag ] != -1 } { 
	set value [lindex $args [expr ([lsearch $args $flag ] + $input_index)]]
    } 
    return $value
}

# Modified from http://wiki.tcl.tk/16939
proc flirt::progress {cur tot} {
    if {$cur % ($tot/10)} { return }
    # set to total width of progress bar
    set total 100
    set percent [expr {100.*$cur/$tot}]
    set val (\ [format "%6.2f%%" $percent]\ )
    set str "[expr {round($percent*$total/100)}]% "
    puts -nonewline $str
}

proc flirt::stats { csum csum2 num } {
    set csum [expr double($csum)]
    set csum2 [expr double($csum2)]
    set avg [format "%.4f" [expr $csum / $num]]
    if { $num <= 1 } { return [ list $avg nan ] }
    set std [format "%.4f" [expr sqrt((($num*$csum2)-pow($csum,2)) / ($num*($num-1)))]]
    return [list $avg $std]
}

# Taken from tcl math package
proc flirt::min { list_in } {
    set min {}
    foreach i $list_in {
     set i [expr {double($i)}]
	if { $min == {} || $i < $min } { set min $i }
    }
    return $min
}
proc flirt::max { list_in } {
    set max {}
    foreach i $list_in {
	set i [expr {double($i)}]
	if { $max == {} || $i > $max } { set max $i }
    }
    return $max
}

proc flirt::dist { com1 com2 } {
    set l [veclength [vecsub $com1 $com2]]
    return $l
}


