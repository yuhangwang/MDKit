# Name: GKAAM.tcl ("Generalized Keep Atoms Away from Membrane-core")
# Purpose: to keep user-defined atoms outside of the space between two membrane layers
# Author: Steven(Yuhang) Wang (stevenwaura@gmail.com)

# Main redesign: allow user to specify any type of atoms to be pushed (or not to be pushed)
#                using the beta field of a PDB file.
# Date: 08/31/2016
# Copyright license: MIT/X11
#
# Example usage:
######################################################################################################
# tclForces                        on
# set betaLabelPDB                 $KAAMrefPDB ;# pdb with beta column labeled with numbers
# set atomLipidUpLayerBV           1   ;# Beta value(BV) for reference atom in membrane upper layer
# set atomLipidDownLayerBV         2   ;# Beta value(BV) for refernce atom in membrane lower layer
# set atomToPushBV                 3   ;# Beta value(BV) for atoms to be pushed
# set printLabelAtom               off ;# to print atoms you labeled?
# set pushForceConstant            0.5 ;# kcal/(mol*Angstrom)
# set targetAtomListUpdateFreq     100 ;# steps
# set lipidAtomListUpdateFreq      5000 ;# steps
# tclForcesScript                  GKAAM.tcl
######################################################################################################


# tolerance for number comparison
set TOLERANCE  [expr 1.0E-7]
##================================================================================
## <<< Step 1: define index list for atoms to be affected >>>
##================================================================================
set lipidAtomUpperLayerIndices {} ;# indices for atoms to be confined in the top/outside the cell
set lipidAtomLowerLayerIndices {} ;# indices for atoms to be confined in the bottom/inside the cell
set atomsCanBePushedIndices {} ;# indices for atoms that can be pushed if the atoms enters the forbidden region
set atomID 1                  ;# in NAMD, atom IDs are 1-based numbers

set inPDB [open $betaLabelPDB]        ;# opening the PDB file

##================================================================================
## <<< Read PDB and find indices for atoms to be confined >>>
##================================================================================
foreach line [split [read $inPDB] \n] {
  set tempRecordType [string range $line 0 3] ;# record type, e.g. "ATOM"
  set tempAtomIndex  [string trim [string range $line  6 10] ] ;# no need to get atom index from PDB, because atom index will be counted by "$atomID"
  set tempAtomName   [string trim [string range $line 13 15] ]
  set tempResName    [string trim [string range $line 17 20] ]
  set tempResidueID  [string trim [string range $line 22 25] ]
  set tempZcoord     [string trim [string range $line 46 53] ]
  set tempSegName    [string trim [string range $line 72 75] ]
  set tempBeta       [string trim [string range $line 60 65] ]
  set tempOccup      [string trim [string range $line 54 59] ]


  ##--------------------------------------------------------------------------------
  if { ([string equal $tempRecordType {ATOM}] || [string equal $tempRecordType {HETA}] ) } { 
    set val_up [expr $tempBeta - $atomLipidUpLayerBV ]
    set val_down [expr $tempBeta - $atomLipidDownLayerBV ]
    set val_push [expr $tempBeta - $atomToPushBV]
    if { [expr abs($val_up)] < $TOLERANCE } {
        lappend lipidAtomUpperLayerIndices $atomID 
        addatom $atomID
    } elseif { [expr abs($val_down)] < $TOLERANCE } {
        lappend lipidAtomLowerLayerIndices $atomID
        addatom $atomID
    } elseif { [expr abs($val_push)] < $TOLERANCE } {
        lappend atomsCanBePushedIndices $atomID
        addatom $atomID
    }
    incr atomID
  }
  ##--------------------------------------------------------------------------------
}
close $inPDB



###################################################################
## <<< Start >>>
###################################################################
print "Starting GKAAM.tcl"

if {([llength $lipidAtomUpperLayerIndices] > 0) && ([llength $lipidAtomLowerLayerIndices] > 0)} {
    set FLAG_PUSH 1
} else {
    print "WARNING: membrane has not been detected"
    print "HINT: you are likely to have specified the wrong residue name for your lipids"
    exit
    set FLAG_PUSH 0
}



# initialize printing counter (independent on step counter)

set pushCounter   $targetAtomListUpdateFreq
set updateCounter $lipidAtomListUpdateFreq
set printCounter  0
set atomToPushUp   {}
set atomToPushDown {}
set activeAtomList [concat $lipidAtomUpperLayerIndices $lipidAtomLowerLayerIndices $atomsCanBePushedIndices]


set zForceUp   [expr $pushForceConstant]
set zForceDown [expr -$pushForceConstant]
set zUpperLayer 0
set zLowerLayer 0

###################################################################
# this procedure is executed at each time step
###################################################################

print ">>> GKAAM.tcl: starting calculating pushing forces..."


proc apply_force {atom_indices force} {
    foreach i $atom_indices {
      set f [list 0.0 0.0 $force ]
      addforce $i $f
    }
}


proc reset_active_atoms {indices} {
    clearconfig

    foreach atom $indices {
        addatom  $atom
    }
}

proc compute_average_z {indices} {
  upvar 1 coords coords
  set avg 0
  foreach i $indices {
    foreach {x y z} $coords($i) {break}
    set avg [expr $avg + $z]
  }

  set n [llength $indices]
  set avg [expr $avg / double($n)]
  return $avg
}


proc calc_membrane_bounds {lipidAtomUpperLayerIndices lipidAtomLowerLayerIndices} {
    upvar 1 coords coords ;# I use upvar because Tcl doesn't allow array in proc arguments
    set zUpperLayer [compute_average_z $lipidAtomUpperLayerIndices]
    set zLowerLayer [compute_average_z $lipidAtomLowerLayerIndices]

    return [list $zUpperLayer $zLowerLayer]
}

proc find_target_atoms {zUpperLayer zLowerLayer atomsCanBePushedIndices} {
    upvar 1 coords coords ;# I use upvar because Tcl doesn't allow array in proc arguments
    ## reset atomToPushUp/Down lists
    set atomToPushUp   {}
    set atomToPushDown {}

    set halfMembraneThickness [expr  ($zUpperLayer - $zLowerLayer)/2.0 ]

    foreach index $atomsCanBePushedIndices {
        foreach {x y z} $coords($index) { break }
        if { $z >= $zLowerLayer && $z <= $zUpperLayer } {
          if { [expr $zUpperLayer - $z] <= $halfMembraneThickness } {
              lappend atomToPushUp $index
          } else {
              lappend atomToPushDown $index
          }
        }
    }
    
    return [list $atomToPushUp $atomToPushDown]
}

proc calcforces {} {
    global zForceUp zForceDown
    global lipidAtomListUpdateFreq targetAtomListUpdateFreq
    global atomsCanBePushedIndices lipidAtomUpperLayerIndices lipidAtomLowerLayerIndices
    global atomToPushUp atomToPushDown
    global zUpperLayer zLowerLayer
    global activeAtomList
    global FLAG_PUSH 


    if {$FLAG_PUSH == 0} { 
      return 
    }

    set step [getstep]

    ##================================================================================
    ## Redefine the list of active atoms (whose coordinates can be requested)
    ##================================================================================    
    if { [expr $step % $lipidAtomListUpdateFreq] == 0 } {
      # reset active atoms
      reset_active_atoms [concat $lipidAtomUpperLayerIndices $lipidAtomLowerLayerIndices $atomsCanBePushedIndices]
      set activeAtomList [concat $lipidAtomUpperLayerIndices $lipidAtomLowerLayerIndices $atomsCanBePushedIndices]
    }

    loadcoords coords
    
    ##================================================================================
    ## <<<< Calculate average z position of the upper & lower membrane layers >>>
    ##================================================================================    
    if { [expr $step % $lipidAtomListUpdateFreq] == 0 } {
      set bounds [calc_membrane_bounds $lipidAtomUpperLayerIndices $lipidAtomLowerLayerIndices] 
      set zUpperLayer [lindex $bounds 0]
      set zLowerLayer [lindex $bounds 1]
    }

    ##================================================================================
    ## Find atoms to be pushed up/down
    ##================================================================================
    if { [expr $step % $targetAtomListUpdateFreq] == 0 } {
      set atomToPush [find_target_atoms $zUpperLayer $zLowerLayer $activeAtomList]
      set atomToPushUp   [lindex $atomToPush 0]
      set atomToPushDown [lindex $atomToPush 1]
      set activeAtomList  [concat $atomToPushUp $atomToPushDown]
    }

    apply_force $atomToPushUp   $zForceUp 
    apply_force $atomToPushDown $zForceDown
  
  return 
}
