mol delete all
set molId [ mol new __FIRSTFILE__ ]
topo clearbonds
mol delrep 1 $molId
mol rep Points 15.000000
mol color Occupancy
mol addrep $molId
mol rep Lines
mol addrep $molId

set M __M__

set sel [atomselect top all]
set NM [ $sel num]
set N [ expr $NM / [expr $M + 1] ]
set k 0
for {set i 0} {$i < $N} {incr i} {
    for {set t 0} {$t < $M} {incr t} {
        topo addbond $k [expr $k+1]
        incr k
    }
    incr k
}

mol modselect 0 $molId x < 9999
mol modselect 1 $molId x < 9999
mol modselect 2 $molId x < 9999

mol selupdate 2 $molId 1
mol colupdate 2 $molId 1

mol selupdate 1 $molId 1
mol colupdate 1 $molId 1


set pdbFiles [ list __PDBFILES__   ]

foreach file $pdbFiles {
    mol addfile $file
}


__PBCSCRIPT__

display resetview


