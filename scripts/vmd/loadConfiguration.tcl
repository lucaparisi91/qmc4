mol delete all
set molId [ mol new /home/luca/source/qmc4/PPA/run/configurations2/particles4.pdb ]
topo clearbonds
mol delrep 1 $molId
mol rep Points 6.000000
mol color Occupancy
mol addrep $molId
mol rep Lines
mol addrep $molId

set M 10
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
display resetview