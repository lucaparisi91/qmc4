mol delete all
set molId [ mol new /home/luca/source/qmc4/PPA/run/configurations/sample49.pdb ]
topo clearbonds
mol delrep 1 $molId
mol rep Points 6.000000
mol color Occupancy
mol addrep $molId
mol rep Lines
mol addrep $molId

set M 100
set sel [atomselect top all]
set NM [ $sel num]
set N [ expr $NM / [expr 100 + 1] ]
set k 0
for {set i 0} {$i < $N} {incr i} {
    for {set t 0} {$t < $M} {incr t} {
        topo addbond $k [expr $k+1]
        incr k
    }
    incr k
}