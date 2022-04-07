mol delete all
set molId [ mol new /home/luca/tmp/N64M20T0.500Seed761r1.000na31.000e-04P0.000/view/sample1.pdb ]
topo clearbonds
mol delrep 1 $molId
mol rep Points 15.000000
mol color Occupancy
mol addrep $molId
mol rep Lines
mol addrep $molId

set M 20

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


set pdbFiles [ list /home/luca/tmp/N64M20T0.500Seed761r1.000na31.000e-04P0.000/view/sample2.pdb /home/luca/tmp/N64M20T0.500Seed761r1.000na31.000e-04P0.000/view/sample3.pdb /home/luca/tmp/N64M20T0.500Seed761r1.000na31.000e-04P0.000/view/sample4.pdb /home/luca/tmp/N64M20T0.500Seed761r1.000na31.000e-04P0.000/view/sample5.pdb /home/luca/tmp/N64M20T0.500Seed761r1.000na31.000e-04P0.000/view/sample6.pdb /home/luca/tmp/N64M20T0.500Seed761r1.000na31.000e-04P0.000/view/sample7.pdb /home/luca/tmp/N64M20T0.500Seed761r1.000na31.000e-04P0.000/view/sample8.pdb   ]

foreach file $pdbFiles {
    mol addfile $file
}




display resetview


