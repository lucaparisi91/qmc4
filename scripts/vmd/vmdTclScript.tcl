mol delete all
set molId [ mol new /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample10.pdb ]
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


set pdbFiles [ list /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample11.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample12.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample13.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample14.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample15.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample16.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample17.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample18.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample19.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample1.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample20.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample21.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample22.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample23.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample24.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample25.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample26.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample27.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample28.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample29.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample2.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample30.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample31.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample32.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample33.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample34.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample35.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample36.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample37.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample38.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample39.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample3.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample40.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample41.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample42.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample43.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample44.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample45.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample46.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample47.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample48.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample49.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample4.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample50.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample51.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample52.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample53.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample54.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample55.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample56.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample57.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample58.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample59.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample5.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample60.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample61.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample62.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample63.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample64.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample65.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample66.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample67.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample68.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample69.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample6.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample70.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample71.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample72.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample7.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample8.pdb /home/luca/data/droplet-finite-temperature/semiCanonical/test/configurations/sample9.pdb   ]

foreach file $pdbFiles {
    mol addfile $file
}



pbc set { 2.889114400454439 2.889114400454439 2.889114400454439 } -all
pbc wrap -all
pbc join connected -all
pbc box_draw
    

display resetview


