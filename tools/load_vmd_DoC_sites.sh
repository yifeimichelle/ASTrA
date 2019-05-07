head -n 9 all.lammpstrj > DoC.lammpstrj
n=$(grep 'BF4\|C1' DoCIndices-all.out | wc -l)
grep 'BF4\|C1' DoCIndices-all.out >> DoC.lammpstrj
sed -i 's/C1/1/' DoC.lammpstrj
sed -i 's/BF4/2/' DoC.lammpstrj
sed -i "4s/.*/${n}/" DoC.lammpstrj
sed -i "9s/.*/ITEM: ATOMS type x y z DoC charge numcarbons avgdist index step/" DoC.lammpstrj

cella=$1
cellb=$2
cellc=$3
alpha=${4:-90}
beta=${5:-90}
gamma=${6:-90}

cat > VMD_visualize_DoC_sites<<EOF
set data [topo readlammpsdata data_symmz.data]
mol delrep 0 top

set sel [atomselect top "type 9"]
set cmol [::TopoTools::selections2mol \$sel]
mol delrep 0 top
mol representation vdw
mol selection {resname anode}
mol color colorid 6
mol addrep top

mol delete \$data

global env
set env(LAMMPSREMAPFIELDS) vx=charge,vy=DoC,vz=avgdist
color scale method BWR
set mol [mol new DoC.lammpstrj type lammpstrj autobonds off]

set nf [molinfo \$mol get numframes]
set sel [atomselect \$mol all]
for {set i 0} {\$i < \$nf} {incr i} {
  \$sel frame \$i
  \$sel set user [\$sel get vx]
  \$sel set vx 0.0
  \$sel set user2 [\$sel get vy]
  \$sel set vy 0.0
  \$sel set user3 [\$sel get vz]
  \$sel set vz 0.0
  #echo [\$sel get user]
  #echo [\$sel get vx]
}
\$sel delete
unset sel
unset env(LAMMPSREMAPFIELDS)
mol delrep 0 top

# set midlist {}
# lappend midlist \$cmol
# lappend midlist \$mol
# set merged [::TopoTools::mergemols \$midlist]

set c [lindex [lindex [pbc get] 0] 2]
set mvlist {}
lappend mvlist 0
lappend mvlist 0
lappend mvlist [expr \$c / 2]
set sel [atomselect \$mol all]
\$sel moveby \$mvlist
set sel [atomselect \$cmol all]
\$sel moveby \$mvlist

set halfz [expr \$c/2]
set sel [atomselect \$cmol "z>\$halfz"]
\$sel set resname "anode"
set sel [atomselect \$cmol "z<\$halfz"]
\$sel set resname "cathode"
set sel [atomselect \$mol "z>\$halfz"]
\$sel set resname "anode"
set sel [atomselect \$mol "z<\$halfz"]
\$sel set resname "cathode"

mol delrep 0 top
mol representation points

mol color User
mol selection {resname anode}
mol addrep top
mol scaleminmax top 0 -0.4 0.4
mol drawframes top 0 {now}

mol color User2
mol selection {resname anode}
mol addrep top
mol scaleminmax top 1 0.1 0.45
mol drawframes top 1 {now}

mol color User3
mol selection {resname anode}
mol addrep top
mol scaleminmax top 2 5 6
mol drawframes top 2 {now}

# set sel [atomselect \$merged "z>\$half"]
# set anode [::TopoTools::selections2mol \$sel]
# set sel [atomselect \$merged "z<\$half"]
# set cathode [::TopoTools::selections2mol \$sel]

set pbclist {}
lappend pbclist $cella
lappend pbclist $cellb
lappend pbclist $cellc
lappend pbclist $alpha
lappend pbclist $beta
lappend pbclist $gamma

pbc set \$pbclist -molid \$mol
pbc set \$pbclist -molid \$cmol
pbc wrap -molid \$mol
pbc wrap -molid \$cmol

EOF