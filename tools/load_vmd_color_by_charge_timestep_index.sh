#!/bin/bash -l

timestep=$1
index=$2
window=${3:-0}

if [ -e ${timestep}_frame ]; then
	frame=$(cat ${timestep}_frame)
else
	n=$(LC_ALL=C grep -n "^${timestep}$" all.lammpstrj)
	n=${n%%:*}
	frame=$(head -n $n all.lammpstrj | LCALL=C grep -c '^ITEM: TIMESTEP$')
	frame=$((frame-1))
	echo $frame > ${timestep}_frame
fi

first=$((frame-window))
last=$((frame+window))
step=1

cat > VMD_color_by_charge_timestep <<EOF
global env
set env(LAMMPSREMAPFIELDS) vx=q
color scale method BWR

topo readlammpsdata data_symmz.data
#mol addfile ztc_electrode_cell.psf
set mol [mol addfile all.lammpstrj type lammpstrj first $first last $last step $step waitfor all]
# insert "mol addfile command here, if needed"
set nf [molinfo \$mol get numframes]
set sel [atomselect \$mol all]
for {set i 0} {\$i < \$nf} {incr i} {
  \$sel frame \$i
  \$sel set user [\$sel get vx]
  \$sel set vx 0.0
  #echo [\$sel get user]
  #echo [\$sel get vx]
}
\$sel delete
unset sel
unset env(LAMMPSREMAPFIELDS)
mol delrep 0 top
mol color User
mol selection {all}
mol addrep top
#mol selupdate 0 top 0
#mol colupdate 0 top 1
mol scaleminmax top 0 -0.010000 0.010000
mol drawframes top 0 {now}

mol representation VDW
mol addrep top
mol scaleminmax top 1 -0.010000 0.010000
mol modselect 1 top "type 9 and within 6.3 of index $index"

mol representation VDW
mol addrep top
mol modselect 2 top "type 3 4 5 and within 7.7 of index $index"
mol modcolor 2 top ColorID 3

mol representation VDW
mol addrep top
mol modselect 3 top "type 1 and within 7.7 of index $index"
mol modcolor 3 top ColorID 7

# solvent
mol representation VDW
mol addrep top
mol modselect 4 top "type 2 6 7 and within 6.5 of index $index"
mol modcolor 4 top ColorID 9
mol modmaterial 4 top Transparent

EOF

vmd -e VMD_color_by_charge_timestep
