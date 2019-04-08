#!/bin/bash -l
grep -n 'Timestep: ' prod-emt-test.xyz > timesteps_prodxyz
sed -i 's/:/ /' timesteps_prodxyz
cat timesteps_prodxyz | awk 'BEGIN{n=0; m=0}; {m=$1-n; n=$1; print $1, $4, m}' | tail -n+2 > line_numbers
