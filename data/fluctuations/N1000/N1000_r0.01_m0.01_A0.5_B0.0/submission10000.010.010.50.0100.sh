#!/bin/sh
for j in {0..100}
do
  ./bin_ddrd_fluctuations -r 0.01 -B 0.0 -N 1000 -m 0.01 -A 0.5 -T 10000 -n $j
done
mv hetero_N1000_r0.01_m0.01_A0.5_B0.0_* hetero/
mv velocity_N1000_r0.01_m0.01_A0.5_B0.0* velocity/
mv profile_N1000_r0.01_m0.01_A0.5_B0.0* profile/
