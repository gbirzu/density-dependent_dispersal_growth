#!/bin/sh
for j in {100..200}
do
  ./bin_ddrd_fluctuations -r 0.01 -B 0.0 -N 10 -m 0.01 -A 3.0 -T 10000 -n $j
done
mv hetero_N10_r0.01_m0.01_A3.0_B0.0_* hetero/
mv velocity_N10_r0.01_m0.01_A3.0_B0.0* velocity/
mv profile_N10_r0.01_m0.01_A3.0_B0.0* profile/
