#!/bin/bash

gcc -o3 turnover.c -lm -o turnover.ce

# Parameters: [mode] [nsims] [mass] [seglength] [branchprob] [het] [pop] ([nseed] [halo] [p] [q] [rho] [clustersize])
# If mode is --snapshots, one must also pass [nseed] [halo] [p] [q] and [rho] [clustersize]
./turnover.ce --snapshots 1 50 0.01 0.02 100 0.5 4 0.0 1.0 1.0 0.25 25 > tmps1 &  # 100-4-50-1.00-1.00-0.00-0.25
./turnover.ce --snapshots 1 50 0.01 0.02 100 0.5 16 0.0 1.0 1.0 0.25 25 > tmps2 & # 100-16-50-1.00-1.00-0.00-0.25
./turnover.ce --snapshots 1 50 0.01 0.02 100 0.5 64 0.0 1.0 1.0 0.0 25 > tmps3 &

./turnover.ce --snapshots 1 50 0.01 0.02 100 0.5 4 0.1 0.5 0.5 0.25 25 > tmps4 &
./turnover.ce --snapshots 1 50 0.01 0.02 100 0.5 16 0.1 0.5 0.5 0.25 25 > tmps5 &
./turnover.ce --snapshots 1 50 0.01 0.02 100 0.5 64 0.1 0.5 0.5 0.25 25 > tmps6 &

#./turnover.ce --simulate 1500 50 0.01 0.02 100 > tmp1 &
