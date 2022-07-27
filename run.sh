#!/bin/bash

gcc -o3 turnover.c -lm -o network.ce

# Parameters: [mode] [nsims] [mass] [seglength] [branchprob] [het] [pop] ([nseed] [halo] [p] [q] [rho] [clustersize])
# If mode is --snapshots, one must also pass [nseed] [halo] [p] [q] and [rho] [clustersize]
./turnover.ce --snapshots 1 50 0.01 0.02 0.5 100 4 0.0 1.0 0.0 0.25 25 > tmps1 &
./turnover.ce --snapshots 1 50 0.01 0.02 0.5 100 16 0.0 1.0 1.0 0.25 25 > tmps2 &
./turnover.ce --snapshots 1 50 0.01 0.02 0.5 100 64 0.0 1.0 1.0 0.0 25 > tmps3 &

./turnover.ce --snapshots 1 50 0.01 0.02 0.5 100 4 0.1 0.5 0.5 0.25 > tmps4 &
./turnover.ce --snapshots 1 50 0.01 0.02 0.5 100 16 0.1 0.5 0.5 0.25 > tmps5 &
./turnover.ce --snapshots 1 50 0.01 0.02 0.5 100 64 0.1 0.5 0.5 0.25 > tmps6 &
