#!/usr/bin/env bash

# example running seir_simulator
./seir_simulator folder=./dec_15_EPI iu=0.92 fu=0.70 rd=15 R0=10 ts=7 force=seasonal T=140 runs=1000

#To do: 
./seir_simulator folder=./ iu=0.92 fu=0.70 rd=15 R0=10 ts=7 force=term T=140 runs=1 

