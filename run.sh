#!/bin/bash
source ~/.bash_profile
moose

#export MOOSE=/Users/andrewstershic/moose_projects/moose/modules/combined/combined-opt
export MOOSE=$HOME/moose_projects/moose/modules/tensor_mechanics/tensor_mechanics-opt

#/sw/bin/python2.7 generate_elasticity.py --style=shertzer 90wt_0bar.ftF2.npy 
#/sw/bin/python2.7 generate_elasticity.py --style=rahmoun 90wt_0bar.ftF4.npy 

#$MOOSE -i finite_strain_elastic_new_test.i
#aprepro -I=90wt_0bar.ftF4.elas test_solid_3d.i test_solid_3d.a
$MOOSE -i test_solid_3d.i

#*** ERROR ***
#ComputeFiniteStrainElasticStress can only be used with elasticity tensor materials that guarantee isotropic tensors.

#(1) convert to elasticity tensor
#note!!! shertzer's F is kanatani's N
#-F2 (shertzer)
#-F2 (rahmoun) # check which FT this assumes
#-F4 (rahmoun)
#-F2/4? (Cowin/Moesen)?  #check which FT this assumes

#(2) figure out how to import into moose (from file?, script-replacement?)
#create input file
#automate create input file

#(3) reconsider overall plan

#(4) contact srdjan
