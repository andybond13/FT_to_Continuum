#!/bin/bash
source ~/.bash_profile
moose

#export MOOSE=/Users/andrewstershic/moose_projects/moose/modules/combined/combined-opt
export MOOSE=/Users/andrewstershic/moose_projects/moose/modules/tensor_mechanics/tensor_mechanics-opt

#/sw/bin/python2.7 generate_elasticity.py --style=shertzer 90wt_0bar.ftF2.npy 
/sw/bin/python2.7 generate_elasticity.py --style=rahmoun 90wt_0bar.ftF4.npy 

#$MOOSE -i finite_strain_elastic_new_test.i
#$MOOSE -i test_solid_3d.i

#create geometry
#automate create geometry

#*** ERROR ***
#ComputeFiniteStrainElasticStress can only be used with elasticity tensor materials that guarantee isotropic tensors.

#read FT
#convert to elasticity tensor
#-F2 (shertzer)
#-F2 (rahmoun)
#-F4 (rahmoun)
#-F2/4? (Cowin/Moesen)?

#figure out how to import into moose (from file?, script-replacement?)
#create input file
#automate create input file

#reconsider overall plan

#contact srdjan
