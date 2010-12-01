#!/bin/sh
Nside=16 #variables must be same as run_ctpp.sh
nsh=10
OmegaL=0.7
H0=73.

In_Dir='Spherical/CTpp_smoothed'
Map_Dir='Spherical'
CTpp='smoothed_ctpp_Nside'${Nside}'_nsh'${nsh}'_Ok'${1}

Out_File='../Output/Likelihood/Spherical/'${CTpp}'_topmarg.out'
#Prints output to file not to screen

./topmarg nelson.ini 2>$Out_File<< EOF
../Output/CTpp_theory/${In_Dir}/${CTpp}
../Output/CTpp_maps/${Map_Dir}/${CTpp}
EOF

# read seed value and if noise from ini file then tail the output file so that lable it fake_map_(noise or nonoise)_seedvalue.fits
#mv fake_test_noise.fits ../../Ctpp_map_Realizations 

