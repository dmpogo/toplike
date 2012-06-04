#!/bin/sh
##########################################################################
#-------------------Ini values-------------------------------------------#
##########################################################################
umask 0003

signal_file='../Data/WMAP/coadd_cleanimap_16.fits'
#signal_file='../Data/WMAP/coadd_map_8.2deg_16.fits'
ring_weight_file='~/Packages/Healpix_2.15a/data/weight_ring_n00016.fits'

space='Spherical'
nside=16
nsh=10
OmegaL=0.7
H0=73.
Ok=${1}

do_rotate='.TRUE.'
find_best_angles='.TRUE.'
# If not find best angles will rotate to orientation below
angles='0.0d0 0.0d0 0.0d0'
lmax=40

mask_file='../Data/WMAP/kp0_8.2deg_thr0.5_16.fits'
mask_file=''

#
# If do_smooth=TRUE  and noise_file is set 
# we smooth noise correlation before adding matrix to CTpp
# If do_smooth=FALSE and noise_file is set 
# the noise correlation matrix is added to CTpp "as is" 
#                                   (presumably should be diagonal )
do_smooth='.TRUE.'
G_fwhm=492.0
noise_file='../Data/WMAP/coadd_noise_8.2deg_16.fits'
#epsil is always added to the diagonal of CTpp. May be used for regularization
epsil='0.0d0'

#make map visualization from CTpp
<<<<<<< HEAD
make_map='.TRUE'
#makes map only does not find max likelihood
make_map_only='.TRUE'
add_map_noise='.FALSE.'
=======
make_map='.TRUE.'
#makes map only does not find max liklihood
make_map_only='.TRUE.'
add_map_noise='.TRUE.'
>>>>>>> 7f622a98736f332550eb381fe77d10f30cd01749

#if set to 0 will use system clock
iseed=1

nice_output='.TRUE.'

#Use nice out put file
output_file='.TRUE.'
nice_out_file='../Output/Likelihood/'${space}'/Results/topmarge_smooth'${G_fwhm}'_nside'${nside}'_Ok'${1}'epsil'${epsil}'.out'
##########################################################################
#--------------------Running script--------------------------------------#
##########################################################################
#Prints output to file not to screen
<<<<<<< HEAD

#Construct file names from parameters given

if [ "$mask_file" == "" ]; then
   mask=""
=======
#Constructs file names from parameters given
if [ "$do_smooth" == ".TRUE."  ]; then
run_out_file='../Output/Likelihood/'${space}'/Complete_run/topmarge_fullrun_smooth'${G_fwhm}'_nside'${nside}'_Ok'${1}'epsil'${epsil}'.out'
beam_file='../Output/CTpp_theory/'${space}'/CTpp_beams/beam_array_gaussian'${G_fwhm}
CTpp='../Output/CTpp_theory/'${space}'/CTpp_smoothed/ctpp_smoothed'${G_fwhm}'_Nside'${nside}'_nsh'${nsh}'_Ok'${1}
if [ "$add_map_noise" == ".TRUE." ]; then
map_out_file='../Output/CTpp_maps/Spherical/map_smoothed'${G_fwhm}'_nside'${nside}'_Ok'${1}'_noise.fits'
>>>>>>> 7f622a98736f332550eb381fe77d10f30cd01749
else
   mask="_kp0"
fi

if [ "$do_smooth" == ".TRUE." ]; then
   smooth="_smoothed${G_fwhm}"
else
   smooth=''
fi

if [ "$add_map_noise" == ".TRUE." ]; then
   noise='_cnoise'
else
   noise='_cnonoise'
fi

map_out_file=../Output/CTpp_maps/${space}/map${mask}${smooth}_nside${nside}_Ok${1}${noise}.fits
run_out_file=../Output/Likelihood/${space}/Complete_run/topmarge_fullrun${smooth}_nside${nside}_Ok${1}epsil${epsil}.out

# Perform preprocessing if needed
if [ "$do_smooth" == ".TRUE."  ]; then
   beam_file='../Output/CTpp_theory/'${space}'/CTpp_beams/beam_array_gaussian'${G_fwhm}
   CTpp='../Output/CTpp_theory/'${space}'/CTpp_smoothed/ctpp_smoothed'${G_fwhm}'_Nside'${nside}'_nsh'${nsh}'_Ok'${1}

#Smooth Ctpp if not already smoothed
   if [ ! -e "$CTpp" ]; then
      echo "Smoothing CTpp"
      if [ -e "$beam_file" ]; then
         save_beam=0
      else
         save_beam=1
      fi

../Spherical/test_process << EOF
../Output/CTpp_theory/${space}/CTpp/ctpp_Nside${nside}_nsh${nsh}_Ok$1
${CTpp}
0
1
${save_beam}
${G_fwhm}
${beam_file}
EOF
echo "Done smoothing CTpp"
   fi

   else
      beam_file=''
      CTpp='../Output/CTpp_theory/'${space}'/CTpp/ctpp_Nside'${nside}'_nsh'${nsh}'_Ok'${1}
   fi

#Main call, if statement to screen
echo 'Starting Ok'${1}
if [ "$2" == -screen ]; then
./topmarg << EOF
${nice_out_file}
Printed to screen
${signal_file}
${noise_file}
${mask_file}
${beam_file}
${ring_weight_file}
${CTpp}
${map_out_file}

${output_file}

${nside}
${nsh}
${OmegaL}
${H0}
${Ok}

${make_map}
${add_map_noise}
${iseed}
${make_map_only}
${nice_output}

${do_rotate}
${angles}
${find_best_angles}
${lmax}
${epsil}
EOF

else
./topmarg 2>$run_out_file<< EOF
${nice_out_file}
${run_out_file}
${signal_file}
${noise_file}
${mask_file}
${beam_file}
${ring_weight_file}
${CTpp}
${map_out_file}

${output_file}

${nside}
${nsh}
${OmegaL}
${H0}
${Ok}

${make_map}
${add_map_noise}
${iseed}
${make_map_only}
${nice_output}

${do_rotate}
${angles}
${find_best_angles}
${lmax}
${epsil}
EOF
fi
#epsilfile="epsil_ok"${Ok}".data"
#eigenfile="epsil_ok"${Ok}".eigen"
#mv epsilplot.data ${epsilfile}
#mv epsil.eigen ${eigenfile}
echo 'DONE'


