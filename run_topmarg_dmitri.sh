#!/bin/sh
##########################################################################
#-------------------Ini values-------------------------------------------#
##########################################################################
umask 0003

space=I
healpix_dir='~/Packages/Healpix_2.15a'
data_basedir='../Data/WMAP'
ctpp_basedir='../Output1/Spherical-'${space}'/CTpp_theory'
output_basedir='../Output1/Spherical-'${space}
extras=''

signal_file=${data_basedir}'/coadd_cleanimap_16.fits'
#signal_file='../Data/WMAP/coadd_map_8.2deg_16.fits'
ring_weight_file=${healpix_dir}'/data/weight_ring_n00016.fits'

nside=16
nsh=10
OmegaL=0.7
H0=73.
Ok=${1}

do_rotate='.TRUE.'
find_best_angles='.TRUE.'
# If false find best angles will rotate to orientation below
angles='0.0d0 0.0d0 0.0d0'
lmax=40

# Cut is determined by the presence of mask_file
mask_file=${data_basedir}'/kp0_8.2deg_thr0.5_16.fits'
#mask_file=''

#
# If do_smooth=TRUE  and noise_file is set 
# we smooth noise correlation before adding matrix to CTpp
# If do_smooth=FALSE and noise_file is set 
# the noise correlation matrix is added to CTpp "as is" 
#                                   (presumably should be diagonal )
do_smooth='.TRUE.'
G_fwhm=492.0
noise_file=${data_basedir}'/coadd_noise_8.2deg_16.fits'
#epsil is always added to the diagonal of CTpp. May be used for regularization
#It must be set to a value, put 0 (or negative) for no regularization
epsil='1.0d-7'

#make map visualization from CTpp
make_map='.FALSE.'
#makes map only does not find max likelihood
make_map_only='.FALSE.'
add_map_noise='.FALSE.'

#if set to 0 will use system clock
iseed=1

##########################################################################
#--------------------Running script--------------------------------------#
##########################################################################
#Prints output to file not to screen

#Construct file names from parameters given

if [ "$mask_file" == "" ]; then
   mask=""
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

#Use nice output file
output_file='.TRUE.'
nice_out_file=${output_basedir}/Likelihood/topmarg_${space}_smooth${G_fwhm}_nside${nside}_nsh${nsh}_Ok${1}_epsil${epsil}${extras}.out
# Debug output
run_out_file=${output_basedir}/Likelihood/topmarg_${space}_fullrun${smooth}_nside${nside}_nsh${nsh}_Ok${1}_epsil${epsil}${extras}.out
# Concise output
short_out_file=${output_basedir}/Likelihood/topmarg_${space}_smooth${G_fwhm}_nside${nside}_nsh${nsh}_allOk_epsil${epsil}${extras}.out
# Map output
map_out_file=${output_basedir}/CTpp_maps/map${mask}${smooth}_nside${nside}_nsh${nsh}_Ok${1}${noise}${extras}.fits

# Perform preprocessing if needed
if [ "$do_smooth" == ".TRUE."  ]; then
# We want smoothed CTpp
   beam_file=${output_basedir}'/CTpp_beams/beam_array_gaussian'${G_fwhm}
   CTpp_unsmoothed=${ctpp_basedir}/CTpp_unsmoothed/ctpp_${space}_Nside${nside}_nsh${nsh}_Ok${1}${extras}
   CTpp=${ctpp_basedir}'/CTpp_smoothed/ctpp_'${space}'_smoothed'${G_fwhm}'_Nside'${nside}'_nsh'${nsh}'_Ok'${1}${extras}

#Smooth Ctpp if not already smoothed
   if [ ! -e "$CTpp" ]; then
      echo "Smoothing CTpp"
      if [ -e "$beam_file" ]; then
         save_beam=0
      else
         save_beam=1
      fi

echo ${extras}
../Spherical/test_process << EOF
${CTpp_unsmoothed}
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
# We want unsmoothed CTpp
      beam_file=''
      CTpp=${ctpp_basedir}'/CTpp_unsmoothed/ctpp_'${space}'_Nside'${nside}'_nsh'${nsh}'_Ok'${1}${extras}
fi

#Main call, if statement to screen
echo 'Starting Ok'${1}
if [ "$2" == -screen ]; then
./topmarg << EOF
${nice_out_file}
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

${do_rotate}
${angles}
${find_best_angles}
${lmax}
${epsil}
EOF

else
./topmarg 2>$run_out_file >> $short_out_file << EOF
${nice_out_file}
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


