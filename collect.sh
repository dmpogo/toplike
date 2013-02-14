#!/bin/bash

space=${1:?"First parameter must define data directory"}

cd ${space}/Likelihood
files=`ls -1 *${2:-}*all*${3:-}*${4:-}`
for file in $files
do
   awk -f ${HOME}/work/cmb/topology/Likelihood/collect.awk $file 
done

cd "$OLDPWD" 
