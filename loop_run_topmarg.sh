#!/bin/sh
inifile=nelson.ini
epsil='1.0E-6'
for i in '-1.0E-1'
#'-1.0E-3' '-6.0E-3' '-8.0E-3' '-1.0E-2' '-1.5E-2' '-2.0E-2' '-2.5E-2' '-3.5E-2' '-4.0E-2' '-5.0E-2' '-6.0E-2' '-7.0E-2' '-8.0E-2' '-9.0E-2' '-1.0E-1' '-1.25E-1' '-1.5E-1' '-1.75E-1' '-2.0E-1'
do
./run_topmarg.sh ${i} ${inifile} ${epsil} 
done