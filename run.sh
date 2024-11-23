#!/bin/bash

EH=1.0
EH_bg=1.0
EG=20
CvG=0.3
CvH=0.0
CvH_bg=0.0
lmd_bg=0.9
mu=1.0

outfile="test.txt"
echo "#lmd CvG EH CvH lmd_bg EH_bg CvH_bg EA mu" > $outfile

outfile="test.txt"
./a.out --EG=${EG} --EH=${EH} --EH_bg=${EH_bg} --CvG=${CvG} --CvH=${CvH} --CvH_bg=${CvH_bg} --lmd_bg=${lmd_bg} --mu=${mu} --output=${outfile}
