#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
hisat2 -p 30 -q -x /mnt/diskRAID/maizev4/genome/genome  $1\_1_dealed.R1  -S $1\_1_dealed.sam
