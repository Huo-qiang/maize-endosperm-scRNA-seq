#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe make 1


 samtools view -bS $1\_1_dealed.sam  -o $1\_1_dealed.bam
