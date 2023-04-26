#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

perl deduplicate.pl $1\_1.R2  $1\_1_dealed.sorted.uniq.bam  |samtools view -bS - -o  $1\_1_dealed.sorted.uniq.dedup.bam
