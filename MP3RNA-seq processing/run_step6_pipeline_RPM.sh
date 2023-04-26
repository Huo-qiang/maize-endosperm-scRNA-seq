#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe make 1

intersectBed -bed -abam $1\_1_dealed.sorted.uniq.dedup.bam -b /mnt/diskRAID/huo/MP3code/1/mazie_V4_genes.bed -wb|perl cal_reads.pl - /mnt/diskRAID/huo/MP3code/1/mazie_V4_genes.bed >1.sorted.uniq.dedup.bam_counts
