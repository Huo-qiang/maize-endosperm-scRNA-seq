#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash

perl deal.R1.pl  $1\_1.R1   >$1\_1_dealed.R1
