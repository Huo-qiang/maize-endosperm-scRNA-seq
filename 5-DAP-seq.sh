#-----------------------------trim--------------------------------------------------------------##
#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
for filename in `ls /mnt/diskRAID/Huo/DAPseq/*.R1.fastq.gz`

do
  # set path names
  dirpaired="/mnt/diskRAID/Huo/DAPseq/trimed"
  dirunpaired="/mnt/diskRAID/Huo/DAPseq/trimout"
  base=$(basename $filename ".R1.fastq.gz")
  
  echo $base

     java -jar /mnt/disk2T/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 14 -phred33 \
     ${base}.R1.fastq.gz \
     ${base}.R2.fastq.gz \
     ${dirpaired}/${base}_R1_paired.fq.gz ${dirunpaired}/${base}_R1_unpaired.fq.gz \
     ${dirpaired}/${base}_R2_paired.fq.gz ${dirunpaired}/${base}_R2_unpaired.fq.gz \
     ILLUMINACLIP:/mnt/disk2T/WangQun/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 \
     LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:0 MINLEN:75
done;
#------------------------------bowtie2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##
#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64

for sample in `ls /mnt/diskRAID/Huo/DAPseq/trimed/*.R1.fastq.gz`

do

   # set path names
   dirfastq="/mnt/diskRAID/Huo/DAPseq/trimed/"
   dirsam="/mnt/diskRAID/Huo/DAPseq/bam/"
   base=$(basename $sample ".R1.fastq.gz")

   # alignment by bowtie
   bowtie2 -p 64 -I 75 -X 1000 --no-discordant --no-mixed -x /mnt/disk2T/WangQun/bowtie2-build/index -1 ${dirfastq}/${base}.R1.fastq.gz -2 ${dirfastq}/${base}.R2.fastq.gz -S ${dirsam}/${base}.sam 2> ${dirsam}/${base}.bowtie2.v4prefixchr.log

   # Convert file from SAM to BAM format  
   samtools view -bS -h ${dirsam}/${base}.sam > ${dirsam}/${base}.bam  
   
   # Sort BAM file  
   samtools sort ${dirsam}/${base}.bam -o ${dirsam}/${base}.sorted.bam
   
   # index the bam files  
   samtools index ${dirsam}/${base}.sorted.bam

   # write bam statistic 
   samtools flagstat ${dirsam}/${base}.sorted.bam > ${dirsam}/${base}.bam-log

   # Remove intermediate files
   rm ${dirsam}/${base}.sam
   rm ${dirsam}/${base}.bam
 
done;
##----------------------------------------------------------------------Duplicates----------------------------------------------------------------------##
#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64

for sample in `ls /mnt/diskRAID/Huo/DAPseq/bam/*.sorted.bam`

do

   # set path names
   dirbam="/mnt/diskRAID/Huo/DAPseq/bam"
   dirdedupbam="/mnt/diskRAID/Huo/DAPseq/dedupbam"
   base=$(basename $sample ".sorted.bam")


java -jar -Xmx50G /mnt/diskRAID/Huo/picard/picard.jar MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
I=${dirbam}/${base}.sorted.bam \
O=${dirdedupbam}/${base}.dedup.bam \
M=${dirdedupbam}/${base}.picard.txt

wait
 
   # Sort BAM file  
   samtools sort ${dirdedupbam}/${base}.dedup.bam -o ${dirdedupbam}/${base}.dedup.sorted.bam
   
   # index the bam files  
   samtools index ${dirdedupbam}/${base}.dedup.sorted.bam

   # write bam statistic 
   samtools flagstat ${dirdedupbam}/${base}.dedup.sorted.bam > ${dirdedupbam}/${base}.dedup.bam-log

   # Remove intermediate files
   rm ${dirdedupbam}/${base}.dedup.bam
 
done;
##----------------------------------------------------------------------GEM----------------------------------------------------------------------##
#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64

for sample in `ls /mnt/diskRAID/Huo/DAPseq/dedupbam/*.dedup.sorted.bam`

do

   # set path names
   dirbam="/mnt/diskRAID/Huo/DAPseq/dedupbam"
   base=$(basename $sample ".dedup.sorted.bam")


java -Xmx20G -jar /mnt/diskRAID/Huo/DAPseq/gem/gem.jar --t 64 --g /mnt/diskRAID/Huo/DAPseq/gem/zeamays.chrom.sizes --d /mnt/diskRAID/Huo/DAPseq/gem/Read_Distribution_default.txt --genome /public3/home/pgv3167/ZM/gem/Zeamaysv4prefixchr --expt ${dirbam}/${base}.dedup.sorted.bam --ctrl /mnt/diskRAID/Huo/DAPseq/dedup/puc57_FKDL192541313.dedup.sorted.bam --f SAM --out /mnt/diskRAID/Huo/DAPseq/gem/results/${base} --k_min 4 --k_max 12 --fold 3 --q 3 --nrf --outNP --outHOMER --outMEME

done;

##---------------------------------------------------------------bw------------------------------------------------------------------------------##
#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64

for sample in `ls /mnt/diskRAID/Huo/DAPseq/dedupbam/*.dedup.sorted.bam`

do

   # set path names
   dirbam="/mnt/diskRAID/huo/bam"
   dirbw="/mnt/diskRAID/huo/bam/bw10"
   base=$(basename $sample ".dedup.sorted.bam")

   bamCoverage -b ${dirbam}/${base}.dedup.sorted.bam -o  ${dirbw}/${base}.bw --binSize 300 --normalizeUsing CPM -p 40   

done;
##--------------------------------------------------------------FRiP----------------------------------------------------------------------------##
for sample in `ls /mnt/diskRAID/Huo/DAPseq/dedupbam/*.dedup.sorted.bam`

do

   # set path names
   dirdedupbam="/mnt/diskRAID/Huo/DAPseq/dedupbam"
   dirgem="/mnt/diskRAID/Huo/DAPseq/gem/results"
   dirFRIP="/mnt/diskRAID/Huo/DAPseq/FRiP"
   base=$(basename $sample ".dedup.sorted.bam")

# total reads
total_reads=$(samtools view -c ${dirdedupbam}/${base}.dedup.sorted.bam)

# reads in peaks
reads_in_peaks=$(bedtools sort -i ${dirgem}/${base}.GEM_rmgreylist.bed \
| bedtools merge -i stdin | bedtools intersect -u -nonamecheck \
-a ${dirdedupbam}/${base}.dedup.sorted.bam -b stdin -ubam | samtools view -c)

# FRiP score
(echo "${base}"; awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}"; echo "FRiP") > ${dirFRIP}/${base}.FRiP

done;
##---------------------------------------------------peak.fa--------------------------------------------------------------------------------------##

for sample in `ls /mnt/diskRAID/Huo/DAPseq/gem/results/narrowpeak/*.GEM_events.narrowPeak`

do
	      dirmeme="/mnt/diskRAID/Huo/DAPseq/gem/results/narrowpeak"
	         dirfa="/mnt/diskRAID/Huo/DAPseq/gem/results/meme/fa"
		    base=$(basename $sample ".GEM_events.narrowPeak")

		    bedtools getfasta -fi mnt/diskRAID/Huo/gem/Zea_mays.B73_RefGen_v4.chr.fa -bed ${dirmeme}/${base}.GEM_events.narrowPeak -fo ${dirfa}/${base}.fa
	    done;
##---------------------------------------------------meme--------------------------------------------------------------------------------------##    
      
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32

for sample in `ls /mnt/diskRAID/Huo/DAPseq/gem/results/meme/fa/*.fa`;

do

	   # set path names
	      dirmeme="/mnt/diskRAID/Huo/DAPseq/gem/results/meme/results"
	         dirfa="/mnt/diskRAID/Huo/DAPseq/gem/results/meme/fa"
		    base=$(basename $sample ".rmblacklist.bed")

		    meme-chip -minw 6 -maxw 12 -ccut 100 -meme-nmotifs 0 -spamo-skip -fimo-skip -oc ${dirmeme}/${base} ${dirfa}/${base}.fa

	    done;
##############################################################################################################################################################
