######################################overlap DAP-seq and GRN#############################################################
for i in {1..50}; do head -n 1000000 GRN${i}.tsv > GRN${i}_extracted.tsv; done
cat GRN*_extracted.tsv > combined.tsv
awk '{print $1"_"$2,$4}' combined.tsv | cut -f2- > modified_combined.tsv
cat modified_combined.tsv | sort | uniq -c  > modified_combined_with_counts.tsv

for filename in `ls /mnt/diskRAID/Huo/GRN/*.tsv`

do
  # set path names
  dirGRN="/mnt/diskRAID/Huo/DAPseq/trimed"
  diroverlap="/mnt/diskRAID/Huo/DAPseq/trimout"
  base=$(basename $filename ".tsv")
grep -F -f ${dirGRN}/${base}.txt /mnt/diskRAID/Huo/DAP-seq/all-peak.txt | sort | uniq > $diroverlap}/${base}-dap.txt
done;
#####################################合并并统计调控关系出现次数###################################################
cat grn1-dap.txt grn2-dap.txt grn3-dap.txt grn4-dap.txt grn5-dap.txt grn6-dap.txt grn7-dap.txt grn8-dap.txt grn9-dap.txt grn10-dap.txt grn11-dap.txt grn12-dap.txt grn13-dap.txt grn14-dap.txt grn15-dap.txt grn16-dap.txt grn17-dap.txt grn18-dap.txt grn19-dap.txt grn20-dap.txt grn21-dap.txt grn22-dap.txt grn23-dap.txt grn24-dap.txt grn25-dap.txt grn26-dap.txt grn27-dap.txt grn28-dap.txt grn29-dap.txt grn30-dap.txt grn31-dap.txt grn32-dap.txt grn33-dap.txt grn34-dap.txt grn35-dap.txt grn36-dap.txt grn37-dap.txt grn38-dap.txt grn39-dap.txt grn40-dap.txt grn41-dap.txt grn42-dap.txt grn43-dap.txt grn44-dap.txt grn45-dap.txt grn46-dap.txt grn47-dap.txt grn48-dap.txt grn49-dap.txt grn50-dap.txt > quanjidap.txt
cat quanjidap.txt | sort | uniq -c >resultdap.txt
#########################################################################################################################################################################
