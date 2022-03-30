#!/bin/sh

###############################################################################
# Define threads
###############################################################################

[ -z "$threads" ] && threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null | awk '{print int($0) - 1}')
[ -z "$threads" ] && threads=1

###############################################################################
# Genome indexing by Subread
###############################################################################

mkdir -p data/mouse_index/subread

time subread-buildindex \
  -o data/mouse_index/subread/subread \
  data/mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa

###############################################################################
# Mapping to genome by Subread
###############################################################################

mkdir -p data/bam

find data/fastq/ -type f |
  sort |
  awk 'NR%2==1 {printf $0" "; next}1' |
  while read -r R1 R2; do
    out_prefix="$(basename ${R1%%_R*})"
    #
    subread-align -t 0 \
      -T "${threads:-1}" \
      -i data/mouse_index/subread/subread \
      -r "$R1" -R "$R2" \
      -o data/bam/tmp.bam
    #
    samtools sort -@ "$threads" data/bam/tmp.bam >data/bam/"$out_prefix".bam
    samtools index -@ "$threads" data/bam/"$out_prefix".bam
    samtools stats -@ "$threads" data/bam/"$out_prefix".bam >data/bam/"$out_prefix"_stats
    rm data/bam/tmp.bam*
  done

###############################################################################
# BigWig files to visualize by IGV
###############################################################################

mkdir -p data/bigwig

for bam in bam/*bam; do
  out_f=$(echo "$bam" |
    sed "s|bam/|bigwig/|" |
    sed "s/.bam$/_bin5_cpm.bw/g")

  bamCoverage -b "$bam" -o "$out_f" -p "$threads" \
    --binSize 5 --normalizeUsing CPM
done
