#!/bin/bash
#SBATCH -p shared,edwards
#SBATCH -c 1
#SBATCH -t 0-12:00
#SBATCH --mem=32000
#SBATCH --mail-type=END
#SBATCH -o samtools_counts_all_BRD2_q3_array_%A_%a.out
#SBATCH -e samtools_counts_all_BRD2_q3_array_%A_%a.err
#SBATCH --array=0-58    # wc -l all_star_dirs.txt

set -euo pipefail

MIN_MAPQ=3

# BRD2 BLAST hits
brd2_csv=/n/netscratch/edwards_lab/Lab/kelsielopez/MHC_test_scott/BRD2.hits.Kelsie.csv

# STAR mapping directories. Mapping all reads to the entire haplotype assemblies
mapdir_AC=/n/netscratch/edwards_lab/Lab/kelsielopez/MHC_test_scott/aphCoe_gene_expression/STAR_mapping
mapdir_AW=/n/netscratch/edwards_lab/Lab/kelsielopez/MHC_test_scott/aphWoo_gene_expression/STAR_mapping

bedtools_path="/n/home03/kelsielopez/bedtools2/bin/bedtools"

cd /n/netscratch/edwards_lab/Lab/kelsielopez/MHC_test_scott

# One entry per haplotype mapping directoriy
dir_list=all_star_dirs.txt

d_rel=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$dir_list" || true)
[[ -z "$d_rel" ]] && exit 0

base=$(basename "$d_rel")        # AC_1713_89780.hap1 etc
species_prefix=${base%%_*}       # AC or AW

case "$species_prefix" in
    AC) mapdir="$mapdir_AC" ;;
    AW) mapdir="$mapdir_AW" ;;
    *)  exit 0 ;;
esac

d="${mapdir}/${d_rel}"

ind=${base%%.*}                  # AC_1713_89780 etc
hap=${base##*.}                  # hap1 / hap2

case "$hap" in
    hap1) hapnum=1 ;;
    hap2) hapnum=2 ;;
    *)    exit 0 ;;
esac

# BRD2 copy ID in CSV, AC_1713_89780.1 etc
copy_tag="${ind}.${hapnum}"

# Convert BLAST query contigs to BAM contig names:
#   CSV: AC_1713_89780#1#h1tg000512l
#   BAM: AC_1713_89780#hap1#h1tg000512l
contig_from="${ind}#${hapnum}#"
contig_to="${ind}#${hap}#"

bed="${ind}.${hap}_BRD2_300bp.bed"

# Build BED file with BRD2 blast hit window for this haplotype
awk -F'\t' -v OFS='\t' \
    -v copy_tag="$copy_tag" \
    -v from="$contig_from" \
    -v to="$contig_to" '
    NR == 1 { next }                  # skip header
    $12 == copy_tag {
        chrom = $4
        gsub(from, to, chrom)
        q1 = $8; q2 = $9
        if (q1 < q2) { start = q1; end = q2 } else { start = q2; end = q1 }
        start = start - 1             # BED start
        name = $12
        print chrom, start, end, name
    }' "$brd2_csv" > "$bed"

[[ ! -s "$bed" ]] && exit 0

# For each tissue BAM for this haplotype:
#   filter (MAPQ, primary, non-supplementary) and count overlaps with BRD2 window
shopt -s nullglob

bams=( "$d"/${ind}_*_Aligned.sortedByCoord.out.bam )
((${#bams[@]} == 0)) && exit 0

for bam in "${bams[@]}"; do
    bam_base=$(basename "$bam")

    # AC_1713_89780_Blood_Aligned.sortedByCoord.out.bam -> Blood
    tissue=$(echo "$bam_base" \
             | sed -E "s/^${ind}_([A-Za-z]+)_Aligned\.sortedByCoord\.out\.bam/\1/")

    [[ -z "$tissue" ]] && continue

    out="${ind}.${hap}_${tissue}_BRD2_q${MIN_MAPQ}_primary_counts.tsv"

    samtools view -b -q "$MIN_MAPQ" -F 0x900 "$bam" \
        | "$bedtools_path" bamtobed -i - \
        | "$bedtools_path" coverage -a "$bed" -b - -counts \
        | awk -v OFS='\t' '{print $4, $NF}' \
        > "$out"
done
