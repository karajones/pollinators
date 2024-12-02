#!/bin/bash

## SETTINGS USED

# 16S
# Chiar16SF–Chiar16SR
# Forward: TARTYCAACATCGRGGTC
# Reverse: CYGTRCDAAGGTAGCATA
# Min length: 175
# Max length: 300


# CO1 (Bombus)
# BombusF-BombusR
# Forward: AGWCAYCCTGGAATATGAA
# Reverse: GTGGRAAAGCTATATCAGG
# Min length: 100
# Max length: 140


# CO1 (Jusino)
# ANML: LCO1490--CO1‐CFMRa
# Forward: GGTCAACAAATCATAAAGATATTGG
# Reverse: GGWACTAATCAATTTCCAAATCC
# Min length: 100
# Max length: 190

### CHANGE THESE PARAMETERS IF NEEDED ###

# PRIMERS
# Use 5'-3' sequence; do NOT reverse-complement the reverse primer!
FPRIMER="GGTCAACAAATCATAAAGATATTGG"
RPRIMER="GGWACTAATCAATTTCCAAATCC"

### DON'T CHANGE ANYTHING BELOW THIS UNLESS YOU KNOW WHAT YOU ARE DOING! :-) ###

# SET STANDARD VARIABLES
# Number of threads available
NPROC=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)

# DEPENDENCY CHECK
echo -n "Checking dependencies... "
for name in fastp seqkit cutadapt rename fastqc multiqc
do
  [[ $(which $name 2>/dev/null) ]] || { 
    echo -en "\n$name cannot be found.";deps=1; 
    }
done
    [[ $deps -ne 1 ]] && echo "OK" || { 
        echo -en "\nPlease install or activate the above dependencies and rerun this script.\n";exit 1; 
}

# Make necessary directories
mkdir -p stats stats/raw stats/trimmed trimmed results


# Move all files from individual directories into a single directory

find . -name "*.fastq.gz" ! -path "./raw/*" -exec cp {} ./raw/ \;

# rename files
for f in ./raw/*_R1_*.fastq.gz
    do mv $f ${f%%_*}_R1.fastq.gz
done
for f in ./raw/*_R2_*.fastq.gz
    do mv $f ${f%%_*}_R2.fastq.gz
done


# run fastqc to generate html reports of raw reads
fastqc -t $NPROC --svg --noextract -o ./stats/raw ./raw/*.fastq.gz
multiqc --interactive -n multiqc_report_raw_reads.html -o ./stats -i "Read quality before processing" ./stats/raw/.


# Deduplicate, error correct, and detect any remaining adapters with fastp (also output quality report)
### note: fastp doesn't recognize non-ACTG characters so primers must be removed after deduplication
for f in ./raw/*_R1.fastq.gz
    do

    base=$(basename ${f%%_*})    
    R1=${base}_R1
    R2=${base}_R2

    echo -en "\nProcessing sample $base ..."

    # cutadapt call to remove primers
    # -e 0 --no-indels: do not allow errors in adapters
    # --discard-untrimmed: remove sequences with no adapter matches
    # -m 100: discard reads shorter than 100 bp (gets rid of adapter contamination and small junk)
    # --max-n 0: remove reads with Ns (not allowed by DADA2)
    # --poly-a: remove polyA/T tails (sequencing artifacts) <-- note: requires newer version of cutadapt!
    # note: when json is supported by multiqc, `add --json=./stats/raw/$R1.cutadapt.json` to call and remove `> ./stats/raw/$R1.cutadapt.log`
    cutadapt -e 0 --no-indels -m 100 --max-n 0 --poly-a --cores=0 --discard-untrimmed -g $FPRIMER -G $RPRIMER -o ./trimmed/$R1.trimmed.fq.gz -p ./trimmed/$R2.trimmed.fq.gz ./raw/$R1.fastq.gz ./raw/$R2.fastq.gz > ./stats/raw/$R1.cutadapt.log

	done

echo "Generating read statistics..."

# run fastqc to generate html reports of trimmed reads
fastqc -t $NPROC --svg --noextract -o ./stats/trimmed ./trimmed/*.trimmed.fq.gz

# compile stats from fastqc and cutadapt
multiqc --interactive -n multiqc_report_trimmed_reads.html -o ./stats -i "Read quality after processing" ./stats/trimmed/.

#seqkit stats -aT ./raw/*.fastq.gz > ./stats/raw.stats.tsv 
seqkit stats -aT ./trimmed/*.trimmed.fq.gz > ./stats/trimmed.stats.tsv

# cleanup


echo "Finished!"


