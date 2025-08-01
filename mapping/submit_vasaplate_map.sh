#!/bin/bash

### input paths (to modify by user)
p2s=/home/projects/dp_immunoth/data/pbmc/vasaseq/VASAseq/mapping   # path to mapping scripts in your computer/HPC
p2trimgalore=/home/projects/dp_immunoth/data/pbmc/vasaseq/trimgalore_env/bin  # path to TrimGalore
p2cutadapt=/home/projects/dp_immunoth/data/pbmc/vasaseq/trimgalore_env/bin                # path to cutadapt
p2bwa=/services/tools/bwa/0.7.17               # path to BWA
p2samtools=/services/tools/samtools/1.20/bin       #  path to samtools
p2star=/services/tools/star/2.7.10b/bin    # path to STAR
p2bedtools=/services/tools/bedtools/2.30.0/bin          # path to bedtools

### check input parameters
if [ $# -ne 7 ]
then
    echo "Please, give:"
    echo "1) library name (prefix of the fastq files, name before _R1.fastq.gz and _R2.fastq.gz)"
    echo "2) genome: MOUSE /  HUMAN / MIXED "
    echo "3) read length (for MOUSE: 59, 74, 96, 246; for HUMAN: 74, 136; for MIXED: 74, 91, 135)"
    echo "4) prefix for output files"
    echo "5) folder for output files"
    echo "6) fastqfile extraction (y/n)"
    echo "7) cellID from filename or readname (f/r)"
    exit
fi

lib=$1
ref=$2
n=$3
out=$4
folder=$5
fqext=$6
cellidori=$7

### check existence of input fastq files
r1=$(ls ${lib}*_R1*.fastq.gz 2>/dev/null)
r2=$(ls ${lib}*_R2*.fastq.gz 2>/dev/null)
echo "Found R1 files: $r1"
echo "Found R2 files: $r2"
if [ ${#r1} == 0 ]
then
    echo "R1 fastq files not found"
    exit 1
fi
if [ ${#r2} == 0 ]
then
    echo "R2 fastq files not found"
    exit 1
fi

### check python version (we want version 3)
v=$(python -c 'import sys; print(".".join(map(str, sys.version_info[:3])))' | awk -F "." '{print $1}')
if [ $v -ne "3" ]
then
    echo "python needs to be 3"
    exit 1
fi

### Create output folder if it doesn't exist
mkdir -p ${folder}

### set references
if [[ $ref == "MOUSE" ]]
then
    riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/rRNA_mouse_Rn45S.fa
    genome=/hpc/hub_oudenaarden/group_references/ensembl/99/mus_musculus/star_v273a_NOMASK_NOERCC_index_$n
    refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed
elif [[ $ref == "HUMAN" ]]
then
    riboref=/home/projects/dp_immunoth/data/pbmc/vasaseq/references_114/unique_rRNA_human.fa
    genome=/home/projects/dp_immunoth/data/pbmc/vasaseq/references_114/star_v2114_index_$n
    refBED=/home/projects/dp_immunoth/data/pbmc/vasaseq/references_114/Homo_sapiens.GRCh38.114.homemade_IntronExonTrna.bed
fi 

if [ ! -d $genome ]
then
    echo "genome not found: $genome"
    exit 1
fi

### Function to print status with timestamp
print_status() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

### extract cell barcodes
if [ $fqext == "y" ]
then
    print_status "Starting barcode extraction..."
    ${p2s}/extractBC.sh ${lib} vasaplate ${p2s} ${folder}
    if [ $? -ne 0 ]; then
        echo "Error in barcode extraction"
        exit 1
    fi
    print_status "Barcode extraction completed successfully"
    echo "Run the script again with fastqfile extraction set to 'n' to continue with mapping"
    exit 0
fi

### Process each barcode file
print_status "Starting processing of barcode files..."
for file in ${folder}/${lib}*cbc.fastq.gz
do
    if [ ! -f "$file" ]; then
        echo "No barcode files found matching pattern: ${folder}/${lib}*cbc.fastq.gz"
        exit 1
    fi
    
    lib_current=${file%_cbc.fastq.gz}
    lib_current=${lib_current/}
    print_status "Processing file: $lib_current"
    
    ### trim
    print_status "Starting trimming for $lib_current..."
    ${p2s}/trim.sh ${lib_current}_cbc.fastq.gz ${folder} ${p2trimgalore} ${p2cutadapt}
    if [ $? -ne 0 ]; then
        echo "Error in trimming for $lib_current"
        exit 1
    fi
    print_status "Trimming completed for $lib_current"
    
    ### Check if trimmed file exists
    trimmed_file="${lib_current}_cbc_trimmed_homoATCG.fq.gz"
    if [ ! -f "$trimmed_file" ]; then
        echo "Trimmed file not found: $trimmed_file"
        exit 1
    fi

    ### ribo-map
    print_status "Starting ribosomal mapping for $lib_current..."
    ${p2s}/ribo-bwamem.sh $riboref ${trimmed_file} ${lib_current}_cbc_trimmed_homoATCG $p2bwa $p2samtools y $p2s
    if [ $? -ne 0 ]; then
        echo "Error in ribosomal mapping for $lib_current"
        exit 1
    fi
    print_status "Ribosomal mapping completed for $lib_current"
    
    ### Check if non-ribosomal file exists
    nonribo_file="${lib_current}_cbc_trimmed_homoATCG.nonRibo.fastq.gz"
    if [ ! -f "$nonribo_file" ]; then
        echo "Non-ribosomal file not found: $nonribo_file"
        exit 1
    fi

    ### map to genome
    print_status "Starting genome mapping for $lib_current..."
    ${p2s}/map_star.sh ${p2star} ${p2samtools} ${genome} ${nonribo_file} ${lib_current}_cbc_trimmed_homoATCG.nonRibo_E114_
    if [ $? -ne 0 ]; then
        echo "Error in genome mapping for $lib_current"
        exit 1
    fi
    print_status "Genome mapping completed for $lib_current"
    
    ### Check if BAM file exists
    bam_file="${lib_current}_cbc_trimmed_homoATCG.nonRibo_E114_Aligned.out.bam"
    if [ ! -f "$bam_file" ]; then
        echo "BAM file not found: $bam_file"
        exit 1
    fi

    ### deal with single mappers
    print_status "Processing single mappers for $lib_current..."
    ${p2s}/deal_with_singlemappers.sh ${bam_file} ${refBED} y ${p2samtools} ${p2bedtools}
    if [ $? -ne 0 ]; then
        echo "Error processing single mappers for $lib_current"
        exit 1
    fi
    print_status "Single mapper processing completed for $lib_current"

    ### deal with multi mappers
    print_status "Processing multi mappers for $lib_current..."
    ${p2s}/deal_with_multimappers.sh ${bam_file} ${refBED} y ${p2samtools} ${p2bedtools}
    if [ $? -ne 0 ]; then
        echo "Error processing multi mappers for $lib_current"
        exit 1
    fi
    print_status "Multi mapper processing completed for $lib_current"
    
    print_status "Completed processing for $lib_current"
done

### Generate count tables
print_status "Starting count table generation..."
python3 ${p2s}/countTables_2pickle_cellsSpliced.py ${folder} ${folder}/${out} vasa $cellidori
if [ $? -ne 0 ]; then
    echo "Error in count table generation (pickle step)"
    exit 1
fi
print_status "Pickle generation completed"

### Check if pickle file exists
pickle_file="${folder}/${out}.pickle.gz"
if [ ! -f "$pickle_file" ]; then
    echo "Pickle file not found: $pickle_file"
    exit 1
fi

print_status "Processing pickle file to generate final count tables..."
python3 ${p2s}/countTables_fromPickle.py ${pickle_file} ${folder}/${out} vasa y
if [ $? -ne 0 ]; then
    echo "Error in final count table generation"
    exit 1
fi

print_status "Pipeline completed successfully!"
print_status "Output files are located in: ${folder}"
print_status "Count tables have prefix: ${folder}/${out}"

### List generated count table files
echo "Generated count table files:"
ls -la ${folder}/${out}*Counts.tsv 2>/dev/null || echo "No count table files found with expected pattern"
