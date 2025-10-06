#!/usr/bin/bash

#alpha version 0.1

#########################################################################
#########################################################################
# Ploymarker flanking region sequence extraction script
# SCRIPT USAGE:
# ./1.2_extract_surrounding_markerseq.sh flanking_region_size
#
# Author: GVarma Saripella
#########################################################################

# for short range 
# Pangenome path:
INDIR=/mnt/c/proj/proj_xxx/data/
OUTDIR=/mnt/c/proj/proj_xxx/work/


# Set Up and Down window size range and Check if input is provided
#updownrange=$1

if [ -z "$1" ]; then
  echo "Error: No input window size provided for Up and Down sequences."
  echo "Usage: $0 <input_up-down-seq-size>"
  echo "Please provide an input sequence size in the range of 50-200."
  echo "Using default number: 50"
  updownrange=50
else
  # Read number from input file
  updownrange=$1
fi

echo "Number used: ${updownrange}"

############

REF=${INDIR}"genome_v3_genomic.fasta"
INFILE=${INDIR}"selected_makrers_xxxx.csv"
INFILE_SNP=${INDIR}"Genotypes_targeted_Bowtie2-Genome_Freebayes-xxxx.csv"
OUTFILE=${OUTDIR}"selected_makrers_xxxx_seq-updown${updownrange}.tsv"

[ -f "$OUTFILE" ] && rm "$OUTFILE"

# Define the header
header="Sno\tChr\tPos\tPos_updown${updownrange}\tSeq_updown${updownrange}"

echo -e "$header" > "${OUTFILE}"

while read p; do
  #echo "$p"
  posNO=$(echo $p | awk -F"," '{print $1}' | sed 's/"//g')
  chrID=$(echo $p | awk -F"," '{print $2}' | sed 's/-/./g'| sed 's/"//g')
  markerpos=$(echo $p | awk -F"," '{print $3}' | grep -o '[0-9]*')

  ##############
  #echo ${posNO}
  #echo ${chrID}
  #echo ${markerpos}
  ##############
  # Define position ranges 
  markerup=$(expr ${markerpos} - ${updownrange} )
  markerdown=$(expr ${markerpos} + ${updownrange})
  markerposrange="${chrID}:${markerup}-${markerdown}"
  #echo ${markerposrange}

  ###############
  # Search for SNPs # 
  chr_markerpos=$(echo ${chrID}"_"${markerpos} | tr '.' '-')
  # Delete monomorphic allele and keep only polymorphic allele
  chr_markerpos_find=$(grep ''${chr_markerpos}'' ${INFILE_SNP} | grep -v '"' | awk -F';' '{ gsub(/\r/, "", $5); if ($5 != "." && length($4) == 1 && length($5) == 1) print "["$4"/"$5"]" }')
  #echo ${chr_markerpos_find}

  ###############

  # Check if the string is empty
  if [ -z "$chr_markerpos_find" ]; then
    continue
    #echo "Skipping pos, missing data!"
  else

    ##############
    # Search in the reference genome 
    # Use sed to replace newlines with empty string except lines starting with '>'
    findstring=$(samtools faidx ${REF} ${markerposrange} | sed -r 's/>/>'${markerpos}'_@_/g')

    # Use awk to Process the FASTA format with tab separation after headers and join all the sequences
    findstring_output=$(echo "$findstring" | awk '/^>/ {if (NR > 1) printf("\t"); printf("%s\t", $0); next} {printf("%s", $0)} END {printf("\n")}' | tr '[:lower:]' '[:upper:]')

    ##############
    # match and check the pos
    findstring_seq=$(echo ${findstring_output} | awk -F" " '{print $2}')
    #echo ${findstring_output}
    #echo ${findstring_seq}

    # Calculate the central position (1-based)
    center=$(expr ${#findstring_seq} / 2)
    
    # Extract the central base (0-based index)
    #central_base=${findstring_seq:$center:1}
    #echo "Central base: ${central_base}"

    ###
    # Construct the new sequence with SNP format at the central base ([G/A])
    findstring_seq_snp="${findstring_seq:0:$center}${chr_markerpos_find}${findstring_seq:$((center + 1))}"
    #echo ${findstring_seq_snp}

    ###############
    # Print the result outfile
    echo -e "$posNO\t$chrID\t$markerpos\t$findstring_output\t$findstring_seq_snp" >> ${OUTFILE}
    #########
  fi
  
  ################
done < ${INFILE}

#########################################################################
