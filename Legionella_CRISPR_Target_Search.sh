#!/bin/bash

# Usage: bash Legionella_CRISPR_Target_Search.sh <SPACER LIST HERE> <ALIGNMENT LENGTH FILTER HERE> <MISMATCH FILTER HERE>
# To find targets of Legionella's CRISPR system using a BLAST-based approach
# The script is commented as much as possible to make my reasoning clearer for others using this script

######################################
# Remove line breaks from FASTA file #
######################################

# This awk command will remove line breaks from the FASTA file, while still leaving in the header as denoted by ">"
	# This will allow us to mask out CRISPR arrays in a downstream step
	# This will minimize "hits" being called to CRISPR arrays in genomes
	
# In this step and every step downstream of it, setting the variable filename="${input%%.*}" allows us to call each file by its unique identifier
	# Removes all the extensions at the end of the file since they're separated by .
	# Leaves only the "basename" of the file so that it can be invoked through out the script and we can add on any extension we want to the end of it
		# This also gives the script more flexibility to analyze any file in a database without having to worry about specific database nomenclature for their individual fasta files

echo "Removing line breaks from FASTA file for downstream masking"

for input in *.fasta; do
	filename="${input%%.*}"
	awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' "$input" > "$filename".formatted.fasta
	done
	
#########################################
# Mask CRISPR arrays from input genomes #
#########################################

# This sed command will find all instances of a specific sequence and "mask" it by replacing it with a series of "N"s
	# Keep in mind these are from Legionella CRISPR systems so they many vary slightly compared with repeats from other organisms

# Breaking it down
	# In general, the format is -e s/pattern1/pattern2/g, which tells the program to replace every instance of pattern 1 with pattern 2
	# In our case, pattern 1 is the repeat for a given CRISPR system (i.e.: GTT.ACTGCCG.ACAGGCAGCTTAGAA.), .* (i.e.: everything, this accounts for spacers between the repeats), and another repeat
		# This tells the program to replace the CRISPR repeats and everything in between them with a row of 32 Ns
	# The downside of this approach is that it shifts the positions of downstream nucleotides so that when we go to map hits them, it'll be a little trickier
		# The positional information from BLAST won't match up with the positional information in the original file
	# However, we'll have the extracted sequence available in the output file so we could always use that and "cntl + F" to find our target region
		# Alternatively, the masked sequence file will be saved so you could always refer back to that because the positions will match the BLAST output

echo "Masking the CRISPR array for type I-C, I-F and II-B CRISPR-Cas systems"

for input in *.formatted.fasta; do
	filename="${input%%.*}"
	#Type I-F repeat
	sed -e 's/GTT.ACTGCCG.ACAGGCAGCTTAGAA..*GTT.ACTGCCG.ACAGGCAGCTTAGAA./NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/g' "$input" > "$filename".temp1.fasta
	#Type I-F repeat (reverse complement)
	sed -e 's/.TTCTAAGCTGCCTGT.CGGCAGT.AAC.*.TTCTAAGCTGCCTGT.CGGCAGT.AAC/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/g' "$filename".temp1.fasta > "$filename".temp2.fasta
	#Type I-C repeat
	sed -e 's/GTCGCGCCCCGTGCGGGCG.G.GGATTGAA.C.*GTCGCGCCCCGTGCGGGCG.G.GGATTGAA.C/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/g' "$filename".temp2.fasta > "$filename".temp3.fasta
	#Type I-C repeat (reverse complement)
	sed -e 's/G.TTCAATCC.C.CGCCCGCACGGGGCGCGAC.*G.TTCAATCC.C.CGCCCGCACGGGGCGCGAC/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/g' "$filename".temp3.fasta > "$filename".temp4.fasta
	#Type II-B repeat
	sed -e 's/CCAATAATCCCTCATCTAAAAAT..AAC.A.TGAAA..*CCAATAATCCCTCATCTAAAAAT..AAC.A.TGAAA./NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/g' "$filename".temp4.fasta > "$filename".temp5.fasta
	#Type II-B repeat (reverse complement)
	sed -e 's/.TTTCA.T.GTT..ATTTTTAGATGAGGGATTATTGG.*.TTTCA.T.GTT..ATTTTTAGATGAGGGATTATTGG/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/g' "$filename".temp5.fasta > "$filename".masked.fasta
	done

#####################################
# BLAST spacers against the genomes #
#####################################

# Make BLAST database for putative targets using all the masked FASTA sequences
	# This pipeline does not need parse_seq_ids flag to run because not using the NCBI parser downstream to look at hits
	
echo "Making BLAST database with masked sequences"

for input in *.masked.fasta; do
	makeblastdb -in $input -dbtype nucl
	done
	
# BLAST spacers against target database in parallel format to speed up process
	# blastn parameters are similar to those used in CRISPRTarget online tool (which this script was inspired by)

# This will align spacers to the database and any potential hits could be potential targets for that spacer

# Also BLASTed for hits on the "plus strand" and the "minus strand" to keep them separate for easier PAM filtering downstream
	# As a result, the rest of the pipeline has every step duplicated in order to keep the two strands separate

# The spacers=${1?Error: no spacer list detected} variable uses the fasta formated spacer list given by the user as the query
	# Will show an error message saying "no spacer list detected" if the spacer list was not specified
	# Things that will cause this to not work properly
		# Spacer list has dashes, underscores or dots in the name
		# Spacer list does not have a .txt file extension
		# Spacer list is not a fasta formated file
		# The program is unable to find the path to the file, ensure it's in the same directory or the entire path to the file is in the terminal window
	
echo "Starting BLAST search on plus strand"

for input in *.masked.fasta; do
	spacers=${1?Error: no spacer list detected}
	find "$input" | parallel -j0 -N 60 'blastn -db {} -query '$spacers' -gapopen 10 -gapextend 2 -reward 1 -penalty -1 -evalue 0.01 -word_size 7 -strand plus -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" -out {}.plus.txt'
	done

echo "Starting BLAST search on minus strand"

for input in *.masked.fasta; do
	spacers=${1?Error: no spacer list detected}
	find "$input" | parallel -j0 -N 60 'blastn -db {} -query '$spacers' -gapopen 10 -gapextend 2 -reward 1 -penalty -1 -evalue 0.01 -word_size 7 -strand minus -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" -out {}.minus.txt'
	done
	
########################
# Score the BLAST hits #
########################

#The awk command here takes the raw BLAST output file, and prints the following columns
	# $1 -> the query ID aka the spacer ID
	# $2 -> the sequence ID aka the putative target
	# $3 -> the percent ID for the match
	# $11 -> the evalue for the match
	# $4 -> the length of the match
	# $13 -> the length of the spacer
	# $13-$4+$5 -> the length of the spacer ($13) - the length of the match ($4) + the number of mismatches between the two ($5) to assign a "score" to the match
		# This adjusts the number of mismatches reported by BLAST to take into account the alignment length between the spacer and its hit
		# For example, an alignment length of 30 with 0 mismatches against a spacer length of 32 would have a "score" of 2
		# This allows us to sort the matches for downstream analyses, and not have a hit for example with 0 mismatches but a short alignment length be ranked higher than a hit with 2 mismatches but an alignment across the whole spacer

# The sort command then sorts the "scored" file in ascending order using column 7 (-k 7; the "score" column), which has a numerical value so we have to specify that with the "n"

echo "Scoring BLAST results"

for input in *.plus.txt; do
	filename="${input%%.*}"
	awk '{print $1"\t"$2"\t"$3"\t"$11"\t"$4"\t"$13"\t"$13-$4+$5}' "$input" > "$filename".plus.scored.txt
	sort -k 7n "$filename".plus.scored.txt > "$filename".plus.scored.sorted.txt
	done
	
for input in *.minus.txt; do
	filename="${input%%.*}"
	awk '{print $1"\t"$2"\t"$3"\t"$11"\t"$4"\t"$13"\t"$13-$4+$5}' "$input" > "$filename".minus.scored.txt
	sort -k 7n "$filename".minus.scored.txt > "$filename".minus.scored.sorted.txt
	done

#######################################
# Convert BLAST results to BED format #
#######################################

# This will let us use BEDTools downstream to pull out the target sequence and flanking sequence on either side for filtering purposes

# In this step, we're also adding flanking sequence to each hit to filter downstream for PAMs and taking into account the alignment length so we can extend properly into the PAM region

# The alignment=${2?Error: no alignment length specified} variable allows the user to specify what alignment length cut-off they would like to use during filtering

# The mismatch=${3?Error: no mismatch number specified} variable allows the user to specify how many mismatches they would like to use during filtering

# Breaking the awk commands down step by step starting with the first awk command:
	# if ($9<$10) print($2"\t"$9-$7-3"\t"$10+$13-$8+3"\t"$1"-"$2"-"$4"-"$13"-"$3"-"$5"-"$11"\t"$4"\t"$5) "$input" > "$filename".temp.plus.output.bed
		# If the start of the subject sequence is smaller than the end of the subject sequence, print
			# The subject ID in column 1
			# The subject start - the spacer start - 3 nt for a PAM in column 2
				# If the spacer alignment starts at let's say position 14, this will extend 14 nucleotides out and then an additional 3 for the PAM to account for the alignment length
			# The subject end + the length of the spacer - the spacer end + 3 nt into the PAM region in column 3
				# Again, this will extend the appropriate number of nucleotides out to finish up the match to account for alignment length and then also add a PAM region to the end of the match
			# The spacer ID, the subject ID, the alignment length, the query length, the %ID, the number of mismatches and the evalue, separated by dashes, in column 4 so we can retain this information at the end of the pipeline
			# The alignment length in column 5 to use for filtering within this loop
			# The number of mismatches in column 6 to use for the second step filtering within this loop
	# The else if step prints the same thing, excepts moves the end position to column 2 if it's smaller than the start position and then puts the start position in column 3
		# BED files need the start position to always be smaller than the end position

# Second awk command:
	# awk '{if ($5 >= '$alignment') print($1"\t"$2"\t"$3"\t"$4"\t"$6)}' "$filename".temp.plus.output.bed > "$filename".temp2.plus.output.bed
		# Take the temp file you just made and filter out any alignment length that is less than the user specified length
	# The else if step prints a "fake" BED entry with the subject name and the coordinates for the two nucleotides of sequence
		# Prevents the script from throwing an error during the extraction step
		# All this will do is extract two nucleotides of sequence, which won't get filtered in the filtering step
		# The main output file will be available, but the filtered files will simply be empty
		# Set at positions 8 and 10 to prevent potential errors from setting the position too close to the end of the sequence

# Third awk command:
	# awk '{if ($5 <= '$mismatch') print($1"\t"$2"\t"$3"\t"$4)}' "$filename".temp2.plus.output.bed > "$filename".temp3.plus.output.bed
		# Take the temp file you just made in the second awk step and filter out anything with more than the user specified number of mismatches between spacer/protospacer
		# This will create a BED file for the final extraction
	# The else if step does the same thing as the second awk command

echo "Converting BLAST results to BED files"

for input in *.masked.fasta.plus.txt; do
	filename="${input%%.*}"
	alignment=${2?Error: no alignment length specified}
	mismatch=${3?Error: no mismatch number specified}
	awk '{if ($9<$10) print($2"\t"$9-$7-3"\t"$10+$13-$8+3"\t"$1"-"$2"-"$4"-"$13"-"$3"-"$5"-"$11"\t"$4"\t"$5);
		else if ($9>$10) print($2"\t"$10+$13-$8+3"\t"$9-$7-3"\t"$1"-"$2"-"$4"-"$13"-"$3"-"$5"-"$11"\t"$4"\t"$5)}' "$input" > "$filename".temp.plus.output.bed
	awk '{if ($5 >= '$alignment') print($1"\t"$2"\t"$3"\t"$4"\t"$6);
		else if ($5 < '$alignment') print($1"\t"$2-$2+8"\t"$2-$2+10"\t"$4"\t"$6)}' "$filename".temp.plus.output.bed > "$filename".temp2.plus.output.bed
	awk '{if ($5 <= '$mismatch') print($1"\t"$2"\t"$3"\t"$4);
		else if ($5 > '$mismatch') print($1"\t"$2-$2+8"\t"$2-$2+10"\t"$4)}' "$filename".temp2.plus.output.bed > "$filename".temp3.plus.output.bed
	done

for input in *.masked.fasta.minus.txt; do
	filename="${input%%.*}"
	alignment=${2?Error: no alignment length specified}
	mismatch=${3?Error: no mismatch number specified}
	awk '{if ($9<$10) print($2"\t"$9+$7+2"\t"$10-$13+$8-4"\t"$1"-"$2"-"$4"-"$13"-"$3"-"$5"-"$11"\t"$4"\t"$5);
		else if ($9>$10) print($2"\t"$10-$13+$8-4"\t"$9+$7+2"\t"$1"-"$2"-"$4"-"$13"-"$3"-"$5"-"$11"\t"$4"\t"$5)}' "$input" > "$filename".temp.minus.output.bed
	awk '{if ($5 >= '$alignment') print($1"\t"$2"\t"$3"\t"$4"\t"$6);
		else if ($5 < '$alignment') print($1"\t"$2-$2+8"\t"$2-$2+10"\t"$4"\t"$6)}' "$filename".temp.minus.output.bed > "$filename".temp2.minus.output.bed
	awk '{if ($5 <= '$mismatch') print($1"\t"$2"\t"$3"\t"$4);
		else if ($5 > '$mismatch') print($1"\t"$2-$2+8"\t"$2-$2+10"\t"$4)}' "$filename".temp2.minus.output.bed > "$filename".temp3.minus.output.bed
	done

# In this step, this changes potential malformed BED entries that would have a negative value as their starting coordinate
	# This happens when the start coordinate is let's say 0, then it is extended for PAM filtering so it's now -3
	# This causes the script to abort

# Now, if the start coordinate is less than 0, it changes the start coordinate to 0
	# If the start coordinate is greater than or equal to 0, it leaves it as it is

for input in *.temp3.plus.output.bed; do
	filename="${input%%.*}"
	awk '{if ($2 < 0) print($1"\t"0"\t"$3"\t"$4"\t"$4"\t"$5);
		else if ($2 >= 0) print($1"\t"$2"\t"$3"\t"$4"\t"$4"\t"$5)}' "$input" > "$filename".plus.output.bed
	done
	
for input in *.temp3.minus.output.bed; do
	filename="${input%%.*}"
	awk '{if ($2 < 0) print($1"\t"0"\t"$3"\t"$4"\t"$4"\t"$5);
		else if ($2 >= 0) print($1"\t"$2"\t"$3"\t"$4"\t"$4"\t"$5)}' "$input" > "$filename".minus.output.bed
	done

#################################################
# Clean folder by creating a BLAST db directory #
#################################################

# Make a directory to move the blast db into since it's not needed for downstream applications

mkdir BLAST-db

mv *.fasta.ndb BLAST-db
mv *.fasta.nhr BLAST-db
mv *.fasta.nin BLAST-db
mv *.fasta.not BLAST-db
mv *.fasta.nsq BLAST-db
mv *.fasta.ntf BLAST-db
mv *.fasta.nto BLAST-db

#####################################
# Extract sequences using BED files #
#####################################

# Extract sequence from original fasta files using the BED file with the extended flank sequences
	# By adding the "name" flag, this will use the "name" column in the BED file as the header for the sequence extracted from the genome
	# This will allow us to use the spacer/protospacer match information as the header
	# Also makes it easier when looking at the PAM results to see if the spacer/protospacer match actually belongs in that bin
		# i.e.: If a I-C hit ends up in the I-F file, it's easy to pick out

# KNOWN ERROR MESSAGES THAT COULD OCCUR AT THIS STEP
	# Error when building the fai index that indexes the genome for sequence extraction
		# This occurs when the fasta file that is being fed into samtools has "ligatures" in the file
		# Honestly I'm not sure what exactly this entails, but it's a formatting problem that has been previously flagged in the samtools GitHub page
			# [fai_load] build FASTA index.
			# [fai_build_core] ignoring duplicate sequence "insert genome sequence header here" at byte offset "insert position here"
		# Upon inspection of the output files, it appears the samples that were flagged with this in the terminal page were processed and whatever duplication was caused by the genome fasta formatting was simply ignored

echo "Extracting sequence using BED files"

for input in *.masked.fasta; do
	filename="${input%%.*}"
	samtools faidx "$input" -o "$filename".fai
	bedtools getfasta -fi "$input" -bed "$filename".plus.output.bed -fo "$filename".plus.output.txt -name
	done
	
for input in *.masked.fasta; do
	filename="${input%%.*}"
	samtools faidx "$input" -o "$filename".fai
	bedtools getfasta -fi "$input" -bed "$filename".minus.output.bed -fo "$filename".minus.output.txt -name
	done

#########################################
# Filtering hits based on PAM sequences #
#########################################

# Grep will only call the "canonical" PAMs, but it's one way of parsing through the data quickly

# We can look manually for "non-canonical" PAMs in the output files

# Since we processed all our data at this point to keep plus and minus strand hits separate, we can use that information to orient our PAMs

# Searches for the PAM on the 5' end of the sequence for the plus strand hits
	# II-B is opposite to the I-C and I-F
		# i.e.: PAM on 3' end for the plus strand
	# II-B has a wild-card because the current canonical PAM is NGG on the 3' end
	# I-F has a wild-card because it's a dinucleotide PAM not a trinucleotide PAM
		# Without the wildcard character, the PAM filtering for I-F doesn't work properly
		
echo "Starting PAM filtering"

for input in *.plus.output.txt; do
	filename="${input%%.*}"
	grep -e '^TTC' "$input" -B 1 > "$filename".plus.temp.PAM.IC.txt
	grep -e '^.CC' "$input" -B 1 > "$filename".plus.temp.PAM.IF.txt
	grep -e '.GG$' "$input" -B 1 > "$filename".plus.temp.PAM.IIB.txt
	done

# Searches for the PAM on the 3' end of the sequence for the minus strand hits
	# II-B is opposite to the I-C and I-F
		# i.e.: PAM on 5' end for the minus strand
	# II-B has a wild-card because the current canonical PAM is CCN on the 5' end
	# I-F has a wild-card because it's a dinucleotide PAM not a trinucleotide PAM
		# Without the wildcard character, the PAM filtering for I-F doesn't work properly

for input in *.minus.output.txt; do
	filename="${input%%.*}"
	grep -e 'GAA$' "$input" -B 1 > "$filename".minus.temp.PAM.IC.txt
	grep -e 'GG.$' "$input" -B 1 > "$filename".minus.temp.PAM.IF.txt
	grep -e '^CC.' "$input" -B 1 > "$filename".minus.temp.PAM.IIB.txt
	done

# Does a second round of filtering to remove mis-sorted sequences from the final output file
	# For example, it's possible that some false positives get sorted into other bins, like a I-F hit in the I-C bin
	# Since the spacer headers contain information that includes what system it comes from, can use this for the second round of sorting

for input in *.plus.output.txt; do
	filename="${input%%.*}"
	grep -e 'IC' "$filename".plus.temp.PAM.IC.txt -A 1 > "$filename".plus.PAM.IC.txt
	grep -e 'IF' "$filename".plus.temp.PAM.IF.txt -A 1 > "$filename".plus.PAM.IF.txt
	grep -e 'IIB' "$filename".plus.temp.PAM.IIB.txt -A 1 > "$filename".plus.PAM.IIB.txt
	done

for input in *.minus.output.txt; do
	filename="${input%%.*}"
	grep -e 'IC' "$filename".minus.temp.PAM.IC.txt -A 1 > "$filename".minus.PAM.IC.txt
	grep -e 'IF' "$filename".minus.temp.PAM.IF.txt -A 1 > "$filename".minus.PAM.IF.txt
	grep -e 'IIB' "$filename".minus.temp.PAM.IIB.txt -A 1 > "$filename".minus.PAM.IIB.txt
	done

#############################################
# Moving outputs into their own directories #
#############################################

# Move the BLAST outputs and the extracted target sequences into a folder for easy retrieval

mkdir BLAST-output

mv *.fasta.minus.txt BLAST-output
mv *.fasta.plus.txt BLAST-output
mv *.output.txt BLAST-output
mv *.sorted.txt BLAST-output
mv *.scored.txt BLAST-output

# Move the filtered PAM hits into a folder for easy retrieval

mkdir PAM-hits

mv *.plus.PAM.*.txt PAM-hits
mv *.minus.PAM.*.txt PAM-hits

# Move the genomes into a folder for easy retrieval

mkdir Genomes

mv *.fasta Genomes

#####################################################
# Removing the temporary files used in the workflow #
#####################################################

# Remove temporary files to tidy the directory
# The 2> /dev/null part at the end of each command prevents the rm output from being printed to the terminal in the event an error occurs because the file doesn't exist

for input in *.fai; do
	filename="${input%%.*}"
	rm "$input" 2> /dev/null
	rm "$filename".temp1.fasta 2> /dev/null
	rm "$filename".temp2.fasta 2> /dev/null
	rm "$filename".temp3.fasta 2> /dev/null
	rm "$filename".temp4.fasta 2> /dev/null
	rm "$filename".temp5.fasta 2> /dev/null
	rm "$filename".minus.output.bed 2> /dev/null
	rm "$filename".plus.output.bed 2> /dev/null
	rm "$filename".temp.minus.output.bed 2> /dev/null
	rm "$filename".temp.plus.output.bed 2> /dev/null
	rm "$filename".temp2.minus.output.bed 2> /dev/null
	rm "$filename".temp2.plus.output.bed 2> /dev/null
	rm "$filename".temp3.minus.output.bed 2> /dev/null
	rm "$filename".temp3.plus.output.bed 2> /dev/null
	rm "$filename".plus.temp.PAM.IC.txt 2> /dev/null
	rm "$filename".plus.temp.PAM.IF.txt 2> /dev/null
	rm "$filename".plus.temp.PAM.IIB.txt 2> /dev/null
	rm "$filename".minus.temp.PAM.IC.txt 2> /dev/null
	rm "$filename".minus.temp.PAM.IF.txt 2> /dev/null
	rm "$filename".minus.temp.PAM.IIB.txt 2> /dev/null
	done
	
