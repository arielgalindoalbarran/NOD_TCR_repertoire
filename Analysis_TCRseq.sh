#!/bin/bash
#SBATCH -J TCRseq
#SBATCH -o TCR_O.txt
#SBATCH -e TCR_E.txt
#SBATCH -t 48:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=10
############################################## Important!:
# -> The fastq file need to be in compressed format,
# the read 1 finish in "R1.fastq.gz",
# and "R2.fastq.gz" for the read 2
# -> Fill or modify the next section 
############################################## Fill or modify this section:
#
Chain=""  #"Alpha" or "Beta"
PREFIX_NAME="" #To anotate in the consecutive analysis, for example number and name of sample: "1.WildType"
#
#Change the paths to the location of bin, Presto and MIXCR
export PATH=/here/your/path/:$PATH #For example:/home/.local/bin/:$PATH
export PATH=/here/your/path/:$PATH #For example:/home/.local/lib/python3.6/site-packages/presto:$PATH
export PATH=/here/your/path/:$PATH #For example:/home/Mixcr/mixcr-2.1.12/:$PATH
#Purge any previous modules and load Python, blast, Java, and R; Modify if is necessary:
module purge
module load system/Python-3.6.3
module load bioinfo/ncbi-blast-2.6.0+
module load system/Java8
module load system/R-3.5.1
#Load the VDJ tools, for example: "java -jar /home/vdjtools/vdjtools-1.1.10.jar"
VDJt="java -jar /pathway/here/to/VDJtools/vdjtools-1.1.10.jar" 
############################################################## End of fill and modifications
####################
WorkDir=$(pwd)
#Get the fastq files
gunzip *.gz
fastq_R1="*R1.fastq"
fastq_R2="*R2.fastq"
echo 'Analysis TCR using pRESTO and MIXCR, by Ariel Galindo'
#choose TCR chain 
if [ "$Chain" == "Alpha" ]; then
TRC="$WorkDir/fasta_files/TRAC.fasta"; TRChain="TRA"; TRV="$WorkDir/fasta_files/TRAV.fasta"; 
else
TRC="$WorkDir/fasta_files/TRBC.fasta"; TRChain="TRB"; TRV="$WorkDir/fasta_files/TRBV.fasta"; 
fi;
set -x;
echo 'Start pRESTO'
echo 'Filter reads with quality <20'
FilterSeq.py quality -s $fastq_R1 -q 20 --outname $PREFIX_NAME-R1 --log FS1.log --nproc 10;
FilterSeq.py quality -s $fastq_R2 -q 20 --outname $PREFIX_NAME-R2 --log FS2.log --nproc 10;
ParseLog.py -l FS1.log FS2.log -f ID QUALITY;
echo 'Mask the PCR primer in constant region and the switch template, TS=Template switch, CR= primer in the constant region.'
MaskPrimers.py score -s $PREFIX_NAME-R1_quality-pass.fastq -p TS.fasta --start 15 --barcode --mode cut --maxerror 0.5 --outname $PREFIX_NAME-R1 --log MP1.log --nproc 10;
MaskPrimers.py score -s $PREFIX_NAME-R2_quality-pass.fastq -p $TRC --start 0 --mode cut --outname $PREFIX_NAME-R2 --log MP2.log --nproc 10;
ParseLog.py -l MP1.log MP2.log -f ID PRIMER BARCODE ERROR;
echo 'Copy UMI from file 1 to file 2'
PairSeq.py -1 $PREFIX_NAME-R1_primers-pass.fastq -2 $PREFIX_NAME-R2_primers-pass.fastq --1f BARCODE --coord illumina;
echo 'Generating UMI consensus reads'  #maxerror=mismatch rate from consensus, prcons= (a) removes individual sequences that do not share a common primer 
BuildConsensus.py -s $PREFIX_NAME-R1_primers-pass_pair-pass.fastq --bf BARCODE --maxerror 0.1 --maxgap 0.5 --outname $PREFIX_NAME-R1 --log BC1.log --nproc 10;
BuildConsensus.py -s $PREFIX_NAME-R2_primers-pass_pair-pass.fastq --bf BARCODE --pf PRIMER --prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outname $PREFIX_NAME-R2 --log BC2.log --nproc 10;
ParseLog.py -l BC1.log BC2.log -f BARCODE SEQCOUNT CONSCOUNT PRCONS PRFREQ ERROR;
echo 'Paired-end assembly of UMI consensus sequences'
echo 'Syncronizing paired-end files'
PairSeq.py -1 $PREFIX_NAME-R1_consensus-pass.fastq -2 $PREFIX_NAME-R2_consensus-pass.fastq --coord presto;
echo 'Assembling UMI consensus mate-pairs'	
AssemblePairs.py sequential -1 $PREFIX_NAME-R2_consensus-pass_pair-pass.fastq -2 $PREFIX_NAME-R1_consensus-pass_pair-pass.fastq -r $TRV --coord presto --rc tail --scanrev --2f CONSCOUNT --1f CONSCOUNT PRCONS --aligner blastn --outname $PREFIX_NAME-C --log AP.log --nproc 10;
ParseLog.py -l AP.log -f ID REFID LENGTH OVERLAP GAP ERROR IDENTITY;
echo 'Combining UMI read group size annotations'
ParseHeaders.py collapse -s $PREFIX_NAME-C_assemble-pass.fastq -f CONSCOUNT --act min;
echo 'Removal of duplicate sequences'
CollapseSeq.py -s $PREFIX_NAME-C_assemble-pass_reheader.fastq -n 20 --inner --uf CREGION --cf CONSCOUNT --act sum --outname $PREFIX_NAME-C;
#
echo 'Filtering to sequences with at least two representative reads'
SplitSeq.py group -s $PREFIX_NAME-C_collapse-unique.fastq -f CONSCOUNT --num 2 --outname $PREFIX_NAME-C_2;
echo 'Creating an annotation table'
ParseHeaders.py table -s $PREFIX_NAME-C_2_atleast-2.fastq -f ID CREGION CONSCOUNT DUPCOUNT;
echo 'Align with TCR using MiXCR'
mixcr align -t 10 --library imgt --species mmu -c $TRChain -r $PREFIX_NAME-align_2.log.txt $PREFIX_NAME-C_2_atleast-2.fastq $PREFIX_NAME-alignments_2.vdjca;
mixcr assemble -t 10 -r $PREFIX_NAME-assemble_2.log.txt $PREFIX_NAME-alignments_2.vdjca $PREFIX_NAME-clones_2.clns;
mixcr exportClones $PREFIX_NAME-clones_2.clns $PREFIX_NAME-clones_2.txt;
echo 'Export to VDJtools'
java -jar /home/agalindo/work/vdjtools/vdjtools-1.1.10.jar Convert -S mixcr $PREFIX_NAME-clones_2.txt vdjtools
#STATISTICS:
printf "Statistics of reads, MIGS, and clonotypes for $PREFIX_NAME \n" > Statistics_$PREFIX_NAME.txt;
printf "Number of start analysis reads:\n" >> Statistics_$PREFIX_NAME.txt;
awk '{s++}END{print s/4}' *_R1.fastq >> Statistics_$PREFIX_NAME.txt;
printf "Number of reads with quality higher than 20:\n" >> Statistics_$PREFIX_NAME.txt;
awk '{s++}END{print s/4}' $PREFIX_NAME-R1_quality-pass.fastq >> Statistics_$PREFIX_NAME.txt;
printf "Reads with Template switch and constant region:\n" >> Statistics_$PREFIX_NAME.txt;
awk '{s++}END{print s/4}' $PREFIX_NAME-R1_primers-pass_pair-pass.fastq >> Statistics_$PREFIX_NAME.txt;
printf "Reads aligned:\n" >> Statistics_$PREFIX_NAME.txt;
grep 'CONSCOUNT' $PREFIX_NAME-R1_consensus-pass_pair-pass.fastq | awk -F '=' '{print $2 }' | awk '{s+=$1} END {printf "%.0f\n", s}' >> Statistics_$PREFIX_NAME.txt;
printf "Number of MIGs Pre- filter:\n" >> Statistics_$PREFIX_NAME.txt;
grep 'CONSCOUNT' $PREFIX_NAME-C_assemble-pass.fastq | awk '{s++}END{print s}' >> Statistics_$PREFIX_NAME.txt;
printf "Number of MIGs Post- filter at least 2 read:\n" >> Statistics_$PREFIX_NAME.txt;
awk '{s++}END{print s/4}' $PREFIX_NAME-C_2_atleast-2.fastq >> Statistics_$PREFIX_NAME.txt;
printf "MIGS Aligned to TCR locus, and used to clonotypes, at least 2 reads:\n" >> Statistics_$PREFIX_NAME.txt;
awk '{print $1 }' vdjtools.$PREFIX_NAME-clones_2.txt | awk '{s+=$1} END {printf "%.0f\n", s}' >> Statistics_$PREFIX_NAME.txt;
printf "Number of clonotypes, at least 2 reads:\n" >> Statistics_$PREFIX_NAME.txt;
awk '{print $1 }' vdjtools.$PREFIX_NAME-clones_2.txt | wc -l >> Statistics_$PREFIX_NAME.txt;
printf "Aligment percentages, at least 2 read:\n" >> Statistics_$PREFIX_NAME.txt;
sed -n 8p $PREFIX_NAME-align_2.log.txt >> Statistics_$PREFIX_NAME.txt;
sed -n 9p $PREFIX_NAME-align_2.log.txt >> Statistics_$PREFIX_NAME.txt;
sed -n 10p $PREFIX_NAME-align_2.log.txt >> Statistics_$PREFIX_NAME.txt;
sed -n 11p $PREFIX_NAME-align_2.log.txt >> Statistics_$PREFIX_NAME.txt;
sed -n 16p $PREFIX_NAME-align_2.log.txt >> Statistics_$PREFIX_NAME.txt;
printf "Good luck with the analysis. Ariel Galindo" >> Statistics_$PREFIX_NAME.txt;
#Generate a table to make histograms:
grep 'CONSCOUNT' $PREFIX_NAME-C_2_atleast-2.fastq | awk -F "|" '{print $3}' | awk -F "=" '{print $2}' |  sort -n -r > $PREFIX_NAME-umiVSreads_atleast-2.csv;
#Move and copy files
gzip *.fastq;
gzip *.log;
gzip *.tab;
gzip *.vdjca; 
mkdir results;
cp vdjtools* results;
mv *.csv results; 
unset fastq_R1;
unset fastq_R2;
unset PREFIX_NAME;
####################################################
#Graphics
####################################################
#### DiversityStats
$VDJt CalcDiversityStats vdjtools.$PREFIX_NAME-clones_2.txt DivStats_$PREFIX_NAME;
#### Spectratype:
Rscript $WorkDir/Rscripts/spectratype.R $Chain $PREFIX_NAME vdjtools.$PREFIX_NAME-clones_2.txt spectratype_$PREFIX_NAME-$Chain.svg;
## circo plot vj 
$VDJt PlotFancyVJUsage vdjtools.$PREFIX_NAME-clones_2.txt CircoPlotVJ_2reads_$PREFIX_NAME;
rm *fancyvj.wt.pdf vj_pairing_plot.r;
Rscript $WorkDir/Rscripts/VJ_Circlopolt.R $Chain CircoPlotVJ_2reads_$PREFIX_NAME.fancyvj.wt.txt CircoPlotVJ_2reads_$PREFIX_Cond2.fancyvj.wt.txt CircoPlotVJ_2reads_$PREFIX_NAME-$Chain.svg CircoPlotVJ_2reads_$PREFIX_Cond2-$Chain.svg;
#Pie chart V segments
Rscript $WorkDir/Rscripts/Piechart_Vsegments.R $Chain vdjtools.$PREFIX_NAME-clones_2.txt vdjtools.$PREFIX_Cond2-clones_2.txt PieCharttV_2reads_$PREFIX_NAME-$Chain.svg PieCharttV_2reads_$PREFIX_Cond2-$Chain.svg
rm Rplots.pdf
#Histograms UMIs vs Clonotypes, and UMIs vs Reads
Rscript $WorkDir/Rscripts/Histograms.R $PREFIX_NAME $PREFIX_NAME-umiVSreads_atleast-2.csv $PREFIX_NAME-umiVSreads_atleast-2.csv vdjtools.$PREFIX_NAME-clones_2.txt vdjtools.$PREFIX_NAME-clones_2.txt;
#Move files:
mkdir Graphics;
mv spectratype_2* Graphics;
mv RarefactionPlot_2* Graphics;
mv DivStats_2* Graphics;
mv vjusage_2* Graphics;
mv V-Heatmap_2reads* Graphics;
mv J-Heatmap_2reads* Graphics;
mv CircoPlotVJ_2* Graphics;
mv PieCharttV_2* Graphics;
mv *UMIvsReads_2read.svg Graphics;
mv *UMIvsClonotype_2read.svg Graphics;
rm *segments.txt;
unset VDJt;
unset PREFIX_NAME;
unset PREFIX_Cond2;
unset WorkDir























