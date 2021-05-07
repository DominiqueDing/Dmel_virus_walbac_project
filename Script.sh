########################
### Trimming process ###
########################
# Trim with trimgalore
./trim_galore --cores 4 \
--quality 5 \
--phred33 \
--fastqc \
--illumina \
--paired \
--length 50 \
--output_dir ./ \
--retain_unpaired \
SLX-17476.RNAA027.H3G27BBXY.s_3.r_1.fq.gz \
SLX-17476.RNAA027.H3G27BBXY.s_3.r_2.fq.gz


#######################
### Mapping process ###
#######################

# Creat genome indices
./STAR --runThreadN 40 \
--runMode genomeGenerate \
--genomeDir Directory_to/genome_indices \
--genomeFastaFiles Directory_to_genome/genome.fasta \
--sjdbGTFfile Directory_to_genome/genome.gtf \
--sjdbOverhang 100

# Map with STAR
./STAR --runThreadN 46 \
--outSAMprimaryFlag AllBestScore \
--outSAMtype BAM SortedByCoordinate \
--genomeDir Directory_to/genome_indices \
--readFilesIn SLX-17476.RNAA027.H3G27BBXY.s_3.r_1.fq.gz SLX-17476.RNAA027.H3G27BBXY.s_3.r_2.fq.gz \
--readFilesCommand zcat \
--outReadsUnmapped Fastx \
--quantMode GeneCounts \
--outFileNamePrefix output_dir/




#######################
### Virus discovery ###
#######################

# Get rid of rRNA reads before trinity assembly
bowtie2 \
--quiet \
--very-sensitive-local \
--phred33  \
-x Directory _to_rRNA_DB/rRNA \
-1 Unmapped.out.mate1 \
-2 Unmapped.out.mate2 \
--threads 12 \
--met-file metrics.txt \
--al-conc-gz blacklist_paired_aligned_fq.gz \
--un-conc-gz blacklist_paired_unaligned_fq.gz \
--al-gz blacklist_unpaired_aligned_fq.gz \
--un-gz blacklist_unpaired_unaligned_fq.gz

# Trinity assembly:
Trinity \
--left blacklist_paired_unaligned_1.fq \
--right blacklist_paired_unaligned_2.fq \
--seqType fq \
--SS_lib_type RF \
--max_memory 20G \
--min_contig_length 100 \
--output output_dir/ \
--CPU 35

# Select candidate novel virus sequences (pipeline) 
TransDecoder.LongOrfs \
-t Trinity.fasta \
-m 30

diamond blastx \
-p 40 \
-d Directory_to_virus_DB/viral.dmnd \
-q Trinity_new.fasta \
-o ncbi_viruses.out \
-f 6 qseqid evalue length stitle


cat ncbi_viruses.out | cut -f 1 | sort | uniq | grep -A 1 --no-group-separator -F -f - Trinity_new.fasta > candidates.fasta

diamond blastx \
-p 40 \
-d Directory_to_nr_DB/nr.dmnd \
-q candidates.fasta \
-o candidates_ncbi_nr.out \
-f 6 qseqid evalue length stitle

grep ">" candidates.fasta | grep -o "^[^ ]*" | sed 's/>//' > candnames.txt
while IFS= read -r candidate
do
	grep -m 1 "${candidate}" candidates_ncbi_nr.out >> topblasthits.txt
done < candnames.txt
grep "[Vv]irus" topblasthits.txt | grep -v "transpos*" | cut -f 1 | sort | uniq > verified_viruses.txt
cat verified_viruses.txt | grep -A 1 --no-group-separator -F -f - Trinity_new.fasta > verified_viruses.fasta



