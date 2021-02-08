echo "##### PIPELINE FOR RNA-Seq & VARIANT CALLING ##################"

echo ""
echo "Checking Quality and Trimming..."
echo ""
echo "Running FastQC..."
fastqc *.fastq.gz

echo "Moving FastQC results"
mkdir fastqc_untrimmed/

mv *.zip fastqc_untrimmed/
mv *.html fastqc_untrimmed/

echo "Finished moving FastQC results"

echo ""
echo "NOW tRIMMING..."
echo "Running Trimmomatic..."
for filename in *.fastq.gz
do
 # first. getting the unique part of the file name alone by removing .fastq.gz
 name=$(basename $filename .fastq.gz)
 
 echo "On sample: $name"
 
 trimmomatic.jar SE ${name}.fastq.gz \
  ${name}.qc.fastq.gz \
  ILLUMINACLIP:TruSeq2-SE.fa:2:0:15 \
  LEADING:15 TRAILING:15 \
  SLIDINGWINDOW:10:20 \
  MINLEN:25
done

echo ""
echo "Running FastQC on trimmed sequences..."
fastqc *.qc.fastq.gz

echo ""
echo "Moving FastQC results of the trimmed reads to 'fastqc_trimmed' folder"
mkdir fastqc_trimmed

mv *.zip fastqc_trimmed
mv *.html fastqc_trimmed

echo "########## RNA-Seq: READS QUANTIFICATION WITH SALMON ##################"
ech ""
echo "Indexing ref. transcriptome for salmon..."
salmon index --index sc_index --transcripts GCA_000146045.2_R64_rna_from_genomic.fna.gz

echo ""
echo "Quantifying with salmon..."
for filename in *.qc.fastq.gz
do
  # Getting the unique part of the file name alone by removing ".qc.fastq.gz"
  name=$(basename $filename .qc.fastq.gz)
  
  salmon quant -i sc_index --libType A \
    -r ${filename} -o  ${name}_quant --seqBias \
    --gcBias --validateMappings
done

echo "############### NOW VARIANT CALLING #######################"
echo ""
echo "Preparing index for bwa..."
bwa index orf_coding.fasta

echo ""
echo "Now mapping with bwa mem..."

echo "indexing the ref. fasta file for samtools variant calling..."
samtools faidx orf_coding.fasta

echo ""
for filename in *.qc.fastq.gz
do
  # get the unique part of the file name by removing .qc.fastq.gz
  name=$(basename $filename .qc.fastq.gz)
  
  echo "Mapping $name ..."
  bwa mem -t 4 orf_coding.fasta $filename > ${name}.sam
  
  echo "Converting sam to BAM..."
  samtools view -S -b ${name}.sam > ${name}.bam
  
  echo "Sorting the BAM file..."
  samtools sort ${name}.bam -o ${name}.sorted.bam
  
  echo "Indexing the BAM file..."
  samtools index ${name}.sorted.bam
  
  echo ""
  echo "Now calling variants..."
  bcftools mpileup -O b -f orf_coding.fasta ${name}.sorted.bam | \
    bcftools call -m -v -o ${name}_variants.vcf
  
  echo ""
  echo "Filtering out variants..."
  vcfutils.pl varFilter ${name}_variants.vcf > ${name}_variants_filtered.vcf
  echo ""
done

echo ""
echo "END!!"
