0. Download RNAseq data

1. Trim data with fastp and align your RNA-seq data to your genome with HISAT2
https://www.reneshbedre.com/blog/hisat2-sequence-aligner.html



           1. Build index on genome assembly

                      sbatch --partition=pibu_el8 --job-name=H1Hisatindex --time=0-21:00:00 --mem-per-cpu=16G --ntasks=12 --cpus-per-task=1 --output=Hap1_HiSat2index.log --error=Hap1_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/01_hap1; module load HISAT2; hisat2-build out_JBAT_hap1.FINAL.fa.mod.MAKER.softmasked hap1_hisat_index -p 12"

                     sbatch --partition=pibu_el8 --job-name=H2Hisatindex --time=0-21:00:00 --mem-per-cpu=16G --ntasks=12 --cpus-per-task=1 --output=Hap2_HiSat2index.log --error=Hap2_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/02_hap2; module load HISAT2; hisat2-build out_JBAT_hap2.FINAL.fa.mod.MAKER.softmasked hap2_hisat_index -p 12"

   2. Map reads to genome


                sbatch --partition=pibu_el8 --job-name=H1HisatMap --time=3-21:00:00 --mem-per-cpu=16G --ntasks=32 --cpus-per-task=1 --output=Hap1_HiSatMap.log --error=Hap1_HiSatMap.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/02_HISAT2_mapping/01_Hap1; module load HISAT2; hisat2 --phred33 --dta -x hap1_hisat_index -1 /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/1A_1_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/2A_1_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/3A_1_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/4A_1_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/5A_1_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/6A_1_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/7A_1_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/8A_1_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/9A_1_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/10A_1_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/11A_1_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/12A_1_trimmed.fastq.gz -2 /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/1A_2_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/2A_2_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/3A_2_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/4A_2_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/5A_2_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/6A_2_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/7A_2_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/8A_2_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/9A_2_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/10A_2_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/11A_2_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/12A_2_trimmed.fastq.gz -S ALLSamples_Hap1.sam -p 36"

                sbatch --partition=pibu_el8 --job-name=H1SAMTOOLS --time=0-21:00:00 --mem-per-cpu=16G --ntasks=32 --cpus-per-task=1 --output=Hap1_SAMTOOLS.log --error=Hap1__SAMTOOLS.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/02_HISAT2_mapping/01_Hap1; module load SAMtools; samtools view --threads 12 -b -o ALLSamples_Hap1.bam ALLSamples_Hap1.sam; samtools sort -m 7G -o ALLSamples_Hap1_sorted.bam -T ALLSamples_Hap1_temp --threads 32 ALLSamples_Hap1.bam"

                
        2b. convert to bam and sort
                
                module load samtools
                samtools view --threads ${PROC} -b -o ${R1_FQ%.*}.bam ${R1_FQ%.*}.sam
                samtools sort -m 7G -o ${R1_FQ%.*}_sorted.bam -T ${R1_FQ%.*}_temp --threads ${PROC} ${R1_FQ%.*}.bam

