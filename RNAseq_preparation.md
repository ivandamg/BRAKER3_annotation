

# 0. Build index on genome assembly

                      sbatch --partition=pibu_el8 --job-name=H1Hisatindex --time=0-21:00:00 --mem-per-cpu=16G --ntasks=12 --cpus-per-task=1 --output=Hap1_HiSat2index.log --error=Hap1_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/01_hap1; module load HISAT2; hisat2-build out_JBAT_hap1.FINAL.fa.mod.MAKER.softmasked hap1_hisat_index -p 12"

                     sbatch --partition=pibu_el8 --job-name=H2Hisatindex --time=0-21:00:00 --mem-per-cpu=16G --ntasks=12 --cpus-per-task=1 --output=Hap2_HiSat2index.log --error=Hap2_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/02_hap2; module load HISAT2; hisat2-build out_JBAT_hap2.FINAL.fa.mod.MAKER.softmasked hap2_hisat_index -p 12"


# 1. own data

## 1. Trim data with fastp and align your RNA-seq data to your genome with HISAT2
https://www.reneshbedre.com/blog/hisat2-sequence-aligner.html


                        cd /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/00_rawreads;
                        for FILE in $(ls *1.fastq.gz); do echo $FILE; sbatch --partition=pshort_el8 --job-name=$(echo $FILE | cut -d'_' -f1)fastp --time=0-04:00:00 --mem-per-cpu=24G --ntasks=1 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1)_fastp.out --error=$(echo $FILE | cut -d'_' -f1)_fastp.error --mail-type=END,FAIL --wrap " cd /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/00_rawreads ; module load FastQC; ~/00_Software/fastp --in1 $FILE --in2 $(echo $FILE | cut -d'_' -f1)_2.fastq.gz --out1 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_1_trimmed.fastq.gz --out2 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_2_trimmed.fastq.gz -h ../02_TrimmedData/$(echo $FILE | cut -d',' -f1)_fastp.html --thread 4; fastqc -t 4 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_1_trimmed.fastq.gz; fastqc -t 4 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_2_trimmed.fastq.gz"; sleep  1; done



## 2. Map reads to genome and compresss

1. Hap1

                sbatch --partition=pibu_el8 --job-name=FrH1Hisatmap --time=3-21:00:00 --mem-per-cpu=16G --ntasks=48 --cpus-per-task=1 --output=FrHap1_HiSat2index.log --error=FrHap1_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/10_RNAseq/00_trimmedReads; module load HISAT2; hisat2 --phred33 --dta -x /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/01_hap1/hap1_hisat_index -1 1A_1_trimmed.fastq.gz,2A_1_trimmed.fastq.gz,3A_1_trimmed.fastq.gz,4A_1_trimmed.fastq.gz,5A_1_trimmed.fastq.gz,6A_1_trimmed.fastq.gz,7A_1_trimmed.fastq.gz,8A_1_trimmed.fastq.gz,9A_1_trimmed.fastq.gz,10A_1_trimmed.fastq.gz,11A_1_trimmed.fastq.gz,12A_1_trimmed.fastq.gz -2 1A_2_trimmed.fastq.gz,2A_2_trimmed.fastq.gz,3A_2_trimmed.fastq.gz,4A_2_trimmed.fastq.gz,5A_2_trimmed.fastq.gz,6A_2_trimmed.fastq.gz,7A_2_trimmed.fastq.gz,8A_2_trimmed.fastq.gz,9A_2_trimmed.fastq.gz,10A_2_trimmed.fastq.gz,11A_2_trimmed.fastq.gz,12A_2_trimmed.fastq.gz -S /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/01_hap1/Ufribourg_12samples_Hap1.sam -p 48"

                sbatch --partition=pibu_el8 --job-name=FrH1SAMTOOLS --time=0-21:00:00 --mem-per-cpu=128G --ntasks=48 --cpus-per-task=1 --output=FrHap1_SAMTOOLS.log --error=FrHap1_SAMTOOLS.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/01_hap1/; module load SAMtools; samtools view --threads 48 -b -o Ufribourg_12samples_Hap1.bam Ufribourg_12samples_Hap1.sam; samtools sort -o Ufribourg_12samples_Hap1_sorted.bam -T Ufribourg_12samples_Hap1_temp --threads 48 Ufribourg_12samples_Hap1.bam"

3. Hap2
   
                sbatch --partition=pibu_el8 --job-name=FrH2Hisatmap --time=3-21:00:00 --mem-per-cpu=16G --ntasks=48 --cpus-per-task=1 --output=FrHap2_HiSat2index.log --error=FrHap2_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData; module load HISAT2; hisat2 --phred33 --dta -x /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/02_hap2/hap2_hisat_index -1 1A_1_trimmed.fastq.gz,2A_1_trimmed.fastq.gz,3A_1_trimmed.fastq.gz,4A_1_trimmed.fastq.gz,5A_1_trimmed.fastq.gz,6A_1_trimmed.fastq.gz,7A_1_trimmed.fastq.gz,8A_1_trimmed.fastq.gz,9A_1_trimmed.fastq.gz,10A_1_trimmed.fastq.gz,11A_1_trimmed.fastq.gz,12A_1_trimmed.fastq.gz -2 1A_2_trimmed.fastq.gz,2A_2_trimmed.fastq.gz,3A_2_trimmed.fastq.gz,4A_2_trimmed.fastq.gz,5A_2_trimmed.fastq.gz,6A_2_trimmed.fastq.gz,7A_2_trimmed.fastq.gz,8A_2_trimmed.fastq.gz,9A_2_trimmed.fastq.gz,10A_2_trimmed.fastq.gz,11A_2_trimmed.fastq.gz,12A_2_trimmed.fastq.gz -S /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/02_hap2/Ufribourg_12samples_Hap2.sam -p 48"

                sbatch --partition=pibu_el8 --job-name=FrH2SAMTOOLS --time=0-21:00:00 --mem-per-cpu=256G --ntasks=48 --cpus-per-task=1 --output=FrHap2_SAMTOOLS.log --error=FrHap2_SAMTOOLS.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/02_hap2/; module load SAMtools; samtools view --threads 48 -b -o Ufribourg_12samples_Hap2.bam Ufribourg_12samples_Hap2.sam; samtools sort -o Ufribourg_12samples_Hap2_sorted.bam -T Ufribourg_12samples_Hap2_temp --threads 48 Ufribourg_12samples_Hap2.bam"



# 2. Public available data

## 1. Download public data with SRAtools

          # Download raw data
          for i in $(echo SRX18026584 SRX18026583 SRX18026582 SRX18026581 SRX18026580 SRX18026579 SRX18026578 SRX18026577 SRX18026576 SRX18026575 SRX18026574 SRX18026573 SRX18026572 SRX18026571 SRX18026570 SRX18026569 SRX18026568 SRX18026567 SRX18026566 SRX18026565 SRX18026564 SRX18026563 SRX18026562 SRX18026561 SRX18026560 SRX18026559 SRX18026558 SRX18026557 SRX18026556 SRX18026555); do echo $i; sbatch --partition=pshort_el8 --job-name=SRAtools1--time=0-04:00:00 --mem-per-cpu=4G --ntasks=1 --cpus-per-task=1 --output=SRAtools1.out --error=SRAtools1.error --mail-type=END,FAIL --wrap "module load SRA-Toolkit;cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/00_rawreads/; prefetch $i"; done

           # mv to parent folder
          mv SRR*/*.sra . 

          # split en paired fastq files
          for i in $(ls SRR*.sra); do echo $i;sbatch --partition=pshort_el8 --job-name=SRAtools2--time=0-04:00:00 --mem-per-cpu=4G --ntasks=1 --cpus-per-task=1 --output=SRAtools2.out --error=SRAtools2.error --mail-type=END,FAIL --wrap "module load SRA-Toolkit;cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/00_rawreads/; fastq-dump --split-files $(echo $i | cut -d'.' -f1)"; done

          # gzip fastq files
          for i in $(ls SRR*.fastq); do echo $i; sbatch --partition=pshort_el8 --job-name=gzip --time=0-04:00:00 --mem-per-cpu=4G --ntasks=1 --cpus-per-task=1 --output=SRAtools1.out --error=SRAtools1.error --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/00_rawreads/; gzip $i"; done

## 2. Trim data with fastp 

                        
                        cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/00_rawreads;
                        for FILE in $(ls *1.fastq.gz); do echo $FILE; sbatch --partition=pshort_el8 --job-name=$(echo $FILE | cut -d'_' -f1)fastp --time=0-04:00:00 --mem-per-cpu=24G --ntasks=1 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1)_fastp.out --error=$(echo $FILE | cut -d'_' -f1)_fastp.error --mail-type=END,FAIL --wrap " cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/00_rawreads ; module load FastQC; ~/00_Software/fastp --in1 $FILE --in2 $(echo $FILE | cut -d'_' -f1)_2.fastq.gz --out1 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_1_trimmed.fastq.gz --out2 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_2_trimmed.fastq.gz -h ../02_TrimmedData/$(echo $FILE | cut -d',' -f1)_fastp.html --thread 4; fastqc -t 4 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_1_trimmed.fastq.gz; fastqc -t 4 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1)_2_trimmed.fastq.gz"; sleep  1; done


## 3. Map reads to genome, sort and compress


1. Hap1

              # run from /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/02_TrimmedData
            sbatch --partition=pibu_el8 --job-name=H1Hisatmap3 --time=7-21:00:00 --mem-per-cpu=64G --ntasks=48 --cpus-per-task=1 --output=Hap1_HiSat2index.log --error=Hap1_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/02_TrimmedData; module load HISAT2; hisat2 --phred33 --dta -x /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/01_hap1/hap1_hisat_index -1 SRR22045430_1_trimmed.fastq.gz,SRR22045431_1_trimmed.fastq.gz,SRR22045432_1_trimmed.fastq.gz,SRR22045433_1_trimmed.fastq.gz,SRR22045434_1_trimmed.fastq.gz,SRR22045435_1_trimmed.fastq.gz,SRR22045436_1_trimmed.fastq.gz,SRR22045437_1_trimmed.fastq.gz,SRR22045438_1_trimmed.fastq.gz,SRR22045439_1_trimmed.fastq.gz,SRR22045440_1_trimmed.fastq.gz,SRR22045441_1_trimmed.fastq.gz,SRR22045442_1_trimmed.fastq.gz,SRR22045443_1_trimmed.fastq.gz,SRR22045444_1_trimmed.fastq.gz,SRR22045445_1_trimmed.fastq.gz,SRR22045446_1_trimmed.fastq.gz,SRR22045447_1_trimmed.fastq.gz,SRR22045448_1_trimmed.fastq.gz,SRR22045449_1_trimmed.fastq.gz,SRR22045450_1_trimmed.fastq.gz,SRR22045451_1_trimmed.fastq.gz,SRR22045452_1_trimmed.fastq.gz,SRR22045453_1_trimmed.fastq.gz,SRR22045454_1_trimmed.fastq.gz,SRR22045455_1_trimmed.fastq.gz,SRR22045456_1_trimmed.fastq.gz,SRR22045457_1_trimmed.fastq.gz,SRR22045458_1_trimmed.fastq.gz -2 SRR22045430_2_trimmed.fastq.gz,SRR22045431_2_trimmed.fastq.gz,SRR22045432_2_trimmed.fastq.gz,SRR22045433_2_trimmed.fastq.gz,SRR22045434_2_trimmed.fastq.gz,SRR22045435_2_trimmed.fastq.gz,SRR22045436_2_trimmed.fastq.gz,SRR22045437_2_trimmed.fastq.gz,SRR22045438_2_trimmed.fastq.gz,SRR22045439_2_trimmed.fastq.gz,SRR22045440_2_trimmed.fastq.gz,SRR22045441_2_trimmed.fastq.gz,SRR22045442_2_trimmed.fastq.gz,SRR22045443_2_trimmed.fastq.gz,SRR22045444_2_trimmed.fastq.gz,SRR22045445_2_trimmed.fastq.gz,SRR22045446_2_trimmed.fastq.gz,SRR22045447_2_trimmed.fastq.gz,SRR22045448_2_trimmed.fastq.gz,SRR22045449_2_trimmed.fastq.gz,SRR22045450_2_trimmed.fastq.gz,SRR22045451_2_trimmed.fastq.gz,SRR22045452_2_trimmed.fastq.gz,SRR22045453_2_trimmed.fastq.gz,SRR22045454_2_trimmed.fastq.gz,SRR22045455_2_trimmed.fastq.gz,SRR22045456_2_trimmed.fastq.gz,SRR22045457_2_trimmed.fastq.gz,SRR22045458_2_trimmed.fastq.gz -S /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/01_hap1/PRJNA863910_30samples_Hap1.sam -p 48"

                sbatch --partition=pibu_el8 --job-name=H1SAMTOOLS --time=0-21:00:00 --mem-per-cpu=64G --ntasks=48 --cpus-per-task=1 --output=Hap1_SAMTOOLS.log --error=Hap1_SAMTOOLS.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/01_hap1/; module load SAMtools; samtools view --threads 48 -b -o PRJNA863910_30samples_Hap1.bam PRJNA863910_30samples_Hap1.sam; samtools sort -o PRJNA863910_30samples_Hap1_sorted.bam -T PRJNA863910_30samples_Hap1_temp --threads 48 PRJNA863910_30samples_Hap1.bam"

2. Hap2

              # run from /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/02_TrimmedData
            sbatch --partition=pibu_el8 --job-name=H2Hisatmap3 --time=3-21:00:00 --mem-per-cpu=64G --ntasks=48 --cpus-per-task=1 --output=Hap2_HiSat2index.log --error=Hap2_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/PRJNA863910_Graft_Tzeela_2022/02_TrimmedData; module load HISAT2; hisat2 --phred33 --dta -x /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/02_hap2/hap2_hisat_index -1 SRR22045430_1_trimmed.fastq.gz,SRR22045431_1_trimmed.fastq.gz,SRR22045432_1_trimmed.fastq.gz,SRR22045433_1_trimmed.fastq.gz,SRR22045434_1_trimmed.fastq.gz,SRR22045435_1_trimmed.fastq.gz,SRR22045436_1_trimmed.fastq.gz,SRR22045437_1_trimmed.fastq.gz,SRR22045438_1_trimmed.fastq.gz,SRR22045439_1_trimmed.fastq.gz,SRR22045440_1_trimmed.fastq.gz,SRR22045441_1_trimmed.fastq.gz,SRR22045442_1_trimmed.fastq.gz,SRR22045443_1_trimmed.fastq.gz,SRR22045444_1_trimmed.fastq.gz,SRR22045445_1_trimmed.fastq.gz,SRR22045446_1_trimmed.fastq.gz,SRR22045447_1_trimmed.fastq.gz,SRR22045448_1_trimmed.fastq.gz,SRR22045449_1_trimmed.fastq.gz,SRR22045450_1_trimmed.fastq.gz,SRR22045451_1_trimmed.fastq.gz,SRR22045452_1_trimmed.fastq.gz,SRR22045453_1_trimmed.fastq.gz,SRR22045454_1_trimmed.fastq.gz,SRR22045455_1_trimmed.fastq.gz,SRR22045456_1_trimmed.fastq.gz,SRR22045457_1_trimmed.fastq.gz,SRR22045458_1_trimmed.fastq.gz -2 SRR22045430_2_trimmed.fastq.gz,SRR22045431_2_trimmed.fastq.gz,SRR22045432_2_trimmed.fastq.gz,SRR22045433_2_trimmed.fastq.gz,SRR22045434_2_trimmed.fastq.gz,SRR22045435_2_trimmed.fastq.gz,SRR22045436_2_trimmed.fastq.gz,SRR22045437_2_trimmed.fastq.gz,SRR22045438_2_trimmed.fastq.gz,SRR22045439_2_trimmed.fastq.gz,SRR22045440_2_trimmed.fastq.gz,SRR22045441_2_trimmed.fastq.gz,SRR22045442_2_trimmed.fastq.gz,SRR22045443_2_trimmed.fastq.gz,SRR22045444_2_trimmed.fastq.gz,SRR22045445_2_trimmed.fastq.gz,SRR22045446_2_trimmed.fastq.gz,SRR22045447_2_trimmed.fastq.gz,SRR22045448_2_trimmed.fastq.gz,SRR22045449_2_trimmed.fastq.gz,SRR22045450_2_trimmed.fastq.gz,SRR22045451_2_trimmed.fastq.gz,SRR22045452_2_trimmed.fastq.gz,SRR22045453_2_trimmed.fastq.gz,SRR22045454_2_trimmed.fastq.gz,SRR22045455_2_trimmed.fastq.gz,SRR22045456_2_trimmed.fastq.gz,SRR22045457_2_trimmed.fastq.gz,SRR22045458_2_trimmed.fastq.gz -S /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/02_hap2/PRJNA863910_30samples_Hap2.sam -p 48"

                sbatch --partition=pibu_el8 --job-name=H2SAMTOOLS --time=0-21:00:00 --mem-per-cpu=64G --ntasks=48 --cpus-per-task=1 --output=Hap2_SAMTOOLS.log --error=Hap2_SAMTOOLS.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/02_hap2/; module load SAMtools; samtools view --threads 48 -b -o PRJNA863910_30samples_Hap2.bam PRJNA863910_30samples_Hap2.sam; samtools sort -o PRJNA863910_30samples_Hap2_sorted.bam -T PRJNA863910_30samples_Hap2_temp --threads 48 PRJNA863910_30samples_Hap2.bam"
