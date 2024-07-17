
# BRAKER3_annotation
Argan_annotation_BRAKER3_container


Script to annotate genomes with braker3


## work in interactive node

1. open interactive node

        srun -p pibu_el8 --mem=4G --cpus-per-task=1 --pty /bin/bash

2. change aptainer tmp dir

        export APPTAINER_TMPDIR=/data/users/imateusgonzalez/Z_Soft/
   
3. build braker3

                singularity build braker3.sif docker://teambraker/braker3:latest
   
4.   export the container image

                export BRAKER_SIF=/data/users/imateusgonzalez/Z_Soft/braker3.sif
     
6.   Run to see the help menu
   
                singularity exec ${BRAKER_SIF} braker.pl

7. copy the Augustus config directory locally to make it writable

               singularity exec --no-home -B $PWD:$PWD braker3.sif cp -R /opt/Augustus/config $PWD
   
8.   remove output directory (eg. test1) if it already exists

              wd=hap1_braker
              if [ -d $wd ]; then rm -r $wd ;fi


9.  Copy all files to annotation in folder: softmasked assembly, aligned sorted bam, and protein from other species.
  
10.  run braker3 in cluster in /data/users/imateusgonzalez/Z_Soft/
    
              sbatch --partition=pibu_el8 --job-name=hap1BRAKER3 --time=12-24:00:00 --mem-per-cpu=50G --ntasks=48 --cpus-per-task=1 --output=hap1BRAKER3.out --error=hap1BRAKER3.error --mail-type=END,FAIL --wrap "export APPTAINER_TMPDIR=/data/users/imateusgonzalez/Z_Soft/; cd /data/users/imateusgonzalez/Z_Soft/; singularity exec --no-home -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl --AUGUSTUS_CONFIG_PATH=${PWD}/config/ --species=Argania --softmasking --gff3 --useexisting --genome=01_Hap1/hap1.fasta.masked --bam=01_Hap1/ALLSamples_Hap1_sorted.bam --prot_seq=01_Hap1/Medicago_3ericales_2sapotaceae_proteins.fasta --workingdir=hap1_braker --GENEMARK_PATH=${ETP}/gmes --threads 48 --busco_lineage embryophyta_odb10 &> hap1.log"

11.  run in interactive with export steps 2 and 4 braker3 only with vitellaria proteins and all RNASeq in cluster in /data/users/imateusgonzalez/Z_Soft/
    
              sbatch --partition=pibu_el8 --job-name=hap1BRAKER3 --time=1-24:00:00 --mem-per-cpu=64G --ntasks=48 --cpus-per-task=1 --output=hap1BRAKER3.out --error=hap1BRAKER3.error --mail-type=END,FAIL --wrap "cd /data/users/imateusgonzalez/Z_Soft/; singularity exec --no-home -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl --AUGUSTUS_CONFIG_PATH=${PWD}/config/ --species=Argania --gff3 --genome=01_Hap1/hap1.fasta.masked --bam=01_Hap1/ALLSamples_Hap1_sorted.bam --prot_seq=01_Hap1/miracle_pep.fasta --workingdir=hap1_braker_v2 --GENEMARK_PATH=${ETP}/gmes --threads 48 --busco_lineage embryophyta_odb10 &> hap1.log"

with protein evidence including EMbryophytaReviewedSwissprot_Medicago_2Sapotaceae_3ericales

                sbatch --partition=pibu_el8 --job-name=hap1BRAKER3_v3 --time=1-24:00:00 --mem-per-cpu=64G --ntasks=48 --cpus-per-task=1 --output=hap1BRAKER3_v3.out --error=hap1BRAKER3_v3.error --mail-type=END,FAIL --wrap "cd /data/users/imateusgonzalez/Z_Soft/; singularity exec --no-home -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl --AUGUSTUS_CONFIG_PATH=${PWD}/config/ --species=Argania --gff3 --genome=01_Hap1/hap1.fasta.masked --bam=01_Hap1/ALLSamples_Hap1_sorted.bam --prot_seq=01_Hap1/2_EmbryophytaSwissProt_Medicago_3ericales_2sapotaceae_proteins.fasta --useexisting --workingdir=hap1_braker_v3 --GENEMARK_PATH=${ETP}/gmes --threads 48 --busco_lineage embryophyta_odb10 &> hap1.log"

12. with 221 Embryophyta proteomes as protein evidence and 12 + 30 RNAseqs 

                sbatch --partition=pibu_el8 --job-name=hap1BRAKER3_v8 --time=1-24:00:00 --mem-per-cpu=64G --ntasks=48 --cpus-per-task=1 --output=hap1BRAKER3_v8.out --error=hap1BRAKER3_v8.error --mail-type=END,FAIL --wrap "cd /data/users/imateusgonzalez/Z_Soft/; singularity exec --no-home -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl --AUGUSTUS_CONFIG_PATH=${PWD}/config/ --species=Argania --gff3 --genome=01_Hap1/hap1.fasta.masked --bam=/data/users/imateusgonzalez/Z_Soft/01_Hap1/ALLSamples_Hap1_sorted.bam,/data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/02_HISAT2_mapping/01_Hap1/PRJNA863910_30samples_Hap1_sorted.bam --prot_seq=/data/users/imateusgonzalez/Z_Soft/01_Hap1/01_Proteomes/221_Uniprot_Embryophyta_Protome.fasta --useexisting --workingdir=hap1_braker_v8 --GENEMARK_PATH=${ETP}/gmes --threads 48 --busco_lineage embryophyta_odb10 &> hap1_v8.log"



13. VF: 221 Embryophyta proteom_ 12 + 30 RNASEQ

                     sbatch --partition=pibu_el8 --job-name=hap2BRAKER3_v9 --time=1-24:00:00 --mem-per-cpu=64G --ntasks=48 --cpus-per-task=1 --output=hap2BRAKER3_v9.out --error=hap2BRAKER3_v9.error --mail-type=END,FAIL --wrap "export APPTAINER_TMPDIR=/data/users/imateusgonzalez/Z_Soft/; export BRAKER_SIF=/data/users/imateusgonzalez/Z_Soft/braker3.sif; cd /data/users/imateusgonzalez/Z_Soft/; singularity exec --no-home -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl --AUGUSTUS_CONFIG_PATH=${PWD}/config/ --species=Argania --gff3 --genome=/data/projects/p782_RNA_seq_Argania_spinosa/30_FinalAssemblyPaper/04_RNAseqMapping/02_hap2/Aspinosa_hap2.fa.masked --bam=/data/projects/p782_RNA_seq_Argania_spinosa/30_FinalAssemblyPaper/04_RNAseqMapping/02_hap2/Ufribourg_12samples_Hap2_sorted.bam,/data/projects/p782_RNA_seq_Argania_spinosa/30_FinalAssemblyPaper/04_RNAseqMapping/02_hap2/PRJNA863910_30samples_Hap2_sorted.bam --prot_seq=/data/projects/p782_RNA_seq_Argania_spinosa/30_FinalAssemblyPaper/05_BRAKER3/Ref_proteomes/221_Embryophyta_proteomes.fasta --useexisting --workingdir=hap2_braker_v9 --GENEMARK_PATH=${ETP}/gmes --threads 48 --busco_lineage eudicots_odb10 &> hap2_v9.log"

14. VF: 221 Embryophyta proteom_ 12 RNASEQ

                sbatch --partition=pibu_el8 --job-name=hap1BRAKER3_v9 --time=1-24:00:00 --mem-per-cpu=64G --ntasks=48 --cpus-per-task=1 --output=hap1BRAKER3_v9.out --error=hap1BRAKER3_v9.error --mail-type=END,FAIL --wrap "export APPTAINER_TMPDIR=/data/users/imateusgonzalez/Z_Soft/; export BRAKER_SIF=/data/users/imateusgonzalez/Z_Soft/braker3.sif; cd /data/users/imateusgonzalez/Z_Soft/; singularity exec --no-home -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl --AUGUSTUS_CONFIG_PATH=${PWD}/config/ --species=Argania --gff3 --genome=/data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/05_BRAKER3/01_hap1/out_JBAT_hap1.FINAL.fa.mod.MAKER.softmasked --bam=/data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/04_RNAseqMapping/01_hap1/Ufribourg_12samples_Hap1_sorted.bam --prot_seq=/data/projects/p782_RNA_seq_Argania_spinosa/50_FinalArgan/05_BRAKER3/Ref_proteomes/221_Embryophyta_proteomes.fasta --useexisting --workingdir=Final_hap1_braker_v9 --GENEMARK_PATH=${ETP}/gmes --threads 48 --busco_lineage eudicots_odb10 &> hap1_v9.log"

14. Evaluate annotation BUSCO

                        sbatch --partition=pibu_el8 --job-name=hap2Busco --time=0-10:00:00 --mem-per-cpu=50G --ntasks=12 --cpus-per-task=1 --output=BuscoHap2.out --error=BuscoHap2.error --mail-type=END,FAIL --wrap "module load BUSCO; cd /data/users/imateusgonzalez/Z_Soft/hap2_braker_v9; busco -o Busco/ -i braker.aa -l eudicots_odb10 --cpu 12 -m protein -f"


