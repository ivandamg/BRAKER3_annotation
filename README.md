# BRAKER3_annotation
Argan_annotation_BRAKER3_container


Script to annotate genomes with braker3


## work in interactive node

1. open interactive node

        srun -p pibu_el8 --mem=50G --cpus-per-task=8 --pty /bin/bash

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
    
              sbatch --partition=pibu_el8 --job-name=hap1BRAKER3 --time=12-24:00:00 --mem-per-cpu=50G --ntasks=36 --cpus-per-task=1 --output=hap1BRAKER3.out --error=hap1BRAKER3.error --mail-type=END,FAIL --wrap "cd /data/users/imateusgonzalez/Z_Soft/; singularity exec --no-home -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl --AUGUSTUS_CONFIG_PATH=${PWD}/config/ --species=Argan --genome=01_Hap1/hap1.fasta.masked --bam=01_Hap1/ALLSamples_Hap1_sorted.bam --prot_seq=01_Hap1/Medicago_3ericales_2sapotaceae_proteins.fasta --workingdir=hap1_braker --GENEMARK_PATH=${ETP}/gmes --threads 36 --busco_lineage embryophyta_odb10 &> hap1.log"







