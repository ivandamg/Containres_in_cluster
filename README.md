# Containers_in_cluster
Use of apptainer in cluster


To use apptainer in cluster, you shoud sent a job or move into an interactive session because there is no access from the front node.

        srun --pty --mem-per-cpu=2G /bin/bash


Follow instructions here: https://docs.pages.bioinformatics.unibe.ch/cluster-docs/documentation/apptainer/ 
and here: https://docs.pages.bioinformatics.unibe.ch/cluster-docs/tutorials/containers/apptainer/ to initialize and use apptainer.

Then in the interactive session

# BRAKER3

mount directories

      apptainer exec --bind "/data" /data/projects/p782_RNA_seq_Argania_spinosa/braker3_latest.sif ls /data






# TE-tools

pull the image


      apptainer shell docker://dfam/tetools:latest


Now inside the container, can work on TE-tools


BuildDatabase -name hap1_DB -engine ncbi hap1.fasta 
RepeatModeler -database hap1_DB -engine ncbi -pa 12


or in cluster

1. Make db

              sbatch --partition=pibu_el8 --job-name=TEtools --time=2-10:00:00 --mem-per-cpu=10G --ntasks=12 --cpus-per-task=1 --output=hap1db.out --error=hap1db.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/01_RepeatMasker/; export APPTAINER_TMPDIR=/data/users/imateusgonzalez; export APPTAINER_CACHEDIR=/data/users/imateusgonzalez; export APPTAINER_BIND=/data ; apptainer run docker://dfam/tetools:latest BuildDatabase -name hap1_DB -engine ncbi hap1.fasta"


2. Repeat modeler

              sbatch --partition=pibu_el8 --job-name=TEtools --time=7-10:00:00 --mem-per-cpu=10G --ntasks=12 --cpus-per-task=1 --output=hap1RepeatMod.out --error=hap1RepeatMod.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/01_RepeatMasker/; export APPTAINER_TMPDIR=/data/users/imateusgonzalez; export APPTAINER_CACHEDIR=/data/users/imateusgonzalez; export APPTAINER_BIND=/data ; apptainer run docker://dfam/tetools:latest RepeatModeler -database hap1_DB -engine ncbi -threads 12 -LTRStruct "

2b. continue Repeat modeler

              sbatch --partition=pibu_el8 --job-name=TEtools --time=7-10:00:00 --mem-per-cpu=10G --ntasks=12 --cpus-per-task=1 --output=hap1RepeatMod.out --error=hap1RepeatMod.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/01_RepeatMasker/; export APPTAINER_TMPDIR=/data/users/imateusgonzalez; export APPTAINER_CACHEDIR=/data/users/imateusgonzalez; export APPTAINER_BIND=/data ; apptainer run docker://dfam/tetools:latest RepeatModeler -database hap1_DB -engine ncbi -threads 12 -recoverDir /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/01_RepeatMasker/RM_133275.MonApr151121242024 -LTRStruct"

if need to dot the same for hap2 the folder is /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/01_RepeatMasker/RM_173348.MonApr151121572024


3. Repeat masker

              sbatch --partition=pibu_el8 --job-name=RMaskHap1 --time=4-10:00:00 --mem-per-cpu=10G --ntasks=12 --cpus-per-task=1 --output=hap1RepeatMask.out --error=hap1RepeatMask.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/01_RepeatMasker/01_Hap1; export APPTAINER_TMPDIR=/data/users/imateusgonzalez; export APPTAINER_CACHEDIR=/data/users/imateusgonzalez; export APPTAINER_BIND=/data ; apptainer run docker://dfam/tetools:latest RepeatMasker hap1.fasta hap1_DB-families.fa -xsmall -pa 12"

                sbatch --partition=pibu_el8 --job-name=RMaskHap2 --time=1-10:00:00 --mem-per-cpu=10G --ntasks=12 --cpus-per-task=1 --output=hap2RepeatMask.out --error=hap2RepeatMask.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/01_RepeatMasker/02_Hap2; export APPTAINER_TMPDIR=/data/users/imateusgonzalez; export APPTAINER_CACHEDIR=/data/users/imateusgonzalez; export APPTAINER_BIND=/data ; apptainer run docker://dfam/tetools:latest RepeatMasker hap2.fasta hap2_DB-families.fa -xsmall -pa 12"


# BRAKER 3
https://bioinformaticsworkbook.org/dataAnalysis/GenomeAnnotation/Intro_to_Braker2.html#gsc.tab=0

mount directories

      apptainer exec --bind "/data" /data/projects/p782_RNA_seq_Argania_spinosa/braker3_latest.sif ls /data

Steps:
1. create repeats library and softmask your genome with RepeatMasker (DONE in previous step)
2. Trim and align your RNA-seq data to your genome
3. Align transcripts, ESTs, or other long transcriptional reads to genome with splice aware aligner
4. Convert Sam to Bam and sort
5. Run Braker3


2. Trim data with fastp and align your RNA-seq data to your genome with HISAT2
https://www.reneshbedre.com/blog/hisat2-sequence-aligner.html



           1. Build index on genome assembly

                      sbatch --partition=pibu_el8 --job-name=H1Hisatindex --time=0-21:00:00 --mem-per-cpu=16G --ntasks=12 --cpus-per-task=1 --output=Hap1_HiSat2index.log --error=Hap1_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/02_HISAT2_mapping/01_Hap1; module load HISAT2; hisat2-build hap1.fasta.masked hap1_hisat_index -p 12"

   2. Map reads to genome

           hisat2 -p 12 -x ${GENOME%.*} -1 ${R1_FQ%%.*}_val_1.fq.gz -2 ${R2_FQ%%.*}_val_2.fq.gz -S ${R1_FQ%.*}.sam


        
5. Run braker

                braker.pl --genome=genome.fa --prot seq=proteins.fa --rnaseq sets ids=SRA ID1,SRA ID2 --rnaseq sets dirs=/path/to/RNASeq/

   



