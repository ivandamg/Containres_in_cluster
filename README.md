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

              sbatch --partition=pibu_el8 --job-name=RMaskHap1 --time=1-10:00:00 --mem-per-cpu=10G --ntasks=36 --cpus-per-task=1 --output=hap1RepeatMask.out --error=hap1RepeatMask.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/01_RepeatMasker/01_Hap1; export APPTAINER_TMPDIR=/data/users/imateusgonzalez; export APPTAINER_CACHEDIR=/data/users/imateusgonzalez; export APPTAINER_BIND=/data ; apptainer run docker://dfam/tetools:latest RepeatMasker -lib hap1_DB-families.fa -xsmall -pa 36 -gff -dir MaskerOutput hap1.fasta"

                sbatch --partition=pibu_el8 --job-name=RMaskHap2 --time=1-10:00:00 --mem-per-cpu=10G --ntasks=36 --cpus-per-task=1 --output=hap2RepeatMask.out --error=hap2RepeatMask.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/01_RepeatMasker/02_Hap2; export APPTAINER_TMPDIR=/data/users/imateusgonzalez; export APPTAINER_CACHEDIR=/data/users/imateusgonzalez; export APPTAINER_BIND=/data ; apptainer run docker://dfam/tetools:latest RepeatMasker -lib hap2_DB-families.fa -xsmall -pa 36 -gff -dir MaskerOutput hap2.fasta"


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
hisat2 --phred33 --dta -x Athaliana -1 ath_seed_sample_R1.fastq,ath_root_sample_R1.fastq -2 ath_seed_sample_R2.fastq,ath_root_sample_R2.fastq -S alignment.sam

                sbatch --partition=pibu_el8 --job-name=H1Hisatmap3 --time=3-21:00:00 --mem-per-cpu=16G --ntasks=12 --cpus-per-task=1 --output=Hap1_HiSat2index.log --error=Hap1_HiSat2index.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/02_HISAT2_mapping/01_Hap1; module load HISAT2; hisat2 --phred33 --dta -x hap1_hisat_index -1 /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/7A_1_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/8A_1_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/9A_1_trimmed.fastq.gz -2 /data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/7A_2_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/8A_2_trimmed.fastq.gz,/data/projects/p782_RNA_seq_Argania_spinosa/21_RNAseqV2/02_TrimmedData/9A_2_trimmed.fastq.gz -S threeControlSamples_Hap1.sam -p 12"

                sbatch --partition=pibu_el8 --job-name=H1SAMTOOLS --time=0-21:00:00 --mem-per-cpu=16G --ntasks=12 --cpus-per-task=1 --output=Hap1_SAMTOOLS.log --error=Hap1__SAMTOOLS.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/02_HISAT2_mapping/01_Hap1; module load SAMtools; samtools view --threads 12 -b -o threeControlSamples_Hap1.bam threeControlSamples_Hap1.sam; samtools sort -m 7G -o threeControlSamples_Hap1_sorted.bam -T threeControlSamples_Hap1_temp --threads 12 threeControlSamples_Hap1.bam"

                sbatch --partition=pibu_el8 --job-name=H1SAMTOOLS --time=0-21:00:00 --mem-per-cpu=16G --ntasks=12 --cpus-per-task=1 --output=Hap1_SAMTOOLS.log --error=Hap1__SAMTOOLS.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/02_HISAT2_mapping/01_Hap1; module load SAMtools; "

        2b. convert to bam and sort
                
                module load samtools
                samtools view --threads ${PROC} -b -o ${R1_FQ%.*}.bam ${R1_FQ%.*}.sam
                samtools sort -m 7G -o ${R1_FQ%.*}_sorted.bam -T ${R1_FQ%.*}_temp --threads ${PROC} ${R1_FQ%.*}.bam


        
5. Run braker


                   export APPTAINER_TMPDIR=/data/users/imateusgonzalez; export APPTAINER_CACHEDIR=/data/users/imateusgonzalez; export APPTAINER_BIND=/data ;apptainer exec /data/projects/p782_RNA_seq_Argania_spinosa/braker3_latest.sif --bind "/data" print_braker3_setup.py

                apptainer exec /data/projects/p782_RNA_seq_Argania_spinosa/braker3_latest.sif braker.pl --genome=/data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/01_RepeatMasker/01_Hap1/hap1.fasta.masked --species=Argan --prot_seq=/data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/03_BRAKER/Vitpa_pep_LATEST.fa,/data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/03_BRAKER/miracle_fruit.pep.fa --bam=/data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/02_HISAT2_mapping/01_Hap1/threeControlSamples_Hap1_sorted.bam

                sbatch --partition=pibu_el8 --job-name=BRAKER3 --time=7-10:00:00 --mem-per-cpu=10G --ntasks=12 --cpus-per-task=1 --output=BRAKER3.out --error=BRAKER3.err --mail-type=END,FAIL --wrap "cd /data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/03_BRAKER; export APPTAINER_TMPDIR=/data/users/imateusgonzalez; export APPTAINER_CACHEDIR=/data/users/imateusgonzalez; export APPTAINER_BIND=/data ; apptainer exec /data/projects/p782_RNA_seq_Argania_spinosa/braker3_latest.sif braker.pl --AUGUSTUS_ab_initio --gff3 --threads=12 --workingdir=/data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/03_BRAKER/01_Hap1 --genome=/data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/01_RepeatMasker/01_Hap1/hap1.fasta.masked --prot_seq=/data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/03_BRAKER/Vitpa_pep_LATEST.fa,/data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/03_BRAKER/miracle_fruit.pep.fa --bam=/data/projects/p782_RNA_seq_Argania_spinosa/21_GenomeAnnotation/02_HISAT2_mapping/01_Hap1/threeControlSamples_Hap1_sorted.bam"
