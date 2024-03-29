## GuppyGermlineDNMs
Scripts for Lin, Y., Darolti, I., van der Bijl, W., Morris, J., Mank, J. E. (2023) Extensive variation in germline *de novo* mutations in *Poecilia reticulata*. Genome Research https://doi.org/10.1101/gr.277936.123

More detailed tutorial to be finished

1. Quality control for raw sequencing reads using [FastQC](https://github.com/s-andrews/FastQC) and [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
   ```
   python 01.QC_multi.py [-t THREADS] \
                        [-fq FASTQ_PATH] \
                        [-o OUTPUT] \
                        [-trim] \
                        [-qc]
   ```

2. Reference genome reconstruction based on [Ensembl reference genome](http://useast.ensembl.org/Poecilia_reticulata/Info/Index) using [LongRanger](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger), [LAST](https://gitlab.com/mcfrith/last/), [ARCS](https://github.com/bcgsc/arcs) and [RagTag](https://github.com/malonge/RagTag)
   ```
   sh 02.refGenomeReconstruct.sh [seq] [ref]
   ```
3. Read alignment using [BWA MEM](https://github.com/lh3/bwa).
   ```
   sh 03.alignment.sh [ref] [prefix]
   ```

4. Genotyping and SNP filtering

   #(1) Genotyping using [BCFtools](https://samtools.github.io/bcftools/howtos/index.html) and [GATK](https://gatk.broadinstitute.org/hc/en-us)

   #(2) SNP filtering using [SNPfiltering.py](./SNPfiltering.py) and [BAMfilter.py](./BAMfilter.py)
   ```
   sh 04.DNM.sh [father] [mother] [prefix] [ref]
   ```
  
5. Genotype phasing and Kinship analysis

   #(1) Genotype phasing using [WhatsHap](https://github.com/whatshap/whatshap)

   #(2) Kinship analysis using [KING](https://www.kingrelatedness.com/)
   ```
   sh 05.phasingKinship.sh [fam] [ref]
   ```

7. Repeat identification using [RepeatModeler2](https://www.repeatmasker.org/RepeatModeler/) and [RepeatMasker](https://www.repeatmasker.org/)
   ```
   sh 06.repeatAnnotation.sh [prefix] [output] [ref] [genomeSize]
   ```

8. Callable genome size calculation 
   ```
   sh 07.CallableGenomeSize.sh [ref] [prefix]
   ```

9. Simulation using [bamSurgeon](https://github.com/adamewing/bamsurgeon).
   ```
   sh 08.simBamSurgeon.sh [ref] [prefix]
   ```


