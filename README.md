## GuppyGermlineDNMs
Scripts for Lin, Y., Darolti, I., van der Bijl, W., Morris, J., Mank, J. E. (2023) Extensive variation in germline *de novo* mutations in *Poecilia reticulata*. *Genome Research* https://doi.org/10.1101/2023.03.22.533860

detailed tutorial to be finished

1. Quality control for raw sequencing reads using [FastQC](https://github.com/s-andrews/FastQC) and [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
   ```
   python 01.QC_multi.py [-t THREADS] \
                        [-fq FASTQ_PATH] \
                        [-o OUTPUT] \
                        [-trim] \
                        [-qc]
   ```

2. Reference genome reconstruction based on [female reference genome](http://uswest.ensembl.org/Poecilia_reticulata/Info/Index) using LongRanger, ARCS
   ```
   sh 02.refGenomeReconstruct.sh [seq] [ref]
   ```
3. Read alignment using [BWA MEM](https://github.com/lh3/bwa).
   ```
   sh 03.alignment.sh [ref] [prefix]
   ```

4. Genotyping and SNP filtering

   #(1) Genotyping using [BCFtools](https://samtools.github.io/bcftools/howtos/index.html) and [GATK](https://gatk.broadinstitute.org/hc/en-us).

   #(2) SNP filtering using [SNPfiltering.sh](./SNPfiltering.py) and [BAMfilter.py](./BAMfilter.py).
  ```
  sh 04.DNM.sh [father] [mother] [prefix] [ref]
  ```
  
  5. Genotype phasing and Kinship analysis
  ```
  sh 05.phasingKinship.sh [fam] [ref]
  ```

6. Repeat identification
  ```
  sh 06.repeatAnnotation.sh [prefix] [output] [ref] [genomeSize]
  ```

7. Callable genome size calculation 
  ```
  sh 07.CallableGenomeSize.sh [ref] [prefix]
  ```

8. Simulation using bamSurgeon.
  ```
  sh 08.simBamSurgeon.sh [ref] [prefix]
  ```


