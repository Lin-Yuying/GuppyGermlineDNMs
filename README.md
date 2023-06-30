## GuppyGermlineDNMs
Scripts for Lin, Y., Darolti, I., van der Bijl, W., Morris, J., Mank, J. E. (2023) Extensive variation in germline *de novo* mutations in *Poecilia reticulata*. *Genome Research* https://doi.org/10.1101/2023.03.22.533860

detailed tutorial to be finished

1. Quality control for raw sequencing reads using [fastQC] and [Trimmomatics]
```
python3 01.QC_multi.py [-h] [-t THREADS] [-fq FASTQ_PATH] [-o OUTPUT] [-trim] [-qc]
```

2. Reference genome reconstruction
```
sh 02.refGenomeReconstruct.sh [seq] [ref]
```
3. Read alignment using BWA mem.
```
sh 03.alignment.sh [ref] [prefix]
```

4. Genotyping and SNP filtering
```
sh 04.DNM.sh [father] [mother] [prefix] [ref]
```

5. Genotype phasing and Kinship analysis
```
sh 05.phasingKinship.sh [fam] [ref]
```

6. Repeat identification
```
sh 05.phasingKinship.sh 
```

7. Callable Genome Size Calculation 
```
code
```

8. Simulation using bamSurgeon.
```
code
```


