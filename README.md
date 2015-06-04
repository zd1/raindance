# Souce code repository for manuscript
**Quantifying Spontaneous Germline Mutations Directly in Human Testis**

# User guide

## Preprocess reads in fastq format

```python
python parse_primer.py  --fastq  --sample sampleName
```

## Read mapping

```
bwa mem -R "@RG\tID:$ID\tSM:$SM\tLB:$LB" $ref $read1 $read2 > $sam
```
## Pileup

```
python parse_primer.py  --pileup  --sample sampleName
```

## Integrate quantifications of samples
```
python parse_primer.py --sum
```

## Aggregrate quantifications of each amplicon across samples, generating files/plots per amplicon
```
python parse_primer.py --agg --amp ampliconName
```

## Call elevated alleles

```
python parse_primer.py --call --amp ampliconName
```

# Meta Data

 |    File Name  | Description |
 | --------------| ------------|
 | meta/samplelist | A list of IDs of all samples |
 | meta/fastq.index | A list of all fastq files |
 | meta/amplicons.csv | Information about amplicons |
 | meta/sample_seq.csv | sequencing production information | 

# Contact
  someone at imm.ox.ac.uk