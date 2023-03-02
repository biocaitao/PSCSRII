# PSCSR-seq V2 analysis pipeline



# the mRNA reads analysis pipeline

Pre-required enviroment 

a. STAR (version >=2.6.1 ) installed 

b. STAR reference database 

#### 1. parse raw reads

e.g. 
```bash
> cd mRNA:
> perl parseseq_mRNA.pl whitelist.txt sample1_R1.fq.gz sample1_R2.fq.gz > processed.fq;
```

#### 2. extract mRNAs
e.g.
```bash
> perl extract_long.pl processed.fq > processed2.fq;
```


#### 3. map reads
e.g.
```bash
> STAR --runThreadN 12 --genomeDir mouseSTARindex --readFilesIn  processed2.fq --outFilterType BySJout --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate  --outFileNamePrefix mRNA
```

#### 4. annotate mRNA reads and count RNA expression
e.g.
```bash
> perl filter_mRNA.pl ~/data/CaiT/database/mmu/protein_exon.gff3 mRNAAligned.sortedByCoord.out.bam > map.anno 
> perl count_mRNA.pl map.anno > map.count;
```

# Please check https://github.com/biocaitao/PSCSR-seq for the small RNA reads analysis
