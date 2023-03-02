perl parseseq_mRNA.pl ~/data/CaiT/sRNA/PSCSR-seq/whitelist.txt M3MAO3-M_R1.fq.gz M3MAO3-M_R2.fq.gz > processed.fq;
perl extract_long.pl processed.fq > processed2.fq;
STAR --runThreadN 12 --genomeDir /ddn/HISEQ2500/genome/Mouse/STAR/index --readFilesIn  processed2.fq --outFilterType BySJout --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate  --outFileNamePrefix mRNA ;
perl stat_reads.pl processed.fq processed2.fq mRNAAligned.sortedByCoord.out.bam > stat_reads.txt &
perl filter_mRNA.pl ~/data/CaiT/database/mmu/protein_exon.gff3 mRNAAligned.sortedByCoord.out.bam > map.anno ;
perl count_mRNA.pl map.anno > map.count;
#perl filter_miRNA.pl ~/data/CaiT/sRNA/PSCSR-seq/ncRNA.gff3 mRNAAligned.sortedByCoord.out.bam.deg >ncRNA.anno ;
#perl count_sRNA.pl ncRNA.anno > ncRNA.count;

perl ~/data/CaiT/sRNA/PSCSRII/map_gene.pl ~/data/CaiT/sRNA/PSCSRII/mouse_lung/exon_sorted.bed  ~/data/CaiT/sRNA/PSCSRII/mouse_lung/intron_sorted.bed ~/data/CaiT/sRNA/PSCSRII/mouse_lung/rRNA.gtf mRNAAligned.sortedByCoord.out.bam.deg > qc.map ;
perl ~/data/CaiT/sRNA/PSCSRII/qc.pl mRNAAligned.sortedByCoord.out.bam.deg qc.map > qc.txt;
