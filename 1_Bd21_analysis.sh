
#1_Bd21_analysis
#Sept 2016
#make total cytosine site count 100bp tiles for the entire Bd21 genome
cd rawdata/genomes/Bd21/
samtools faidx Bd21Control_SNPincorp_sgr1_genome.fa
cut -f1,2 Bd21Control_SNPincorp_sgr1_genome.fa.fai > sizes.Bd21.genome
cd ../../../

mkdir 1_Bd21_analysis_output
cd 1_Bd21_analysis_output


coverage2cytosine --genome_folder ../rawdata/genomes/Bd21/ -CX ../0_initial_mapping_output/Bd21_wgbspipeline_2016-09-20-14-09-31/5_output_files/Bd21_CpG.bed.bismark.cov -o allcytosines.txt

grep -P "\tCG\t" allcytosines.txt | awk '{print $1 "\t" $2 "\t" $2 "\t" $6 "\t" $7 "\t" $3}' | sort -k1,1 -k2,2n > Bd21_CGsites.bed
grep -P "\tCHG\t" allcytosines.txt | awk '{print $1 "\t" $2 "\t" $2 "\t" $6 "\t" $7 "\t" $3}' | sort -k1,1 -k2,2n> Bd21_CHGsites.bed
grep -P "\tCHH\t" allcytosines.txt | awk '{print $1 "\t" $2 "\t" $2 "\t" $6 "\t" $7 "\t" $3}' | sort -k1,1 -k2,2n> Bd21_CHHsites.bed

bedtools makewindows -g ../rawdata/genomes/Bd21/sizes.Bd21.genome -w 100 | awk '{print $1 "\t" $2 "\t" $3-1}' | sort -k1,1 -k2,2n > Bd21_allwindows.bed

bedtools intersect -wa -wb -a Bd21_allwindows.bed -b Bd21_CGsites.bed | bedtools groupby -g 1,2,3, -c 4 -o count > CGsites_100bp.bed
bedtools intersect -wa -wb -a Bd21_allwindows.bed -b Bd21_CHGsites.bed | bedtools groupby -g 1,2,3, -c 4 -o count > CHGsites_100bp.bed
bedtools intersect -wa -wb -a Bd21_allwindows.bed -b Bd21_CHHsites.bed | bedtools groupby -g 1,2,3, -c 4 -o count > CHHsites_100bp.bed

#make an annotation file where we take every genomic window to its nearest gene and transposon

cd ../rawdata/annotations/


gzip -d Bdistachyon_192_gene.gff3.gz
gzip -d MIPS_Bd_Transposons_v2.2_16-07-2009.gff3.gz

###############################
R

getAttributeField <- function (x, field, attrsep = ";") {
     s = strsplit(x, split = attrsep, fixed = TRUE)
     sapply(s, function(atts) {
         a = strsplit(atts, split = "=", fixed = TRUE)
         m = match(field, sapply(a, "[", 1))
         if (!is.na(m)) {
             rv = a[[m]][2]
         }
         else {
             rv = as.character(NA)
         }
         return(rv)
     })
}

gffRead <- function(gffFile, nrows = -1) {
     cat("Reading ", gffFile, ": ", sep="")
     gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
     header=FALSE, comment.char="#", nrows = nrows,
     colClasses=c("character", "character", "character", "integer",  
"integer",
     "character", "character", "character", "character"))
     colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
     cat("found", nrow(gff), "rows with classes:",
         paste(sapply(gff, class), collapse=", "), "\n")
     stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
     return(gff)
}

gene_gff <- gffRead('Bdistachyon_192_gene.gff3')
gene_gff.sub <- subset(gene_gff, gene_gff$feature=='gene')
gene_gff.sub$geneID <- getAttributeField(gene_gff.sub$attributes, 'Name')

#gff is 1-based, if we want a bed output, we want 0-based

gene_gff.sub$start <- gene_gff.sub$start - 1

#typo in gff file for gene Bradi1g02575 as it is certainly not 100Mb in size...
gene_gff.sub[which((gene_gff.sub$end - gene_gff.sub$start) > 100000),5] <- 1734141

gene_out <- gene_gff.sub[,c(1,4,5,10,6,7,3)]
write.table(gene_out,'Bdistachyon_192_gene.bed',sep='\t',row.names=F,col.names=F,quote=F)
	
	
	
te_gff <- gffRead('MIPS_Bd_Transposons_v2.2_16-07-2009.gff3')

#remove 'fragments' and keep only full elements
te_gff.sub <- subset(te_gff, te_gff$feature=='transposable_element')
te_gff.sub$class <- getAttributeField(te_gff.sub$attributes, 'class')
te_gff.sub$Name <- getAttributeField(te_gff.sub$attributes, 'Name')

#rename chromosomes to match gene data
te_gff.sub$chrom <- ifelse(te_gff.sub$seqname=='chr01_pseudomolecule','Bd1',
	ifelse(te_gff.sub$seqname=='chr02_pseudomolecule','Bd2',
		ifelse(te_gff.sub$seqname=='chr03_pseudomolecule','Bd3',
			ifelse(te_gff.sub$seqname=='chr04_pseudomolecule','Bd4',
				ifelse(te_gff.sub$seqname=='chr05_pseudomolecule','Bd5',te_gff.sub$seqname)))))
	
	
te_out <- te_gff.sub[,c(12,4,5,11,10,7,3)]
write.table(te_out,'Bdistachyon_v2.2_MIPS_transposons.bed',sep='\t',row.names=F,col.names=F,quote=F)

quit()
n

sort -k1,1 -k2,2n Bdistachyon_192_gene.bed > temp
mv temp Bdistachyon_192_gene.bed

sort -k1,1 -k2,2n Bdistachyon_v2.2_MIPS_transposons.bed > temp
mv temp Bdistachyon_v2.2_MIPS_transposons.bed

cd ../../1_Bd21_analysis_output

bedtools closest -D "ref" -a Bd21_allwindows.bed -b ../rawdata/annotations/Bdistachyon_192_gene.bed > Bd21_allwindows.nearestgene.bed
sort -k1,1 -k2,2n Bd21_allwindows.nearestgene.bed > test
mv test Bd21_allwindows.nearestgene.bed 

bedtools closest -D "ref" -a Bd21_allwindows.nearestgene.bed -b ../rawdata/annotations/Bdistachyon_v2.2_MIPS_transposons.bed > Bd21_allwindows.nearestgene.nearestte.bed

bedtools makewindows -g ../rawdata/genomes/Bd21/sizes.Bd21.genome -w 500000 | sort -k1,1 -k2,2n > Bd21_500k.windows.bed

bedtools intersect -wa -wb -a Bd21_500k.windows.bed -b ../rawdata/annotations/Bdistachyon_192_gene.bed | bedtools groupby -g 1,2,3, -c 4 -o count > gene.density.bd21.txt
bedtools intersect -wa -wb -a Bd21_500k.windows.bed -b ../rawdata/annotations/Bdistachyon_v2.2_MIPS_transposons.bed | bedtools groupby -g 1,2,3, -c 4 -o count > repeat.density.bd21.txt



############################
R
############################

all.cg=read.delim('CGsites_100bp.bed',head=F)
all.chg=read.delim('CHGsites_100bp.bed',head=F)
all.chh=read.delim('CHHsites_100bp.bed',head=F)

t1=merge(all.cg,all.chg,by=c('V1','V2','V3'),all=T)
all.windows=merge(t1,all.chh,by=c('V1','V2','V3'),all=T)
colnames(all.windows)=c('V1','V2','V3','CGsites','CHGsites','CHHsites')

cg.possible=table(is.na(all.windows$CGsites))[1]
chg.possible=table(is.na(all.windows$CHGsites))[1]
chh.possible=table(is.na(all.windows$CHHsites))[1]

#read in Bd21 mapped tile data
cg=read.delim('../0_initial_mapping_output/100bp_tiles/Bd21_CpG_100bp.wig',head=F,skip=1)
chg=read.delim('../0_initial_mapping_output/100bp_tiles/Bd21_CHG_100bp.wig',head=F,skip=1)
chh=read.delim('../0_initial_mapping_output/100bp_tiles/Bd21_CHH_100bp.wig',head=F,skip=1)

#get rid of scaffolds
cg=cg[grep('^Bd',cg$V1),]
chg=chg[grep('^Bd',chg$V1),]
chh=chh[grep('^Bd',chh$V1),]

cg=cg[with(cg,order(cg[,1],cg[,2])),]
chg=chg[with(chg,order(chg[,1],chg[,2])),]
chh=chh[with(chh,order(chh[,1],chh[,2])),]

all=merge(all.windows,cg,by=c('V1','V2','V3'),all=T)
all=merge(all,chg,by=c('V1','V2','V3'),all=T)
colnames(all)=c('V1','V2','V3','CGsites','CHGsites','CHHsites','cg_prop','cg_met','cg_unmet','cg_total','cg_site','chg_prop','chg_met','chg_unmet','chg_total','chg_site')
all=merge(all,chh,by=c('V1','V2','V3'),all=T)
colnames(all)=c('V1','V2','V3','CGsites','CHGsites','CHHsites','cg_prop','cg_met','cg_unmet','cg_total','cg_site','chg_prop','chg_met','chg_unmet','chg_total','chg_site','chh_prop','chh_met','chh_unmet','chh_total','chh_site')

cg.coverage=table(is.na(all$cg_total)==T)[1]
chg.coverage=table(is.na(all$chg_total)==T)[1]
chh.coverage=table(is.na(all$chh_total)==T)[1]

pdf('Supplemental_Fig_S1.pdf',width=10,height=5)
	par(mfrow=c(1,3))
	plot(density(all$cg_prop,na.rm=T),main='CG methylation',xlab='Methylation %',col='#990000',lwd=3)
	plot(density(all$chg_prop,na.rm=T),main='CHG methylation',xlab='Methylation %',col='#003399',lwd=3)
	plot(density(all$chh_prop,na.rm=T),main='CHH methylation',xlab='Methylation %',col='#336600',lwd=3,xlim=c(0,20))
dev.off()

cg.class=ifelse(all$cg_prop>50,1,0)
chg.class=ifelse(all$chg_prop>30,1,0)
chh.class=ifelse(all$chh_prop>10,1,0)

all=cbind(all,cg.class,chg.class,chh.class)

window_class=
ifelse(all$cg.class==1 & all$chg.class==1 & all$chh.class==1,'all_met',
ifelse(all$cg.class==1 & all$chg.class==0 & all$chh.class==0,'cg_only',
ifelse(all$cg.class==0 & all$chg.class==1 & all$chh.class==0,'chg_only',
ifelse(all$cg.class==0 & all$chg.class==0 & all$chh.class==1,'chh_only',
ifelse(all$cg.class==1 & all$chg.class==1 & all$chh.class==0,'cg_chg',
ifelse(all$cg.class==1 & all$chg.class==0 & all$chh.class==1,'cg_chh',
ifelse(all$cg.class==0 & all$chg.class==1 & all$chh.class==1,'chg_chh',
ifelse(all$cg.class==0 & all$chg.class==0 & all$chh.class==0,'no_met',
'NA'))))))))

all=cbind(all,window_class)

pdf('Supplemental_Fig_S2.pdf',width=10,height=5)
	barplot(table(all$window_class,exclude=NULL),las=2)
dev.off()


anno=read.delim('Bd21_allwindows.nearestgene.nearestte.bed',head=F)

anno.class=
  ifelse(anno$V11==0 & anno$V19==0, 'TE-gene boundary',
  ifelse(anno$V11==0,'gene',
  ifelse(anno$V19==0,as.character(anno$V16),
  'intergenic')))

anno2=cbind(anno[,1:3],anno.class)


#fix the NA name for transposable element fragments
anno2$anno.class=ifelse(is.na(anno2$anno.class)==T,'TE_frag',as.character(anno2$anno.class))
arf=all
all=merge(all,anno2,by=c('V1','V2','V3'))
write.table(all,'Bd21_total_tile_annotation.txt',sep='\t',row.names=F,quote=F)


pdf('Supplemental_Fig_S3.pdf',width=10,height=5)
	barplot(table(all$anno.class),las=2,main='Number of tiles across annotation classes')
dev.off()

prop.met.class=prop.table(table(all$anno.class,all$window_class),margin=2)

prop.anno.class=prop.table(table(all$window_class,all$anno.class),margin=2)
prop.anno.class=prop.anno.class[,c(7,8,14,12,1,2,3,4,5,6,9,10,11)]
pdf('fig1c.pdf',width=10,height=5)
	par(mar=c(5,4,4,7),xpd=T)
	barplot(prop.anno.class,col=c('black','#8B008B','#CD6600','#990000','#008B8B','#003399','#336600','grey'),las=2,horiz=T)
	legend("topright",inset=c(-.16,0),fill=c('black','#8B008B','#CD6600','#990000','#008B8B','#003399','#336600','grey'),legend=rownames(prop.anno.class))
dev.off()

plot.table=matrix(NA,nrow=4,ncol=3)
colnames(plot.table)=c('CG','CHG','CHH')
rownames(plot.table)=c('Total tiles with sites','Total tiles with coverage >0','Tiles unmethylated','Tiles methylated')
plot.table[1,]=c(cg.possible,chg.possible,chh.possible)
plot.table[2,]=c(cg.coverage,chg.coverage,chh.coverage)
plot.table[3:4,1]=table(cg.class)
plot.table[3:4,2]=table(chg.class)
plot.table[3:4,3]=table(chh.class)

pdf('fig1a.pdf',width=10,height=5)
	barplot(plot.table[1:2,],beside=T,horiz=T,las=1)
	barplot(plot.table[3:4,],beside=F,horiz=T,las=1,col=c('#336600','grey','#003399','grey','#990000','grey'))
dev.off()


library(fields)
library(ggplot2)

pdf('fig1b.pdf',width=10,height=1)

cg.chr1=subset(cg,cg$V1=='Bd1')
x=stats.bin(cg.chr1$V2,cg.chr1$V4,300)
arf=cbind.data.frame(x$centers,x$stats['mean',])
colnames(arf)=c('pos','met')
ggplot(arf,aes(pos,met)) + geom_bar(stat='identity',fill='#990000',col='#990000') + theme_classic()

chg.chr1=subset(chg,chg$V1=='Bd1')
x=stats.bin(chg.chr1$V2,chg.chr1$V4,300)
arf=cbind.data.frame(x$centers,x$stats['mean',])
colnames(arf)=c('pos','met')
ggplot(arf,aes(pos,met)) + geom_bar(stat='identity',fill='#003399',col='#003399') + theme_classic() + ylim(c(0,100))

chh.chr1=subset(chh,chh$V1=='Bd1')
x=stats.bin(chh.chr1$V2,chh.chr1$V4,300)
arf=cbind.data.frame(x$centers,x$stats['mean',])
colnames(arf)=c('pos','met')
ggplot(arf,aes(pos,met)) + geom_bar(stat='identity',fill='#336600',col='#336600') + theme_classic() + ylim(c(0,5))

gene=read.delim('gene.density.bd21.txt',head=F)
gene.chr1=subset(gene,gene$V1=='Bd1')
gene.chr1$class=rep('a',nrow(gene.chr1))
ggplot(gene.chr1,aes(V2,class)) + geom_tile(aes(fill=V4)) + scale_fill_continuous(low='black',high='yellow') + theme_classic()

te=read.delim('repeat.density.bd21.txt',head=F)
te.chr1=subset(te,te$V1=='Bd1')
te.chr1$class=rep('a',nrow(te.chr1))
ggplot(te.chr1,aes(V2,class)) + geom_tile(aes(fill=V4)) + scale_fill_continuous(low='black',high='red') + theme_classic()

dev.off()


quit()
n
####################################
cd 0_initial_mapping_output/100bp_tiles
../../bed_to_rel_dist.sh -wig ../../rawdata/annotations/Bdistachyon_192_gene.bed Bd21 gene
mv Bd21_gene* ../../1_Bd21_analysis_output/

../../bed_to_rel_dist.sh -wig ../../rawdata/annotations/Bdistachyon_v2.2_MIPS_transposons.bed Bd21 TE
mv Bd21_TE* ../../1_Bd21_analysis_output/

cd ../../1_Bd21_analysis_output/

mv Bd21_gene_methylation.pdf fig1d_gene.pdf
mv Bd21_TE_methylation.pdf fig1d_TE.pdf

#