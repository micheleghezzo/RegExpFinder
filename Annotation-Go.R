library(GenomicFeatures)
library(org.Hs.eg.db)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ReactomePA)
library(ggplot2)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peak <- readPeakFile("position.bed")

peakAnno <- annotatePeak(
  peak, 
  tssRegion=c(-3000, 3000), 
  assignGenomicAnnotation = TRUE, 
  addFlankGeneInfo = FALSE,                                                                            
  TxDb=txdb, 
  annoDb="org.Hs.eg.db"
  )

annotation <- peakAnno@anno@elementMetadata@listData[["annotation"]]
length <- peak@ranges@width
gene <- peakAnno@anno@elementMetadata@listData[["SYMBOL"]]
m = data.frame(gene, length, annotation)
write.csv(m, "output.csv", row.names = FALSE)

print(paste("Length in output: ", toString(mean(length))," ± ", toString(sd(length))))

Promoter <- gene[which(peakAnno@detailGenomicAnnotation[["Promoter"]] == TRUE)]
length <- peak@ranges@width[which(peakAnno@detailGenomicAnnotation[["Promoter"]] == TRUE)]
m = data.frame(Promoter, length)
write.csv(m, "Promoter.csv", row.names = FALSE)

print(paste("Length in Promoter: ", toString(mean(length))," ± ", toString(sd(length))))

Intergenic <- gene[which(peakAnno@detailGenomicAnnotation[["Intergenic"]] == TRUE)]
length <- peak@ranges@width[which(peakAnno@detailGenomicAnnotation[["Intergenic"]] == TRUE)]
m = data.frame(Intergenic, length)
write.csv(m, "Intergenic.csv", row.names = FALSE)

print(paste("Length in Intergenic: ", toString(mean(length))," ± ", toString(sd(length))))

fiveUTR <- gene[which(peakAnno@detailGenomicAnnotation[["fiveUTR"]] == TRUE)]
length <- peak@ranges@width[which(peakAnno@detailGenomicAnnotation[["fiveUTR"]] == TRUE)]
m = data.frame(fiveUTR, length)
write.csv(m, "fiveUTR.csv", row.names = FALSE)

print(paste("Length in fiveUTR: ", toString(mean(length))," ± ", toString(sd(length))))

threeUTR <- gene[which(peakAnno@detailGenomicAnnotation[["threeUTR"]] == TRUE)]
length <- peak@ranges@width[which(peakAnno@detailGenomicAnnotation[["threeUTR"]] == TRUE)]
m = data.frame(threeUTR, length)
write.csv(m, "threeUTR.csv", row.names = FALSE)

print(paste("Length in threeUTR: ", toString(mean(length))," ± ", toString(sd(length))))

Exon <- gene[which(peakAnno@detailGenomicAnnotation[["Exon"]] == TRUE)]
length <- peak@ranges@width[which(peakAnno@detailGenomicAnnotation[["Exon"]] == TRUE)]
m = data.frame(Exon, length)
write.csv(m, "Exon.csv", row.names = FALSE)

print(paste("Length in Exon: ", toString(mean(length))," ± ", toString(sd(length))))

Intron <- gene[which(peakAnno@detailGenomicAnnotation[["Intron"]] == TRUE)]
length <- peak@ranges@width[which(peakAnno@detailGenomicAnnotation[["Intron"]] == TRUE)]
m = data.frame(Intron, length)
write.csv(m, "Intron.csv", row.names = FALSE)

print(paste("Length in Intron: ", toString(mean(length))," ± ", toString(sd(length))))

downstream <- gene[which(peakAnno@detailGenomicAnnotation[["downstream"]] == TRUE)]
length <- peak@ranges@width[which(peakAnno@detailGenomicAnnotation[["downstream"]] == TRUE)]
m = data.frame(downstream, length)
write.csv(m, "downstream.csv", row.names = FALSE)

print(paste("Length in downstream: ", toString(mean(length))," ± ", toString(sd(length))))

distal_intergenic <- gene[which(peakAnno@detailGenomicAnnotation[["distal_intergenic"]] == TRUE)]
length <- peak@ranges@width[which(peakAnno@detailGenomicAnnotation[["distal_intergenic"]] == TRUE)]
m = data.frame(distal_intergenic, length)
write.csv(m, "distal_intergenic.csv", row.names = FALSE)

print(paste("Length in distal_intergenic: ", toString(mean(length))," ± ", toString(sd(length))))

png(filename='AnnoPie.png', width=3500, height=2000, res=600)
plotAnnoPie(peakAnno)
dev.off()

geneId <- peakAnno@anno@elementMetadata@listData[["geneId"]]
go <- enrichPathway(geneId)
p <- dotplot(go)

ggsave(
  file='GO.png',
  scale = 1,
  dpi = 600,
  limitsize = TRUE,
)

