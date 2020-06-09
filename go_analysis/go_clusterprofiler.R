

## clusterProfiler
## 2019-02-24


library(AnnotationHub)
hub <- AnnotationHub()

# query(hub, "Cricetulus")

## ID conversion
library(org.Hs.eg.db)
library(clusterProfiler)
library(org.Mm.eg.db)
data(geneList, package = "DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

ego <- enrichGO(gene = gene,
                universe = names(geneList),
                OrgDb = org.Hs.eg.db,
                ont = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05)
head(ego)

ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego2))

ego3 <- enrichGO(gene         = gene.df$SYMBOL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego3))





































