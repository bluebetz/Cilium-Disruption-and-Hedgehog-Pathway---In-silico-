library(WGCNA)
library(dplyr)
library(limma)
library(edgeR)

#By: Gisele Rocha Edited: Gabriela
#01/06/2025
# ----------------------------
# 1. METADADOS
# ----------------------------
clinical.temp.df <- data.frame(
  ID = dataPrep2$barcode,
  tipo = dataPrep2$shortLetterCode,
  estadio = dataPrep2$ajcc_clinical_stage
)

clinical.temp.df <- clinical.temp.df %>%
  mutate(ID2 = gsub("-", ".", ID),
         estadio = gsub("[ABC]", "", estadio),
         estadio = trimws(estadio)) %>%
  mutate(Grupo_Clinico = case_when(
    tipo == "NT" ~ "Normal",
    estadio %in% c("Stage I", "Stage II", "Stage III") ~ "Inicial",
    estadio %in% c("Stage IV") ~ "Avancado",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Grupo_Clinico)) %>%
  mutate(Progresso = case_when(
    Grupo_Clinico == "Normal" ~ 0,
    Grupo_Clinico == "Inicial" ~ 1,
    Grupo_Clinico == "Avancado" ~ 2
  ))

rownames(clinical.temp.df) <- clinical.temp.df$ID2


# 2. MATRIZ DE EXPRESSÃO (corrigir nomes e VOOM)
colnames(matexpr2.out.2) <- gsub("-", ".", colnames(matexpr2.out.2))
exp.mat <- matexpr2.out.2[, colnames(matexpr2.out.2) %in% clinical.temp.df$ID2]
exp.mat <- exp.mat[, clinical.temp.df$ID2]  # Reordenar para manter a ordem dos metadados

# Voom
d0 <- DGEList(exp.mat)
design <- model.matrix(~0 + clinical.temp.df$Grupo_Clinico)
colnames(design) <- c("Normal", "Inicial", "Avancado")

y <- voom(d0, design, plot = TRUE)
y2 <- y$E  # Matriz normalizada para WGCNA

# 3. VALIDAÇÃO DE DADOS (goodSamplesGenes)
gsg <- goodSamplesGenes(y2, verbose = 3)
if (!gsg$allOK) {
  y2 <- y2[gsg$goodGenes, gsg$goodSamples]
}

# 4. SELEÇÃO DOS GENES MAIS VARIÁVEIS
n <- round(nrow(y2) * 0.25)
var_genes <- apply(y2, 1, sd)
df <- y2[order(var_genes, decreasing = TRUE)[1:n], ]
mat.df <- t(df)

## 5. SOFT THRESHOLD E TOPOLOGIA
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(mat.df, powerVector = powers, verbose = 5)

# Aqui: Plotar dois gráficos lado a lado
sizeGrWindow(12, 6)        
par(mfrow = c(1, 2))        #

cex1 = 0.9              

# Primeiro Gráfico: Scale-free topology model fit
par(mar = c(4, 4, 3, 1) + 0.1)  # Margens
plot(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale Independence",
     cex.lab = 1.2,
     cex.main = 1.4,
     cex.axis = 1.2)
text(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, 
     cex = cex1, 
     col = "red")
abline(h = 0.8, col = "blue", lty = 2)  # Linha de referência R² = 0.8

# Segundo Gráfico: Mean Connectivity
plot(sft$fitIndices[, 1], 
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean Connectivity",
     cex.lab = 1.2,
     cex.main = 1.4,
     cex.axis = 1.2)
text(sft$fitIndices[, 1], 
     sft$fitIndices[, 5],
     labels = powers, 
     cex = cex1, 
     col = "red")

# Escolher manualmente o SoftPower
softPower <- 4
adjacency <- adjacency(mat.df, type = "signed hybrid", power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Salvar TOM e dissTOM
save(y2, TOM, dissTOM, file = "y2_TOMdiss_TCGA_HNSC.RData")

# 6. CLUSTERIZAÇÃO DE GENES
geneTree2 <- hclust(as.dist(dissTOM), method = "average")
minModuleSize <- 100
dynamicMods <- cutreeDynamic(dendro = geneTree2, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)

# 7. FUSÃO DE MÓDULOS (Merge)
MEList <- moduleEigengenes(mat.df, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes")
MEDissThres <- 0.15
abline(h = MEDissThres, col = "red")

merge <- mergeCloseModules(mat.df, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

plotDendroAndColors(geneTree2, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleColors <- mergedColors
MEs <- mergedMEs

# Salvar módulos e árvore de clusterização
save(moduleColors, MEs, geneTree2, file = "SALVEmergedDynamics_TCGA_HNSC.RData")


# 8. ANÁLISE DE ASSOCIAÇÃO COM DADOS CLÍNICOS

meta.categorias <- data.frame(
  Normal = ifelse(clinical.temp.df$Grupo_Clinico == "Normal", 1, 0),
  Inicial = ifelse(clinical.temp.df$Grupo_Clinico == "Inicial", 1, 0),
  Avancado = ifelse(clinical.temp.df$Grupo_Clinico == "Avancado", 1, 0)
)
rownames(meta.categorias) <- clinical.temp.df$ID2

nGenes <- ncol(mat.df)
nSamples <- nrow(mat.df)

MEs0 <- moduleEigengenes(mat.df, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

moduleTraitCor <- cor(MEs, meta.categorias, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

output_dir <- "/Users/macbook/Desktop/TCGA/"  # ajuste o caminho
file_name <- paste0(output_dir, "ModuleTrait_Heatmap_TCGA_HNSC.png")

# Abrir dispositivo gráfico
png(file = file_name, width = 6000, height = 6000, res = 1000)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(meta.categorias),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.lab.y = 0.5,
               cex.lab.x = 0.5,
               cex.text = 0.5,
               zlim = c(-1, 1),
               colorBarWidth = 0.03)

dev.off()

# 9. ANÁLISE INTRAMODULAR
prog.df <- data.frame(Progresso = clinical.temp.df$Progresso)
rownames(prog.df) <- clinical.temp.df$ID2

modNames <- substring(names(MEs), 3)

# Correlação dos genes com os módulos
geneModuleMembership <- as.data.frame(cor(mat.df, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

# Correlação dos genes com a variável de interesse
geneTraitSignificance <- as.data.frame(cor(mat.df, prog.df, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- "GS.Progresso"
names(GSPvalue) <- "p.GS.Progresso"

# Dataframe inicial
geneInfo0 <- data.frame(
  Gene = colnames(mat.df),
  moduleColor = moduleColors,
  geneTraitSignificance,
  GSPvalue
)

# Adicionar MM e p.MM
for (mod in 1:ncol(geneModuleMembership)) {
  geneInfo0 <- cbind(geneInfo0, geneModuleMembership[, mod], MMPvalue[, mod])
  names(geneInfo0)[(ncol(geneInfo0)-1):ncol(geneInfo0)] <- c(
    paste("MM", modNames[mod], sep = "."),
    paste("p.MM", modNames[mod], sep = ".")
  )
}

# Ordenação
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Progresso))
geneInfo <- geneInfo0[geneOrder, ]

# 10. CRIAR DATAFRAMES DOS MÓDULOS
turquoise.df <- geneInfo[geneInfo$moduleColor == "turquoise", ]
yellow.df    <- geneInfo[geneInfo$moduleColor == "yellow", ]
green.df     <- geneInfo[geneInfo$moduleColor == "green", ]

# 11. EXPORTAR MÓDULOS DE INTERESSE
dir.out <- "/Users/macbook/Desktop/TCGA/"

if (!dir.exists(dir.out)) dir.create(dir.out, recursive = TRUE)

# Exportar CSVs
write.csv(turquoise.df, file = file.path(dir.out, "turquoise_module_genes.csv"), row.names = FALSE)
write.csv(yellow.df, file = file.path(dir.out, "yellow_module_genes.csv"), row.names = FALSE)
write.csv(green.df, file = file.path(dir.out, "green_module_genes.csv"), row.names = FALSE)

# ----------------------------
# 12. SALVAR OBJETOS FINAIS
# ----------------------------
save(turquoise.df, yellow.df, green.df, file = "Modulos.df_gabi.RData")

print(getwd())

