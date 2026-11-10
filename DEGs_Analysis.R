### Template DEGs: Oral Cancer
### Script: TCGABiolinks + DESeq2 + Limma + ComplexHeatmap
### Autora: Gisele Rocha, Roberto dos Santos, Gabriela Tofalini, Viviane Silva 
##  and Clarissa Gurgel

# ----------------------------
# 1. CARREGAR PACOTES
# ----------------------------
library(BiocManager)
library(TCGAbiolinks)
library(DESeq2)
library(hgu133plus2.db)
library(limma)
library(VennDiagram)
library(RColorBrewer)
library(WGCNA)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(edgeR)
library(circlize)
library(grid)
library(org.Hs.eg.db)

# ----------------------------
# 2. DOWNLOAD E PREPARO DOS DADOS
# ----------------------------
# Diretórios
dir.dt <- "/Users/macbook/Desktop/Manuscript/"

# Carregar amostras
samples.Tum <- read.table(paste0(dir.dt, "tumoral.txt"), header = TRUE)
samples.norm <- read.table(paste0(dir.dt, "normal.txt"), header = TRUE)

# Query de dados
CancerProject <- "TCGA-HNSC"
queryDown <- GDCquery(
  project = CancerProject,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = c(samples.norm$SampleID, samples.Tum$SampleID)
)

# Download e preparo
GDCdownload(query = queryDown)
dataPrep <- GDCprepare(query = queryDown)

# Separar NT e TP samples
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataPrep), typesample = "NT")
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataPrep), typesample = "TP")

# Nova Query combinando NT + TP
queryDown2 <- GDCquery(
  project = CancerProject,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = c(samplesTP, samplesNT)
)

GDCdownload(query = queryDown2)
dataPrep2 <- GDCprepare(query = queryDown2)

# Obter matriz de expressão
HSCMatrix <- assay(dataPrep2, "stranded_first")
rownames(HSCMatrix) <- substr(rownames(HSCMatrix), 1, 15)

# ----------------------------
# 3. PROCESSAMENTO DA MATRIZ
# ----------------------------
# Anotar genes
Symbols_1 <- mapIds(hgu133plus2.db, keys = rownames(HSCMatrix),
                    column = "SYMBOL", keytype = "ENSEMBL")

# Substituir ENSEMBL por símbolos
rownames(HSCMatrix) <- Symbols_1

# Remover genes sem anotação
gene_ID_na <- is.na(Symbols_1)
gene_ID_blank <- Symbols_1 == ""
gene_ID_remove <- gene_ID_na | gene_ID_blank
expMat <- HSCMatrix[!gene_ID_remove, ]

# Remover pseudogenes e unknown
genetp <- mapIds(hgu133plus2.db, keys = rownames(expMat),
                 column = "GENETYPE", keytype = "SYMBOL")
remove_gene_type <- genetp %in% c("pseudo", "unknown")
expMat <- expMat[!remove_gene_type, ]


# ----------------------------
# 4. COLLAPSE ROWS (EXATAMENTE COMO PEDIU)
# ----------------------------
#### trabalhar com a matriz expMat deve conter 27088 linhas e 376 colunas
teste <- unique(rownames(expMat))

##collapse rows, para os genes com mais de um spot faremos a selecão do spot com maior variabilidade
collapse.object = collapseRows(datET = expMat, rowGroup = teste, rowID = teste, method = "maxRowVariance")
#collapse.object = collapseRows(datET = expressao, rowGroup = gene_ID, rowID = probe_ID, method = "Average")

dat1Collapsed = data.frame(collapse.object$group2row, collapse.object$datETcollapsed)
dim(dat1Collapsed)


# Matriz sem duplicação
expmat <- as.matrix(dat1Collapsed[, -c(1, 2)])

# ----------------------------
# 5. FILTRO DE EXPRESSÃO
# ----------------------------
cutoff <- median(expmat)
is_expressed <- expmat > cutoff
keep <- rowSums(is_expressed) > 10
matexpr2.out.2 <- expmat[keep, ]

# ----------------------------
# 6. PREPARAR METADADOS
# ----------------------------
metadata2 <- data.frame(
  Barcode = gsub("-", ".", dataPrep2$barcode),
  Status = dataPrep2$sample_type,
  Sample = dataPrep2$sample,
  ajcc_clinical_stage = dataPrep2$ajcc_clinical_stage
)

# Corrigir para normais
metadata2$ajcc_clinical_stage[metadata2$Status == "Solid Tissue Normal"] <- NA

# Classificar Grupo Clínico
metadata2$GrupoEstagio <- ifelse(
  metadata2$Status == "Solid Tissue Normal", "Normal",
  ifelse(metadata2$ajcc_clinical_stage %in% c("Stage I", "Stage II", "Stage III"), "Inicial",
         ifelse(metadata2$ajcc_clinical_stage %in% c("Stage IVA", "Stage IVB", "Stage IVC"), "Avancado", NA))
)
metadata2 <- metadata2[!is.na(metadata2$GrupoEstagio), ]

# Ajustar matriz de expressão para amostras válidas
colnames(matexpr2.out.2) <- gsub("-", ".", colnames(matexpr2.out.2))
matexpr2.out.2 <- matexpr2.out.2[, metadata2$Barcode]

# Verificação de alinhamento
stopifnot(all(colnames(matexpr2.out.2) == metadata2$Barcode))

# ----------------------------
# 7. ANÁLISE DE DEGs (LIMMA + VOOM)
# ----------------------------
# Preparar objeto DGEList
d0 <- DGEList(matexpr2.out.2)

# Design matrix para status
design <- model.matrix(~0 + metadata2$Status)  # Sem intercepto
colnames(design) <- c("Normal", "Tumour")

# Normalizar com voom
y <- voom(d0, design, plot = TRUE)

# Ajustar o modelo linear
fit <- lmFit(y, design)

# Definir contraste: Tumour - Normal
contrasts <- makeContrasts(Tumour - Normal, levels = design)

# Ajustar contrastes e calcular estatísticas
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# Obter resultados
full_results_TCGA <- topTable(fit2, number = Inf, p.value = 0.05, lfc = 1.5)
dim(full_results_TCGA)

rdata.dir <- "/Users/macbook/Desktop/TCGA/"
save(fit2, matexpr2.out.2, full_results_TCGA, metadata2, file  = paste0(rdata.dir, "TCGA_DISSERT.RData"))

# ----------------------------
# 8. FILTRAR DEGs e NORMALIZAR
# ----------------------------
p_cutoff <- 0.05
fc_cutoff <- 1.5
DEGs_COCE <- full_results_TCGA %>%
  filter(abs(logFC) > fc_cutoff & adj.P.Val < p_cutoff)

DEGs_COCE$Symbol <- rownames(DEGs_COCE)

# Selecionar genes em comum
genes_comuns <- intersect(rownames(DEGs_COCE), rownames(matexpr2.out.2))
heat.matx <- matexpr2.out.2[genes_comuns, ]

# Normalizar (Z-score por gene)
matz <- t(scale(t(heat.matx)))

# Ordenar colunas e metadata
metadata2 <- metadata2[match(colnames(matz), metadata2$Barcode), ]
grupo_estagio_vector <- factor(metadata2$GrupoEstagio, levels = c("Normal", "Inicial", "Avancado"))

ordem <- order(grupo_estagio_vector)
matz <- matz[, ordem]
metadata2 <- metadata2[ordem, ]
grupo_estagio_vector <- grupo_estagio_vector[ordem]

# ----------------------------
# 9. HEATMAP COMPACTO
# ----------------------------
# Cores dos grupos
colores_estagios <- list(Estagio = c(
  "Normal" = "#A6D854",
  "Inicial" = "#FFD92F",
  "Avancado" = "#E78AC3"
))

col_ha <- HeatmapAnnotation(
  Estagio = grupo_estagio_vector,
  col = colores_estagios,
  annotation_name_side = "right"
)

# Escala de cor Z-score
cores_heat <- colorRamp2(
  c(-2, 0, 2),
  c("#313695", "#FFFFFF", "#A50026")  # Azul escuro → Branco → Vermelho escuro
)

# Criar Heatmap
ht <- Heatmap(
  matz,
  name = "z-score",
  col = cores_heat,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_split = grupo_estagio_vector,
  top_annotation = col_ha,
  use_raster = TRUE,
  raster_device = "png",
  raster_quality = 2,
  heatmap_legend_param = list(
    title = "z-score",
    at = c(-2, -1, 0, 1, 2),
    legend_direction = "horizontal",
    title_position = "topcenter",
    legend_width = unit(6, "cm")
  )
)

# Plotar
draw(ht)

# ----------------------------
# 10. SALVAR HEATMAP
# ----------------------------
png("/Users/macbook/Desktop/TCGA/HEATMAP_PROCESSOSCLINICOS.png",
    width = 1600, height = 1800, res = 300)
draw(ht, heatmap_legend_side = "top", padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()
