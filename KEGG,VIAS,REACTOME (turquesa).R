# Pacotes necessários
library(readxl)
library(dplyr)
library(ggplot2)
library(forcats)
library(janitor)
library(ggtext)

# Função para ler uma aba do Excel
ler_comunidade <- function(sheet_name) {
  read_excel("/Users/macbook/Desktop/GO todas as comunidades.xlsx", sheet = sheet_name) %>%
    clean_names() %>%
    mutate(comunidade = sheet_name)
}

# Lê a aba da comunidade TURQUOISE
dados <- ler_comunidade("Turquoise")

# Lista dos termos que você deseja mostrar
termos_selecionados <- c(
  # KEGG Pathway
  "Chemical carcinogenesis",
  "Drug metabolism – cytochrome P450",
  "Metabolism of xenobiotics by cytochrome P450",
  "Glutathione metabolism",
  "Arachidonic acid metabolism",
  # GO Biological Processes
  "postsynaptic specialization assembly",
  "odontogenesis",
  "odontogenesis of dentin-containing tooth",
  "positive regulation of mitotic nuclear division",
  "acute inflammatory response",
  # REACTOME
  "AGE-RAGE signaling pathway in diabetic complications",
  "Interleukin-4 and Interleukin-13 signaling",
  "Interleukin-10 signaling",
  "Binding and Uptake of Ligands by Scavenger Receptors",
  "Neurexins and neuroligins"
)

# Filtra os dados para incluir apenas os termos desejados
dados_filtrados <- dados %>%
  filter(descricao %in% termos_selecionados) %>%
  mutate(
    categoria = case_when(
      descricao %in% c(
        "Chemical carcinogenesis",
        "Drug metabolism – cytochrome P450",
        "Metabolism of xenobiotics by cytochrome P450",
        "Glutathione metabolism",
        "Arachidonic acid metabolism"
      ) ~ "KEGG Pathway",
      descricao %in% c(
        "postsynaptic specialization assembly",
        "odontogenesis",
        "odontogenesis of dentin-containing tooth",
        "positive regulation of mitotic nuclear division",
        "acute inflammatory response"
      ) ~ "GO Biological Processes",
      TRUE ~ "REACTOME"
    ),
    descricao = fct_rev(factor(descricao, levels = termos_selecionados)),
    comunidade = "Turquoise"
  )

# Cores por categoria
cores_categoria <- c(
  "KEGG Pathway" = "red",
  "GO Biological Processes" = "blue",
  "REACTOME" = "darkgreen"
)

# Gráfico
g_turquoise <- ggplot(dados_filtrados, aes(x = contagem, y = descricao, fill = categoria)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = cores_categoria) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Top 15 Enriched Terms – Module Turquoise",
    x = "Gene Count",
    y = NULL,
    fill = "Category"
  ) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

# Exibir gráfico
print(g_turquoise)

# Salvar imagem (opcional)
ggsave("/Users/macbook/Desktop/turquoise_module_enrichment.png", g_turquoise, width = 10, height = 7, dpi = 300)
