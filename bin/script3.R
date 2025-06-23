# script3.R

# Exemplo simples de análise em R
cat("Executando Análise 3 em R...\n")

# Criar pasta de saída se não existir
output_dir <- file.path("outputs", "WORK3")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Pasta de saída criada em:", output_dir, "\n")
}

# Carregar pacote necessário
library(ggplot2)

# Criar e salvar um gráfico de exemplo
data(mpg)
p <- ggplot(mpg, aes(x=displ, y=hwy)) +
  geom_point() +
  labs(title="Análise 3: Exemplo de Gráfico em R")

# Caminho completo do arquivo
output_file <- file.path(output_dir, "grafico_analise3.png")
ggsave(output_file, plot = p)

cat("Gráfico salvo em:", output_file, "\n")
cat("Análise 3 concluída.\n")
