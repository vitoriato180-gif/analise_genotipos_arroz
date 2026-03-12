

rm(list = ls())


# ============================================================
# 1. PACOTES
# ============================================================
packages <- c(
  "lme4", "lmerTest", "GA", "dplyr", "caret", "glmnet",
  "scales", "doParallel", "readxl", "DHARMa", "lubridate",
  "forcats", "parallel"
)

for(pkg in packages){
  if(!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# ============================================================
# 2. CARREGAMENTO E PREPARAÇÃO DOS DADOS
# ============================================================

setwd("C:\\Users\\elisv\\Downloads\\novo")  # trocar aqui

tryCatch({
  final_data <- read_excel("dados_dezembro.xlsx")
}, error = function(e) {
  stop("Erro ao ler o arquivo. Verifique se 'dados_dezembro.xlsx' está na pasta correta.")
})

# Filtros de Local e Estado
final_data <- final_data %>%
  filter(
    !LOCATION %in% c("SANTO ANTÔNIO DE GOIÁS", "SÃO MATEUS DO MARANHÃO", "CODO", "PARAGOMINAS"),
    ST != "TO"
  )

# Criação da variável Solo Dominante
final_data <- final_data %>%
  mutate(
    Simb_dom = case_when(
      grepl("Latossolo", Simb) ~ "Latossolo",
      grepl("Argissolo", Simb) ~ "Argissolo",
      grepl("Neossolo", Simb) ~ "Neossolo",
      grepl("Plintossolo", Simb) ~ "Plintossolo",
      TRUE ~ "Outros"
    )
  )

# Agregação (Médias por Ensaio e Genótipo)
final_data <- final_data %>%
  mutate(across(where(is.character), as.factor)) %>%
  group_by(TRIAL, GEN) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    across(where(is.factor), ~ first(.x)),
    .groups = "drop"
  )

# Remoção de colunas desnecessárias
remover <- c(
  "TRIAL","DATE","ST","REGIAO","TYPE","DESIGN",
  "LATITUDE","LONGITUDE","Ano_Semeadura",
  "PI_dias","PM_dias","DTF","REP", "N_REP","LOCATION",
  "Simb" 
)

final_data <- final_data %>%
  filter(!is.na(GY), GY > 0) %>%
  select(-any_of(remover)) %>%
  select(where(~ !is.numeric(.x) || var(.x, na.rm = TRUE) != 0)) 

# ============================================================
# 3. DEFINIÇÃO DE CANDIDATOS E CORREÇÃO DE MULTICOL.
# ============================================================
cols_nao_candidatas <- c("GY", "GEN", "LOC", "Simb_dom") 
fixed_candidates <- setdiff(names(final_data), cols_nao_candidatas)

# Removemos 'Classe_AD' pois ela briga com 'AD_UM'
if("Classe_AD" %in% fixed_candidates && "AD_UM" %in% fixed_candidates) {
  message("--> Removendo 'Classe_AD' da lista para evitar multicolinearidade.")
  fixed_candidates <- setdiff(fixed_candidates, "Classe_AD")
}

fixed_candidates <- fixed_candidates[fixed_candidates %in% names(final_data)]

# ============================================================
# 4. FUNÇÕES AUXILIARES (FITNESS E LIMPEZA)
# ============================================================

# 4.1 Função de Fitness para o GA (Minimizar RMSE)
fitness_glmm_Gaussiana_rmse <- function(cromossomo, random_effect){
  vars <- fixed_candidates[cromossomo == 1]
  if(length(vars) == 0) return(-1e6) # Penalidade se vazio
  
  vars_final <- unique(c("Simb_dom", vars))
  form_str <- paste("GY ~", paste(vars_final, collapse = " + "), "+", random_effect)
  
  # REML=FALSE para seleção (comparável via ML)
  fit <- tryCatch(
    lmer(as.formula(form_str), data = final_data, REML = FALSE),
    error = function(e) NULL, warning = function(w) NULL
  )
  
  if(is.null(fit)) return(-1e6)
  
  y_hat <- predict(fit, re.form = NULL)
  rmse <- sqrt(mean((final_data$GY - y_hat)^2))
  return(-rmse) 
}
set.seed(2025)
# 4.2 Função  para rodar o GA
rodar_GA_cenario <- function(nome_cenario, random_effect, cluster_obj){
  message(paste("\n>>> Rodando GA para:", nome_cenario))
  
  resultado_GA <- ga(
    type = "binary",
    fitness = function(x) fitness_glmm_Gaussiana_rmse(x, random_effect),
    nBits = length(fixed_candidates),
    popSize = 50,      
    maxiter = 40,      
    run = 10,
    pmutation = 0.1,
    elitism = 2,
    parallel = cluster_obj,
    monitor = FALSE,
    seed = 2025
  )
  
  sol <- resultado_GA@solution
  if(is.matrix(sol)) sol <- sol[1, ]
  vars_sel <- fixed_candidates[sol == 1]
  
  return(vars_sel)
}

# 4.3 Função de Limpeza
# Remove variáveis uma a uma até todas serem significativas (p < 0.05)
limpar_modelo_misto <- function(modelo_inicial, alpha = 0.05) {
  modelo_atual <- modelo_inicial
  
  message("   Iniciando limpeza estatística...")
  
  while(TRUE) {
    # Recalcula summary com o modelo atual
    summ <- summary(modelo_atual)$coefficients
    
    # Filtra linhas para ignorar Intercepto
    vars_analise <- rownames(summ)[rownames(summ) != "(Intercept)"]
    
    # Se só sobrou intercepto, para
    if(length(vars_analise) == 0) break
    
    # Extrai p-valores
    p_vals <- summ[vars_analise, "Pr(>|t|)", drop = FALSE]
    
    # Acha o pior p-valor
    max_p <- max(p_vals)
    
    # Verifica se o pior ainda é ruim (maior que alpha)
    if(max_p > alpha) {
      # Identifica o nome da linha 
      term_to_remove <- rownames(p_vals)[which.max(p_vals)]
      
      # Tenta descobrir o nome da variável original
      # (Se for fator, o nome no summary é diferente do nome da coluna)
      var_to_remove <- term_to_remove
      
      # Assim: se o termo contém o nome de uma coluna candidata, remove a coluna
      for(col in c(fixed_candidates, "Simb_dom")) {
        if(grepl(col, term_to_remove)) {
          var_to_remove <- col
          break
        }
      }
      
      message(paste("   Removendo:", var_to_remove, "(p-valor max:", round(max_p, 4), ")"))
      
      # Atualiza fórmula: . ~ . - variavel
      nova_form <- update(formula(modelo_atual), paste(". ~ . -", var_to_remove))
      modelo_atual <- lmer(nova_form, data = final_data, REML = TRUE)
      
    } else {
      message("   Limpeza concluída! Todas as variáveis remanescentes são significativas.")
      break
    }
  }
  return(modelo_atual)
}

# ============================================================
# 5. EXECUÇÃO DOS CENÁRIOS
# ============================================================

# Configura Cluster
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl, varlist = c("final_data", "fixed_candidates", "fitness_glmm_Gaussiana_rmse", "lmer"))

# --- CENÁRIO 1: (1 | GEN:LOC) ---
vars_1 <- rodar_GA_cenario("Cenário 1", "(1 | GEN:LOC)", cl)
form_1 <- as.formula(paste("GY ~", paste(unique(c("Simb_dom", vars_1)), collapse = " + "), "+ (1 | GEN:LOC)"))
mod_1_full <- lmer(form_1, data = final_data, REML = TRUE)
mod_1_final <- limpar_modelo_misto(mod_1_full) # LIMPEZA AQUI

# --- CENÁRIO 2: (1 | GEN:Simb_dom) ---
vars_2 <- rodar_GA_cenario("Cenário 2", "(1 | GEN:Simb_dom)", cl)
form_2 <- as.formula(paste("GY ~", paste(unique(c("Simb_dom", vars_2)), collapse = " + "), "+ (1 | GEN:Simb_dom)"))
mod_2_full <- lmer(form_2, data = final_data, REML = TRUE)
mod_2_final <- limpar_modelo_misto(mod_2_full) # LIMPEZA AQUI

# --- CENÁRIO 3: (1 | GEN:LOC) + (1 | GEN:Simb_dom) ---
vars_3 <- rodar_GA_cenario("Cenário 3", "(1 | GEN:LOC) + (1 | GEN:Simb_dom)", cl)
form_3 <- as.formula(paste("GY ~", paste(unique(c("Simb_dom", vars_3)), collapse = " + "), "+ (1 | GEN:LOC) + (1 | GEN:Simb_dom)"))
mod_3_full <- lmer(form_3, data = final_data, REML = TRUE)
mod_3_final <- limpar_modelo_misto(mod_3_full) # LIMPEZA AQUI
# (1 | 1:GEN)
stopCluster(cl)

# ============================================================
# 6. RESULTADOS FINAIS E DIAGNÓSTICOS
# ============================================================

print(summary(mod_1_final))
print(summary(mod_2_final))
print(summary(mod_3_final))

cat("\nCOMPARACÃO DE AIC:\n")
print(AIC(mod_1_final, mod_2_final, mod_3_final))

# ============================================================
# 7. VISUALIZAÇÃO DHARMa (Focado no Modelo 1)
# ============================================================

res_dharma <- simulateResiduals(mod_1_final, n = 1000)

plot(res_dharma)

# Testes formais adicionais
cat("\nTeste de Uniformidade (KS test):\n")
print(testUniformity(res_dharma))

cat("\nTeste de Dispersão:\n")
print(testDispersion(res_dharma))

cat("\nTeste de Outliers:\n")
print(testOutliers(res_dharma))

# ============================================================
# DIAGNÓSTICOS DHARMa - MODELOS 2 E 3
# ============================================================

res_dharma_2 <- simulateResiduals(mod_2_final, n = 1000)
plot(res_dharma_2)

# Testes Formais Modelo 2
print(testUniformity(res_dharma_2)) # Verifica normalidade/uniformidade
print(testDispersion(res_dharma_2)) # Verifica variância


# --- MODELO 3: (1 | GEN:LOC) + (1 | GEN:Simb_dom) ---

res_dharma_3 <- simulateResiduals(mod_3_final, n = 1000)
plot(res_dharma_3)

# Testes Formais Modelo 3
print(testUniformity(res_dharma_3))
print(testDispersion(res_dharma_3))


# ============================================================
# CENÁRIO 4: ESTRUTURA COMPLETA (GEN + GEN:LOC)
# ============================================================

# Vamos aproveitar as variáveis que o GA selecionou para o Cenário 1 (Local),
# pois elas explicam o ambiente, e adicionar o efeito principal de GEN.


# Fórmula: Variáveis do Modelo 1 + GEN (efeito principal) + GEN:LOC (interação)
# Nota: Se 'vars_1' não estiver na memória, substitua pelos nomes das colunas ou rode o GA anterior.
form_4 <- as.formula(paste(
  "GY ~", paste(unique(c("Simb_dom", vars_1)), collapse = " + "), 
  "+ (1 | GEN) + (1 | GEN:LOC)"
))

# Ajusta o modelo completo
mod_4_full <- lmer(form_4, data = final_data, REML = TRUE)

# Aplica a limpeza (Backward Elimination) para garantir significância dos fixos
mod_4_final <- limpar_modelo_misto(mod_4_full)

#ESTUDAR RANDIM EFFECTS
# ============================================================
# COMPARANDO COM O VENCEDOR ANTERIOR
# ============================================================


print(summary(mod_4_final))

# Verifica quanto da variância ficou para GEN e quanto para GEN:LOC
print(VarCorr(mod_4_final))


# Comparando o Modelo 1 (só interação) com o Modelo 4 (Gen + Interação)
aic_tab <- AIC(mod_1_final, mod_4_final)
print(aic_tab)

# Quem ganhou?
melhor_aic <- rownames(aic_tab)[which.min(aic_tab$AIC)]
message(paste("O melhor modelo pelo critério de AIC é:", melhor_aic))

# ============================================================
# DIAGNÓSTICO DHARMa DO NOVO MODELO
# ============================================================
res_dharma_4 <- simulateResiduals(mod_4_final, n = 1000)
plot(res_dharma_4)

# ============================================================
# EXPANDINDO MODELOS 2 E 3 (ADICIONANDO EFEITO PRINCIPAL DE GENÓTIPO)
# ============================================================

# ------------------------------------------------------------
# 1. NOVO MODELO 2: (1 | GEN) + (1 | GEN:Simb_dom)
# ------------------------------------------------------------

# Monta a fórmula com vars_2 (do GA do cenário de solo)
form_2_new <- as.formula(paste(
  "GY ~", paste(unique(c("Simb_dom", vars_2)), collapse = " + "), 
  "+ (1 | GEN) + (1 | GEN:Simb_dom)"
))

# Roda e Limpa
mod_2_new_full <- lmer(form_2_new, data = final_data, REML = TRUE)
mod_2_new_final <- limpar_modelo_misto(mod_2_new_full)


# ------------------------------------------------------------
# 2. NOVO MODELO 3: (1 | GEN) + (1 | GEN:LOC) + (1 | GEN:Simb_dom)
# ------------------------------------------------------------

# Monta a fórmula com vars_3 (do GA do cenário misto)
form_3_new <- as.formula(paste(
  "GY ~", paste(unique(c("Simb_dom", vars_3)), collapse = " + "), 
  "+ (1 | GEN) + (1 | GEN:LOC) + (1 | GEN:Simb_dom)"
))

# Roda e Limpa
mod_3_new_full <- lmer(form_3_new, data = final_data, REML = TRUE)
mod_3_new_final <- limpar_modelo_misto(mod_3_new_full)


# ============================================================
# 3. RESULTADOS E COMPARAÇÕES (QUEM GANHOU?)
# ============================================================


# Comparação Modelo 2 (Só Interação) vs Modelo 2 New (Gen + Interação)
cat("\n--- Batalha do Cenário 2 (Solo) ---\n")
print(AIC(mod_2_final, mod_2_new_final))

# Comparação Modelo 3 (Só Interação) vs Modelo 3 New (Gen + Interação)
cat("\n--- Batalha do Cenário 3 (Completo) ---\n")
print(AIC(mod_3_final, mod_3_new_final))

# ============================================================
# 4. VERIFICAR VARIÂNCIAS (Herdabilidade implícita)
# ============================================================
message("\n--- Variâncias do Novo Modelo 3 ---")

print(VarCorr(mod_3_new_final))

# ============================================================
# 5. DIAGNÓSTICO DHARMa DOS NOVOS MODELOS
# ============================================================
plot(simulateResiduals(mod_2_new_final))
plot(simulateResiduals(mod_3_new_final))

setwd("C:\\Users\\elisv\\Downloads\\ultima_22dez")  # trocar aqui


# ============================================================
# SALVANDO OS GRÁFICOS (Alta Resolução - 300 DPI)
# ============================================================
#sem interações corretas 
# 1. Salvar Modelo 1
png(filename = "DHARMa_Modelo_1.png", 
    units = "in", width = 12, height = 6, res = 300)
plot(res_dharma) # ou res_dharma_1 dependendo do nome que você usou
dev.off() 

# 2. Salvar Modelo 2
png(filename = "DHARMa_Modelo_2.png", 
    units = "in", width = 12, height = 6, res = 300)
plot(res_dharma_2)
dev.off()

# 3. Salvar Modelo 3
png(filename = "DHARMa_Modelo_3.png", 
    units = "in", width = 12, height = 6, res = 300)
plot(res_dharma_3)
dev.off()







# ============================================================
# SALVANDO DIAGNÓSTICOS DHARMa (CENÁRIOS COM GENÓTIPO)
# ============================================================

# Confirma se o pacote está carregado
if(!requireNamespace("DHARMa", quietly = TRUE)) library(DHARMa)



# ------------------------------------------------------------
# 1. CENÁRIO 4: (1 | GEN) + (1 | GEN:LOC)
# ------------------------------------------------------------
if(exists("mod_4_final")) {
  message("Salving DHARMa para Cenário 4...")
  
  # Simula os resíduos
  res_4 <- simulateResiduals(mod_4_final, n = 1000)
  
  # Abre o arquivo PNG
  png(filename = "DHARMa_Cenario_4_Gen_Loc.png", 
      units = "in", width = 12, height = 6, res = 300)
  
  # Plota
  plot(res_4)
  
  # Fecha o arquivo (Salva)
  dev.off()
} else {
  warning("Modelo 'mod_4_final' não encontrado na memória.")
}

# ------------------------------------------------------------
# 2. NOVO MODELO 2: (1 | GEN) + (1 | GEN:Simb_dom)
# ------------------------------------------------------------
if(exists("mod_2_new_final")) {
  message("Salvando DHARMa para Novo Modelo 2...")
  
  res_2n <- simulateResiduals(mod_2_new_final, n = 1000)
  
  png(filename = "DHARMa_Modelo_2_Gen_Solo.png", 
      units = "in", width = 12, height = 6, res = 300)
  
  plot(res_2n)
  
  dev.off()
} else {
  warning("Modelo 'mod_2_new_final' não encontrado na memória.")
}

# ------------------------------------------------------------
# 3. NOVO MODELO 3: (1 | GEN) + (1 | GEN:LOC) + (1 | GEN:Simb_dom)
# ------------------------------------------------------------
if(exists("mod_3_new_final")) {
  message("Salvando DHARMa para Novo Modelo 3 (Completo)...")
  
  res_3n <- simulateResiduals(mod_3_new_final, n = 1000)
  
  png(filename = "DHARMa_Modelo_3_Completo.png", 
      units = "in", width = 12, height = 6, res = 300)
  
  plot(res_3n)
  
  dev.off()
} else {
  warning("Modelo 'mod_3_new_final' não encontrado na memória.")
}






########IMPORTANCIA RELATIVA POR FASE###################




# ==============================================================================
# 0. CARREGAR PACOTES E DEFINIR A FUNÇÃO DE GRÁFICOS
# ==============================================================================
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(lme4)

# Função reutilizável para gerar os gráficos
gerar_analise_importancia <- function(modelo_final, nome_cenario, nome_arquivo_base) {
  
  # Verificação de segurança: O modelo existe?
  if(!exists(deparse(substitute(modelo_final))) && is.null(modelo_final)) {
    warning(paste("O modelo", nome_cenario, "não foi encontrado ou está vazio."))
    return(NULL)
  }
  
  # 1. Extrair Coeficientes
  coefs <- fixef(modelo_final)
  importancia_raw <- abs(coefs[names(coefs) != "(Intercept)"])
  
  if(length(importancia_raw) == 0) {
    message(paste("O modelo", nome_cenario, "não possui efeitos fixos para plotar."))
    return(NULL)
  }
  
  # 2. Calcular %
  total_impacto <- sum(importancia_raw)
  df_imp <- data.frame(
    Variavel = names(importancia_raw),
    Valor_Abs = importancia_raw,
    Porcentagem = (importancia_raw / total_impacto) * 100
  )
  
  # 3. Categorizar Fases (Ajuste os nomes conforme suas colunas reais)
  df_imp <- df_imp %>%
    mutate(Fase = case_when(
      str_detect(Variavel, "veg_") ~ "Vegetativo",
      str_detect(Variavel, "repro_") ~ "Reprodutivo",
      str_detect(Variavel, "gf_") ~ "Enchimento de Grãos",
      str_detect(Variavel, "Simb|Solo|Latossolo|Neossolo") ~ "Solo/Ambiente",
      str_detect(Variavel, "AD_UM") ~ "Solo/Ambiente",
      TRUE ~ "Outros"
    )) %>%
    arrange(desc(Porcentagem))
  
  # 4. Resumo Agrupado
  resumo_fases <- df_imp %>%
    group_by(Fase) %>%
    summarise(Total_Pct = sum(Porcentagem)) %>%
    mutate(Label = paste0(round(Total_Pct, 1), "%"))
  
  # --- PLOT 1: POR FASE ---
  p1 <- ggplot(resumo_fases, aes(x = reorder(Fase, Total_Pct), y = Total_Pct, fill = Fase)) +
    geom_col(width = 0.6, color = "black", alpha = 0.8) +
    geom_text(aes(label = Label), hjust = -0.1, fontface = "bold") +
    scale_y_continuous(limits = c(0, max(resumo_fases$Total_Pct) + 15)) +
    labs(
      title = paste("Importância por Fase -", nome_cenario),
      y = "Importância Relativa (%)", x = NULL
    ) +
    coord_flip() +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
  
  # --- PLOT 2: TOP VARIÁVEIS ---
  p2 <- ggplot(head(df_imp, 15), aes(x = reorder(Variavel, Porcentagem), y = Porcentagem, color = Fase)) +
    geom_segment(aes(xend = Variavel, yend = 0), size = 1.1) +
    geom_point(size = 4) +
    geom_text(aes(label = round(Porcentagem, 1)), vjust = -0.7, color = "black", size = 3) +
    labs(
      title = paste("Top Variáveis -", nome_cenario),
      y = "Contribuição (%)", x = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  
  # Salvar Arquivos
  ggsave(paste0(nome_arquivo_base, "_Fases.png"), p1, width = 8, height = 5, dpi = 300)
  ggsave(paste0(nome_arquivo_base, "_Variaveis.png"), p2, width = 10, height = 6, dpi = 300)
  
  return(list(resumo = resumo_fases, tabela_completa = df_imp))
}


# ==============================================================================
# 1. ANÁLISE DO CENÁRIO 4: (1 | GEN) + (1 | GEN:LOC)
# ==============================================================================

# Certifique-se que mod_4_final está criado na memória
if(exists("mod_4_final")) {
  res_4 <- gerar_analise_importancia(
    mod_4_final, 
    nome_cenario = "(Gen + Gen:Loc)", 
    nome_arquivo_base = "Grafico_Cenario4"
  )
  print(res_4$resumo)
} else {
  message("ATENÇÃO: mod_4_final não encontrado. Rode a etapa de modelagem anterior.")
}


# ==============================================================================
# 2. ANÁLISE DO NOVO MODELO 2: (1 | GEN) + (1 | GEN:Simb_dom)
# ==============================================================================

# Certifique-se que mod_2_new_final está criado na memória
if(exists("mod_2_new_final")) {
  res_2_new <- gerar_analise_importancia(
    mod_2_new_final, 
    nome_cenario = "(Gen + Gen:Solo)", 
    nome_arquivo_base = "Grafico_Modelo2_New"
  )
  print(res_2_new$resumo)
} else {
  message("ATENÇÃO: mod_2_new_final não encontrado.")
}


# ==============================================================================
# 3. ANÁLISE DO NOVO MODELO 3: (1 | GEN) + (1 | GEN:LOC) + (1 | GEN:Simb_dom)
# ==============================================================================

# Certifique-se que mod_3_new_final está criado na memória
if(exists("mod_3_new_final")) {
  res_3_new <- gerar_analise_importancia(
    mod_3_new_final, 
    nome_cenario = "Modelo 3", 
    nome_arquivo_base = "Grafico_Modelo3_New"
  )
  print(res_3_new$resumo)
} else {
  message("ATENÇÃO: mod_3_new_final não encontrado.")
}








# =======================================================
# SCRIPT DE CORREÇÃO: REMOÇÃO DE NA + GRÁFICO
# =======================================================

library(dplyr)
library(lme4)
library(ggplot2)
library(tidytext)




dados_limpos <- final_data %>%
  # Seleciona as colunas que vamos usar (ajuste se tiver AD_UM ou não)
  select(GY, GEN, LOC, Simb_dom, any_of("AD_UM")) %>% 
  # O comando MÁGICO: Remove qualquer linha que tenha buraco (NA)
  na.omit()

message(paste("Linhas originais:", nrow(final_data)))
message(paste("Linhas limpas (usadas no modelo):", nrow(dados_limpos)))

# =======================================================
# PASSO 2: RECRIAR O MODELO COM DADOS LIMPOS
# =======================================================

# Definindo fórmula segura
formula_segura <- "GY ~ Simb_dom + (1 | GEN) + (1 | GEN:LOC) + (1 | GEN:Simb_dom)"
if("AD_UM" %in% names(dados_limpos)) formula_segura <- paste(formula_segura, "+ AD_UM")

# Rodando o modelo na tabela 'dados_limpos'
mod_3_new_final <- lmer(as.formula(formula_segura), data = dados_limpos, REML = TRUE)

# =======================================================
# PASSO 3: GERAR O GRÁFICO 
# =======================================================

# Agora usamos 'dados_limpos' que tem 369 linhas, igualzinho à predição!
dados_grafico <- dados_limpos %>%
  mutate(Predito = predict(mod_3_new_final))

# Calcular Médias por Ambiente
media_ambientes <- dados_grafico %>%
  group_by(Simb_dom) %>%
  summarise(Media_Ambiente = mean(Predito, na.rm = TRUE), .groups = "drop")

# Calcular RPIPY
dados_rpipy <- dados_grafico %>%
  group_by(GEN, Simb_dom) %>%
  summarise(Prod_Genotipo = mean(Predito, na.rm = TRUE), .groups = "drop") %>%
  left_join(media_ambientes, by = "Simb_dom") %>%
  mutate(
    RPIPY = ((Prod_Genotipo - Media_Ambiente) / Media_Ambiente) * 100
  )

# =======================================================
# PASSO 4: PLOTAGEM COM LEGENDA (CORRIGIDO)
# =======================================================

plt <- ggplot(dados_rpipy, aes(x = reorder_within(GEN, RPIPY, Simb_dom), 
                               y = RPIPY, 
                               color = RPIPY)) +
  
  # Linha de referência (Zero)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Haste (Lollipop)
  geom_segment(aes(xend = reorder_within(GEN, RPIPY, Simb_dom), 
                   yend = 0), size = 0.8, alpha = 0.7) +
  
  # Pontos
  geom_point(size = 4) +
  
  # Facetas
  facet_wrap(~ Simb_dom, scales = "free_y", ncol = 2) + 
  
  # Ajustes de Eixos
  coord_flip() +
  scale_x_reordered() +
  
  # --- AQUI ESTÁ A MUDANÇA DA LEGENDA ---
  # Removemos o guide="none" e adicionamos um título para a barra
  scale_color_gradient2(
    low = "#B2182B", 
    mid = "grey95", 
    high = "#2166AC", 
    midpoint = 0, 
    name = "Performance (%)" # Título da legenda
  ) +
  
  # Textos
  labs(
    title = "Performance Relativa por Tipo de Solo",
    subtitle = "Ranking de Genótipos",
    y = "Diferença em relação à média do solo (%)",
    x = NULL
  ) +
  
  # Tema
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    
    # --- POSICIONAMENTO DA LEGENDA ---
    legend.position = "bottom",       # Legenda embaixo
    legend.justification = "center",  # Centralizada
    legend.key.width = unit(1.5, "cm") # Barra de cor mais larguinha para ler melhor
  )

# Exibir e Salvar
print(plt)
ggsave("Grafico_RPIPY_Final_Com_Legenda.png", plot = plt, width = 12, height = 8, dpi = 300)





# =======================================================
# GRÁFICO RPIPY POR LOCAL (LOC)
# =======================================================

library(dplyr)
library(ggplot2)
library(tidytext)


if(!exists("dados_limpos")) stop("Rode o script de limpeza de NAs anterior primeiro!")

dados_grafico_loc <- dados_limpos %>%
  mutate(Predito = predict(mod_3_new_final))

# 2. Calcular Médias por LOCAL
media_locais <- dados_grafico_loc %>%
  group_by(LOC) %>%
  summarise(Media_Local = mean(Predito, na.rm = TRUE), .groups = "drop")

# 3. Calcular RPIPY por Genótipo dentro de cada Local
dados_rpipy_loc <- dados_grafico_loc %>%
  group_by(GEN, LOC) %>%
  summarise(Prod_Genotipo = mean(Predito, na.rm = TRUE), .groups = "drop") %>%
  left_join(media_locais, by = "LOC") %>%
  mutate(
    RPIPY = ((Prod_Genotipo - Media_Local) / Media_Local) * 100
  )

# 4. Plotagem (Ajustada para muitos locais)
plt_loc <- ggplot(dados_rpipy_loc, aes(x = reorder_within(GEN, RPIPY, LOC), 
                                       y = RPIPY, 
                                       color = RPIPY)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_segment(aes(xend = reorder_within(GEN, RPIPY, LOC), 
                   yend = 0), size = 0.8, alpha = 0.7) +
  geom_point(size = 3) + # Pontos um pouco menores para caber
  
  # Facetas por LOCAL (Usei ncol=3 para distribuir melhor se forem muitos)
  facet_wrap(~ LOC, scales = "free_y", ncol = 3) + 
  
  coord_flip() +
  scale_x_reordered() +
  
  scale_color_gradient2(
    low = "#B2182B", mid = "grey95", high = "#2166AC", 
    midpoint = 0, name = "Performance (%)"
  ) +
  
  labs(
    title = "Performance Relativa por Local",
    subtitle = "Adaptação específica dos genótipos em cada município/ensaio",
    y = "Diferença em relação à média do local (%)",
    x = NULL
  ) +
  
  theme_bw(base_size = 12) + # Fonte um pouco menor
  theme(
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm")
  )

# Exibir
print(plt_loc)

# Salvar (Aumentei a altura 'height' porque locais ocupam muito espaço vertical)
ggsave("Grafico_RPIPY_Locais.png", plot = plt_loc, width = 14, height = 12, dpi = 300)




library(tidytext)



dados_para_grafico <- final_data %>%
  na.omit()


gerar_grafico_pb <- function(dados, modelo, nome_arquivo, titulo_grafico) {
  
  message(paste(">>> Processando:", titulo_grafico))
  
  # 1. Predição (Agora vai funcionar pois 'dados' tem as colunas de clima)
  dados_pred <- dados %>%
    mutate(Predito = predict(modelo, newdata = dados))
  
  # 2. RPIPY por Local
  media_locais <- dados_pred %>%
    group_by(LOC) %>%
    summarise(Media_Local = mean(Predito, na.rm = TRUE), .groups = "drop")
  
  dados_rpipy <- dados_pred %>%
    group_by(GEN, LOC) %>%
    summarise(Prod_Genotipo = mean(Predito, na.rm = TRUE), .groups = "drop") %>%
    left_join(media_locais, by = "LOC") %>%
    mutate(
      RPIPY = ((Prod_Genotipo - Media_Local) / Media_Local) * 100
    )
  
  # 3. Plotagem P&B
  plt <- ggplot(dados_rpipy, aes(x = reorder_within(GEN, RPIPY, LOC), 
                                 y = RPIPY)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.4) +
    geom_segment(aes(xend = reorder_within(GEN, RPIPY, LOC), yend = 0), 
                 color = "gray20", size = 0.6) +
    geom_point(color = "black", size = 2.5) +
    facet_wrap(~ LOC, scales = "free_y", ncol = 3) + 
    coord_flip() +
    scale_x_reordered() +
    labs(
      title = titulo_grafico,
      y = "RPIPY (%)",
      x = NULL
    ) +
    theme_bw(base_family = "serif") +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black", size = 9),
      strip.background = element_rect(fill = "white", color = "black", size = 0.5),
      strip.text = element_text(face = "bold", size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray90", linetype = "dotted"),
      panel.border = element_rect(color = "black", fill = NA)
    )
  
  # Salvar
  ggsave(nome_arquivo, plot = plt, width = 10, height = 15, dpi = 300)
  message(paste(">>> Salvo com sucesso:", nome_arquivo))
}

# =======================================================
# GERANDO OS GRÁFICOS
# =======================================================


if(exists("mod_4_final")) {
  gerar_grafico_pb(
    dados = dados_para_grafico,  # Usando o dataset completo!
    modelo = mod_4_final,
    nome_arquivo = "Grafico_RPIPY_Mod1_Local_PB.png",
    titulo_grafico = "Model 1: Genotype + Local Interaction"
  )
} else {
  message("AVISO: mod_4_final não encontrado. Rode a modelagem do Cenário 4 antes.")
}

# --- GRÁFICO 2: MODELO 2 (Interação Solo) ---
# Se o modelo mod_2_new_final existir (que é o Novo Modelo 2), usamos ele.
if(exists("mod_2_new_final")) {
  gerar_grafico_pb(
    dados = dados_para_grafico, # Usando o dataset completo!
    modelo = mod_2_new_final,
    nome_arquivo = "Grafico_RPIPY_Mod2_Solo_PB.png",
    titulo_grafico = "Model 2: Genotype + Soil Interaction"
  )
} else {
  message("AVISO: mod_2_new_final não encontrado. Rode a modelagem do Modelo 2 antes.")
}




# =======================================================
# GRÁFICO RPIPY POR LOCAL (TODOS) 
# =======================================================


dados_limpos_loc <- final_data %>%
  select(GY, GEN, LOC, Simb_dom, any_of("AD_UM")) %>% 
  na.omit()


if(!exists("mod_3_new_final")) {
  message("Recalibrando modelo rápido...")
  form <- "GY ~ Simb_dom + (1 | GEN) + (1 | GEN:LOC) + (1 | GEN:Simb_dom)"
  if("AD_UM" %in% names(dados_limpos_loc)) form <- paste(form, "+ AD_UM")
  mod_3_new_final <- lmer(as.formula(form), data = dados_limpos_loc, REML = TRUE)
}


# Predição
dados_grafico_loc <- dados_limpos_loc %>%
  mutate(Predito = predict(mod_3_new_final, newdata = dados_limpos_loc))

# Média do Local
media_locais <- dados_grafico_loc %>%
  group_by(LOC) %>%
  summarise(Media_Local = mean(Predito, na.rm = TRUE), .groups = "drop")

# RPIPY do Genótipo no Local
dados_rpipy_loc <- dados_grafico_loc %>%
  group_by(GEN, LOC) %>%
  summarise(Prod_Genotipo = mean(Predito, na.rm = TRUE), .groups = "drop") %>%
  left_join(media_locais, by = "LOC") %>%
  mutate(
    RPIPY = ((Prod_Genotipo - Media_Local) / Media_Local) * 100
  )



plt_bw <- ggplot(dados_rpipy_loc, aes(x = reorder_within(GEN, RPIPY, LOC), 
                                      y = RPIPY)) +
  
  # 1. Linha Zero (Referência)
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.4) +
  
  # 2. Haste do Lollipop (Cinza escuro para não pesar)
  geom_segment(aes(xend = reorder_within(GEN, RPIPY, LOC), yend = 0), 
               color = "gray20", size = 0.6) +
  
  # 3. Pontos (Preto Sólido)
  geom_point(color = "black", size = 2.5) +
  
  # 4. Facetas por LOCAL (Aqui mostramos TODOS)
  # ncol = 3 ou 4 dependendo de quantos locais vc tem. 3 costuma ser seguro.
  facet_wrap(~ LOC, scales = "free_y", ncol = 3) + 
  
  # 5. Ajustes Técnicos
  coord_flip() +
  scale_x_reordered() +
  
  # 6. Textos
  labs(
    title = "Relative Performance across Locations", 
    y = "RPIPY (%)",
    x = NULL
  ) +
  
  # 7. TEMA PRETO E BRANCO CLÁSSICO
  theme_bw(base_family = "serif") + 
  theme(
    # Textos pretos
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 9),
    
    # Títulos dos painéis (Caixa branca com borda preta fina)
    strip.background = element_rect(fill = "white", color = "black", size = 0.5),
    strip.text = element_text(face = "bold", size = 9),
    
    # Limpeza visual
    panel.grid.major.y = element_blank(), # Remove linhas horizontais
    panel.grid.major.x = element_line(color = "gray90", linetype = "dotted"),
    panel.border = element_rect(color = "black", fill = NA)
  )

# Exibir
print(plt_bw)

ggsave("Grafico_RPIPY_Locais_PB.png", plot = plt_bw, width = 10, height = 15, dpi = 300)










