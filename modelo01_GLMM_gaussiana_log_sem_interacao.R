rm(list = ls())

inicio <- Sys.time()

# Ajuste o diretório
# setwd("/Users/mariavitoriadiastorres/Downloads")
# setwd("C:\\Users\\elisv\\Downloads\\novo")

# =========================
# 0. Pacotes
# =========================
packages <- c(
  "lme4", "lmerTest", "GA", "dplyr", "caret", "glmnet",
  "scales", "doParallel", "readxl", "DHARMa", "lubridate",
  "forcats"
)

for(pkg in packages){
  if(!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# =========================
# 1. Carregar dados
# =========================
# Certifique-se que o arquivo está no diretório correto
tryCatch({
  final_data <- read_excel("dados_dezembro.xlsx")
}, error = function(e) {
  stop("Erro ao ler o arquivo. Verifique o caminho e o nome do arquivo.")
})

# =========================
# 2. Filtros iniciais
# =========================
final_data <- final_data %>%
  filter(
    !LOCATION %in% c(
      "SANTO ANTÔNIO DE GOIÁS",
      "SÃO MATEUS DO MARANHÃO",
      "CODO",
      "PARAGOMINAS"
    ),
    ST != "TO"
  )

# =========================
# 3. Criar SOLO DOMINANTE
# =========================
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

# =========================
# 4. Cheque de segurança
# =========================
check_simb <- final_data %>%
  group_by(TRIAL, GEN) %>%
  summarise(
    n_simb = n_distinct(Simb_dom),
    .groups = "drop"
  ) %>%
  filter(n_simb > 1)

if(nrow(check_simb) > 0) warning("Atenção: Existem Trials/Gen com mais de um solo dominante!")

# =========================
# 5. AGREGAÇÃO
# =========================

# Converter character para factor ANTES de agregar
final_data <- final_data %>%
  mutate(across(where(is.character), as.factor))

final_data_agg <- final_data %>%
  group_by(TRIAL, GEN) %>%
  summarise(
    # Numéricas → média
    across(where(is.numeric), mean, na.rm = TRUE),
    
    # Fatores → primeiro 
    across(where(is.factor), ~ first(.x)),
    
    .groups = "drop"
  )

final_data <- final_data_agg

# =========================
# 6. Limpezas finais
# =========================
remover <- c(
  "TRIAL","DATE","ST","REGIAO","LOC","TYPE","DESIGN",
  "LOCATION","LATITUDE","LONGITUDE","Ano_Semeadura",
  "PI_dias","PM_dias","DTF"
)

final_data <- final_data %>%
  filter(
    !is.na(GY),
    GY > 0,
    !is.na(Classe_AD),
    !is.na(AD_UM)
  ) %>%
  select(-any_of(remover)) %>%
  # Remove colunas com variância zero (apenas para numéricas)
  select(where(~ !is.numeric(.x) || var(.x, na.rm = TRUE) != 0)) %>%
  # Remove colunas constantes (para todos os tipos)
  select(where(~ n_distinct(.) > 1))

# =========================
# 7. Garantir fatores e estrutura
# =========================
# Convertendo o que for texto para fattor
final_data <- final_data %>%
  mutate(across(where(is.character), as.factor))

# =========================
# Preparação para o GA
# =========================

random_effect <- "GEN"

# Definição das candidatas: APENAS NUMÉRICAS e excluindo a resposta/efeito aleatório
fixed_candidates <- final_data %>%
  select(where(is.numeric)) %>%     # <--- OBRIGATÓRIO: Só numéricas entram no scale e no GA
  select(-GY) %>%                   # Remove a variável resposta
  names()

# Padronizar apenas as candidatas numéricas
for(col in fixed_candidates){
  final_data[[col]] <- as.numeric(scale(final_data[[col]]))
}

# =========================
# Função FITNESS (GLMM)
# =========================
fitness_glmm_Gaussiana_rmse <- function(cromossomo){
  
  # Seleciona variáveis baseada no cromossomo binário
  vars <- fixed_candidates[cromossomo == 1]
  
  # Simb_dom entra forçado no modelo
  vars_final <- unique(c("Simb_dom", vars))
  
  # Cria a fórmula
  form_str <- paste(
    "GY ~", paste(vars_final, collapse = " + "),
    "+ (1 |", random_effect, ")"
  )
  form <- as.formula(form_str)
  
  fit <- tryCatch(
    glmer(
      form,
      data = final_data,
      family = gaussian(link = "log"),
      nAGQ = 0, # Mais rápido para GA
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
    ),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  
  if(is.null(fit)) return(-1e6)
  
  y_hat <- predict(fit, type = "response", re.form = NULL) # re.form=NULL usa efeitos aleatórios
  
  if(any(!is.finite(y_hat))) return(-1e6)
  
  rmse <- sqrt(mean((final_data$GY - y_hat)^2))
  
  return(-rmse) 
}

# Teste rápido antes de rodar
# fitness_glmm_Gaussiana_rmse(rep(1, length(fixed_candidates)))

# =========================
# 8. Rodar GA (Paralelo) 
# =========================

# 1. Configurar Cluster
n_cores <- detectCores(logical = FALSE) - 1
cl <- makeCluster(ifelse(n_cores < 1, 1, n_cores))
registerDoParallel(cl)


# Exporta os objetos de dados que a função fitness usa
clusterExport(cl, varlist = c("final_data", "fixed_candidates", "random_effect"))

# Carrega o pacote lme4 em cada nó do cluster (para o glmer funcionar)
clusterEvalQ(cl, library(lme4))

# ---------------------------------------------------------

set.seed(2025)

resultado_GA <- ga(
  type = "binary",
  fitness = fitness_glmm_Gaussiana_rmse,
  nBits = length(fixed_candidates),
  names = fixed_candidates, 
  popSize = 60,
  maxiter = 50,
  run = 15,
  pmutation = 0.1,
  elitism = 2,
  parallel = cl  # Agora o cluster já tem os dados carregados
)

stopCluster(cl) # Finaliza o cluster
# Extração da melhor solução (tratando caso haja empate)
solution_vector <- resultado_GA@solution
if(is.matrix(solution_vector)) solution_vector <- solution_vector[1, ]

variaveis_selecionadas <- fixed_candidates[solution_vector == 1]
message("Variáveis selecionadas pelo GA: ", paste(variaveis_selecionadas, collapse = ", "))

# =========================
# 9. Função de ajuste final
# =========================
ajustar_glmm <- function(vars_fixas){
  
  form_str <- paste(
    "GY ~ Simb_dom",
    if(length(vars_fixas) > 0) paste("+", paste(vars_fixas, collapse = " + ")),
    "+ (1 |", random_effect, ")"
  )
  
  glmer(
    as.formula(form_str),
    data = final_data,
    family = gaussian(link = "log"),
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  )
}

# =========================
# 10. Backward 
# =========================
alpha <- 0.05
vars_atuais <- variaveis_selecionadas

# Nota: Uso do GLM simples para backward rápido das fixas para captar a estrutura de efeitos fixos, 
# GLM para velocidade.
repeat {
  
  # Fórmula atual
  form_bk <- as.formula(paste("GY ~", paste(c("Simb_dom", vars_atuais), collapse = " + ")))
  
  modelo_glm <- glm(
    form_bk,
    data = final_data,
    family = gaussian(link = "log")
  )
  
  coefs <- summary(modelo_glm)$coefficients
  
  # Filtra apenas as variáveis candidatas (ignora intercepto e Simb_dom fixo)
  pvals <- coefs[rownames(coefs) %in% vars_atuais, "Pr(>|t|)"]
  
  if(length(pvals) == 0) break
  
  max_p_val <- max(pvals, na.rm = TRUE)
  
  if(max_p_val <= alpha) break
  
  var_remover <- names(which.max(pvals))
  
  message("Removendo: ", var_remover, " | p = ", round(max_p_val, 4))
  vars_atuais <- setdiff(vars_atuais, var_remover)
}

message("Variáveis Finais: ", paste(vars_atuais, collapse = ", "))

# =========================
# 11. Modelo FINAL e Diagnóstico
# =========================

modelo_final <- ajustar_glmm(vars_atuais)
print(summary(modelo_final))

rmse_final <- sqrt(mean((final_data$GY - predict(modelo_final, type = "response"))^2))
message("RMSE Final: ", rmse_final)

# DHARMa
sim <- simulateResiduals(modelo_final)
plot(sim)
