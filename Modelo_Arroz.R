rm(list = ls())
gc()

# ============================================================
# REPRODUTIBILIDADE TOTAL E CONFIGURAÇÕES
# ============================================================
set.seed(2025)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

Sys.setenv(
  OMP_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  VECLIB_MAXIMUM_THREADS = 1,
  NUMEXPR_NUM_THREADS = 1
)
options(mc.cores = 1)

# ============================================================
# 1. PACOTES
# ============================================================
packages <- c("lme4", "lmerTest", "GA", "dplyr", "readxl", "DHARMa", 
              "stringr", "MuMIn", "tidyr", "caret", "ggplot2", "scales")

#install.packages(packages)

for(pkg in packages){
  if(!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

set.seed(2025)

# ============================================================
# 2. DADOS E PREPARAÇÃO
# ============================================================

# Carregando a base de dados
message("Por favor, selecione o arquivo de dados 'dados_dezembro.xlsx'")
final_data <- read_excel(file.choose())

# Vetor com as bases climáticas desejadas
vars_keep <- c(
  "T2M", "T2M_MAX", "T2M_MIN", "PRECTOTCORR",
  "GWETROOT", "RH2M", "ALLSKY_SFC_PAR_TOT",
  "EVPTRNS", "WS2M", "CDD0"
)

# Filtros
clim_keep <- grepl(paste(vars_keep, collapse = "|"), names(final_data))
clim_all <- grepl("^(veg_|repro_|gf_)", names(final_data))
final_data <- final_data[, !clim_all | clim_keep]

final_data <- final_data %>%
  group_by(TRIAL, GEN) %>%
  filter(n() == 4) %>%
  ungroup() %>%
  mutate(across(where(is.character), as.factor)) %>%
  group_by(TRIAL, GEN) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    across(where(is.factor), ~ first(.x)),
    .groups = "drop"
  ) %>%
  filter(!is.na(GY), GY > 0) %>% 
  filter(!is.na(AD_UM))



colnames(final_data)

remover <- c(
  "TRIAL","DATE","ST","REGIAO","TYPE","DESIGN",
  "LATITUDE","LONGITUDE","Ano_Semeadura",
  "PI_dias","PM_dias","DTF","REP","N_REP",
  "LOCATION","Simb","AD_UM_cat","Classe_AD"
)

final_data <- final_data %>%
  mutate(LOC = paste(LOC, ST, sep = "_")) %>%
  select(-any_of(remover))

final_data$LOC <- factor(final_data$LOC)
final_data$GEN <- factor(final_data$GEN)

# ============================================================
# 2.5 REMOÇÃO DE VARIÂNCIA ZERO
# ============================================================
vars_numericas <- names(final_data)[sapply(final_data, is.numeric)]
vars_numericas <- setdiff(vars_numericas, c("GY"))

variancias <- sapply(final_data[vars_numericas], var, na.rm = TRUE)
vars_zero_var <- names(variancias[variancias == 0 | is.na(variancias)])

if(length(vars_zero_var) > 0) {
  message("Removendo variáveis com variância zero: ", paste(vars_zero_var, collapse = ", "))
  final_data <- final_data %>% select(-all_of(vars_zero_var))
}

# ============================================================
# 3. O MODELO CHEIO 
# ============================================================
todas_climaticas <- names(final_data)[sapply(final_data, is.numeric)]

todas_climaticas <- setdiff(todas_climaticas, c("GY", "AD_UM"))

formula_cheia <- paste(
  "GY ~ LOC + AD_UM +", 
  paste(todas_climaticas, collapse = " + "), 
  "+ (1|GEN) + (1|GEN:LOC)"
)

modelo_cheio <- tryCatch({
  lmer(as.formula(formula_cheia), data = final_data, REML = TRUE, 
       control = lmerControl(optimizer = "bobyqa"))
}, error = function(e) {
  message("O modelo cheio quebrou: ", e)
  return(NULL)
})

if(!is.null(modelo_cheio)){
  print(summary(modelo_cheio))
  message("Modelo cheio ajustado sem duplicidade de AD_UM. Pressione Enter para ir para o GA...")
  readline()
}


summary(modelo_cheio)
capture.output(
  print(summary(modelo_cheio)),
  file = "/Users/mariavitoriadiastorres/Downloads/summary_modelo_cheio_test.txt"
)

res2 <- simulateResiduals(modelo_cheio)
plot(res2)
# ============================================================
# 4. FUNÇÃO FITNESS PARA O GA 
# ============================================================
set.seed(1234)

fixed_candidates <- final_data %>% select(-GY, -GEN, -LOC, -AD_UM)

fitness_ga <- function(cromossomo){
  
  vars_selecionadas <- names(fixed_candidates)[cromossomo == 1]
  num_vars <- length(vars_selecionadas)
  
  if(num_vars == 0){
    form_str <- "GY ~ LOC + AD_UM + (1|GEN) + (1|GEN:LOC)"
  } else {
    form_str <- paste("GY ~ LOC + AD_UM +", paste(vars_selecionadas, collapse = " + "), "+ (1|GEN) + (1|GEN:LOC)")
  }
  
  fit <- tryCatch(
    lmer(as.formula(form_str), data = final_data, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )
  
  if(is.null(fit)) return(-1e6)
  
  # GA maximiza, então retornamos -RMSE para minimizar o erro
  rmse <- sqrt(mean(residuals(fit)^2))
  return(-rmse)
}

# ============================================================
# 5. EXECUÇÃO DO GA
# ============================================================
ga_res <- ga(
  type = "binary",
  fitness = fitness_ga,
  nBits = ncol(fixed_candidates),
  popSize = 50,
  maxiter = 50,
  run = 15,
  seed = 2025,
  parallel = FALSE,
  monitor = TRUE
)

sol <- ga_res@solution
if(is.matrix(sol)) sol <- sol[1,]

vars_ga <- names(fixed_candidates)[sol == 1]
print("Variáveis selecionadas naturalmente pelo GA:")
print(vars_ga)

# ============================================================
# 6. MODELO BASE DO GA
# ============================================================
if(length(vars_ga) == 0){
  form_final <- "GY ~ LOC + AD_UM + (1|GEN) + (1|GEN:LOC)"
} else {
  form_final <- paste("GY ~ LOC + AD_UM +", paste(vars_ga, collapse = " + "), "+ (1|GEN) + (1|GEN:LOC)")
}

mod_ga <- lmer(
  as.formula(form_final),
  data = final_data,
  REML = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)
# está dando erro aqui

summary(mod_ga)
capture.output(
  print(summary(mod_ga)),
  file = "/Users/mariavitoriadiastorres/Downloads/summary_mod_ga.txt"
)

res3 <- simulateResiduals(mod_ga)
plot(res3)

# ============================================================
# 7. FUNÇÕES DE LIMPEZA BACKWARD 
# ============================================================
limpar_niveis_pos_backward <- function(modelo, data, fator, alpha = 0.05, metodo = "bonferroni"){
  repeat {
    summ <- summary(modelo)$coefficients
    termos <- rownames(summ)[grepl(paste0("^", fator), rownames(summ))]
    if(length(termos) == 0) break
    col_p <- grep("Pr\\(>\\|", colnames(summ), value = TRUE)
    if(length(col_p) == 0) break
    
    p_vals <- summ[termos, col_p]
    p_adj  <- p.adjust(p_vals, method = metodo)
    nao_sig <- names(p_adj[p_adj > alpha])
    
    if(length(nao_sig) == 0) break
    
    niveis <- sub(paste0(fator), "", nao_sig)
    ref <- levels(data[[fator]])[1]
    
    data[[fator]] <- as.character(data[[fator]])
    data[[fator]][data[[fator]] %in% niveis] <- ref
    data[[fator]] <- droplevels(factor(data[[fator]]))
    
    f_atual <- formula(modelo)
    if(nlevels(data[[fator]]) < 2){
      f_atual <- update(f_atual, as.formula(paste(". ~ . -", fator)))
      modelo <- lmer(f_atual, data = data, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
      break
    }
    modelo <- lmer(f_atual, data = data, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
  }
  return(list(modelo = modelo, data = data))
}

limpar_modelo_misto <- function(modelo_inicial, data, alpha = 0.05) {
  modelo_atual <- modelo_inicial
  fixas_obrigatorias <- c("LOC", "AD_UM", "(Intercept)") 
  
  repeat {
    summ <- summary(modelo_atual)$coefficients
    vars_analise <- rownames(summ)
    vars_analise <- vars_analise[!sapply(vars_analise, function(x) any(sapply(fixas_obrigatorias, grepl, x)))]
    
    if(length(vars_analise) == 0) break
    
    p_vals <- summ[vars_analise, "Pr(>|t|)"]
    max_p  <- max(p_vals)
    
    if(max_p > alpha) {
      var_to_remove <- names(which.max(p_vals))
      nova_form <- update(formula(modelo_atual), as.formula(paste(". ~ . -", var_to_remove)))
      modelo_atual <- lmer(nova_form, data = data, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
    } else {
      break
    }
  }
  return(modelo_atual)
}

# ============================================================
# 8. APLICAÇÃO DA LIMPEZA NO MODELO GA
# ============================================================
dados_c1 <- final_data   

# Agrupa níveis não significativos de LOC
res_1 <- limpar_niveis_pos_backward(mod_ga, dados_c1, "LOC", 0.05)
mod_ga_tmp <- res_1$modelo
dados_c1   <- res_1$data

nivel_referencia <- levels(dados_c1$LOC)[1]
levels(dados_c1$LOC)[levels(dados_c1$LOC) == nivel_referencia] <- "OTHERS"

form_atual <- formula(mod_ga_tmp)
mod_ga_tmp <- lmer(form_atual, data = dados_c1, REML = TRUE, 
                   control = lmerControl(optimizer = "bobyqa"))

# Limpa variáveis numéricas climáticas
mod_misto <- limpar_modelo_misto(mod_ga_tmp, dados_c1, 0.05)

summary(mod_misto)
ranef(mod_misto)

# ============================================================
# 9. AVALIAÇÃO DE RESÍDUOS
# ============================================================
res <- simulateResiduals(mod_misto)
plot(res)
