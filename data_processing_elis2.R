
data_processing_elis2 <- function(data) {
  library(dplyr)
  data1 <- data %>%
    mutate(
      DATE = as.Date(DATE),
      Ano_Semeadura = lubridate::year(DATE)
    ) %>%
    # Filtra últimos 10 anos e apenas VCU
    filter(Ano_Semeadura > (max(Ano_Semeadura, na.rm = TRUE) - 10),
           TYPE == "VCU") %>%
    # Cria colunas de PI, PM e DTF em dias e datas
    mutate(
      PI_dias = DTF - 30,
      PM_dias = DTF + 30,
      PI_data = DATE + PI_dias,
      PM_data = DATE + PM_dias,
      DTF_data = DATE + DTF,
      .after = DATE
    ) %>%
    # Remove estados indesejados e NAs
    filter(!is.na(ST),
           !ST %in% c("AC", "SP", "AM","RR"),
           !is.na(DTF)) %>%
    # Remove genótipos que começam com H1, H5 ou AB
    filter(!grepl("^(H1|H5|AB)", GEN)) %>%
    # Remove genótipos desconhecidos
    filter(!GEN %in% c("UNKNOWN_01", "UNKNOWN")) %>%
    # Adiciona a coluna de Região
    mutate(
      REGIAO = case_when(
        ST %in% c("GO", "MT", "MS", "DF") ~ "Centro-Oeste",
        ST %in% c("MA", "PI", "CE", "RN", "PB", "PE", "AL", "SE", "BA") ~ "Nordeste",
        ST %in% c("AC", "AM", "RO", "RR", "PA", "AP", "TO") ~ "Norte",
        ST %in% c("SP", "RJ", "ES", "MG") ~ "Sudeste",
        ST %in% c("PR", "SC", "RS") ~ "Sul",
        TRUE ~ NA_character_
      )
    )
  
  return(data1)
}

# 
# 
# dados_arroz <- read.csv("C:\\Users\\elisv\\Downloads\\novo\\dados_arroz.csv")
# 
# 
# dados_arroz <- data_processing_elis2(dados_arroz)
# 
# 
# unique(dados_arroz$ST)
# unique(dados_arroz$GEN)
# unique(dados_arroz$LOCATION)
# unique(dados_arroz$TYPE)
# unique(dados_arroz$TRIAL)
# 
# table(dados_arroz$ST)
# table(dados_arroz$GEN)
# table(dados_arroz$LOCATION)
# table(dados_arroz$TYPE)
# table(dados_arroz$TRIAL)
# 
# 
# 

































