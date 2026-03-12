
# scritp para fazer a junção dos dados de solo com dados de arroz:

rm(list=ls()) 

library(dplyr)

# lendo os dados
dados_arroz_inicial <- read.csv("C:/Users/elisv/Downloads/novo/dados_arroz.csv")

source("C:\\Users\\elisv\\Downloads\\novo\\data_processing_elis2.R")

dados_arroz_inicial <- data_processing_elis2(dados_arroz_inicial)

dados_modelos <- readxl::read_excel("C:\\Users\\elisv\\Downloads\\novo\\dados_clima_solo_v03_dezembro.xlsx")


climate_results <- dados_arroz_inicial 
dados_com_solo <- dados_modelos

dim(climate_results)
dim(dados_com_solo)


# Organizando o formato de data
climate_results <- climate_results %>%
  mutate(
    DATE = as.Date(DATE),      # Converte DATE para formato Date
    DTF_data = as.Date(DTF_data),      # Converte DTF_data para data só pra ter certeza.
    PI_data = as.Date(PI_data),
    PM_data = as.Date(PM_data)   
  )                                    # #pode remover os dias

str(climate_results$DATE)

# Organizando o formato de data
dados_com_solo$DATE <- as.Date(dados_com_solo$DATE) #
dados_com_solo$PI <- as.Date(dados_com_solo$PI) #
dados_com_solo$PM <- as.Date(dados_com_solo$PM)
dados_com_solo$PF <- as.Date(dados_com_solo$PF)

dados_com_solo <- dados_com_solo %>%
  rename(
    PI_data = PI, #mudar o nome das colunas existente
    PM_data = PM,
    DTF_data = PF,
    LATITUDE = lat,
    LONGITUDE = lon)

names(dados_com_solo)

names(climate_results)

# conferindo a estrutura
str(dados_com_solo$DATE)


# juntando os dados de solo e clima
data_completo <- climate_results %>%
  left_join(dados_com_solo,
            by = c("LATITUDE", "LONGITUDE", "DATE", "DTF_data",
                   "PI_data",
                   "PM_data"))
names(data_completo)

vars_solo <- c(
  "AD_UM", "Classe_AD",
  "code_state", "abbrev_state", "name_state",
  "code_region", "name_region", "id",
  "Simb", "AD_UM_cat"
)

data_completo <- data_completo %>%
  relocate(all_of(vars_solo), .after = REGIAO)

names(data_completo)

# FIM

View(dlookr::diagnose(data_completo))


writexl::write_xlsx(data_completo,"C:\\Users\\elisv\\Downloads\\novo\\dados_completos_dezembro_bruto.xlsx")




# organização final dos dados:


dados_completos <- readxl::read_excel("C:\\Users\\elisv\\Downloads\\novo\\dados_completos_dezembro_bruto.xlsx")

# View(dlookr::diagnose(dados_completos))

library(dplyr) 

final_data <- dados_completos


remover_solo <- c("code_state","abbrev_state","name_state","code_region","name_region","id") # variáveis extras.

remover_NA <- c("BLO", "LBL", "LOD", "PBL", "GDS", "LSC", "BSP") # Variáveis removidas por alta proporção de valores ausentes

remover_sem_sentido <- c("SYST", "PLOT", "YEAR")

remover_variaveis_nao_climaticas <- c(
  "X",
  "PI_data",
  "PM_data",
  "DTF_data",
  "MEAN",
  "H2",
  "CV",
  "PHT",
  "codigo_ibge"
)

# variaveis_nao_climaticas <- c(    # # TODAS AS NÃO CLIMÁTICAS = 31
#   "X",
#   "TRIAL",
#   "DATE",
#   "PI_dias",
#   "PM_dias",
#   "PI_data",
#   "PM_data",
#   "DTF_data",
#   "ST",
#   "LOCATION",
#   "LOC",
#   "TYPE",
#   "DESIGN",
#   "N_REP",
#   "MEAN",
#   "H2",
#   "CV",
#   "REP",
#   "GEN",
#   "GY",
#   "PHT",
#   "codigo_ibge",
#   "LONGITUDE",
#   "LATITUDE",
#   "DTF",
#   "Ano_Semeadura",
#   "REGIAO",
#   "AD_UM",
#   "Classe_AD",
#   "Simb",
#   "AD_UM_cat"
# )

final_data <- final_data %>% select(-all_of(remover_solo)) %>% select(-all_of(remover_NA)) %>% 
  select(-all_of(remover_sem_sentido)) %>% select(-all_of(remover_variaveis_nao_climaticas)) %>% 
  relocate( REP, .after = TRIAL) %>% 
  relocate(GEN, .after = TRIAL) %>% 
  relocate(DTF, .after = PI_dias) 

names(final_data)

library(dplyr)
library(stringr)

final_data <- final_data %>% #removendo AEd que falhou ao encontrar terreno
  filter(!str_detect(Simb, "AEd"))  # remove os valores de AD_UM e AD_UM_cat que tinha NA. (junto com PHT(X) e GY)

final_data$Simb <- str_replace( 
  final_data$Simb,
  " \\+ -?\\d+\\.\\d+ \\+ -?\\d+\\.\\d+$",
  ""
)

View(dlookr::diagnose(final_data))


# writexl::write_xlsx(final_data,"C:\\Users\\elisv\\Downloads\\novo\\dados_dezembro.xlsx")

dados <- readxl::read_excel("C:\\Users\\elisv\\Downloads\\novo\\dados_dezembro.xlsx")

# names(dados)


table(dados$Simb)

table(dados$AD_UM)

table(dados$Classe_AD)

table(dados$AD_UM_cat)

unique(final_data$ST)
unique(final_data$GEN)
unique(final_data$LOCATION)
unique(final_data$TYPE)
unique(final_data$TRIAL)

table(final_data$ST)
table(final_data$GEN)
table(final_data$LOCATION)
table(final_data$TYPE)
table(final_data$TRIAL)

