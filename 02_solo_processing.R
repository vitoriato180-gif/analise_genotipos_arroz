
# dados de solo com o script

rm(list=ls())

# Bibliotecas Utilizadas
library(terra)
library(geobr)
library(ggplot2)
library(sf)
library(readxl)
library(purrr)
library(dplyr)

# Inserir o Shape File do AD_Brasil
#shp <- vect("C:\\Users\\gabri\\Downloads\\AD_Brasil_Tabela_Reduzida\\AD_Brasil.shp")
shp <- vect("C:\\Users\\elisv\\Downloads\\ic_embrapa\\g_dados_de_solo\\AD_Brasil_Tabela_Reduzida\\AD_Brasil.shp")

class(shp)
crs(shp)
# Mapa do Brasil inteiro com os dados do SHP (muito pesado por sinal.)
#shp

# Nomes das variáveis
#names(shp)
#head(as.data.frame(shp))

# Pegar os Shape File dos estados do Brasil
estados <- geobr::read_state(year = 2020, simplified = TRUE)

# Colocar os estados que seus dados estão
#tocantins <- estados[estados$name_state %in% c("Tocantins", "Pará", "Goiás"), ]
tocantins <- estados[estados$name_state %in% c("Mato Grosso", "Goiás", "Pará", "Rondônia",
                                               "Piauí", "Tocantins", "Maranhão"), ]

estados <- vect(estados)  # converte sf → SpatVector
tocantins <- vect(tocantins)

# Cruzar o shape file do Brasil com o da AD_Brasil para ficar mais leve o tratamento
shp_tocantins <- terra::intersect(shp, tocantins)
#saveRDS(shp_tocantins, "shp_tocantins2_inicial_RDS")

# shp_tocantins <- shp[shp$NM_ESTADO == "Tocantins"]

# Se vier lista, junta
if (inherits(shp_tocantins, "list")) {
  shp_tocantins <- do.call(rbind, shp_tocantins)
}

class(shp_tocantins)
class(shp_tocantins_sf)


# Salvar shapefile
shp_tocantins_sf <- st_as_sf(shp_tocantins)
# st_write(shp_tocantins_sf, "shp_tocantins.shp", delete_layer = TRUE)


st_crs(shp)
st_crs(tocantins)

# Plotar para conferir se o curzamento deu certo
plot(tocantins, border = "red", lwd = 2)
plot(shp_tocantins, add = TRUE, col = "lightblue")

shp_tocantins

# Converter para data.frame compatível com ggplot2
df_tocantins <- st_as_sf(shp_tocantins) 

# saveRDS(df_tocantins, "df_tocantins_RDS")


# Plotando com ggplot
ggplot(df_tocantins) +
  geom_sf(aes(fill = C1_Simb), color = NA) +
  labs(title = "Mapa do Tocantins colorido por C1_Simb") +
  theme_minimal()

head(df_tocantins)

# Inserir os dados que está utilizando na modelagem
dados <- read.csv("C:\\Users\\elisv\\Downloads\\novo\\climate_data_dezembro.csv")

head(dados)

# Transformando a LATITUDE e LONGITUDE em uma coordenada geográfica
dados_sf <- st_as_sf(
  dados,
  coords = c("LONGITUDE", "LATITUDE"), # define colunas
  crs = 4674 # ou 2000
)

# Converter os polígonos para 2D
df_tocantins <- st_transform(df_tocantins, 4674)

# Converter também os pontos (se ainda não tiver feito)
dados_sf <- st_transform(dados_sf, 4674)

# Conferir se estão iguais
st_crs(df_tocantins)
st_crs(dados_sf)

# Nomes das variáveis que estão disponíveis
names(df_tocantins)

# Agora fazer o join com as variáveis que deseja inserir no modelo
dados_join <- st_join(dados_sf, df_tocantins)

intersec <- st_intersects(dados_sf, df_tocantins)
# dados_sf$Classe_Solo <- map_chr(intersec, ~ if(length(.x) > 0) df_tocantins$C1_Ord[.x[1]] else NA )

dados_sf$Classe_Solo <- map_chr(
  intersec,
  ~ if (length(.x) > 0) df_tocantins$C1_Ord[.x[1]] else NA
)


# -----------------------------
# 4. Separar pontos problemáticos e não problemáticos
# -----------------------------
dados_edif <- dados_sf %>% filter(Classe_Solo == "Área Edificada")
dados_ok   <- dados_sf %>% filter(Classe_Solo != "Área Edificada" | is.na(Classe_Solo))


# Função robusta para perturbar um ponto dentro do mesmo polígono (evitando "Área Edificada")
perturbar_coord <- function(ponto, polygons, classe_col = "C1_Ord",
                            classe_avoid = "Área Edificada",
                            max_sample_attempts = 200, max_jitter_attempts = 200) {
  # ponto: sfc_POINT (ou sf row geometry)
  # polygons: sf polygons com coluna classe_col (ex: "C1_Ord")
  # retorna: sfc_POINT (novo) ou o ponto original em caso de falha
  
  # garantir que são sfc com mesmo CRS
  stopifnot(inherits(ponto, c("sfc_POINT", "sfc")))
  stopifnot(inherits(polygons, "sf"))
  
  # ver qual polígono contém o ponto (usar primeiro se houver vários)
  hits <- st_intersects(ponto, polygons)[[1]]
  if (length(hits) == 0) {
    # ponto não está em nenhum polígono -> não perturbamos
    return(ponto)
  }
  muni_poly <- polygons[hits[1], , drop = FALSE]
  
  # 1) Tentar amostrar pontos aleatórios *dentro do mesmo polígono* (st_sample)
  # st_sample devolve sfc; se poly inválido pode dar erro, então tryCatch
  for (i in seq_len(max_sample_attempts)) {
    p_new <- tryCatch({
      st_sample(muni_poly, size = 1)        # sfc_POINT
    }, error = function(e) NULL)
    if (is.null(p_new) || length(p_new) == 0) next
    
    # verificar classe do polígono onde o ponto caiu (pode cair em sub-polígono)
    hits2 <- st_intersects(p_new, polygons)[[1]]
    if (length(hits2) == 0) next
    classe <- polygons[[classe_col]][hits2[1]]
    if (!is.na(classe) && classe != classe_avoid) {
      # garantir CRS igual ao ponto e retornar como sfc_POINT
      return(st_sfc(st_geometry(p_new)[[1]], crs = st_crs(ponto)))
    }
  }
  
  # 2) Se falhar, tentar jitter ao redor do ponto original até ficar dentro do polígono e não ser Área Edificada
  coords <- st_coordinates(ponto)
  lon <- coords[1]; lat <- coords[2]
  
  for (i in seq_len(max_jitter_attempts)) {
    lon_n <- lon + runif(1, -0.01, 0.01)
    lat_n <- lat + runif(1, -0.01, 0.01)
    p_new <- st_sfc(st_point(c(lon_n, lat_n)), crs = st_crs(ponto))
    inside <- st_intersects(p_new, muni_poly)[[1]]
    if (length(inside) > 0) {
      classe <- muni_poly[[classe_col]][inside[1]]
      if (!is.na(classe) && classe != classe_avoid) {
        return(p_new)
      }
    }
  }
  
  # 3) fallback: não conseguiu, retorna o ponto original
  return(ponto)
}
##########################



# Supondo que você tenha uma função que gera ruído pequeno
perturbar_coord <- function(x, y, max_offset = 0.001) {
  x + runif(length(x), -max_offset, max_offset)
}

dados_ajustados <- dados_edif %>%
  mutate(
    lon = st_coordinates(geometry)[, 1],
    lat = st_coordinates(geometry)[, 2],
    lon_pert = perturbar_coord(lon),
    lat_pert = perturbar_coord(lat),
    geometry = st_as_sf(data.frame(lon_pert, lat_pert), coords = c("lon_pert", "lat_pert"), crs = st_crs(dados_edif))$geometry
  ) %>%
  st_as_sf()







# -----------------------------
# 6. Aplicar a perturbação
# -----------------------------
dados_ajustados <- dados_edif %>%
  mutate(
    geometry = st_sfc(
      purrr::map(geometry, ~ st_geometry(perturbar_coord(.x, df_tocantins))[[1]]),
      crs = st_crs(dados_edif)
    )
  )


########################
dados_ajustados <- dados_edif %>%
  mutate(
    geometry = st_sfc(
      purrr::map(geometry, ~ {
        ponto_sf <- st_sfc(.x, crs = st_crs(dados_edif))
        st_geometry(perturbar_coord(ponto_sf, df_tocantins))[[1]]
      }),
      crs = st_crs(dados_edif)
    )
  )
########################
# -----------------------------
# 7. Reunir todos os pontos
# -----------------------------
dados_final <- bind_rows(dados_ok, dados_ajustados)
dados_final <- st_as_sf(dados_final)

# Após gerar 'dados_final' com os pontos ajustados
dados_final <- st_as_sf(dados_sf)

# Fazer join espacial para adicionar todas as colunas do shapefile
dados_final_completo <- st_join(dados_final, df_tocantins, join = st_intersects, left = TRUE)

# Conferir as colunas
names(dados_final_completo)

# -----------------------------
# 8. Conferir quantos pontos ainda caem em Área Edificada
# -----------------------------
inter_final <- st_intersects(dados_final, df_tocantins)
classe_final <- map_chr(inter_final, ~ if(length(.x) > 0) df_tocantins$C1_Ord[.x[1]] else NA)
cat("Número de pontos ainda em Área Edificada:", sum(classe_final == "Área Edificada", na.rm = TRUE), "\n")



# -----------------------------
# 9. Salvar shapefile
# -----------------------------
st_write(dados_final_completo, "dados_final_corrigido_completo.shp", delete_layer = TRUE)

# Conferência de variáveis específicas para o modelo
#unique(dados_final_completo$C1_Simb)
#unique(df_tocantins$C1_Ord)



# Exportar a base de dados

writexl::write_xlsx(dados_final_completo, "C:\\Users\\elisv\\Downloads\\novo\\20251006_dados_clima_solo.xlsx")


dados_final_completo <- dados_final_completo %>%
  mutate(id = dplyr::row_number())




library(stringr)


df_aux <- dados_final_completo %>%
  mutate(
    id = dplyr::row_number(),
    Simb = pmap_chr(dplyr::select(., matches("^C[1-5]_Simb$")),
                    ~ stringr::str_c(na.omit(c(...)), collapse = " + "))
  ) %>%
  dplyr::select(-matches("^C[1-5]_(Simb)$"))


# ----------------------------
# 3. Filtrar observações válidas
# ----------------------------

df_kmeans <- df_aux %>%
  filter(!is.na(AD_UM))

set.seed(123)
kmeans_result <- kmeans(df_kmeans$AD_UM, centers = 3, nstart = 25) ## AD_UM foi categorizada em três níveis (Baixa, Média, Alta)
df_kmeans$AD_UM_cluster <- kmeans_result$cluster

cluster_means <- df_kmeans %>%
  group_by(AD_UM_cluster) %>%
  summarise(media_AD_UM = mean(AD_UM, na.rm = TRUE)) %>%
  arrange(media_AD_UM)

nomes_categorias <- c("Baixa", "Média", "Alta")
nomes_ordenados <- setNames(nomes_categorias, cluster_means$AD_UM_cluster)
df_kmeans$AD_UM_cat <- nomes_ordenados[as.character(df_kmeans$AD_UM_cluster)]

# Agora o join funciona perfeitamente
df <- df_aux %>%
  left_join(df_kmeans %>% st_drop_geometry() %>% select(id, AD_UM_cat), by = "id")

# ----------------------------
# 9. Visualização — Boxplot
# ----------------------------
ggplot(df %>% filter(!is.na(AD_UM_cat)),
       aes(x = AD_UM_cat, y = AD_UM, fill = AD_UM_cat)) +
  geom_boxplot(alpha = 0.8, outlier.color = "black") +
  scale_fill_manual(values = c("Baixa" = "#66c2a5",
                               "Média" = "#fc8d62",
                               "Alta"  = "#8da0cb")) +
  labs(
    title = "Distribuição de AD_UM por Categoria (K-Means)",
    x = "Categoria de AD_UM",
    y = "AD_UM"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

colnames(df)

df <- df %>%
  dplyr::mutate(
    Simb = str_replace_all(
      Simb,
      c(
        "FFc"  = "Plintossolo Pétrico Concressionário",
        "RQo"  = "Neossolo Quartzarenico",
        "CXbd" = "Cambissolo",
        "PVAd" = "Argissolo Vermelho Amarelo",
        "PAd"  = "Argissolo Amarelo",
        "LVAd" = "Latossolo Vermelho Amarelo",
        "LVd"  = "Latossolo Vermelho",
        "LVw"  = "Latossolo Vermelho",
        "RYbd" = "Neossolo Flúvico",
        "RLd"  = "Neossolo Litólico"
      )
    )
  )

colnames(df_aux)

# Vetor com os nomes das colunas a remover
cols_remover <- c(
  "Classe_Solo", "OBJECTID", "COD", "Linha",
  # Bloco C1
  "C1_Leg", "C1_Class", "C1_Ord", "C1_Subord", "C1_Ggrup", "C1_Sgrup",
  "C1_Text1", "C1_Text2", "C1_Text3", "C1_Text4", "C1_Horiz",
  "C1_Erosao", "C1_Pedr", "C1_Roch", "C1_Relevo",
  "C1_ADtotal", "C1_Psolo", "C1_ADProp",
  # Bloco C2
  "C2_Leg", "C2_Class", "C2_Ord", "C2_Subord", "C2_Ggrup", "C2_Sgrup",
  "C2_Text1", "C2_Text2", "C2_Text3", "C2_Text4", "C2_Horiz",
  "C2_Erosao", "C2_Pedr", "C2_Roch", "C2_Relevo",
  "C2_ADtotal", "C2_Psolo", "C2_ADProp",
  # Bloco C3
  "C3_Leg", "C3_Class", "C3_Ord", "C3_Subord", "C3_Ggrup", "C3_Sgrup",
  "C3_Text1", "C3_Text2", "C3_Text3", "C3_Text4", "C3_Horiz",
  "C3_Erosao", "C3_Pedr", "C3_Roch", "C3_Relevo",
  "C3_ADtotal", "C3_Psolo", "C3_ADProp",
  # Bloco C4
  "C4_Leg", "C4_Class", "C4_Ord", "C4_Subord", "C4_Ggrup", "C4_Sgrup",
  "C4_Text1", "C4_Text2", "C4_Text3", "C4_Text4", "C4_Horiz",
  "C4_Erosao", "C4_Pedr", "C4_Roch", "C4_Relevo",
  "C4_ADtotal", "C4_Psolo", "C4_ADProp",
  # Bloco C5
  "C5_Leg", "C5_Class", "C5_Ord", "C5_Subord", "C5_Ggrup", "C5_Sgrup",
  "C5_Text1", "C5_Text2", "C5_Text3", "C5_Text4", "C5_Horiz",
  "C5_Erosao", "C5_Pedr", "C5_Roch", "C5_Relevo",
  "C5_ADtotal", "C5_Psolo", "C5_ADProp",
  # Outros
  "A_Sp_km2", "A_Mp_km2"
)

# Remover as colunas
df <- df %>% dplyr::select(-any_of(cols_remover))

colnames(df)

unique(df$Simb)

df <- df %>%
  mutate(
    Simb = sapply(Simb, function(x) {
      solos <- str_split(x, "\\s*\\+\\s*")[[1]]      # separa pelos "+"
      solos_unicos <- unique(solos)                  # remove duplicatas
      paste(solos_unicos, collapse = " + ")          # junta de volta
    })
  )
unique(df$Simb)

colnames(df)


df_export <- df %>%
  mutate(
    lon = st_coordinates(geometry)[,1],
    lat = st_coordinates(geometry)[,2]
  ) %>%
  st_drop_geometry()

writexl::write_xlsx(df_export,
                    "C:\\Users\\elisv\\Downloads\\novo\\dados_clima_solo_v03_dezembro.xlsx")



