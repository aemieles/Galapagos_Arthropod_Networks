################ Analisis Grok #########
setwd("~/Articles/INFORME ZONAS AGRICOLAS/Analisis25")
# Cargar librerías
library(tidyverse)
library(vegan)
library(lme4)

# Leer los datos (ajusta la ruta)
data <- read.csv("abu_brut.csv", check.names = FALSE)

# 1. Abundancia total
abundancia_total <- data %>%
  group_by(Isla, Cultivo) %>%
  summarise(Abundancia_Total = sum(Abundancia, na.rm = TRUE)) %>%
  ungroup()
print("Abundancia total por isla y cultivo:")
print(abundancia_total)

# 2. Diversidad Shannon
# Matriz especie-por-muestra (Isla-Cultivo-Trampa como muestras)
matriz <- data %>%
  group_by(Isla, Cultivo, Trampa, Especie) %>%
  summarise(Abundancia = sum(Abundancia, na.rm = TRUE)) %>%
  pivot_wider(names_from = Especie, values_from = Abundancia, values_fill = 0) %>%
  unite("Muestra", Isla, Cultivo, Trampa, sep = "_")
matriz_sp <- matriz %>% select(-Muestra)
shannon <- diversity(matriz_sp, index = "shannon")
diversidad <- data.frame(Muestra = matriz$Muestra, Shannon = shannon)
print("Diversidad Shannon por muestra (Isla-Cultivo-Trampa):")
print(diversidad)

# 3. Especies clave
# Más abundantes (top 5 global)
top_abundantes <- data %>%
  group_by(Especie) %>%
  summarise(Total = sum(Abundancia, na.rm = TRUE)) %>%
  arrange(desc(Total)) %>%
  slice(1:5)
print("Top 5 especies más abundantes:")
print(top_abundantes)

# Más comunes (frecuencia por isla y cultivo)
frecuencia <- data %>%
  group_by(Isla, Cultivo, Especie) %>%
  summarise(Presencia = sum(Abundancia > 0, na.rm = TRUE)) %>%
  group_by(Especie) %>%
  summarise(Frecuencia = sum(Presencia > 0)) %>%
  arrange(desc(Frecuencia)) %>%
  slice(1:5)
print("Top 5 especies más comunes (en más combinaciones Isla-Cultivo):")
print(frecuencia)

# PERMANOVA para composición
permanova <- adonis2(matriz_sp ~ Isla * Cultivo, data = matriz %>% separate(Muestra, c("Isla", "Cultivo", "Trampa"), sep = "_"), 
                     method = "bray")
print("PERMANOVA para composición:")
print(permanova)

# Guardar resultados
write.csv(abundancia_total, "abundancia_total.csv")
write.csv(diversidad, "diversidad_shannon.csv")

########################################

# Crear matriz especie-por-muestra
matriz <- data %>%
  group_by(Isla, Cultivo, Trampa, Especie) %>%
  summarise(Abundancia = sum(Abundancia, na.rm = TRUE)) %>%
  pivot_wider(names_from = Especie, values_from = Abundancia, values_fill = 0) %>%
  unite("Muestra", Isla, Cultivo, Trampa, sep = "_")

# Separar matriz
matriz_sp <- matriz %>% select(-Muestra)
meta <- matriz %>% select(Muestra) %>% separate(Muestra, c("Isla", "Cultivo", "Trampa"), sep = "_")

# Diagnóstico
print("Dimensiones de la matriz original (filas = muestras, columnas = especies):")
dim(matriz_sp)
print("Primeras 5 filas y 5 columnas de la matriz original:")
print(matriz_sp[1:5, 1:5])
print("Número de muestras con solo ceros en la matriz original:")
sum(rowSums(matriz_sp) == 0)

# Filtrado (frecuencia >= 2, abundancia total >= 3)
frecuencia <- colSums(matriz_sp > 0)
abundancia_total <- colSums(matriz_sp)
especies_a_mantener <- names(frecuencia)[frecuencia >= 2 & abundancia_total >= 3]
matriz_sp_filtrada <- matriz_sp[, especies_a_mantener, drop = FALSE]

# Diagnóstico después del filtrado
print("Dimensiones de la matriz filtrada:")
dim(matriz_sp_filtrada)
print("Primeras 5 filas y 5 columnas de la matriz filtrada:")
print(matriz_sp_filtrada[1:5, 1:5])
print("Número de muestras con solo ceros después del filtrado:")
sum(rowSums(matriz_sp_filtrada) == 0)

################

# Crear matriz especie-por-muestra (solo Isla y Cultivo)
matriz <- data %>%
  group_by(Isla, Cultivo, Especie) %>%
  summarise(Abundancia = sum(Abundancia, na.rm = TRUE)) %>%
  pivot_wider(names_from = Especie, values_from = Abundancia, values_fill = 0) %>%
  unite("Muestra", Isla, Cultivo, sep = "_")

# Separar matriz
matriz_sp <- matriz %>% select(-Muestra)
meta <- matriz %>% select(Muestra) %>% separate(Muestra, c("Isla", "Cultivo"), sep = "_")

# Filtrado (frecuencia >= 2, abundancia total >= 3)
frecuencia <- colSums(matriz_sp > 0)
abundancia_total <- colSums(matriz_sp)
especies_a_mantener <- names(frecuencia)[frecuencia >= 2 & abundancia_total >= 3]
matriz_sp_filtrada <- matriz_sp[, especies_a_mantener, drop = FALSE]

# Transformar abundancias
matriz_sp_trans <- log1p(matriz_sp_filtrada)

# NMDS
set.seed(123)
nmds <- metaMDS(matriz_sp_trans, distance = "bray", k = 2)

# Graficar con colores por Isla
plot(nmds, type = "n", main = "NMDS por Isla-Cultivo (filtrado, log-transformado)")
points(nmds, display = "sites", col = as.factor(meta$Isla), pch = 19, cex = 1.2)
text(nmds, display = "sites", labels = meta$Muestra, cex = 0.7, adj = 1)
legend("topright", legend = levels(as.factor(meta$Isla)), col = 1:length(levels(as.factor(meta$Isla))), pch = 19, title = "Isla")


#Definir colores y símbolos para cada isla (amigable para daltónicos)
colores <- c("Cristobal" = "#0072B2",  # Azul
             "Floreana" = "#D55E00",   # Naranja
             "Isabela" = "#CC79A7")    # Morado
simbolos <- c("Cristobal" = 16,  # Círculo
              "Floreana" = 17,   # Triángulo
              "Isabela" = 15)    # Cuadrado
##Graficar con colores, símbolos y elipses
plot(nmds, type = "n", main = "NMDS por Isla-Cultivo (filtrado, log-transformado)")
points(nmds, display = "sites", col = colores[meta$Isla], pch = simbolos[meta$Isla], cex = 1.2)
text(nmds, display = "sites", labels = meta$Muestra, cex = 0.7, adj = 1)

#Agregar elipses por isla
ord <- ordiellipse(nmds, groups = meta$Isla, display = "sites", kind = "ehull", col = colores, lwd = 2)

#Agregar leyenda
legend("topright", legend = levels(as.factor(meta$Isla)), col = colores, pch = simbolos, title = "Isla", pt.cex = 1.2)

# Diagnóstico adicional: cuántas muestras se graficaron
print("Número de muestras graficadas en el NMDS:")
nrow(nmds$points)

###### Permanova ###
permanova <- adonis2(matriz_sp_trans ~ Isla * Cultivo, data = meta, method = "bray")
print("Resultados de la PERMANOVA:")
print(permanova)

# Guardar la PERMANOVA simplificada en una tabla
permanova_simplificada_df <- as.data.frame(permanova_simplificada$aov.tab)
permanova_df <- as.data.frame(permanova)
write.csv(permanova_df, "permanova_resultados.csv", row.names = TRUE)
permanova_simplificada <- adonis2(matriz_sp_trans ~ Isla + Cultivo, data = meta, method = "bray")
print("Resultados de la PERMANOVA simplificada (Isla + Cultivo):")
print(permanova_simplificada)
permanova_simplificada_df <- as.data.frame(permanova_simplificada)
write.csv(permanova_simplificada_df, "permanova_simplificada_resultados.csv", row.names = TRUE)

########Indices#############
#Calcular el índice de Shannon por Isla y Cultivo
shannon <- diversity(matriz_sp, index = "shannon")
#Crear un data frame con los resultados
diversidad <- data.frame(
  Muestra = matriz$Muestra,
  Isla = meta$Isla,
  Cultivo = meta$Cultivo,
  Shannon = shannon
)
#Mostrar los resultados
print("Índice de Shannon por Isla y Cultivo:")
print(diversidad)
#Guardar los resultados en un CSV
write.csv(diversidad, "diversidad_shannon.csv", row.names = FALSE)
#Resumen de la diversidad por Isla
diversidad_por_isla <- diversidad %>%
  group_by(Isla) %>%
  summarise(Shannon_Media = mean(Shannon), Shannon_SD = sd(Shannon))
print("Diversidad promedio por Isla:")
print(diversidad_por_isla)
#Resumen de la diversidad por Cultivo
diversidad_por_cultivo <- diversidad %>%
  group_by(Cultivo) %>%
  summarise(Shannon_Media = mean(Shannon), Shannon_SD = sd(Shannon))
print("Diversidad promedio por Cultivo:")
print(diversidad_por_cultivo)

#Test de Kruskal-Wallis para comparar diversidad entre Isla
kruskal_isla <- kruskal.test(Shannon ~ Isla, data = diversidad)
print("Test de Kruskal-Wallis para Isla:")
print(kruskal_isla)

#Test de Kruskal-Wallis para comparar diversidad entre Cultivo
kruskal_cultivo <- kruskal.test(Shannon ~ Cultivo, data = diversidad)
print("Test de Kruskal-Wallis para Cultivo:")
print(kruskal_cultivo)

#Guardar los resultados en un archivo de texto
sink("kruskal_diversidad.txt")
print("Test de Kruskal-Wallis para Isla:")
print(kruskal_isla)
print("Test de Kruskal-Wallis para Cultivo:")
print(kruskal_cultivo)
sink()

# Especies más abundantes (top 5 global)
top_abundantes <- data %>%
  group_by(Especie) %>%
  summarise(Total = sum(Abundancia, na.rm = TRUE)) %>%
  arrange(desc(Total)) %>%
  slice(1:5)
print("Top 5 especies más abundantes (global):")
print(top_abundantes)

# Especies más abundantes por Isla
top_abundantes_isla <- data %>%
  group_by(Isla, Especie) %>%
  summarise(Total = sum(Abundancia, na.rm = TRUE)) %>%
  arrange(Isla, desc(Total)) %>%
  group_by(Isla) %>%
  slice(1:5)
print("Top 5 especies más abundantes por Isla:")
print(top_abundantes_isla)

# Especies más comunes (frecuencia por Isla y Cultivo)
frecuencia <- data %>%
  group_by(Isla, Cultivo, Especie) %>%
  summarise(Presencia = sum(Abundancia > 0, na.rm = TRUE)) %>%
  group_by(Especie) %>%
  summarise(Frecuencia = sum(Presencia > 0)) %>%
  arrange(desc(Frecuencia)) %>%
  slice(1:5)
print("Top 5 especies más comunes (en más combinaciones Isla-Cultivo):")
print(frecuencia)

# Guardar los resultados en CSV
write.csv(top_abundantes, "top_abundantes.csv", row.names = FALSE)
write.csv(top_abundantes_isla, "top_abundantes_isla.csv", row.names = FALSE)
write.csv(frecuencia, "especies_comunes.csv", row.names = FALSE)

## Cargar el paquete lme4 (instálalo si no lo tienes: install.packages("lme4"))
library(lme4)

# Preparar datos: sumar abundancia total por Isla, Cultivo, Trampa y Sector
datos_glmm <- data %>%
  group_by(Isla, Cultivo, Trampa, Sector) %>%
  summarise(Abundancia_Total = sum(Abundancia, na.rm = TRUE))

# Ajustar el GLMM
glmm <- glmer(Abundancia_Total ~ Isla + Cultivo + (1 | Trampa) + (1 | Sector),
              data = datos_glmm, family = poisson)

# Mostrar los resultados
print("Resultados del GLMM:")
summary(glmm)

# Guardar los resultados en un archivo de texto
sink("glmm_resultados.txt")
print("Resultados del GLMM:")
summary(glmm)
sink()

# Instalar lme4 si no está instalado
if (!requireNamespace("lme4", quietly = TRUE)) {
  install.packages("lme4")
}

# Inspeccionar los datos
print("Primeras filas de datos_glmm:")
head(datos_glmm)

print("Resumen de datos_glmm:")
summary(datos_glmm)

# Verificar si hay valores NA o problemáticos en Abundancia_Total
print("Número de NA en Abundancia_Total:")
sum(is.na(datos_glmm$Abundancia_Total))

# Verificar los niveles de las variables categóricas
print("Niveles de Isla:")
levels(factor(datos_glmm$Isla))

print("Niveles de Cultivo:")
levels(factor(datos_glmm$Cultivo))

# Verificar si glmer.nb está disponible
if (exists("glmer.nb")) {
  print("La función glmer.nb está disponible.")
} else {
  print("La función glmer.nb NO está disponible. Actualiza lme4 o revisa tu instalación.")
}
install.packages("lme4", dependencies = TRUE)
library(lme4)

# Instalar y cargar MASS si no está instalado
if (!requireNamespace("MASS", quietly = TRUE)) {
  install.packages("MASS")
}
library(MASS)

# Ajustar un GLM simplificado (sin efectos aleatorios, solo Isla)
glm_nb_simple <- glm.nb(Abundancia_Total ~ Isla, data = datos_glmm)

# Mostrar resultados
print("Resultados del GLM simplificado (Negativa Binomial):")
summary(glm_nb_simple)

# Intentar ajustar el modelo y capturar errores/advertencias
result <- try(glmm_nb <- glmer.nb(Abundancia_Total ~ Isla + Cultivo + (1 | Trampa) + (1 | Sector),
                                  data = datos_glmm,
                                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
              silent = FALSE)

# Mostrar el resultado del intento
if (inherits(result, "try-error")) {
  print("Error al ajustar el modelo:")
  print(result)
} else {
  print("Modelo ajustado correctamente. Resultados:")
  summary(glmm_nb)
}

# Instalar y cargar glmmTMB
if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  install.packages("glmmTMB")
}
library(glmmTMB)

# Ajustar el modelo con glmmTMB
glmm_tmb <- glmmTMB(Abundancia_Total ~ Isla + Cultivo + (1 | Trampa) + (1 | Sector),
                    data = datos_glmm,
                    family = nbinom2)

# Mostrar los resultados
print("Resultados del GLMM con glmmTMB (Negativa Binomial):")
summary(glmm_tmb)

# Guardar los resultados
sink("glmm_tmb_resultados.txt")
print("Resultados del GLMM con glmmTMB (Negativa Binomial):")
summary(glmm_tmb)
sink()

# Instalar y cargar los paquetes necesarios
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

# Cargar los paquetes
library(dplyr)
library(ggplot2)

# Abundancia de las especies clave por Isla
top_especies <- data %>%
  filter(Especie %in% c("Wasmannia_auropunctata", "Hypogastruridae_sp._1", "Solenopsis_geminata")) %>%
  group_by(Isla, Especie) %>%
  summarise(Abundancia = sum(Abundancia, na.rm = TRUE))

# Crear el gráfico
ggplot(top_especies, aes(x = Isla, y = Abundancia, fill = Especie)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Abundancia de Especies Clave por Isla", y = "Abundancia", x = "Isla")

# Guardar el gráfico (opcional)
ggsave("abundancia_especies_clave_por_isla.png", width = 8, height = 6)

# Abundancia de las especies clave por Cultivo
top_especies_cultivo <- data %>%
  filter(Especie %in% c("Wasmannia_auropunctata", "Hypogastruridae_sp._1", "Solenopsis_geminata")) %>%
  group_by(Cultivo, Especie) %>%
  summarise(Abundancia = sum(Abundancia, na.rm = TRUE))

# Crear el gráfico
ggplot(top_especies_cultivo, aes(x = Cultivo, y = Abundancia, fill = Especie)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Abundancia de Especies Clave por Cultivo", y = "Abundancia", x = "Cultivo") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotar etiquetas del eje X

# Guardar el gráfico
ggsave("abundancia_especies_clave_por_cultivo.png", width = 10, height = 6)

# Instalar y cargar paquetes
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}

#####################

setwd("~/Articles/INFORME ZONAS AGRICOLAS/Analisis25/Red")
install.packages("bipartite")
library(bipartite)

isab_funct  = read.csv("isafun.csv", row.names=1)
cris_funct  = read.csv("crifun.csv", row.names=1)
flor_funct  = read.csv("flofun.csv", row.names=1)

# Specialization H2'
networklevel(isab_funct, index="H2")
networklevel(cris_funct, index="H2")
networklevel(flor_funct, index="H2")

# Species Strength
specieslevel(isab_funct, index="species strength")
specieslevel(cris_funct, index="species strength")
specieslevel(flor_funct, index="species strength")


# Robustness
networklevel(isab_funct, index="robustness")
networklevel(cris_funct, index="robustness")
networklevel(flor_funct, index="robustness")

# Connectance
conectancia_manual_i <- sum(isab_funct > 0) / (nrow(isab_funct) * ncol(isab_funct))
print(conectancia_manual_i)
conectancia_manual_c <- sum(cris_funct > 0) / (nrow(cris_funct) * ncol(cris_funct))
print(conectancia_manual_c)
conectancia_manual_f <- sum(flor_funct > 0) / (nrow(flor_funct) * ncol(flor_funct))
print(conectancia_manual_f)

# Plots
isab_funct  = read.csv("isafun.csv", row.names=1, check.names=FALSE)
plotweb(
  isab_funct,  # Usa la matriz SIN ordenar (elimina sortweb)
  arrow = "down",  # Flechas hacia abajo
  col.interaction = "#00008B",  # Color de conexiones
  bor.col.interaction = NA,  # Sin bordes en las conexiones
  text.rot = 90,  # Rotación de etiquetas
  col.high = "darkorange2",  # Color nivel alto (ej: plantas)
  col.low = "#A80A20",  # Color nivel bajo (ej: polinizadores)
  y.lim = c(0, 2.5),  # Ajuste de espacio vertical
  method = "normal"  # Método de visualización
)

cris_funct  = read.csv("crifun.csv", row.names=1, check.names=FALSE)
plotweb(
  cris_funct,  # Usa la matriz SIN ordenar (elimina sortweb)
  arrow = "down",  # Flechas hacia abajo
  col.interaction = "#00008B",  # Color de conexiones
  bor.col.interaction = NA,  # Sin bordes en las conexiones
  text.rot = 90,  # Rotación de etiquetas
  col.high = "darkorange2",  # Color nivel alto (ej: plantas)
  col.low = "#A80A20",  # Color nivel bajo (ej: polinizadores)
  y.lim = c(0, 2.5),  # Ajuste de espacio vertical
  method = "normal"  # Método de visualización
)

# Plots
flor_funct  = read.csv("flofun.csv", row.names=1, check.names=FALSE)
plotweb(
  flor_funct,  # Usa la matriz SIN ordenar (elimina sortweb)
  arrow = "down",  # Flechas hacia abajo
  col.interaction = "#00008B",  # Color de conexiones
  bor.col.interaction = NA,  # Sin bordes en las conexiones
  text.rot = 90,  # Rotación de etiquetas
  col.high = "darkorange2",  # Color nivel alto (ej: plantas)
  col.low = "#BD1100",  # Color nivel bajo (ej: polinizadores)
  y.lim = c(0, 2.5),  # Ajuste de espacio vertical
  method = "normal"  # Método de visualización
)
