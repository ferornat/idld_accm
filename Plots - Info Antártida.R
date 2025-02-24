################################
### Instalamos las librerías ###
################################

# Genéricas
install.packages("remotes")
install.packages("tictoc")
install.packages("repr") # Ésta es para el tamaño de los gráficos

# Trabajamos con la librería de arquetipos. Éstas son algunos links de relevancia:
install.packages("archetypes") # https://cran.r-project.org/web/packages/archetypes/vignettes/archetypes.pdf

# Ésta es para hacer las simulaciones con los arquetipos
install.packages("microbenchmark")

# Instalamos la librería del algoritmo al cual le tenemos que ganar
remotes::install_github("lfernandezpiana/idld")

# Instalamos las librerías que vamos a usar para implementar los cálculos de las distancias de Fréchet y todo el trabajo con la info de HYSPLIT
remotes::install_github("lfernandezpiana/tcd")
install.packages("SimilarityMeasures") # Esto es para calcular la distancia de Fréchet
install.packages("trend")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("mapproj")

# Esta librería es para el clustering de los arquetipos
# install.packages("ktaucenters") # Esto fue para una prueba con clusters robustos.
install.packages("stats") # Esta es para hacer hclust ward.d2

install.packages("reticulate") # Esta es para interpelar a la distancia de Fréchet implementada en Python
# py_install("similaritymeasures") # Vamos a usar esta en realidad más adelante,
# porque la implementación que está en R no funciona bien siempre

install.packages("ggdendro") # Esto es para el dendograma

install.packages("dotenv") # Esta es para hiddear mis direcciones locales en github
install.packages("sessioninfo")

##############################
### Cargamos las librerías ###
##############################

library("dplyr")
library("tictoc")
library("archetypes") # https://cran.r-project.org/web/packages/archetypes/vignettes/archetypes.pdf
library("ggplot2")
library("repr")

library("microbenchmark") # Ésto es para hacer las iteraciones de los arquetipos

# La librería de Lucas y Marcela está albergada en el siguiente repositorio. Allí están tanto los algoritmos de cluster como los de proyecciones de profundidades, bootstrap y demás:
# https://github.com/lfernandezpiana/idld/tree/master/R
library("idld") # Librería del algoritmo a ganar

# Cargamos las librerías que vamos a usar para implementar los cálculos de las distancias de Fréchet y todo el trabajo con la info de HYSPLIT
library("tcd")
library("trend")
library("tidyverse")
library("rnaturalearth")
library("rnaturalearthdata")
library("mapproj")
library("SimilarityMeasures")

# library("ktaucenters") # Esta librería es para el clustering de los arquetipos
library("stats") # Ésta es para el k-means normal, pero también para el Ward.D2
library("cluster") # Ésta es para hacer el PAM
library("reticulate")

# library("ggdendro") # Esto es por si quisieramos hacer un dendograma

library("dotenv")
library("sessioninfo") # Vamos a chequear versiones de paquetes utilizados

# Vemos las versiones de las librerías que estamos usando
# session_info(pkgs = c("dplyr","tictoc","archetypes","ggplot2",
#                       "repr","microbenchmark","idld","tcd",
#                       "trend","tidyverse","rnaturalearth","rnaturalearthdata",
#                       "mapproj","SimilarityMeasures","stats","reticulate",
#                       "dotenv","sessioninfo", "ggdendro"),
#              to_file = TRUE)


#############################
### Proceso de clustering ###
#############################


# Vamos a usar la versión de python de SimilarityMeasures
# https://github.com/cjekel/similarity_measures/blob/master/frechet_distance_recursion_vs_dp.ipynb
# https://github.com/cjekel/similarity_measures/issues/6#issuecomment-544350039. Este link habla de potenciales problemas de convergencia
py_install("similaritymeasures")
frechet <- import("similaritymeasures")

# Creo que lo conveniente es trabajar el formato de la base de datos antes. Asumimos entonces en nuestra función que los datos van a ser
# con TANTOS REGISTROS COMO OBSERVACIONES, con TANTAS COLUMNAS COMO DIMENSIONES TENGA DICHA OBSERVACIÓN.
# En nuestro caso particular tendremos los datos de invierno de un año en particular, donde va a tener 368 filas (pues hay 368 observaciones),
# con 242 columnas (las 121 longitudes observadas en esos cinco días seguidas de las 121 latitudes observadas en los mismos momentos).

#############################
### Estructuramos la base ###
#############################

# Levantamos los datos. Usamos "winter_traj_ver2.rds". Esto tiene la data de "Byers_2005_2016" con las profundidades ya computadas.
# load_dot_env(file = "config.env") # PC 1: Ubuntu (es necesario estás situado dentro del repositorio como directorio de trabajo)
load_dot_env(file = "config_mac.env") # PC 2: Mac (es necesario estás situado dentro del repositorio como directorio de trabajo)
data_path <- Sys.getenv("DATA_PATH") # Esto está definido en DATA_PATH
winter <- readRDS(data_path)


#######################
### Elegimos un año ###
#######################

# En particular, hacemos primero un ejercicio con el invierno de 2007
anio = 3 # 3 es el 2007, 6 es el 2010, 7 es 2011

curvas = winter[[anio]]$curves
length(curvas) # Son 368 curvas
nrow(curvas[[1]]) # Para cada curva tenemos 121 ubicaciones

# Reestructuramos para poder usar el procedimiento
lon_list <- list()
lat_list <- list()

# Iterar sobre cada curva
for (i in 1:length(curvas)) {
  # Convertir lon y lat en matrices, si no lo son ya
  lon_list[[i]] <- curvas[[i]]$lon
  lat_list[[i]] <- curvas[[i]]$lat
}

# Combinar las listas en dataframes
lon_df <- do.call(rbind, lon_list)
lat_df <- do.call(rbind, lat_list)

# Crear el dataframe estructurado y renombrar las columnas
df_struc <- data.frame(lon_df, lat_df) 
colnames(df_struc) <- c(paste0("lon_", 1:ncol(lon_df)), paste0("lat_", 1:ncol(lat_df)))

# Hay que verlo como matriz y ponerlo en formato "double"
df_struc <- as.matrix(df_struc) # Tenemos las 368
storage.mode(df_struc) <- "double"

# Nos quedamos con las curvas más profundas
q_max = which(winter[[anio]]$depth > quantile(winter[[anio]]$depth,0.60)) # 40% más profundo
df_struc <- df_struc[q_max,] 
nrow(df_struc) # Nos quedan 147 curvas profundas

##########################################################################################
### Aquí sería el camino correcto. Es decir, primero la instancia de los core-clusters ###
### y luego el cálculo de los arquetipos.                                              ###
##########################################################################################

#########################################
### Calculamos la matriz de distancia ###
#########################################

n_deep_curves <- nrow(df_struc)

# Hay que poner las matrices en el formato original de trayectorias para poder calcular las distancias
trayectorias_orig <- replicate(n_deep_curves, list(), simplify = FALSE)

# Suponemos que ambas tienen tantas longitudes como latitudes, lo cual siempre va a ser cierto en estos casos
for(cant in 1:n_deep_curves) {
  trayectorias_orig[[cant]] <- df_struc[cant, ]
}

for (i in 1:n_deep_curves) {
  # Extraemos las listas
  current_list <- trayectorias_orig[[i]]
  
  # Obtenemos longitud de las mismas
  n <- length(current_list)
  
  # Partimos en dos la cantidad de columnas
  lon_values <- as.numeric(current_list[1:(n/2)])  
  lat_values <- as.numeric(current_list[(n/2 + 1):n]) 
  
  # Creamos las trayectorias como las teníamos inicialmente
  trayectorias_orig[[i]] <- data.frame(lon = round(lon_values,3), lat = round(lat_values,3)) # Esto es para que nos quede en la misma precisión que lo que venía en el .RDS
}

# Esto nos queda como lo que levantábamos inicialmente del .RDS
trayectorias_orig[[1]]
head(df_struc) # Chequeamos que las curvas nos queden bien

# Inicializamos la matriz de distancias de Fréchet vacía. ACA TUVIMOS QUE MODIFICAR LA DISTANCIA DE FRÉCHET
# Vamos a primero generar todas las distancias entre las curvas pero utilizando Fréchet.

frechet_tot_dist_matrix <- matrix(NA, n_deep_curves, n_deep_curves)

tic()
# Vamos a hacer las distancias, pero únicamente las que están debajo de la diagonal
for (j in 2:n_deep_curves) {
  for (i in 1:(j - 1)) {
    # Hacemos los pares de trayectorias
    traj_i <- as.matrix(trayectorias_orig[[i]])
    traj_j <- as.matrix(trayectorias_orig[[j]])
    
    # Computamos las distancias
    dist_ij <- frechet$frechet_dist(traj_i, traj_j)
    frechet_tot_dist_matrix[j, i] <- dist_ij
  }
}
toc() # Ésto tarda 61.462 segundos


# Lo convertimos en distancia para poder hacer los clusters
frechet_tot_dist_matrix[is.na(frechet_tot_dist_matrix)] <- 0
frechet_tot_dist_matrix <- as.dist(frechet_tot_dist_matrix)

# Armamos la proyección con las referencias para el mapa
wm <- map_data("world")

byers_map = ggplot() +
  geom_polygon(
    data = wm, aes(x = long, y = lat, group = group),
    fill = "white", colour = "black", alpha = 0.8, size=0.3
  ) +
  scale_y_continuous(
    limits = c(-90, -30), breaks = seq(-45, -90, by = -10),
    labels = NULL, expand = c(0, 0)
  ) +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks=element_blank()) +
  coord_map("ortho", orientation = c(-90, -60, 0), xlim=c(-180,-10), ylim=c(-30,-90)) +
  labs(x = NULL, y = NULL)

mapas_path_deep <- Sys.getenv("MAPAS_PATH") # Esto está definido en DATA_PATH

# Primer plot: Año 2007 - Trayectorias por profundidad

tic()
for (k in anio:anio) { # Es 1:length(winter) en realidad, e hicimos 1:1 para que haga el año 2007
  titulo_w <- paste0("Trayectorias por profundidades")
  byers_map_w <- byers_map + ggtitle(titulo_w)
  
  q_max <- which(winter[[k]]$depth > quantile(winter[[k]]$depth, 0.60)) # El n% más profundo
  q_resto <- which(winter[[k]]$depth <= quantile(winter[[k]]$depth, 0.60)) # El resto

  # Curvas menos profundas
  for (j in q_resto) {
    curve <- winter[[k]]$curves[[j]]
    curve$group <- "Menos profundas"
    byers_map_w <- byers_map_w + geom_path(data = curve, aes(x = lon, y = lat, color = group), alpha = 0.5)
  }
  
  # Curvas más profundas
  for (i in q_max) {
    curve <- winter[[k]]$curves[[i]]
    curve$group <- "Más profundas"
    byers_map_w <- byers_map_w + geom_path(data = curve, aes(x = lon, y = lat, color = group), alpha = 0.5)
  }
  
  # Apply the custom color scale and legend
  byers_map_w <- byers_map_w + 
    scale_color_manual(values = c(
      "Menos profundas" = "lightgrey", 
      "Más profundas" = "darkgrey" 
    )) +
    theme(
      plot.title = element_text(size = 40, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 30),
      legend.position = "right"
    )
  
  # Save the plot
  setwd(mapas_path_deep)
  png(filename = paste0(2005 + k - 1, "_depth.png"), width = 800, height = 800)
  print(byers_map_w)
  dev.off()
}
toc()

#####################################
### Calculamos ahora los clusters ###
#####################################

# Este va a ser nuestro primer plot. Vamos a ver, dentro de la información más profunda, cuantos clusters
# nos quedarían.

cantidad_clusters <- 2
pam_result <- pam(frechet_tot_dist_matrix, cantidad_clusters)
clusters <- pam_result$clustering
length(pam_result$clustering) # Son las 147 curvas


# Segundo plot: Año 2007 - Clusters en profundidad

tic()
for (k in anio:anio) { # Es 1:length(winter) en realidad, e hicimos 1:1 para que haga el año 2007
  titulo_w <- paste0("Core-Clusters")
  byers_map_w <- byers_map + ggtitle(titulo_w)
  
  q_max <- which(winter[[k]]$depth > quantile(winter[[k]]$depth, 0.60)) # El n% más profundo
  q_resto <- which(winter[[k]]$depth <= quantile(winter[[k]]$depth, 0.60)) # El resto
  n_clust <- length(clusters)
  
  # Curvas menos profundas
  for (j in q_resto) {
    curve <- winter[[k]]$curves[[j]]
    curve$group <- "Menos profundas"
    byers_map_w <- byers_map_w + geom_path(data = curve, aes(x = lon, y = lat, color = group), alpha = 0.5)
  }
  
  # Curvas más profundas
  for (i in q_max) {
    curve <- winter[[k]]$curves[[i]]
    curve$group <- "Más profundas"
    byers_map_w <- byers_map_w + geom_path(data = curve, aes(x = lon, y = lat, color = group), alpha = 0.5)
  }
  
  # Clusters 
  for (l in 1:n_clust) {
    curve <- trayectorias_orig[[l]]
    cluster_label <- paste0("Cluster ", clusters[l])  # "Cluster 1" o "Cluster 2"
    curve$group <- factor(cluster_label)  # Asegura que se trata como factor
    
    byers_map_w <- byers_map_w + geom_path(data = curve, aes(x = lon, y = lat, color = group), alpha = 0.5)
  }
  
  # Apply the custom color scale and legend
  byers_map_w <- byers_map_w + 
    scale_color_manual(values = c(
      "Menos profundas" = "lightgrey", 
      "Cluster 1" = "darkgreen", 
      "Cluster 2" = "darkred"
    )) +
    theme(
      plot.title = element_text(size = 40, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 30),
      legend.position = "right"
    )
  
  # Save the plot
  setwd(mapas_path_deep)
  png(filename = paste0(2005 + k - 1, "_corecluster_depth.png"), width = 800, height = 800)
  print(byers_map_w)
  dev.off()
}
toc()

# Tercer plot: Año 2007 - Arquetipos de clusters en profundidad

########################################################
### Ahora vamos a sacar 3 arquetipos para cada grupo ###
########################################################

n_arquetipos = 3

# Estrategia para evitar la singularidad en los resultados del cálculo de los arquetipos:
# - Arrancamos con 6 arquetipos.
# - Si funciona, continuamos con el procedimiento.
# - Si no funciona, volvemos a intentar con otra semilla, a ver si se soluciona la iteración
# - Y así sucesivamente n-veces. Si en ninguna de ellas se hallan los arquetipos solicitados, se vuelve a iniciar el procedimiento de la búsqueda de arquetipos con 4 - 1 arquetipos. Se repite todo lo previo, y de no hallarse resultados, se vuelve a inicializar la búsqueda de arquetipos con 4 - 2 arquetipos. Y así, 
# hasta que en el peor de los casos llegue a un arquetipo, total esto siempre va a ser calculable.


# ESTO LO TENEMOS QUE HACER. EN CASO DE QUE LA MATRIZ NO SEA SINGULAR LE TENEMOS QUE AGREGAR UN POCO DE 
# RUIDO PORQUE SI NO SE ROMPE LA FUNCIÓN DE SVD() DE ADENTRO DEL PAQUETE.

retry_archetypes <- function(df, n_arch, seed_arch, max_iter, epsilon = 1e-6) {
  for (i in 1:max_iter) {
    set.seed(seed_arch + i)
    
    # Probamos primero sin ruido
    arch <- tryCatch({
      archetypes(data = df, k = n_arch, verbose = FALSE)
    }, error = function(e) {
      return(NULL)  # NULL si hay un error
    })
    
    if (!is.null(arch)) {
      return(parameters(arch))  # Trae los arquetipos si es exitoso
    }
    
    # Si falla, agrega un ruido
    message("Error encountered, trying with added noise...")
    
    df_noisy <- df + matrix(rnorm(length(df), mean = 0, sd = epsilon), nrow = nrow(df))
    arch <- tryCatch({
      archetypes(data = df_noisy, k = n_arch, verbose = FALSE)
    }, error = function(e) {
      return(NULL)  # NULL si hay un error
    })
    
    if (!is.null(arch)) {
      return(parameters(arch))  # Trae los arquetipos si es exitoso
    }
    
    # Si acá falla, prueba la siguiente iteración
    message("Failed to compute archetypes after adding noise.")
  }
  
  return(NULL)  # Return NULL if no successful result after max_iter
}


# Hacemos un loop de arquetipos para, en caso de que a pesar de las iteraciones y el ruido, no encuentre la cantidad deseada,
# y entonces baje a buscar una cantidad una unidad menor.

# Nos quedamos con los sujetos del grupo 1 y 2 por separado
grupo1 <- trayectorias_orig[clusters == 1]
grupo2 <- trayectorias_orig[clusters == 2]
length(grupo1) + length(grupo2) # Esto me debe dar los 147 con los que venía trabajando

# Obtenemos los arquetipos para cada grupo. Para esto vamos a tener que retransformar las columnas
# y demás

# Grupo 1
lon_list_g1 <- list()
lat_list_g1 <- list()

# Iterar sobre cada curva
for (i in 1:length(grupo1)) {
  # Convertir lon y lat en matrices, si no lo son ya
  lon_list_g1[[i]] <- grupo1[[i]]$lon
  lat_list_g1[[i]] <- grupo1[[i]]$lat
}

# Combinar las listas en dataframes
lon_df_g1 <- do.call(rbind, lon_list_g1)
lat_df_g1 <- do.call(rbind, lat_list_g1)

# Crear el dataframe estructurado y renombrar las columnas
df_struc_g1 <- data.frame(lon_df_g1, lat_df_g1) 
colnames(df_struc_g1) <- c(paste0("lon_", 1:ncol(lon_df_g1)), paste0("lat_", 1:ncol(lat_df_g1)))

# Hay que verlo como matriz y ponerlo en formato "double"
df_struc_g1 <- as.matrix(df_struc_g1) # Tenemos las 368
storage.mode(df_struc_g1) <- "double"

# Obtenemos los arquetipos 
arch_g1 <- retry_archetypes(df_struc_g1, n_arch = 3, seed_arch = 123, max_iter = 10, epsilon = 1e-6)

if (is.null(arch_g1)) {
  message("Failed to compute archetypes with n_arch = 3. Retrying with reduced number of archetypes.")
  
  # Retry with reduced number of archetypes
  arch_g1 <- retry_archetypes(df_struc_g1, n_arch = 2, seed_arch = 123, max_iter = 10, epsilon = 1e-6)
  
  if (is.null(arch_g1)) {
    stop("Failed to compute archetypes with reduced number of archetypes as well.")
  }
}

# Finalmente tenemos los arquetipos para el grupo 1. En este caso obtuvimos los 3 arquetipos
arch_g1

# Grupo 1
lon_list_g2 <- list()
lat_list_g2 <- list()

# Iterar sobre cada curva
for (i in 1:length(grupo2)) {
  # Convertir lon y lat en matrices, si no lo son ya
  lon_list_g2[[i]] <- grupo2[[i]]$lon
  lat_list_g2[[i]] <- grupo2[[i]]$lat
}

# Combinar las listas en dataframes
lon_df_g2 <- do.call(rbind, lon_list_g2)
lat_df_g2 <- do.call(rbind, lat_list_g2)

# Crear el dataframe estructurado y renombrar las columnas
df_struc_g2 <- data.frame(lon_df_g2, lat_df_g2) 
colnames(df_struc_g2) <- c(paste0("lon_", 1:ncol(lon_df_g2)), paste0("lat_", 1:ncol(lat_df_g2)))

# Hay que verlo como matriz y ponerlo en formato "double"
df_struc_g2 <- as.matrix(df_struc_g2) # Tenemos las 368
storage.mode(df_struc_g2) <- "double"

# Obtenemos los arquetipos 
arch_g2 <- retry_archetypes(df_struc_g2, n_arch = 3, seed_arch = 123, max_iter = 10, epsilon = 1e-6)

if (is.null(arch_g2)) {
  message("Failed to compute archetypes with n_arch = 3. Retrying with reduced number of archetypes.")
  
  # Retry with reduced number of archetypes
  arch_g2 <- retry_archetypes(df_struc_g2, n_arch = 2, seed_arch = 123, max_iter = 10, epsilon = 1e-6)
  
  if (is.null(arch_g2)) {
    stop("Failed to compute archetypes with reduced number of archetypes as well.")
  }
}

# Finalmente tenemos los arquetipos para este grupo. En este caso pudimos obtener 2 arquetipos
arch_g2

# Agrupamos los arquetipos y lo estructuramos
arch_g1g2 <- rbind(arch_g1, arch_g2) 
n_arch_g1g2_curves <- nrow(arch_g1g2)

# Hay que poner las matrices en el formato original de trayectorias para poder calcular las distancias
trayectorias_g1g2 <- replicate(n_arch_g1g2_curves, list(), simplify = FALSE)

# Suponemos que ambas tienen tantas longitudes como latitudes, lo cual siempre va a ser cierto en estos casos
for(cant in 1:n_arch_g1g2_curves) {
  trayectorias_g1g2[[cant]] <- arch_g1g2[cant, ]
}

for (i in 1:n_arch_g1g2_curves) {
  # Extraemos las listas
  current_list_g1g2 <- trayectorias_g1g2[[i]]
  
  # Obtenemos longitud de las mismas
  n_g1g2 <- length(current_list_g1g2)
  
  # Partimos en dos la cantidad de columnas
  lon_values_g1g2 <- as.numeric(current_list_g1g2[1:(n_g1g2/2)])  
  lat_values_g1g2 <- as.numeric(current_list_g1g2[(n_g1g2/2 + 1):n_g1g2]) 
  
  # Creamos las trayectorias como las teníamos inicialmente
  trayectorias_g1g2[[i]] <- data.frame(lon = round(lon_values_g1g2,3), lat = round(lat_values_g1g2,3)) # Esto es para que nos quede en la misma precisión que lo que venía en el .RDS
}

# Imprimimos las curvas

tic()
for (k in anio:anio) { # Es 1:length(winter) en realidad, e hicimos 1:1 para que haga el año 2007
  titulo_w <- paste0("Arquetipos de Core-Clusters")
  byers_map_w <- byers_map + ggtitle(titulo_w)
  
  q_max <- which(winter[[k]]$depth > quantile(winter[[k]]$depth, 0.60)) # El n% más profundo
  q_resto <- which(winter[[k]]$depth <= quantile(winter[[k]]$depth, 0.60)) # El resto
  clusters <- c(rep(1,nrow(arch_g1)), rep(2,nrow(arch_g2))) # Esto viene de arriba, de cuando obtuvimos la base con las trayectorias
  n_arc <- length(trayectorias_g1g2)
  
  # Curvas menos profundas
  for (j in q_resto) {
    curve <- winter[[k]]$curves[[j]]
    curve$group <- "Menos profundas"
    byers_map_w <- byers_map_w + geom_path(data = curve, aes(x = lon, y = lat, color = group), alpha = 0.5)
  }
  
  # Curvas más profundas
  for (i in q_max) {
    curve <- winter[[k]]$curves[[i]]
    curve$group <- "Más profundas"
    byers_map_w <- byers_map_w + geom_path(data = curve, aes(x = lon, y = lat, color = group), alpha = 0.5)
  }
  
  # Arquetipos 
  for (l in 1:n_arc) {
    curve <- trayectorias_g1g2[[l]]
    cluster_label <- paste0("Cluster ", clusters[l])  # "Cluster 1" o "Cluster 2"
    curve$group <- factor(cluster_label)  # Asegura que se trata como factor
    
    byers_map_w <- byers_map_w + geom_path(data = curve, aes(x = lon, y = lat, color = group), alpha = 0.5)
  }
  
  # Apply the custom color scale and legend
  byers_map_w <- byers_map_w + 
    scale_color_manual(values = c(
      "Menos profundas" = "lightgrey",
      "Más profundas" = "darkgrey",
      "Cluster 1" = "darkgreen", 
      "Cluster 2" = "darkred"
    )) +
    theme(
      plot.title = element_text(size = 40, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 30),
      legend.position = "right"
    )
  
  # Save the plot
  setwd(mapas_path_deep)
  png(filename = paste0(2005 + k - 1, "_arc_corecluster_depth.png"), width = 800, height = 800)
  print(byers_map_w)
  dev.off()
}
toc()


# Cuarto plot: Año 2007 - Clusters finales

# Ahora debiéramos computarle la distancia de Fréchet a todas las curvas, y clusterizarlas en función de a cual se encuentran
# más cerca

gc()

tic()
# Inicializamos el vector donde vamos a guardar el arquetipo más cercano a la curva
closest_archetype <- numeric(length(winter[[anio]]$curves))
clust_archetypes <- c(rep(1,nrow(arch_g1)), rep(2,nrow(arch_g2))) # Esto viene de arriba, de cuando obtuvimos la base con las trayectorias

# Hacemos el Loop para todas las trayectorias del año
for (j in 1:length(winter[[anio]]$curves)) {
  # Extraemos la observación
  current_curve <- as.matrix(winter[[anio]]$curves[[j]])
  
  # Hacemos la listas a donde va a ir cada distancia respecto a cada arquetipo vis-à-vis
  distances <- numeric(length(trayectorias_g1g2))
  
  # Loopeamos a través de todos las trayectorias arquetipo que tengamos previamente calculadas
  for (i in 1:length(trayectorias_g1g2)) {
    # Exctraemos la i-ésimo trayectoria arquetipo
    archetype_curve <- as.matrix(trayectorias_g1g2[[i]])
    
    # Computamos la distancia de Fréchet
    # distances[i] <- Frechet(current_curve, archetype_curve)[[1]] # Esto era con la implementación vieja
    distances[i] <- frechet$frechet_dist(current_curve, archetype_curve)
  }
  
  # Esto era ANTES
  # # Asignamos el arquetipo más cercano para cada curva
  # closest_archetype[j] <- which.min(distances)
  
  # Esto es AHORA
  # Buscamos el índica del más cercano
  # distances[distances < 0] <- NA # Esto lo pusimos de barrera en "Prueba interna.R", pero porque usábamos la librería de R
  closest_archetype_index <- which.min(distances)
  
  # Asignamos su cluster
  closest_archetype[j] <- clust_archetypes[closest_archetype_index]
}
toc()

# Imprimimos todas las curvas

tic()
for (k in anio:anio) { # Es 1:length(winter) en realidad, e hicimos 1:1 para que haga el año 
  byers_map_w = byers_map
  titulo_w = paste0("Asignación Final")
  byers_map_w = byers_map_w + ggtitle(titulo_w)
  
  for (i in 1:length(winter[[k]]$curves)) {
    curve <- winter[[k]]$curves[[i]]
    # Asignamos el color en función del cluster
    curve$group <- ifelse(closest_archetype[i] == 1, "Cluster 1", "Cluster 2")
    byers_map_w <- byers_map_w + geom_path(data = curve, aes(x = lon, y = lat, color = group), alpha = 0.5)
  }
  
  # Apply the custom color scale and legend
  byers_map_w <- byers_map_w + 
    scale_color_manual(values = c(
      "Cluster 1" = "darkgreen", 
      "Cluster 2" = "darkred"
    )) +
    theme(
      plot.title = element_text(size = 40, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 30),
      legend.position = "right"
    )  
  # Save the plot
  setwd(mapas_path_deep)
  png(filename = paste0(2005 + k - 1, "_asignacion_cluster_final.png"), width = 800, height = 800)
  print(byers_map_w)
  dev.off()
}
toc()

#############################################################################
### Esto es para generar los gráficos para los clusters de todos los años ###
#############################################################################

mapas_path_deep <- Sys.getenv("MAPAS_PATH") # Esto está definido en DATA_PATH

tic()
for(datos_inv in 1:length(winter)){
  # Estructuramos la data
  anio = datos_inv # 3 es el 2007, 6 es el 2010, 7 es 2011
  
  curvas = winter[[anio]]$curves
  
  # Reestructuramos para poder usar el procedimiento
  lon_list <- list()
  lat_list <- list()
  
  # Iterar sobre cada curva
  for (i in 1:length(curvas)) {
    # Convertir lon y lat en matrices, si no lo son ya
    lon_list[[i]] <- curvas[[i]]$lon
    lat_list[[i]] <- curvas[[i]]$lat
  }
  
  # Combinar las listas en dataframes
  lon_df <- do.call(rbind, lon_list)
  lat_df <- do.call(rbind, lat_list)
  
  # Crear el dataframe estructurado y renombrar las columnas
  df_struc <- data.frame(lon_df, lat_df) 
  colnames(df_struc) <- c(paste0("lon_", 1:ncol(lon_df)), paste0("lat_", 1:ncol(lat_df)))
  
  # Hay que verlo como matriz y ponerlo en formato "double"
  df_struc <- as.matrix(df_struc) # Tenemos las 368
  storage.mode(df_struc) <- "double"
  
  # Nos quedamos con las curvas más profundas
  q_max = which(winter[[anio]]$depth > quantile(winter[[anio]]$depth,0.60)) # 40% más profundo
  df_struc <- df_struc[q_max,] 
  
  #########################################
  ### Calculamos la matriz de distancia ###
  #########################################
  
  n_deep_curves <- nrow(df_struc)
  
  # Hay que poner las matrices en el formato original de trayectorias para poder calcular las distancias
  trayectorias_orig <- replicate(n_deep_curves, list(), simplify = FALSE)
  
  # Suponemos que ambas tienen tantas longitudes como latitudes, lo cual siempre va a ser cierto en estos casos
  for(cant in 1:n_deep_curves) {
    trayectorias_orig[[cant]] <- df_struc[cant, ]
  }
  
  for (i in 1:n_deep_curves) {
    # Extraemos las listas
    current_list <- trayectorias_orig[[i]]
    
    # Obtenemos longitud de las mismas
    n <- length(current_list)
    
    # Partimos en dos la cantidad de columnas
    lon_values <- as.numeric(current_list[1:(n/2)])  
    lat_values <- as.numeric(current_list[(n/2 + 1):n]) 
    
    # Creamos las trayectorias como las teníamos inicialmente
    trayectorias_orig[[i]] <- data.frame(lon = round(lon_values,3), lat = round(lat_values,3)) # Esto es para que nos quede en la misma precisión que lo que venía en el .RDS
  }
  
  # Esto nos queda como lo que levantábamos inicialmente del .RDS
  
  # Inicializamos la matriz de distancias de Fréchet vacía. 
  
  frechet_tot_dist_matrix <- matrix(NA, n_deep_curves, n_deep_curves)
  
  # Vamos a hacer las distancias, pero únicamente las que están debajo de la diagonal
  for (j in 2:n_deep_curves) {
    for (i in 1:(j - 1)) {
      # Hacemos los pares de trayectorias
      traj_i <- as.matrix(trayectorias_orig[[i]])
      traj_j <- as.matrix(trayectorias_orig[[j]])
      
      # Computamos las distancias
      dist_ij <- frechet$frechet_dist(traj_i, traj_j)
      frechet_tot_dist_matrix[j, i] <- dist_ij
    }
  }
  
  
  # Lo convertimos en distancia para poder hacer los clusters
  frechet_tot_dist_matrix[is.na(frechet_tot_dist_matrix)] <- 0
  frechet_tot_dist_matrix <- as.dist(frechet_tot_dist_matrix)
  
  #####################################
  ### Calculamos ahora los clusters ###
  #####################################
  
  # Vamos a ir con dos clusters como vimos en el ejemplo
  
  cantidad_clusters <- 2
  pam_result <- pam(frechet_tot_dist_matrix, cantidad_clusters)
  clusters <- pam_result$clustering
  
  # Obtenemos los arquetipos. 
  n_arquetipos = 3
  
  retry_archetypes <- function(df, n_arch, seed_arch, max_iter, epsilon = 1e-6) {
    for (i in 1:max_iter) {
      set.seed(seed_arch + i)
      
      # Probamos primero sin ruido
      arch <- tryCatch({
        archetypes(data = df, k = n_arch, verbose = FALSE)
      }, error = function(e) {
        return(NULL)  # NULL si hay un error
      })
      
      if (!is.null(arch)) {
        return(parameters(arch))  # Trae los arquetipos si es exitoso
      }
      
      # Si falla, agrega un ruido
      message("Error encountered, trying with added noise...")
      
      df_noisy <- df + matrix(rnorm(length(df), mean = 0, sd = epsilon), nrow = nrow(df))
      arch <- tryCatch({
        archetypes(data = df_noisy, k = n_arch, verbose = FALSE)
      }, error = function(e) {
        return(NULL)  # NULL si hay un error
      })
      
      if (!is.null(arch)) {
        return(parameters(arch))  # Trae los arquetipos si es exitoso
      }
      
      # Si acá falla, prueba la siguiente iteración
      message("Failed to compute archetypes after adding noise.")
    }
    
    return(NULL)  # Return NULL if no successful result after max_iter
  }
  
  # Nos quedamos con los sujetos del grupo 1 y 2 por separado
  grupo1 <- trayectorias_orig[clusters == 1]
  grupo2 <- trayectorias_orig[clusters == 2]
  length(grupo1) + length(grupo2) # Esto me debe dar los 147 con los que venía trabajando
  
  # Obtenemos los arquetipos para cada grupo. Para esto vamos a tener que retransformar las columnas
  # y demás
  
  # Grupo 1
  lon_list_g1 <- list()
  lat_list_g1 <- list()
  
  # Iterar sobre cada curva
  for (i in 1:length(grupo1)) {
    # Convertir lon y lat en matrices, si no lo son ya
    lon_list_g1[[i]] <- grupo1[[i]]$lon
    lat_list_g1[[i]] <- grupo1[[i]]$lat
  }
  
  # Combinar las listas en dataframes
  lon_df_g1 <- do.call(rbind, lon_list_g1)
  lat_df_g1 <- do.call(rbind, lat_list_g1)
  
  # Crear el dataframe estructurado y renombrar las columnas
  df_struc_g1 <- data.frame(lon_df_g1, lat_df_g1) 
  colnames(df_struc_g1) <- c(paste0("lon_", 1:ncol(lon_df_g1)), paste0("lat_", 1:ncol(lat_df_g1)))
  
  # Hay que verlo como matriz y ponerlo en formato "double"
  df_struc_g1 <- as.matrix(df_struc_g1) # Tenemos las 368
  storage.mode(df_struc_g1) <- "double"
  
  # Obtenemos los arquetipos 
  arch_g1 <- retry_archetypes(df_struc_g1, n_arch = 3, seed_arch = 123, max_iter = 10, epsilon = 1e-6)
  
  if (is.null(arch_g1)) {
    message("Failed to compute archetypes with n_arch = 3. Retrying with reduced number of archetypes.")
    
    # Retry with reduced number of archetypes
    arch_g1 <- retry_archetypes(df_struc_g1, n_arch = 2, seed_arch = 123, max_iter = 10, epsilon = 1e-6)
    
    if (is.null(arch_g1)) {
      stop("Failed to compute archetypes with reduced number of archetypes as well.")
    }
  }
  
  # Finalmente tenemos los arquetipos para el grupo 1. En este caso obtuvimos los 3 arquetipos
  
  # Grupo 1
  lon_list_g2 <- list()
  lat_list_g2 <- list()
  
  # Iterar sobre cada curva
  for (i in 1:length(grupo2)) {
    # Convertir lon y lat en matrices, si no lo son ya
    lon_list_g2[[i]] <- grupo2[[i]]$lon
    lat_list_g2[[i]] <- grupo2[[i]]$lat
  }
  
  # Combinar las listas en dataframes
  lon_df_g2 <- do.call(rbind, lon_list_g2)
  lat_df_g2 <- do.call(rbind, lat_list_g2)
  
  # Crear el dataframe estructurado y renombrar las columnas
  df_struc_g2 <- data.frame(lon_df_g2, lat_df_g2) 
  colnames(df_struc_g2) <- c(paste0("lon_", 1:ncol(lon_df_g2)), paste0("lat_", 1:ncol(lat_df_g2)))
  
  # Hay que verlo como matriz y ponerlo en formato "double"
  df_struc_g2 <- as.matrix(df_struc_g2) # Tenemos las 368
  storage.mode(df_struc_g2) <- "double"
  
  # Obtenemos los arquetipos 
  arch_g2 <- retry_archetypes(df_struc_g2, n_arch = 3, seed_arch = 123, max_iter = 10, epsilon = 1e-6)
  
  if (is.null(arch_g2)) {
    message("Failed to compute archetypes with n_arch = 3. Retrying with reduced number of archetypes.")
    
    # Retry with reduced number of archetypes
    arch_g2 <- retry_archetypes(df_struc_g2, n_arch = 2, seed_arch = 123, max_iter = 10, epsilon = 1e-6)
    
    if (is.null(arch_g2)) {
      stop("Failed to compute archetypes with reduced number of archetypes as well.")
    }
  }
  
  # Finalmente tenemos los arquetipos para este grupo. En este caso pudimos obtener 2 arquetipos
  
  # Agrupamos los arquetipos y lo estructuramos
  arch_g1g2 <- rbind(arch_g1, arch_g2) 
  n_arch_g1g2_curves <- nrow(arch_g1g2)
  
  # Hay que poner las matrices en el formato original de trayectorias para poder calcular las distancias
  trayectorias_g1g2 <- replicate(n_arch_g1g2_curves, list(), simplify = FALSE)
  
  # Suponemos que ambas tienen tantas longitudes como latitudes, lo cual siempre va a ser cierto en estos casos
  for(cant in 1:n_arch_g1g2_curves) {
    trayectorias_g1g2[[cant]] <- arch_g1g2[cant, ]
  }
  
  for (i in 1:n_arch_g1g2_curves) {
    # Extraemos las listas
    current_list_g1g2 <- trayectorias_g1g2[[i]]
    
    # Obtenemos longitud de las mismas
    n_g1g2 <- length(current_list_g1g2)
    
    # Partimos en dos la cantidad de columnas
    lon_values_g1g2 <- as.numeric(current_list_g1g2[1:(n_g1g2/2)])  
    lat_values_g1g2 <- as.numeric(current_list_g1g2[(n_g1g2/2 + 1):n_g1g2]) 
    
    # Creamos las trayectorias como las teníamos inicialmente
    trayectorias_g1g2[[i]] <- data.frame(lon = round(lon_values_g1g2,3), lat = round(lat_values_g1g2,3)) # Esto es para que nos quede en la misma precisión que lo que venía en el .RDS
  }

  # Inicializamos el vector donde vamos a guardar el arquetipo más cercano a la curva
  closest_archetype <- numeric(length(winter[[anio]]$curves))
  clust_archetypes <- c(rep(1,nrow(arch_g1)), rep(2,nrow(arch_g2))) # Esto viene de arriba, de cuando obtuvimos la base con las trayectorias
  
  # Hacemos el Loop para todas las trayectorias del año
  for (j in 1:length(winter[[anio]]$curves)) {
    # Extraemos la observación
    current_curve <- as.matrix(winter[[anio]]$curves[[j]])
    
    # Hacemos la listas a donde va a ir cada distancia respecto a cada arquetipo vis-à-vis
    distances <- numeric(length(trayectorias_g1g2))
    
    # Loopeamos a través de todos las trayectorias arquetipo que tengamos previamente calculadas
    for (i in 1:length(trayectorias_g1g2)) {
      # Exctraemos la i-ésimo trayectoria arquetipo
      archetype_curve <- as.matrix(trayectorias_g1g2[[i]])
      
      # Computamos la distancia de Fréchet
      # distances[i] <- Frechet(current_curve, archetype_curve)[[1]] # Esto era con la implementación vieja
      distances[i] <- frechet$frechet_dist(current_curve, archetype_curve)
    }
  
    closest_archetype_index <- which.min(distances)
    
    # Asignamos su cluster
    closest_archetype[j] <- clust_archetypes[closest_archetype_index]
  }
  
  # Mapa del mundo
  wm <- map_data("world")
  
  byers_map = ggplot() +
    geom_polygon(
      data = wm, aes(x = long, y = lat, group = group),
      fill = "white", colour = "black", alpha = 0.8, size=0.3
    ) +
    scale_y_continuous(
      limits = c(-90, -30), breaks = seq(-45, -90, by = -10),
      labels = NULL, expand = c(0, 0)
    ) +
    scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
    theme(panel.background = element_rect(fill = "white"),
          axis.ticks=element_blank()) +
    coord_map("ortho", orientation = c(-90, -60, 0), xlim=c(-180,-10), ylim=c(-30,-90)) +
    labs(x = NULL, y = NULL)
  
  # Armamos el plot del año
  for (k in anio:anio) {  
    byers_map_w = byers_map
    # titulo_w = paste0(2005 + k - 1) # Queda mejor si le ponemos el año en Latex
    byers_map_w = byers_map_w # + ggtitle(titulo_w)
    
    # Add geom_path with color mapping based on closest_archetype
    for (i in 1:length(winter[[k]]$curves)) {
      byers_map_w = byers_map_w + 
        geom_path(data = winter[[k]]$curves[[i]], 
                  aes(x = lon, y = lat), 
                  color = ifelse(closest_archetype[i] == 1, "darkred", 
                                 ifelse(closest_archetype[i] == 2, "darkgreen", "grey")), 
                  alpha = 0.5)
    }
    
    # Apply title and theme
    byers_map_w = byers_map_w + 
      theme(
        plot.title = element_text(size = 40, face = "bold"),
        legend.text = element_text(size = 30)
      )
    
    # Save the plot
    setwd(mapas_path_deep)
    png(filename = paste0(2005 + k - 1, "_clusters_final.png"), width = 800, height = 800)
    print(byers_map_w)
    dev.off()
  }
}
toc()