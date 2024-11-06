################################
### Instalamos las librerías ###
################################

# Genéricas
install.packages("remotes")
install.packages("tictoc")
install.packages("repr") # Ésta es para el tamaño de los gráficos

# Para evaluar clasificación contra los clusters reales
install.packages("clue") # Ésta es para la función de reasignado de etiquetas usando la asignación húngara y hacer CCR
install.packages("fossil") # Esto es para usar el RI
install.packages("mclust") # Estp es para usar el ARI

# Trabajamos con la librería de arquetipos. Éstas son algunos links de relevancia:
install.packages("archetypes") # https://cran.r-project.org/web/packages/archetypes/vignettes/archetypes.pdf
install.packages("pracma")
install.packages("mclust")

# Ésta es para hacer las simulaciones con los arquetipos
install.packages("microbenchmark")

# Instalamos la librería del algoritmo al cual le tenemos que ganar
remotes::install_github("lfernandezpiana/idld")

# Instalamos las librerías que vamos a usar para implementar los cálculos de las distancias de Fréchet y todo el trabajo con la info de HYSPLIT
remotes::install_github("lfernandezpiana/tcd")
install.packages("SimilarityMeasures") # Esto es para calcular la distancia de Fréchet
install.packages("trend")
# install.packages("maptools") # Esta era en la primera versión de "byers_map"
# install.packages("ggmap")    # Esta era en la primera versión de "byers_map"
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("mapproj")

# Esta librería es para el clustering de los arquetipos
# install.packages("ktaucenters") 
install.packages("stats") # Esta es para hacer K-means

install.packages("reticulate")
# py_install("similaritymeasures") # Vamos a usar esta en realidad más adelante,
# porque la implementación que está en R no funciona bien siempre

##############################
### Cargamos las librerías ###
##############################

library("mvnfast")
library('dplyr')
library("tictoc")
library("cluster")
library("archetypes") # https://cran.r-project.org/web/packages/archetypes/vignettes/archetypes.pdf
library("ggplot2")
library("repr")
library("clue") # Esto es para hacer CCR
library("fossil") # Esto es para hacer RI
library("mclust") # Esto es para hacer el ARI

library("microbenchmark") # Ésto es para hacer las iteraciones de los arquetipos

library("mclust")
library("mvtnorm") # Esto es para poder usar "rmrnorm".

require('geometry')

# La librería de Lucas y Marcela está albergada en el siguiente repositorio. Allí están tanto los algoritmos de cluster como los de proyecciones de profundidades, bootstrap y demás:
# https://github.com/lfernandezpiana/idld/tree/master/R
library("idld") # Librería del algoritmo a ganar

# Cargamos las librerías que vamos a usar para implementar los cálculos de las distancias de Fréchet y todo el trabajo con la info de HYSPLIT
library("tcd")
library("trend")
library("tidyverse")
# library("maptools") # Esta era en la primera versión de "byers_map"
# library("ggmap")    # Esta era en la primera versión de "byers_map"
library("rnaturalearth")
library("rnaturalearthdata")
library("mapproj")
library("SimilarityMeasures")

# library("ktaucenters") # Esta librería es para el clustering de los arquetipos
library("stats") # Ésta es para el k-means normal, pero también para el Ward.D2
library("reticulate")

#######################
### cluster_alloc() ###
#######################

# Lo que vamos a hacer ahora es hacer una comparativa entre la estrategia de clustering actual y la propuesta con arquetipos, 
# comparando sus tiempos de cómputos y errores de clasificación para algún dataset simulado.
#
# Definimos las funciones que utilizamos en ambos casos.
#
# Ésta es la función que proponemos con la estrategia del cambio de semilla y la reducción de arquetipos.
#
# RECORDAMOS: Nosotros habíamos empezado a armar una función que pusiera una serie detrás de la otra, pero al final usamos la función genérica
# de cluster_alloc() que trabajábamos previamente en las simulaciones.

# Definimos las distancias. Calculamos la distancia a cada punto y asignamos la mínima

distancia_euclidea_fer <- function(point, hull) {
  hull_distances = apply(hull, 1, function(hull_point) {
    sqrt(sum((point - hull_point)^2))
  })
  return(min(hull_distances))
}

distancia_mahalanobis_fer <- function(point, hull_points, cov_matrix) {    # Le cambiamos el nombre porque si no nos rompe de arriba
  cluster_distances = apply(hull_points, 1, function(hull_point) {
    diff = point - hull_point
    sqrt(t(diff) %*% solve(cov_matrix) %*% diff)
  })
  return(min(cluster_distances))
}

distancia_frechet_fer <-function(){}

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

# Levantamos los datos
data_path = "/home/pepino/Desktop/Tesis Maestría Ciencia de Datos - Udesa - 2024/Paper Antártida/winter_traj_ver2.rds"
winter <- readRDS(data_path)

# En particular, hacemos primero un ejercicio con el invierno de 2025
anio = 7 # 3 es el 2007, 6 es el 2010, 7 es 2011
# curvas_2005 = winter[[1]]$curves

curvas_2011 = winter[[anio]]$curves
length(curvas_2011)
nrow(curvas_2011[[1]])

# Reestructuramos para poder usar el procedimiento
lon_list <- list()
lat_list <- list()

# Iterar sobre cada curva
for (i in 1:length(curvas_2011)) {
  # Convertir lon y lat en matrices, si no lo son ya
  lon_list[[i]] <- curvas_2011[[i]]$lon
  lat_list[[i]] <- curvas_2011[[i]]$lat
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
df_struc <- df_struc[q_max,] # Nos quedan 92
nrow(df_struc)

# ESTO LO TENEMOS QUE HACER. EN CASO DE QUE LA MATRIZ NO SEA SINGULAR LE TENEMOS QUE AGREGAR UN POCO DE 
# RUIDO PORQUE SI NO SE ROMPE LA FUNCIÓN DE SVD() DE ADENTRO DEL PAQUETE.

# Definimos un ruido
epsilon <- 1e-6
n_arquetipos <- 6 # Esto puede romper

# Tratamos de hacer archetipos. Si rompe así como viene, vamos a agregar un ruido
arquetipos_prueba <- tryCatch({
  # Intento sin ruido
  arquetipos_prueba <-archetypes(df_struc, k = n_arquetipos)
}, error = function(e) {
  message("Error encountered: ", e$message)
  message("Lo hacemos agregando un ruido")
  
  # Intento con ruido
  df_struc_noisy <- df_struc + matrix(rnorm(length(df_struc), mean = 0, sd = epsilon), nrow = nrow(df_struc))
  arquetipos_prueba <- archetypes(df_struc_noisy, k = n_arquetipos) 
})

# OBS: Ojo que acá si le ponemos muchos puede rompernos y quedarnos de nuevo un problema de triangulación.
# En realidad el problema que estamos resolviendo acá es que queda con problemas de singularidad al momento de comenzar,
# no es que comienza y luego no llega a converger (eso fue lo que nos pasó con la simulación que estábamos trabajando antes)

# Aquí tenemos que clusterizar los arquetipos
# https://github.com/cran/ktaucenters

# tic()
# clusters_arquetipos <- ktaucenters(arquetipos_prueba$archetypes, K = 3, nstart = 1234)
# toc()

# tic()
# clusters_arquetipos <- kmeans(arquetipos_prueba$archetypes, centers = 3, nstart = 1234)
# toc()

# length(clusters_arquetipos$cluster) # Te los clusteriza a los 16, aunque nos vamos a quedar con los centros
# arquecentros <- clusters_arquetipos$centers
# colnames(arquecentros) <- colnames(arquetipos_prueba$archetypes) # Actualizamos los nombres de las columnas


# Volvemos a armar los arquetipos como trayectorias
arquetipos <- arquetipos_prueba$archetypes
trayectorias_arquetipicas <- replicate(n_arquetipos, list(), simplify = FALSE)

# Suponemos que ambas tienen tantas longitudes como latitudes, lo cual siempre va a ser cierto en estos casos
for(cant in 1:n_arquetipos) {
 trayectorias_arquetipicas[[cant]] <- arquetipos[cant, ]
}

for (i in 1:n_arquetipos) {
 # Extraemos las listas
 current_list <- trayectorias_arquetipicas[[i]]
  
 # Obtenemos longitud de las mismas
 n <- length(current_list)
  
# Partimos en dos la cantidad de columnas
 lon_values <- as.numeric(current_list[1:(n/2)])  
 lat_values <- as.numeric(current_list[(n/2 + 1):n]) 
  
# Creamos las trayectorias como las teníamos inicialmente
 trayectorias_arquetipicas[[i]] <- data.frame(lon = round(lon_values,3), lat = round(lat_values,3)) # Esto es para que nos quede en la misma precisión que lo que venía en el .RDS
 }

# Esto nos queda como lo que levantábamos inicialmente del .RDS
trayectorias_arquetipicas

# ANTES DE CLUSTERIZAR VAMOS A NECESITAR ARMAR LA MATRIZ DE DISTANCIAS DE FRÉCHET
class(trayectorias_arquetipicas)
length(trayectorias_arquetipicas)

# El cálculo de la matriz puede llevar algún tiempo
# https://cran.r-project.org/web/packages/SimilarityMeasures/SimilarityMeasures.pdf

# Inicializamos la matriz de distancias de Fréchet vacía. ACA TUVIMOS QUE MODIFICAR LA DISTANCIA DE FRÉCHET
frechet_dist_matrix <- matrix(NA, n_arquetipos, n_arquetipos)

tic()
# Vamos a hacer las distancias, pero únicamente las que están debajo de la diagonal
for (j in 2:n_arquetipos) {
  for (i in 1:(j - 1)) {
    # Hacemos los pares de trayectorias
    traj_i <- as.matrix(trayectorias_arquetipicas[[i]])
    traj_j <- as.matrix(trayectorias_arquetipicas[[j]])
    
    # Computamos las distancias
    dist_ij <- frechet$frechet_dist(traj_i, traj_j)
    frechet_dist_matrix[j, i] <- dist_ij
  }
}
toc()

# 253.447 sec elapsed para 2007, con 10 arquetipos
# 55.999 para 2007, con 6 arquetipos
# Con esta implementación en Python vuela. Bajé de 20 segundos a 0.198 !!!!!!!!!!!!!!!!

# Ya verificamos que nos da los mismos valores en los caso que se puede calcular con la función
# que veníamos utilizando
# frechet_dist_matrix # Valor de librería nueva
# Frechet(as.matrix(trayectorias_arquetipicas[[1]]), as.matrix(trayectorias_arquetipicas[[5]])) # Valor de librería vieja

# Lo convertimos en distancia
frechet_dist_matrix[is.na(frechet_dist_matrix)] <- 0
frechet_dist_matrix <- as.dist(frechet_dist_matrix)

dendograma = hclust(frechet_dist_matrix, method = "ward.D2", members = NULL)
plot(dendograma, main = "Dendrograma de ward.D2 sobre la matriz de distancia de Fréchet", hang = -1)

# Cortamos el dendograma para hacer los clusters
clusters <- cutree(dendograma, k = 2)
print(clusters)
rect.hclust(dendograma, k = 2, border = "red")

### Ploteamos a ver si tiene sentido

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

byers_map

length(winter) # Esto es todo invierno, que son 12 años
length(curvas_2011) # Esto es el caso que tomamos, que tiene 368 trayectorias


# Prueba Manual
titulo_w = paste0("Probamos - para borrar")

byers_map_w = byers_map + ggtitle(titulo_w)
n_arc = length(trayectorias_arquetipicas)
clust_arc = clusters # Van a ser los colores de los clusters del arquetipo
clust_arc

byers_map_w = byers_map_w + geom_path(data = trayectorias_arquetipicas[[1]], aes(x=lon, y=lat), color="darkblue", alpha = 0.5)
byers_map_w = byers_map_w + geom_path(data = trayectorias_arquetipicas[[2]], aes(x=lon, y=lat), color="darkblue", alpha = 0.5)
byers_map_w = byers_map_w + geom_path(data = trayectorias_arquetipicas[[3]], aes(x=lon, y=lat), color="darkblue", alpha = 0.5)
byers_map_w = byers_map_w + geom_path(data = trayectorias_arquetipicas[[4]], aes(x=lon, y=lat), color="darkred", alpha = 0.5)
byers_map_w = byers_map_w + geom_path(data = trayectorias_arquetipicas[[5]], aes(x=lon, y=lat), color="darkblue", alpha = 0.5)
byers_map_w = byers_map_w + geom_path(data = trayectorias_arquetipicas[[6]], aes(x=lon, y=lat), color="darkblue", alpha = 0.5)

byers_map_w = byers_map_w + theme(plot.title = element_text(size = 30, face = "bold"))
print(byers_map_w)

# Continuamos con los gráficos de rigor

# Ploteamos
mapas_path_deep = "/home/pepino/Desktop/Tesis Maestría Ciencia de Datos - Udesa - 2024/Paper Antártida/Archivos de Lucas/mapas/2011/3) Pruebas con clusters/6 archetipos - 2 clusters - 40% deeper - CORREGIDO"

# Esto es para el plot de los datos totales, los datos más profundos y los arquetipos
tic()
for (k in anio:anio) { # Es 1:length(winter) en realidad, e hicimos 1:1 para que haga el año 2005
  titulo_w = paste0("Prueba de arquetipos al 40% - ", 2005+k-1)
  byers_map_w = byers_map + ggtitle(titulo_w)
  q_max = which(winter[[k]]$depth > quantile(winter[[k]]$depth,0.60)) # el n% más profundo
  q_resto = which(winter[[k]]$depth <= quantile(winter[[k]]$depth,0.60)) # El resto
  n_arc = length(trayectorias_arquetipicas)
  clust_arc = clusters # Van a ser los colores de los clusters del arquetipo
  
  # Curvas menos profundas
  for (j in q_resto) 
    byers_map_w = byers_map_w + geom_path(data = winter[[k]]$curves[[j]], aes(x=lon, y=lat), color="lightgrey", alpha = 0.5)
  
  # Curvas más profundas
  for (i in q_max)
    byers_map_w = byers_map_w + geom_path(data = winter[[k]]$curves[[i]], aes(x=lon, y=lat), color = "darkred", alpha = 0.5)
  
  # Esto es para el plot con todos los arquetipos con el mismo color
  for (i in 1:n_arc)
   byers_map_w = byers_map_w + geom_path(data = trayectorias_arquetipicas[[i]], aes(x=lon, y=lat), color="darkblue", alpha = 0.5)
  
  
  byers_map_w = byers_map_w + theme(plot.title = element_text(size = 30, face = "bold"))
  
  setwd(mapas_path_deep)
  png(filename=paste0(2005+k-1,"_depth.png"), width=800, height=800)
  print(byers_map_w)
  dev.off()
}
toc()

# Esto es para ver los colores de los clusters del arquetipo
tic()
for (k in anio:anio) { # Es 1:length(winter) en realidad, e hicimos 1:1 para que haga el año 2005
  titulo_w = paste0("Prueba de arquetipos al 40% - ", 2005+k-1)
  byers_map_w = byers_map + ggtitle(titulo_w)
  q_max = which(winter[[k]]$depth > quantile(winter[[k]]$depth,0.60)) # el n% más profundo
  q_resto = which(winter[[k]]$depth <= quantile(winter[[k]]$depth,0.60)) # El resto
  n_arc = length(trayectorias_arquetipicas)
  clust_arc = clusters # Van a ser los colores de los clusters del arquetipo
  
  # Curvas menos profundas
  for (j in q_resto) 
    byers_map_w = byers_map_w + geom_path(data = winter[[k]]$curves[[j]], aes(x=lon, y=lat), color="lightgrey", alpha = 0.5)
  
  # Esto es para tener a cada arquetipo con los colores del cluster
  for (i in 1:n_arc)
   byers_map_w = byers_map_w + geom_path(data = trayectorias_arquetipicas[[i]], aes(x=lon, y=lat), color = clust_arc[i] , alpha = 0.5)
  
  byers_map_w = byers_map_w + theme(plot.title = element_text(size = 30, face = "bold"))
  
  setwd(mapas_path_deep)
  png(filename=paste0(2005+k-1,"_arc_para_cluster.png"), width=800, height=800)
  print(byers_map_w)
  dev.off()
}
toc()


# Ahora debiéramos computarle la distancia de Fréchet a todas las curvas, y clusterizarlas en función de a cual se encuentran
# más cerca


# prueba1 <- as.matrix(winter[[1]]$curves[[29]])
# prueba2 <- as.matrix(winter[[1]]$curves[[30]])
# Esto es con "SimilarityMeasures"
# tic()
# frechet_distance <- Frechet(prueba1, prueba2)[[1]] # Tarda 18.924 sec.
# toc()
# print(frechet_distance) # 5.207116 de distancia. ¿En qué unidad está?
# frechet_distance[[1]]

gc()

tic()
# Inicializamos el vector donde vamos a guardar el arquetipo más cercano a la curva
closest_archetype <- numeric(length(winter[[anio]]$curves))
clust_archetypes <- clusters # Esto viene de arriba, de hacer un cluster jerárquico sobre la matriz de distancia de Fréchet

# Hacemos el Loop para todas las trayectorias del año
for (j in 1:length(winter[[anio]]$curves)) {
  # Extraemos la observación
  current_curve <- as.matrix(winter[[anio]]$curves[[j]])
  
  # Hacemos la listas a donde va a ir cada distancia respecto a cada arquetipo vis-à-vis
  distances <- numeric(length(trayectorias_arquetipicas))
  
  # Loopeamos a través de todos las trayectorias arquetipo que tengamos previamente calculadas
  for (i in 1:length(trayectorias_arquetipicas)) {
    # Exctraemos la i-ésimo trayectoria arquetipo
    archetype_curve <- as.matrix(trayectorias_arquetipicas[[i]])
    
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


# Para 3 arquetipos tardó 7022.415 sec elapsed. Osea, casi dos horas
# Para 2 arquetipos tardó 5197.599 sec elapsed
# Para 6 arquetipos tardó 10486 sec, osea 2.91 horas

# CON LA NUEVA IMPLEMENTACIÓN DE PYTHON TARDÓ 21.908 SEGUNDOS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Ploteamos
mapas_path_deep = "/home/pepino/Desktop/Tesis Maestría Ciencia de Datos - Udesa - 2024/Paper Antártida/Archivos de Lucas/mapas/2011/3) Pruebas con clusters/6 archetipos - 2 clusters - 40% deeper - CORREGIDO"

tic()
for (k in anio:anio) { # Es 1:length(winter) en realidad, e hicimos 1:1 para que haga el año 2005
  byers_map_w = byers_map
  titulo_w = paste0("Prueba de clusters al 40% - ", 2005+k-1)
  byers_map_w = byers_map_w + ggtitle(titulo_w)
  for (i in 1:length(winter[[k]]$curves))
    byers_map_w = byers_map_w + geom_path(data=winter[[k]]$curves[[i]], aes(x = lon, y = lat, ), color = as.factor(closest_archetype[i]),  alpha=0.5)
  byers_map_w = byers_map_w + theme(plot.title = element_text(size = 30, face = "bold"))
  
  setwd(mapas_path_deep)
  png(filename=paste0(2005+k-1,"_cluster.png"), width=800, height=800)
  print(byers_map_w)
  dev.off()
}
toc()



# Tomamos 10 arquetipos cómo mucho
# Tomamos las distancia de Fréchet entre esos 10 arquetipos
# Sobre esta matriz de distancia armamos un cluster con Ward 2 (es parecido a K-medias)
# Nos quedamos con 2 clusters, cada uno formado por pocos arquetipos. Estos van a ser el core cluster del procedimiento general.









cluster_alloc_fer <- function(X, ld, alpha_q = seq(0.5, 0.9, 0.1), K, distance = "euclidean", n_arch = 3, seed_arch = 9102, max_iter = 10) {
  d <- dim(X)
  alpha_length <- length(alpha_q)
  cluster <- matrix(NA, nrow = d[1], ncol = alpha_length)
  archetypes <- vector("list", length = K * alpha_length)
  
  retry_archetypes <- function(df, n_arch, seed_arch, max_iter) {
    for (i in 1:max_iter) {
      set.seed(seed_arch + i)
      tryCatch({
        arch <- archetypes(data = df, k = n_arch, verbose = FALSE)
        return(parameters(arch))  # Return the archetypes if successful
      }, warning = function(w) {
        if (grepl("alphas > maxKappa", w$message)) {
          return(NULL)
        } else {
          return(NULL)
        }
      }, error = function(e) {
        return(NULL)
      })
    }
    return(NULL)  # Return NULL if no successful result after max_iter
  }
  
  for (j in 1:alpha_length) {
    region_index <- which(ld > quantile(ld, alpha_q[j]))
    others_index <- which(ld <= quantile(ld, alpha_q[j]))
    
    if (length(region_index) <= K + 1) {
      next
    }
    
    # Clustering for the main region
    region_cluster <- kmeans(X[region_index, ], centers = K)$cluster
    cluster[region_index, j] <- region_cluster
    
    # Compute archetypes for each cluster
    for (k in 1:K) {
      df_cluster <- data.frame(X[which(cluster[, j] == k), ])
      arch <- retry_archetypes(df_cluster, n_arch, seed_arch, max_iter)
      
      if (is.null(arch)) {
        message("Failed to compute archetypes for cluster ", k, " at alpha_q[", j, "] with n_arch = ", n_arch)
        # Try with one less archetype
        arch <- retry_archetypes(df_cluster, n_arch - 1, seed_arch, max_iter)
        
        if (is.null(arch)) {
          stop("Failed to compute archetypes with reduced number as well.")
        }
      }
      archetypes[[k + (j - 1) * K]] <- arch
    }
    
    # Assign clusters to remaining observations based on convex hull distances or archetypes
    for (i in others_index) {
      min_dist <- Inf
      closest_cluster <- NA
      
      for (k in 1:K) {
        archetype_point <- archetypes[[k + (j - 1) * K]]
        
        if (distance == "euclidean") {
          dist_to_archetype <- distancia_euclidea_fer(X[i, ], archetype_point)
        } else if (distance == "mahalanobis") {
          cluster_points <- X[region_index[region_cluster == k], ]
          cov_matrix <- var(cluster_points)
          dist_to_archetype <- distancia_mahalanobis_fer(X[i, ], archetype_point, cov_matrix)
        } else {
          stop("Invalid distance parameter")
        }
        
        if (!is.na(dist_to_archetype) && dist_to_archetype < min_dist) {
          min_dist <- dist_to_archetype
          closest_cluster <- k
        }
      }
      
      cluster[i, j] <- closest_cluster
    }
  }
  
  return(list(cluster = cluster, archetypes = archetypes))
}




