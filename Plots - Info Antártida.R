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
library("reticulate")

library("ggdendro")

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
load_dot_env(file = "config.env") # PC 1: Ubuntu (es necesario estás situado dentro del repositorio como directorio de trabajo)
# load_dot_env(file = "config_mac.env") # PC 2: Mac (es necesario estás situado dentro del repositorio como directorio de trabajo)
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

###########################
### Buscamos arquetipos ###
###########################

n_arquetipos = 6

# Estrategia para evitar la singularidad en los resultados del cálculo de los arquetipos:
# - Arrancamos con 6 arquetipos.
# - Si funciona, continuamos con el procedimiento.
# - Si no funciona, volvemos a intentar con otra semilla, a ver si se soluciona la iteración
# - Y así sucesivamente n-veces. Si en ninguna de ellas se hallan los arquetipos solicitados, se vuelve a iniciar el procedimiento de la búsqueda de arquetipos con 4 - 1 arquetipos. Se repite todo lo previo, y de no hallarse resultados, se vuelve a inicializar la búsqueda de arquetipos con 4 - 2 arquetipos. Y así, 
# hasta que en el peor de los casos llegue a un arquetipo, total esto siempre va a ser calculable.

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

arch <- retry_archetypes(df_struc, n_arch = 6, seed_arch = 123, max_iter = 10, epsilon = 1e-6)

if (is.null(arch)) {
  message("Failed to compute archetypes with n_arch = 6. Retrying with reduced number of archetypes.")
  
  # Retry with reduced number of archetypes
  arch <- retry_archetypes(df_struc, n_arch = 5, seed_arch = 123, max_iter = 10, epsilon = 1e-6)
  
  if (is.null(arch)) {
    stop("Failed to compute archetypes with reduced number of archetypes as well.")
  }
}

# Store the result
arquetipos_prueba <- arch


#################################
### Clustering de  arquetipos ###
#################################

# Volvemos a armar los arquetipos como trayectorias
arquetipos <- arquetipos_prueba
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

# Lo convertimos en distancia
frechet_dist_matrix[is.na(frechet_dist_matrix)] <- 0
frechet_dist_matrix <- as.dist(frechet_dist_matrix)

# Perform hierarchical clustering
dendograma <- hclust(frechet_dist_matrix, method = "ward.D2", members = NULL)

par(mar = c(5, 5, 4, 2))  # Adjust the margins for better space
plot(dendograma, 
     main = "Fréchet-Clusters de Arquetipos", 
     hang = -1, 
     xlab = "Arquetipos", 
     ylab = "Alturas", 
     sub = "",
     axes = FALSE,           # Suppress default axes
     cex.main = 2,           # Title size
     cex.lab = 1.5,          # Axis labels size
     cex.axis = 1.2,         # Axis tick labels size
     cex = 1.5)              # Leaf label values size

# Add custom axis if needed (optional)
axis(side = 2)  # Add y-axis (Height)

# Primer plot: Año 2007 - Dendograma
rect.hclust(dendograma, k = 2, border = "red") # Modifico el tamaño al momento de bajarlo. Va a ser 800 x 800

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

###########################
### Plots de las curvas ###
###########################

mapas_path_deep <- Sys.getenv("MAPAS_PATH") # Esto está definido en DATA_PATH


# Segundo plot: Año 2007 - Totalidad de los datos, datos más profundos y sus arquetipos
tic()
for (k in anio:anio) { # Es 1:length(winter) en realidad, e hicimos 1:1 para que haga el año 2007
  titulo_w <- paste0("Trayectorias, Profundidades y Arquetipos")
  byers_map_w <- byers_map + ggtitle(titulo_w)
  
  q_max <- which(winter[[k]]$depth > quantile(winter[[k]]$depth, 0.60)) # El n% más profundo
  q_resto <- which(winter[[k]]$depth <= quantile(winter[[k]]$depth, 0.60)) # El resto
  n_arc <- length(trayectorias_arquetipicas)
  
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
  
  # Trayectorias arquetípicas
  for (l in 1:n_arc) {
    curve <- trayectorias_arquetipicas[[l]]
    curve$group <- "Arquetipos"
    byers_map_w <- byers_map_w + geom_path(data = curve, aes(x = lon, y = lat, color = group), alpha = 0.5)
  }
  
  # Apply the custom color scale and legend
  byers_map_w <- byers_map_w + 
    scale_color_manual(values = c(
      "Menos profundas" = "lightgrey", 
      "Más profundas" = "red", 
      "Arquetipos" = "black"
    )) +
    theme(
      plot.title = element_text(size = 40, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 30),
      legend.position = "right"
    )
  
  # Save the plot
  setwd(mapas_path_deep)
  png(filename = paste0(2005 + k - 1, "_depth_arc.png"), width = 800, height = 800)
  print(byers_map_w)
  dev.off()
}
toc()

# Tercer plot: Año 2007 - Totalidad de los datos y los arquetipos con sus respectivos clusters
tic()
for (k in anio:anio) { # Es 1:length(winter) en realidad, e hicimos 1:1 para que haga el año 2007
  titulo_w <- paste0("Arquetipos por clusters")
  byers_map_w <- byers_map + ggtitle(titulo_w)
  n_arc <- length(trayectorias_arquetipicas)
  clusters = cutree(dendograma, k = 2)
  renom_clusters = paste0("Cluster ", clusters)  # Assuming 'clusters' are your cluster identifiers
  
  # Curvas totales
  for (j in 1:length(winter[[k]]$curves)) {
    curve <- winter[[k]]$curves[[j]]
    curve$group <- "Curvas totales"  # Ensure 'group' is a column in the data
    byers_map_w <- byers_map_w + geom_path(data = curve, aes(x = lon, y = lat, color = group), alpha = 0.5)
  }
  
  # Archetypes with cluster names in the legend
  for (l in 1:n_arc) {
    curve <- trayectorias_arquetipicas[[l]]
    curve$group <- "Arquetipo"  # Ensure 'group' is a column in the data
    curve$cluster <- renom_clusters[l]  # Assign cluster names instead of numeric values
    
    # Add paths for each archetype, colored by its cluster name
    byers_map_w <- byers_map_w + 
      geom_path(data = curve, aes(x = lon, y = lat, color = as.factor(cluster)), alpha = 0.7)
  }
  
  # Apply custom color scale and legend for archetypes based on cluster names
  byers_map_w <- byers_map_w + 
    scale_color_manual(values = c(
      "Curvas totales" = "lightgrey", 
      "Arquetipo" = "black",  # Set a default color for archetypes if needed
      "Cluster 1" = "darkred", 
      "Cluster 2" = "darkgreen"
    )) +
    theme(
      plot.title = element_text(size = 40, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 30),
      legend.position = "right"
    )
  
  # Save the plot
  setwd(mapas_path_deep)
  png(filename = paste0(2005 + k - 1, "_total_arc.png"), width = 800, height = 800)
  print(byers_map_w)
  dev.off()
}
toc()


# Ahora debiéramos computarle la distancia de Fréchet a todas las curvas, y clusterizarlas en función de a cual se encuentran
# más cerca

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


# Cuarto plot: Año 2007 - Asignación de clusters final

tic()
for (k in anio:anio) {  # Es 1:length(winter) en realidad, e hicimos 1:1 para que haga el año 2007
  byers_map_w = byers_map
  titulo_w = paste0("Asignación Final")
  byers_map_w = byers_map_w + ggtitle(titulo_w)
  
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
  png(filename = paste0(2005 + k - 1, "_asignacion_cluster_final.png"), width = 800, height = 800)
  print(byers_map_w)
  dev.off()
}
toc()

