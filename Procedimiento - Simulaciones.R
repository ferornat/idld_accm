################################
### Instalamos las librerías ###
################################

# Para evaluar clasificación contra los clusters reales
install.packages("clue") # Ésta es para la función de reasignado de etiquetas usando la asignación húngara y hacer CCR
install.packages("fossil") # Esto es para usar el RI
install.packages("mclust") # Estp es para usar el ARI

# Trabajamos con la librería de arquetipos. Éstas son algunos links de relevancia:
# https://cran.r-project.org/web/packages/archetypal/vignettes/archetypal.html
# https://cran.r-project.org/web/packages/archetypal/archetypal.pdf
# https://cran.r-project.org/web/packages/archetypal/index.html
# install.packages("archetypal")
install.packages("archetypes") # https://cran.r-project.org/web/packages/archetypes/vignettes/archetypes.pdf

# install.packages("grDevices") # Esta sólo la usamos cuando recién probábamos arquetipos, pero ya no para la parte de
# clustering que después nos interesa
install.packages("pracma")
install.packages("mvnfast")
install.packages("mvtnorm")

install.packages("microbenchmark") # Ésta es para hacer las simulaciones con los arquetipos

# Instalamos la librería del algoritmo al cual le tenemos que ganar
remotes::install_github("lfernandezpiana/idld")

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

# library("grDevices")# Esta sólo la usamos cuando recién probábamos arquetipos, pero ya no para la parte de
# clustering que después nos interesa
library("mclust")
library("mvtnorm") # Esto es para poder usar "rmrnorm".

# La librería geometry es la que me va a permitir generar las coberturas convexas para los arquetipos que genere.
# https://www.rdocumentation.org/packages/geometry/versions/0.4.7/topics/convhulln # Documentación de convhull() en la librería 'geometry'.
# http://www.qhull.org/ # Documentación de la librería 'Qhull', a la cual 'geometry' le monta la interfase.
require('geometry') # Es para definir la función de distancia


# La librería de Lucas y Marcela está albergada en el siguiente repositorio. Allí están tanto los algoritmos de cluster como los de proyecciones de profundidades, bootstrap y demás:
# https://github.com/lfernandezpiana/idld/tree/master/R
library("idld") # Librería del algoritmo a ganar

library("sessioninfo") # Vamos a chequear versiones de paquetes utilizados

# Vemos las versiones de las librerías que estamos usando
# session_info(pkgs = c("idld","microbenchmark","mvtnorm","mvnfast",
#                       "mclust","pracma","archetypes","fossil",
#                       "clue"),
#              to_file = TRUE)

##################
### Arquetipos ###
##################

# Hacemos algunas pruebas para ver el tiempo de cómputo

### Base en 2D

p1 = c(1,2)
p2 = c(3,5)
p3 = c(7,3)
dp = rbind(p1,p2,p3)

set.seed(9102)
pts=t(sapply(1:100, function(i,dp){
  cc=runif(3)
  cc=cc/sum(cc)
  colSums(dp*cc)
},dp))

df = data.frame(pts)
colnames(df)=c("x","y")

tic()
set.seed(1234)
borrar_arch = archetypes(df, 3, verbose = TRUE)
toc()

suppressWarnings(plot(dp, colour = 'black', pch = 19))
points(df)
points(parameters(borrar_arch), col = 'red', pch = 19)

# Lo hacemos en GGPLOT
dp_df <- as.data.frame(dp)
colnames(dp_df) <- c("x", "y")

archetypes_df <- as.data.frame(parameters(borrar_arch))
colnames(archetypes_df) <- c("x", "y")

dp_df$type <- "Originales"  
archetypes_df$type <- "Estimados"

legend_data <- rbind(dp_df, archetypes_df)

ggplot() +
  geom_point(data = legend_data, aes(x = x, y = y, color = type), shape = 19, size = 3) +
  geom_point(data = df, aes(x = x, y = y), alpha = 0.7) +
  scale_color_manual(values = c("Arquetipos originales" = "blue",
                                "Arquetipos estimados" = "red")) +
  labs(title = "Obtención de Arquetipos en 2D",
       x = "X-axis",
       y = "Y-axis",
       color = "Arquetipos") +  
  theme_minimal()

### Múltiples arquetipos

# Función para muestra aleatoria
generate_data <- function(n) {
  X <- rnorm(n, mean = 0, sd = 1)
  Y <- rnorm(n, mean = 0, sd = 1)
  data.frame(X = X, Y = Y)
}

borrar = generate_data(700)

# Vemos tiempo de cómputo. Va a tardar ~0.4 seg
tic()
set.seed(1234)
borrar_arch = archetypes(borrar, 3, verbose = TRUE)
toc()

# Probamos lo que es la iteración para múltiples cantidades de arquetipos, pero usando una pseudo-matriz inversa de Monroe-Penrose
# Este igual puede romper. Va a tardar ~75 seg
tic()
borrar_arch = stepArchetypes(data = borrar,
                             k = 1:10,
                             verbose = F,
                             family = archetypesFamily("original",
                                                       zalphasfn = archetypes:::ginv.zalphasfn),
                             nrep = 10)
toc()

screeplot(borrar_arch)

### ¿Qué pasa si probamos con muchas dimensiones?

generate_data_bis <- function(n, dimensions) {
  data <- matrix(rnorm(n * dimensions), ncol = dimensions)
  colnames(data) <- paste0("V", 1:dimensions)
  as.data.frame(data)
}

borrar = generate_data_bis(700,20)

# No tarda un tiempo excesivo, siempre que no rompa por no poder invertir las matrices.
tic()
set.seed(1234)
borrar_arch = archetypes(borrar,
                         20,
                         verbose = TRUE,
                         family = archetypesFamily("robust"))
toc()

tic()
set.seed(1234)
borrar_arch = archetypes(borrar, 3, verbose = TRUE)
toc()

# Aquí directamente hicimos la implementación con la librería "archetypes", pero en realidad tuvimos muchas pruebas con la librería "archetypal",
# la cuál es más robusta para obtener resultados, aunque es excesivamente lenta en comparación con la primera. Además, esta última trae muchos
# atributos y resultados que, a nuestros fines prácticos, no nos interesan ya que únicamente necesitamos las representaciones arquetípicas.

###########################################################################################
### Comenzamos simulaciones para poner a prueba dimensiones y cantidad de observaciones ###
###########################################################################################

# Definimos:
# - Tamaños de las muestras
# - # de arquetipos a calcular
# - Semilla para garantizar reproductibilidad

sample_sizes <- c(300, 500, 1000, 5000, 10000)
archetypes_list <- c(2,3,4,5,6,7,8,9,10,15,20)
seed <- 123

# Función para medir el tiempo con manejo de advertencias y errores
measure_time <- function(df, num_archetypes, seed) {
  time_result <- tryCatch({
    times <- microbenchmark(
      {
        set.seed(seed)
        # archetypal(df = df, kappas = num_archetypes, verbose = FALSE, rseed = seed, save_history = FALSE) # "Verbose" es para que te haga print de los resultados
        archetypes(data = df, k = num_archetypes, verbose = FALSE) # Esta es la alternativa que estábamos probando
        # archetypes(data = df, k = num_archetypes, verbose = FALSE, family = archetypesFamily("original", zalphasfn = archetypes:::ginv.zalphasfn)) # Usamos la pseudoinversa de Moore-Penrose
      },
      times = 1
    )
    median(times$time) / 1e9  # De nanosegundos a segundos
  }, warning = function(w) {
    message("Warning encountered: ", w)
    NA  # Devuelve NA si se encuentra una advertencia
  }, error = function(e) {
    message("Error encountered: ", e)
    NA  # Devuelve NA si se encuentra un error
  })
  
  return(time_result)
}

### Inicializamos el loop

tic()
results <- data.frame()

# Combinamos la muestra con la cantidad de arquetipos
for (n in sample_sizes) {
  df <- generate_data(n)
  for (arc in archetypes_list) {
    cat("// Muestra: ", n, " - ", "Cantidad de Arquetipos: ", arc, "Errores: ")
    time_taken <- measure_time(df, arc, seed)
    results <- rbind(results, data.frame(SampleSize = n, Archetypes = arc, Time = time_taken))
  }
}
# Esto con la librería archetypal tardaba ~5065.692 segundos. Esto fue con 7 variantes en cantidades de arquetipos.
# Ahora me tarda ~233.801 segundos, pero con muchos errores de matrices singulares en el medio. Esto fue con 11 variantes en cantidades de arquetipos.
# En el caso de usar la pseudo-inversa de Moore-Penrose, me tarda ~238.879 segundos (familia de arquetipos "original")

# (https://www.researchgate.net/post/What_is_the_most_efficient_way_for_computing_Archetypal_Analysis)
toc()

# Trabajamos un poco sobre los datos de resultados
arch_computables <- results %>%
  group_by(SampleSize) %>%
  filter(!is.na(Time)) %>%
  summarize("Máximo de Arquetipos Computables" = max(Archetypes, na.rm = TRUE))

# Esto es usando la matriz inversa . Usando la pseudo inversa de Moore-Penrose la verdad es que no cambiaba tanto, y tardaba incluso un poco más de 
# tiempo

# Graficamos los resultados
options(warn = -1) # Silenciamos además los warnings
options(repr.plot.width = 25, repr.plot.height = 8)

ggplot(results, aes(x = Archetypes, y = Time, color = as.factor(SampleSize))) +
  geom_line(size = 1.5) +  # Increase line thickness
  geom_point(size = 3) +   # Increase point size
  labs(title = "Tiempo de Cómputo para Diferentes Cantidades de Arquetipos",
       x = "Cantidad de Arquetipos",
       y = "Tiempo de Cómputo (segundos)",
       color = "Tamaño de muestra") +
  coord_cartesian(xlim = c(2, 8)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# Vemos arquetipos computables
ggplot(arch_computables, aes(x = arch_computables$`Máximo de Arquetipos Computables`, y = factor(SampleSize))) +
  geom_bar(stat = "identity", aes(fill = factor(SampleSize))) +
  labs(x = "Arquetipos Computables", y = "Tamaño de la muestra") +
  theme_minimal() +
  scale_x_continuous(breaks = 0:8) +  # Define los valores del eje X
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"  # Elimina la leyenda
  )

## Ahora vamos a modificar la cantidad de dimensiones, dejando la cantidad de observaciones fija en 300

# Funciones para generar la información y además para medir el tiempo
generate_data_bis <- function(n, dimensions) {
  data <- matrix(rnorm(n * dimensions), ncol = dimensions)
  colnames(data) <- paste0("V", 1:dimensions)
  as.data.frame(data)
}

# Definimos:
# - Tamaños de las muestras
# - # de arquetipos a calcular
# - Cantidad de dimensiones
# - Semilla para garantizar reproductibilidad

sample_size <- 300
archetypes_list <- c(2,3,4,5,6,7,8,9,10,15,20)
dimensions_list <- c(2, 3, 4, 5, 10, 20,5,100)
seed <- 123

# Inicializamos el loop
tic()
set.seed(seed)
results <- data.frame()

# Cruzamos las dimensiones con la cantidad de arquetipos
for (d in dimensions_list) {
  df <- generate_data_bis(sample_size, d)
  for (arc in archetypes_list) {
    cat("// Dimensión: ", d, " - ", "Cantidad de Arquetipos: ", arc, "Errores: ")
    time_taken <- measure_time(df, arc, seed)
    results <- rbind(results, data.frame(Dimensions = d, Archetypes = arc, Time = time_taken))
  }
}
# Esto con la librería archetypal tardaba ~3391.382 segundos . Esto fue con 7 variantes en cantidades de arquetipos.
# Ahora me tarda ~88.573 segundos, pero con muchos errores de matrices singulares en el medio. Esto fue con 11 variantes en cantidades de arquetipos.
# (https://www.researchgate.net/post/What_is_the_most_efficient_way_for_computing_Archetypal_Analysis)
toc()

arch_computables <- results %>%
  group_by(Dimensions) %>%
  filter(!is.na(Time)) %>%
  summarize("Máximo de Arquetipos Computables" = max(Archetypes, na.rm = TRUE))

# Graficamos los resultados
options(warn = -1) # Silenciamos además los warnings
options(repr.plot.width = 25, repr.plot.height = 8)

ggplot(results, aes(x = Archetypes, y = Time, color = as.factor(Dimensions))) +
  geom_line() +
  geom_point() +
  labs(title = "Tiempo de Cómputo para Diferentes Cantidades de Arquetipos",
       x = "Cantidad de Arquetipos",
       y = "Tiempo de Cómputo (segundos)",
       color = "Cantidad de dimensiones") +
  theme_minimal() +
  scale_x_continuous(breaks = 0:17) +  # Define los valores del eje X
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# Vemos arquetipos computables
ggplot(arch_computables, aes(x = arch_computables$`Máximo de Arquetipos Computables`, y = factor(Dimensions))) +
  geom_bar(stat = "identity", aes(fill = factor(Dimensions))) +
  labs(x = "Arquetipos Computables", y = "Dimensiones en la muestra") +
  theme_minimal() +
  scale_x_continuous(breaks = 0:17) +  # Define los valores del eje X
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"  # Elimina la leyenda
  )

# Vemos que en muchos casos, al tratarse de un procedimiento que detrás está computando inversa de matrices, cuando llegamos a 
# elevar mucho la cantidad de arquetipos en un contexto donde las observaciones son muy parecidas terminamos recayendo en una 
# matriz prácticamente singular. Es por es eso que pensamos, en función de los gráficos de arriba en una heurística que trate 
# de ser estable en cuando varía tanto el tamaño de la muestra como las dimensiones.

# Función para muestra aleatoria en dimensiones arbitrarias
generate_data_bis_2 <- function(n, dims) {
  data <- as.data.frame(matrix(rnorm(n * dims), ncol = dims))
  colnames(data) <- paste0("V", 1:dims)
  return(data)
}

# Definimos:
# - Tamaño de la muestra (fijo)
# - # de dimensiones
# - Semilla para garantizar reproductibilidad

sample_size <- 300
dimensions_list <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 40, 50, 60) # Hasta 100 no hay chances de que corra
seed <- 123

#########################################
### Probamos estrategia dimensión + 1 ###
#########################################

# Inicializamos el loop
tic()
set.seed(seed)
results <- data.frame()

# Combinamos la cantidad de dimensiones y arquetipos
for (dims in dimensions_list) {
  df <- generate_data_bis_2(sample_size, dims)
  archetypes <- dims + 1  # "x" dimesiones  y "x+1" arquetipos
  cat("// Dimensión: ", dims, " - ", "Cantidad de Arquetipos: ", archetypes, "Errores: ")
  time_taken <- measure_time(df, archetypes, seed)
  results <- rbind(results, data.frame(Dimensions = dims, Archetypes = archetypes, Time = time_taken))
}
# Esto antes tardaba ~321.224 segundos, usando la librería "archetypal"
# Ahora tarda ~8.61 segundos.
toc()

# Graficamos los resultados
options(warn = -1) # Silenciamos además los warnings
options(repr.plot.width = 25, repr.plot.height = 8)

ggplot(results, aes(x = Dimensions, y = Time)) +
  geom_line(size = 1.5) +  # Increase line thickness
  geom_point(size = 3) +   # Increase point size
  labs(title = "Tiempo de Cómputo para Diferentes Cantidades de Dimensiones y Arquetipos",
       x = "Cantidad de Dimensiones (d), con cantidad Arquetipos fija (d+1)",
       y = "Tiempo de Cómputo (segundos)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

#######################################################
### Probamos estrategia de 3 arquetipos por default ###
#######################################################

# Inicializamos el loop
tic()
set.seed(seed)
results <- data.frame()

# Combinamos la cantidad de dimensiones y arquetipos
for (dims in dimensions_list) {
  df <- generate_data_bis_2(sample_size, dims)
  archetypes <- 3  # Esto lo sacamos de haber observado ambos gráficos arriba por separado
  cat("// Dimensión: ", dims, " - ", "Cantidad de Arquetipos: ", archetypes, "Errores: ")
  time_taken <- measure_time(df, archetypes, seed)
  results <- rbind(results, data.frame(Dimensions = dims, Archetypes = archetypes, Time = time_taken))
}
# Tarda ~18 segundos.
toc()

# Graficamos los resultados
options(warn = -1) # Silenciamos además los warnings
options(repr.plot.width = 25, repr.plot.height = 8)

ggplot(results, aes(x = Dimensions, y = Time)) +
  geom_line(size = 1.5) +  # Increase line thickness
  geom_point(size = 3) +   # Increase point size
  labs(title = "Tiempo de Cómputo para Diferentes Cantidades de Dimensiones y Arquetipos",
       x = "Cantidad de Dimensiones (d), con cantidad Arquetipos fija (3)",
       y = "Tiempo de Cómputo (segundos)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

#####################################
### Probamos una nueva estrategia ###
#####################################

# Arrancamos con 4 arquetipos.
# 
# Si funciona, continuamos con el procedimiento.
# 
# Si no funciona, volvemos a intentar con otra semilla. Y así sucesivamente n-veces. Si en ninguna de ellas se hallan los 
# arquetipos solicitados, se vuelve a iniciar el procedimiento de la búsqueda de arquetipos con 4 - 1 arquetipos. Se repite 
# todo lo previo, y de no hallarse resultados, se vuelve a inicializar la búsqueda de arquetipos con 4 - 2 arquetipos. Y 
# así, hasta que se llegue a un arquetipo.

measure_time_propuesta <- function(df, num_archetypes, initial_seed, max_retries = 5, epsilon = 1e-6) {
  final_time <- NA
  all_seeds <- list()  # To store seeds used for each attempt
  num_archetypes_orig <- num_archetypes
  
  while (num_archetypes > 0) {
    seeds_used <- c()
    time_result <- NA
    retries <- 0
    
    while (retries < max_retries) {
      set.seed(initial_seed + retries)
      seeds_used <- c(seeds_used, initial_seed + retries)  # Store the seed
      
      # Try fitting archetypes with retry logic and noise handling
      archetypes_result <- NULL
      for (i in 0:max_retries) {
        set.seed(initial_seed + retries + i)  # Adjust seed for noise attempts
        
        # Try without noise
        archetypes_result <- tryCatch({
          archetypes(data = df, k = num_archetypes, verbose = FALSE)
        }, error = function(e) {
          NULL  # Return NULL if an error occurs
        })
        
        if (!is.null(archetypes_result)) break  # Break if successful
        
        # If failed, add noise and retry
        message("Attempt ", i + 1, ": Adding noise and retrying...")
        df_noisy <- df + matrix(rnorm(length(df), mean = 0, sd = epsilon), nrow = nrow(df))
        
        archetypes_result <- tryCatch({
          archetypes(data = df_noisy, k = num_archetypes, verbose = FALSE)
        }, error = function(e) {
          NULL  # Return NULL if an error occurs
        })
        
        if (!is.null(archetypes_result)) break  # Break if successful
      }
      
      if (!is.null(archetypes_result)) {
        # Measure time only if successful
        time_result <- tryCatch({
          times <- microbenchmark(
            {
              archetypes(data = df, k = num_archetypes, verbose = FALSE)
            },
            times = 1
          )
          median(times$time) / 1e9
        }, error = function(e) {
          NA  # Return NA if an error occurs
        })
        break  # Break retry loop if time measurement succeeds
      }
      
      retries <- retries + 1
    }
    
    all_seeds[[as.character(num_archetypes)]] <- seeds_used  # Store the seeds used for this num_archetypes
    
    if (!is.na(time_result)) {
      final_time <- time_result
      break  # Break the num_archetypes loop if a result is achieved
    }
    
    num_archetypes <- num_archetypes - 1  # Reduce the number of archetypes and try again
  }
  
  return(list(time = final_time, seeds = all_seeds, archetypes = num_archetypes, num_archetypes_deseados = num_archetypes_orig))
}

# Chequeamos dimensiones
paste0("Recordamos tamaños: ")
paste0(sample_sizes)
paste0("Recordamos arquetipos: ")
paste0(archetypes_list)


# Ponemos algo más chico para ver que funcione
sample_sizes <- c(300,500)
archetypes_list <- c(2, 3, 4, 5, 10)
seed <- 123


# Inicializamos el loop
tic()
results <- data.frame()

# Combinamos la muestra con la cantidad de arquetipos
for (n in sample_sizes) {
  df <- generate_data(n)
  for (arc in archetypes_list) {
    cat("// Muestra: ", n, " - ", "Cantidad de Arquetipos: ", arc, "Errores: ")
    time_taken <- measure_time_propuesta(df, arc, seed, max_retries = 10)
    
    results <- rbind(results, data.frame(sample_size = n, arc_obtenidos = time_taken$archetypes, arc_deseados = time_taken$num_archetypes_deseados, tiempo = time_taken$time))
  }
}
toc()

# Vemos los resultados
options(warn = -1) # Silenciamos además los warnings
options(repr.plot.width = 25, repr.plot.height = 8)

ggplot(results, aes(x = arc_deseados, y = arc_obtenidos)) +
  geom_line(size = 1.5) +  # Increase line thickness
  geom_point(size = 3) +   # Increase point size
  labs(title = "Arquetipos Deseados vs Obtenidos",
       x = "Cantidad de Arquetipos Deseados",
       y = "Cantidad de Arquetipos Obtenidos") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

############################################
### Función de clustering: cluster_alloc ###
############################################

# **cluster_alloc()**
# Lo que vamos a hacer ahora es hacer una comparativa entre la estrategia de clustering actual y la propuesta con arquetipos, 
# comparando sus tiempos de cómputos y errores de clasificación para algún dataset simulado.

# library(cluster)
# library(geometry)

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

cluster_alloc_fer <- function(X, ld, alpha_q = seq(0.5, 0.9, 0.1), K, distance = "euclidean", n_arch = 3, seed_arch = 9102, max_iter = 10) {
  d <- dim(X)
  alpha_length <- length(alpha_q)
  cluster <- matrix(NA, nrow = d[1], ncol = alpha_length)
  archetypes <- vector("list", length = K * alpha_length)
  
  # Vamos a utilizar la versión final de retry_archetypes(), que incluye iteraciones y ruido en la inicialización
  
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

##########################################
### cluster_alloc() de Lucas y Marcela ###
##########################################

# Estas son las funciones de distancia que establecieron Lucas y Marcela

#### MAHALANOBIS DIST
mahalanobis_dist = function(x, data, sigma=var(data)) {
  if (any(class(data)=="numeric")) {
    #print("Sigma can't be approximated if nrow(data)=1")
    return(NA)
  }
  else {
    n = nrow(as.matrix(data))
    cond_number = rcond(sigma)
    if (cond_number < 0.01) {
      #message = paste0(
      #                paste("Condition number is:", cond_number),
      #                paste(". Matrix is ill-conditioned, return NA")
      #                )
      #print(message)
      return(NA)
    } else {
      sigma_inv = solve(sigma)
      md = numeric(n)
      for (i in 1:n) md[i] = (x-data[i,])%*%sigma_inv%*%as.matrix(x-data[i,])
      return(md)
    }
  }
}

#### EUCLIDEAN DIST
euclidean_dist = function(x, data) {
  if (any(class(data)=="numeric")) data = matrix(data, nrow=1, ncol=length(data))
  n = nrow(data)
  ed = numeric(n)
  for (i in 1:n) ed[i] = (x-data[i,])%*%as.matrix(x-data[i,])
  return(ed)
}

# La función propiamente dicha

#################################
#### CLUSTER ALLOCATION VER2 ####
#################################

# source("distances.R") # Esto lo tenemos definido acá arriba

cluster_alloc = function(X, ld, alpha_q=seq(0.5,0.9,0.1), K, distance="euclidean") {
  #INPUT
  #X: data that will be assigned local depth. Matrix each row is an observation.
  #ld: local depth calculated for the data.
  #alpha_q: a grid of values between 0.5 and 0.9 to select the percentage of deepest points for
  #                the core region. For example, alpha_quantile = 0.9 implies that the core region
  #                contains the ten percent of deepest points.
  #K: number of clusters.
  #distance: euclidean or mahalanobis.
  #verbose: if TRUE prints the progress
  #OUTPUT
  #clusters: clusters allocation.
  
  d = dim(X)
  alpha_length = length(alpha_q)
  cluster = matrix(NA, nrow=d[1], ncol=alpha_length)
  
  for (j in 1:alpha_length) {
    # Step 1: core
    region_index = which(ld>quantile(ld,alpha_q[j]))
    others_index = which(ld<=quantile(ld,alpha_q[j]))
    if (length(region_index)<=K+1) {
      #print("Warning: the number of core points is less than the number of clusters")
      next
    }
    region_cluster = kmeans(X[region_index,], centers=K)$cluster
    cluster[region_index,j] = region_cluster
    
    # Aquí habría que calcular los arquetipos
    
    # Step 2: distance. # Esta distancia sería luego al arquetipo
    if (distance=="euclidean") {
      for (i in others_index) {
        dis = numeric(K)
        for (k in 1:K) dis[k] = min(euclidean_dist(x=X[i,], data=X[region_index[region_cluster==k],]))
        cluster[i,j] = which.min(dis)
      }
    }
    else if(distance=="mahalanobis") {
      for (i in others_index) {
        dis = numeric(K)
        for (k in 1:K) {
          md = mahalanobis_dist(x=X[i,], data=X[region_index[region_cluster==k],])
          if (any(is.na(md))) {
            cluster[i,j] = NA
          } else {
            dis[k] = min(md)
            cluster[i,j] = which.min(dis)
          }
        }
      }
    }
    else {
      print("Distance parameter invalid")
    }
  }
  return(cluster)
}

# Tenemos en particular tres ejemplos:
# - **Ejemplo 5**: Tibshirani retomado por salibian con contaminacion en la variable informativa, esférico, sin ruido.
# - **Ejemplo 8**: Tibshirani retomado por salibian con contaminacion en la variable informativa, elipsoidal, con ruido sin Paindavaine.
# - **Ejemplo 10**: Tibshirani retomado por Salibian con contaminacion en la variable noinformativa, elipsoidal, sin ruido.
# 
# Tomamos el caso 5, porque es el que veníamos trabajando. Pero da lo mismo cualquiera de ellos.

## Ejemplo 5
# Ejemplo tibshirani retomado por salibian con contaminacion en la variable informativa,
# esférico, sin ruido

n = 300 # tamaño muestral
p = 3 # nro de variables no informativas
mu = 3
k = 3
sig = diag(rep(1,p),p)
NN = 500

classgood=c(rep(1,n/3),rep(2,n/3),rep(3,n/3))

set.seed(1234)
X1 = rmvnorm(n/3,c((-1)*mu,(-1)*mu,rep(0,p-2)),sigma=sig)
X2 = rmvnorm(n/3,rep(0,p),sigma=sig)
X3 = rmvnorm(n/3,c(mu,mu,rep(0,p-2)),sigma=sig)

X = rbind(X1,X2,X3)
X[1:5,1]=rep(25,5)+runif(5)*0.01 ######agrego outliers en la variable informativa

##

# ## Ejemplo 8
# # Ejemplo Tibshirani retomado por salibian
# # con contaminacion en la variable informativa, elipsoidal, con ruido
# # sin Paindavaine

# n = 300
# p=5
# mu=3
# k=3
# sig=diag(c(4,1/3,rep(1,p-2)),p)
# NN=500

# classgood=c(rep(1,n/3),rep(2,n/3),rep(3,n/3))

# set.seed(1234)
# X1 = rmvnorm(n/3,c((-1)*mu,(-1)*mu,rep(0,p-2)),sigma=sig)
# X2 = rmvnorm(n/3,rep(0,p),sigma=sig)
# X3 = rmvnorm(n/3,c(mu,mu,rep(0,p-2)),sigma=sig)

# X = rbind(X1,X2,X3)
# X[1:5,1]=25+runif(5)*0.01

##

# ## Ejemplo 10

# # Ejemplo Tibshirani retomado por Salibian con
# # contaminacion en la variable noinformativa, elipsoidal, sin ruido

# n = 300
# p=3
# mu=3
# k=3
# sig=diag(c(4,1/3,rep(1,p-2)),p)
# NN=500

# classgood=c(rep(1,n/3),rep(2,n/3),rep(3,n/3))

# set.seed(1234)
# X1 = rmvnorm(n/3,c((-1)*mu,(-1)*mu,rep(0,p-2)),sigma=sig)
# X2 = rmvnorm(n/3,rep(0,p),sigma=sig)
# X3 = rmvnorm(n/3,c(mu,mu,rep(0,p-2)),sigma=sig)

# X = rbind(X1,X2,X3)
# X[1:5,p]=25+runif(5)*0.01
# #plot(X,col=c(rep(1,100),rep(2,100),rep(3,100)))

## Hacemos la comparativa de errores de cómputo y clasificación.

# Fijamos parámetros y profundidades iniciales, que no las vamos a cambiar
# Esto va a ser insumo de ambos criterios de clasificación

tic()
beta = 0.3
m = 500
alpha_q = c(0.5)
K = 3
distance = "euclidean"

ld_lucmar = idld_m(X,X,beta,m,verbose=FALSE)
toc()

# Algoritmo original de Lucas y Marcela:
tic()
cluster_original_lucmar = cluster_alloc(X = X, ld = ld_lucmar, alpha_q = alpha_q, K = K, distance = distance)
toc()

plot(X[, 1], X[, 2], col = cluster_original_lucmar[,1], pch=19, main = "Data clusterizada con algoritmo cluster_alloc() de Lucas y Marcela")
colors = rainbow(K)

## Vemos métricas de acierto
# Asignaciones de clúster para el primer alpha_q
pred_clusters_lucmar = cluster_original_lucmar[, 1]

# Matriz de confusión
conf_matrix_lucmar = table(classgood, pred_clusters_lucmar)
conf_matrix_lucmar
(94+91+100)/sum(conf_matrix_lucmar) # 0.95


# Esto se hace dificil de automatizar, puesto que tenemos que probar qué algoritmo funciona bien para reasignar las etiquetas. 
# El algoritmo húngaro, por ejemplo, no funciona bien, porque en vez de darme el 95% de arriba me da 3%.

# # Crear una matriz de confusión
# cm <- table(classgood, pred_clusters_lucmar)
# 
# # Usar el algoritmo húngaro para encontrar la mejor correspondencia de clusters
# assignment <- solve_LSAP(cm, maximum = TRUE)
# 
# # Asignar los clusters según la correspondencia encontrada
# aligned_cluster_labels <- sapply(pred_clusters_lucmar, function(x) assignment[x])
# 
# # Calcular el CCR
# ccr <- sum(aligned_cluster_labels == classgood) / length(classgood)
# print(paste("Correct Classification Rate (CCR):", ccr))
# 
# Limpiamos esto que entonces no vamos a usar
# rm(cm, assignment, aligned_cluster_labels, ccr)

### Hacemos entonces el Rand Index (RI) y el Adjusted Rand Index (ARI).

# Calcular el índice de Rand (RI)
ri <- rand.index(classgood, pred_clusters_lucmar)
print(paste("Rand Index (RI):", ri))

# Calcular el índice Rand ajustado (ARI)
ari <- adjustedRandIndex(classgood, pred_clusters_lucmar)
print(paste("Adjusted Rand Index (ARI):", ari))

## Algoritmo de Fer:
tic()
cluster_original_fer = cluster_alloc_fer(X = X, ld = ld_lucmar, alpha_q = alpha_q, K = K, distance = distance, n_arch = 4, seed_arch = 9102)
# Esto antes, con la versión 1, tardaba ~36.181 segundos
# Ahora tarda ~0.194 segundos
toc()

plot(X[, 1], X[, 2], col = cluster_original_fer$cluster[,1], pch = 19, main = "Data clusterizada con algoritmo cluster_alloc_fer()")
colors = rainbow(K)

# Vemos entonces el RI y el ARI.

# Calcular el índice de Rand (RI)
ri <- rand.index(classgood, cluster_original_fer$cluster)
print(paste("Rand Index (RI):", ri))

# Calcular el índice Rand ajustado (ARI)
ari <- adjustedRandIndex(classgood, cluster_original_fer$cluster)
print(paste("Adjusted Rand Index (ARI):", ari))

#################################################################################################################################
### Ahora queremos ver entonces, cómo evolucionan tanto los tiempos de ejecución como las métricas de ajuste de los clusters. ###
#################################################################################################################################

# Armamos los datasets que vamos a querer comparar
# Hacemos el ejercicio con la estructura de muestra del ejemplo 5

# Función para generar el conjunto de datos
generate_data_ejemplo_5 <- function(n, p, mu, sig) {
  X1 <- rmvnorm(n/3, c((-1) * mu, (-1) * mu, rep(0, p - 2)), sigma = sig)
  X2 <- rmvnorm(n/3, rep(0, p), sigma = sig)
  X3 <- rmvnorm(n/3, c(mu, mu, rep(0, p - 2)), sigma = sig)
  
  X <- rbind(X1, X2, X3)
  X[1:5, 1] <- rep(25, 5) + runif(5) * 0.01 # Agregar outliers en la variable informativa
  
  # Crear las etiquetas de clase
  classgood <- c(rep(1, n/3), rep(2, n/3), rep(3, n/3))
  
  # Devolver los datos y las etiquetas como una lista
  return(list(data = X, class_labels = classgood))
}

# Parámetros
p <- 3
mu <- 3
sig <- diag(rep(1, p), p)
semilla <- 1234

# Vamos a hacer estas muestras
# sample_sizes <- c(300, 900, 1200, 3000, 6000, 9000, 15000, 18000) # Esto explota, + 6hs y no termina
sample_sizes <- c(300, 900, 3000, 9000, 18000) # Esto me mata la ejecución en local, y no me permite que me asignen un equipo en el Colab
# sample_sizes <- c(300, 900, 3000, 6000, 10000)

# Crear la lista de conjuntos de datos
datasets <- lapply(sample_sizes, function(n) {
  set.seed(semilla) # Para reproducibilidad
  generate_data_ejemplo_5(n, p, mu, sig)
})

# Acceso al output:
# datasets[[1]]$data
# datasets[[1]]$class_labels

tic()
# Inicialiamos un df vacío que vamos a ir llenando a medida que vaya avanzando el bucle
results <- data.frame(
  Dataset = character(),
  Method = character(),
  RI = numeric(),
  ARI = numeric(),
  ExecutionTime = numeric(),
  stringsAsFactors = FALSE
)

# Hacemos el bucle de cada base simulada
for (i in seq_along(datasets)) {
  dataset <- datasets[[i]]$data
  classgood <- datasets[[i]]$class_labels
  
  dataset_name <- paste("Dataset ", dim(dataset)[1]) # Le asignamos como nombre su cantidad de observaciones
  
  # Establecemos parámetros fijos para cada profundidad y clustering
  beta_loop = 0.3
  m_loop = 500
  alpha_q_loop = c(0.5)
  K_loop = 3
  distance = "euclidean"
  cantidad_de_arquetipos = 3
  semilla_arquetipos = 9102
  
  # Calculamos la profundidad
  ld_lucmar_loop = idld_m(dataset, dataset, beta_loop, m_loop, verbose = FALSE)
  
  
  # Corremos cluster_alloc_fer() para medir su tiempo y resultados
  start_time_fer <- Sys.time()
  cluster_original_fer <- cluster_alloc_fer(X = dataset, ld = ld_lucmar_loop, alpha_q = alpha_q_loop, K = K_loop, distance = distance, n_arch = cantidad_de_arquetipos, seed_arch = semilla_arquetipos)
  end_time_fer <- Sys.time()
  execution_time_fer <- as.numeric(difftime(end_time_fer, start_time_fer, units = "secs"))
  
  # Corremos cluster_alloc() para medir su tiempo y resultados
  start_time_lucmar <- Sys.time()
  cluster_original_lucmar <- cluster_alloc(X = dataset, ld = ld_lucmar_loop, alpha_q = alpha_q_loop, K = K_loop, distance = distance)
  end_time_lucmar <- Sys.time()
  execution_time_lucmar <- as.numeric(difftime(end_time_lucmar, start_time_lucmar, units = "secs"))
  
  # Computamos el RI y el ARI
  ri_fer <- rand.index(classgood, cluster_original_fer$cluster)
  ari_fer <- adjustedRandIndex(classgood, cluster_original_fer$cluster)
  
  ri_lucmar <- rand.index(classgood, cluster_original_lucmar[,1])
  ari_lucmar <- adjustedRandIndex(classgood, cluster_original_lucmar[,1])
  
  # Guardamos resultados
  results <- rbind(results, data.frame(
    Dataset = dataset_name,
    Method = "fer",
    RI = ri_fer,
    ARI = ari_fer,
    ExecutionTime = execution_time_fer
  ))
  
  results <- rbind(results, data.frame(
    Dataset = dataset_name,
    Method = "lucmar",
    RI = ri_lucmar,
    ARI = ari_lucmar,
    ExecutionTime = execution_time_lucmar
  ))
}

# Convertir la columna Dataset a factor con el orden deseado
results$Muestra <- factor(results$Dataset, levels = unique(results$Dataset))

# Esto tardaba > 20_000 segundos (~ 6 horas). En mi computadora local tarda ~2.7236 horas
toc()

results

# Dejamos el print para no tener que correr todo de nuevo.

# Dataset Method        RI       ARI ExecutionTime        Muestra
# 1    Dataset  300    fer 0.9311260 0.8448292    0.09111452   Dataset  300
# 2    Dataset  300 lucmar 0.9395987 0.8638254    0.07613993   Dataset  300
# 3    Dataset  900    fer 0.9391843 0.8631642    0.27689409   Dataset  900
# 4    Dataset  900 lucmar 0.9459177 0.8782419    0.66091466   Dataset  900
# 5   Dataset  3000    fer 0.9340465 0.8518516    0.57539916  Dataset  3000
# 6   Dataset  3000 lucmar 0.9448214 0.8759981    7.75903654  Dataset  3000
# 7   Dataset  9000    fer 0.9288728 0.8402979    0.97625661  Dataset  9000
# 8   Dataset  9000 lucmar 0.9544236 0.8975275   69.22658205  Dataset  9000
# 9  Dataset  18000    fer 0.9099782 0.7978114    6.28738475 Dataset  18000
# 10 Dataset  18000 lucmar 0.9566967 0.9026260  297.80853844 Dataset  18000

# Modificamos para poder hacer los plots
results$Muestra <- factor(results$Dataset, levels = unique(results$Dataset))
names(results)[names(results) == 'Method'] <- 'Estrategia'

###########################
### Comparación gráfica ###
###########################

# Comparación de RI
options(warn = -1) # Silenciamos además los warnings
options(repr.plot.width = 25, repr.plot.height = 8)

# Crear el gráfico de barras con etiquetas del mismo color que las barras
ggplot(results, aes(x = Muestra, y = RI, fill = Estrategia)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = round(RI, 3), color = Estrategia),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3.5) +  # Ajusta vjust para la posición vertical y size para el tamaño del texto
  labs(title = "Comparación de Rand Indexes", y = "Rand Index") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Elimina la cuadrícula principal
    panel.grid.minor = element_blank(),  # Elimina la cuadrícula menor
    axis.line = element_line(color = "black"),  # Añade líneas de eje
    text = element_text(size = 16),  # Aumenta el tamaño del texto en general
    axis.title = element_text(size = 18),  # Aumenta el tamaño del título de los ejes
    axis.text = element_text(size = 14),  # Aumenta el tamaño del texto de los ejes
    legend.text = element_text(size = 14),  # Aumenta el tamaño del texto de la leyenda
    legend.title = element_text(size = 16),  # Aumenta el tamaño del título de la leyenda
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  # Aumenta el tamaño del título del gráfico
  ) +
  # scale_color_manual(values = c("fer" = "blue", "lucmar" = "red")) +  # Ajusta los colores de las etiquetas
  coord_cartesian(ylim = c(0.85, 1.0))  # Establece el límite del eje y en 0.70 a 1.0

# Comparación de ARI
options(warn = -1) # Silenciamos además los warnings
options(repr.plot.width = 25, repr.plot.height = 8)

# Crear el gráfico de barras con etiquetas del mismo color que las barras
ggplot(results, aes(x = Muestra, y = RI, fill = Estrategia)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = round(ARI, 3), color = Estrategia),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3.5) +  # Ajusta vjust para la posición vertical y size para el tamaño del texto
  labs(title = "Comparación de Adjusted Rand Indexes", y = "Adjusted Rand Index") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Elimina la cuadrícula principal
    panel.grid.minor = element_blank(),  # Elimina la cuadrícula menor
    axis.line = element_line(color = "black"),  # Añade líneas de eje
    text = element_text(size = 16),  # Aumenta el tamaño del texto en general
    axis.title = element_text(size = 18),  # Aumenta el tamaño del título de los ejes
    axis.text = element_text(size = 14),  # Aumenta el tamaño del texto de los ejes
    legend.text = element_text(size = 14),  # Aumenta el tamaño del texto de la leyenda
    legend.title = element_text(size = 16),  # Aumenta el tamaño del título de la leyenda
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  # Aumenta el tamaño del título del gráfico
  ) +
  coord_cartesian(ylim = c(0.70, 1.0))  # Establece el límite del eje y en 0.70 a 1.0


# Comparamos los tiempos de ejecución
options(warn = -1) # Silenciamos además los warnings
options(repr.plot.width = 25, repr.plot.height = 8)

ggplot(results, aes(x = Muestra, y = ExecutionTime, color = Estrategia, group = Estrategia)) +
  geom_line(linewidth = 1.2) +  # Ajusta el grosor de las líneas
  geom_point(size = 3) +       # Ajusta el tamaño de los puntos
  labs(title = "Comparación del Tiempo de Ejecución",
       x = "Tamaño del Dataset",
       y = "Tiempo de Ejecución (segundos)") +
  theme_minimal() +           # Aplica un tema minimalista
  theme(
    panel.grid.major = element_blank(),  # Elimina la cuadrícula principal
    panel.grid.minor = element_blank(),  # Elimina la cuadrícula menor
    axis.line = element_line(color = "black"),  # Añade líneas de eje
    text = element_text(size = 16),  # Aumenta el tamaño del texto en general
    axis.title = element_text(size = 18),  # Aumenta el tamaño del título de los ejes
    axis.text = element_text(size = 14),  # Aumenta el tamaño del texto de los ejes
    legend.text = element_text(size = 14),  # Aumenta el tamaño del texto de la leyenda
    legend.title = element_text(size = 16),  # Aumenta el tamaño del título de la leyenda
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  # Aumenta el tamaño del título del gráfico
  )