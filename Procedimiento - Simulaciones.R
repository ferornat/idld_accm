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

# Vemos que en muchos casos, al tratarse de un procesimiento que detrás está computando inversa de matrices, cuando llegamos a 
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