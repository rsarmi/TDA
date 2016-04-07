# SIMULA EL KEPPLER MAPPER

# STEPS:
#     1. CREAR JSON

# -------  PARAMETROS ---------------------------------------------
library(jsonlite)
set.seed(110104)
#


# ------- SIMULACION DE UNA MATRIZ SIMETRICA CON PESOS 1 o 2 ------
nodes.n <- 3
#adj.matrix <- matrix(rbinom(n=nodes.n^2, size=5, p=.1), nrow=nodes.n, ncol=nodes.n)
adj.matrix<-as.matrix(read.csv("/home/pedrohserrano/TDA/adj_matrix.csv"))
##AQI LA MATRIZ
diag(adj.matrix) <- 0
adj.matrix <- adj.matrix + t(adj.matrix)
head(adj.matrix)

#nodes.size <- runif(nodes.n, 4, 10)
nodes.size<-100*as.matrix(read.csv("/home/pedrohserrano/TDA/summary_cluster.csv"))
nodes.tooltips <- paste("Grupo:", 1:nodes.n)
nodes.names <- 1:nodes.n
nodes.color <- as.character(1:nodes.n)
#

# ------- AHORA TENEMOS QUE CREAR UN JSON DE ESO -----------------------------

aux_mat <- data.frame()
for(i in 1:nodes.n) for(j in 1:nodes.n) if(adj.matrix[i, j]!=0) aux_mat <- rbind(aux_mat, data.frame(source=i-1, target=j-1, value=adj.matrix[i, j]))
linksJSON <- toJSON(aux_mat)
nodesJSON <- toJSON(data.frame(color=nodes.color, group=nodes.size, name=nodes.names, tooltip=nodes.tooltips))
graphJSON <- sprintf("{\"nodes\": %s, \"links\": %s}", nodesJSON, linksJSON)

head(graphJSON)


# ------------  CREAMOS EL HTML ----------------------------------------------------------
htmlFile <- readLines("/home/pedrohserrano/TDA/ManifoldLearning/www/index.html")
#htmlFile <- readLines("www/index.html")
graph_def_line <- which(grepl("graph =", htmlFile))
#htmlFile[graph_def_line] <- sprintf("graph = %s;", graphJSON)
htmlFile[graph_def_line] <- sprintf("graph = %s;", graphJSON)
#writeLines(htmlFile, "www/index.html")
writeLines(htmlFile, "/home/pedrohserrano/TDA/ManifoldLearning/www/index.html")

browseURL("file:///home/pedrohserrano/TDA/ManifoldLearning/www/index.html")
