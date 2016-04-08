# ------------------------------------------------------
# ------------------------------------------------------
# LOAD DOG AND CATA SAMPLE DATA ------------------------
# ------------------------------------------------------
# ------------------------------------------------------

library("doParallel")
library("jpeg")
library("rARPACK")
library("Matrix")
library("bigmemory")

set.seed(110104)
sample.size <- 2500
nfactors <- 1500
pixWidth <- 120
pixHeight <- 90
filter <- "kpca" # OPCIONES: kpca, pca, centralizador

n_int <- 20  # number of intervals for subdividing
p <- 0.2    # proportion of each interval that should overlap with the next
eps <- 0.7  # max clustering value


# ------------------------------------------------------
# LEE Y ESCRIBE ---------------------------------------
# ------------------------------------------------------
baseStr <- "data/catsDogs/"
fileList <- paste0(baseStr, list.files(baseStr))
fileList <- fileList[sample.int(length(fileList), size=sample.size)]
catDummy <- grepl("cat\\.", fileList)
# cl <- makeCluster(4)
# registerDoParallel(cl)
# dta <- parLapply(cl, as.list(fileList), function(fileName){
#   dta <- readJPEG(fileName)
#   dta <- (1/3)*(dta[ , , 1] + dta[ , , 2] + dta[ , , 3])
#   rows <- seq(1, nrow(dta), length.out=pixHeight)
#   cols <- seq(1, ncol(dta), length.out=pixWidth)
#   dta <- dta[rows, cols]
#   return(as.numeric(dta))
# })
# stopCluster(cl)
# dta <- do.call("rbind", dta)
# dta <- as.big.matrix(dta) # es la mejor forma de guardar los datos
# saveRDS(dta, "data/sampleDta/catsDogs.RDS")
dta <- readRDS("data/sampleDta/catsDogs.RDS")
dta <- as.big.matrix(dta) # no space
# ------------------------------------------------------
# ------------------------------------------------------


if(filter=="kpca"){
# ------------------------------------------------------
# KERNEL PCA -------------------------------------------
# ------------------------------------------------------
sigma0 <- sqrt(ncol(dta))*sd(dta[1, ]) # mis heurísticas
# para poder hacer esto paralelo necesito paralelizar el C++ o compilar mi función en un paquete
Rcpp::sourceCpp('kernelMatrix.cpp')
# # lapply(as.list(1:nrow(dta)), function(i){
# #   vec <- rbfdotDist(dta[i, ], dta, sigma0)
# #   saveRDS(vec, sprintf("data/sampleDta/row%i.RDS", i))
# # })
# # cl <- makeCluster(4)
# # registerDoParallel(cl)
# # kmatrix <- foreach(i=1:nrow(dta), .combine=rbind) %dopar% readRDS(sprintf("data/sampleDta/row%i.RDS", i))
# # stopCluster(cl)
# system.time(
#   kmatrix <- rbfdotMatrix(as.matrix(dta), sigma0)
#   ) # claramente C++ es la opcion (de horas a 12m; paralelizando es a 2m)
# normalizer <- forceSymmetric(matrix(1/nrow(kmatrix), nrow=nrow(kmatrix), ncol=nrow(kmatrix)))
# kmatrix <- kmatrix - normalizer%*%kmatrix - kmatrix%*%normalizer + normalizer%*%kmatrix%*%normalizer
# rm(normalizer)
# saveRDS(kmatrix, "data/sampleDta/kmatrix.RDS")
kmatrix <- readRDS("data/sampleDta/kmatrix.RDS")
# ------------------------------------------------------
# ANALYZE COMPONENTS -----------------------------------
# ------------------------------------------------------
system.time(res <- eigs(forceSymmetric(as.matrix(kmatrix)), nfactors))
res$vectors <- apply(res$vectors, 2, function(x) x/sd(x))
nind <- nrow(dta)
ncomps <- min(which(res$values/sum(res$values) <= 1/nind)) 
plot((res$values/sum(res$values))[1:(2*ncomps)], type="l")
abline(h=1/nind, col="red")
### Compute projections
prkcomp <- kmatrix %*% res$vectors
plot(prkcomp[ ,1], prkcomp[ ,2], type="p", col=factor(catDummy, labels=c("dog", "cat")), 
     pch=20, main="Primeras componentes principales", xlab="componente 1", ylab="componente 2") 
legend("topright", c("dog", "cat"), col=c("red", "black"), pch=20)
filterVar <- prkcomp[ ,1]
# comp1.order <- order(prkcomp[ ,1]) # orden.
# par(mfrow=c(3,3))
# Algunos ejemplos
# plot(pixmapGrey(dta[comp1.order[1], ], ncol=pixWidth))
# plot(pixmapGrey(dta[comp1.order[2], ], ncol=pixWidth))
# plot(pixmapGrey(dta[com p1.order[3], ], ncol=pixWidth))
# plot(pixmapGrey(dta[comp1.order[sample.size/2], ], ncol=pixWidth))
# plot(pixmapGrey(dta[comp1.order[sample.size/2+1], ], ncol=pixWidth))
# plot(pixmapGrey(dta[comp1.order[sample.size/2+2], ], ncol=pixWidth))
# plot(pixmapGrey(dta[comp1.order[sample.size-2], ], ncol=pixWidth))
# plot(pixmapGrey(dta[comp1.order[sample.size-1], ], ncol=pixWidth))
# plot(pixmapGrey(dta[comp1.order[sample.size], ], ncol=pixWidth))
# ------------------------------------------------------
# ------------------------------------------------------
}


if(filter=="svd"){
# ------------------------------------------------------
# SVD --------------------------------------------------
# ------------------------------------------------------
system.time(res <- svds(apply(as.matrix(dta), 2, function(x) x-mean(x)), 2))
filterVar <- res$u[ ,1]
# ------------------------------------------------------
# ------------------------------------------------------ 
}

if(filter=="pca"){
# ------------------------------------------------------
# PCA --------------------------------------------------
# ------------------------------------------------------
system.time({
  temp.mat <- apply(as.matrix(dta), 2, function(x) x-mean(x))
  temp.mat <- (1/(nrow(dta)-1))*t(as.matrix(dta))%*%as.matrix(dta) # this takes time....
  res <- eigs(temp.mat, 2)
  rm(temp.mat)
  comps <- as.matrix(dta)%*% res$vectors
  plot(comps[ ,1], comps[ ,2], col=catDummy+2, xlab="CP1", ylab="CP2", main="Componentes principales")
  filterVar <- comps[ ,1]
})
# ------------------------------------------------------
# ------------------------------------------------------ 
}

if(filter=="centrality"){
# ------------------------------------------------------
# CENTRALITY --------------------------------------------------
# ------------------------------------------------------
centerOfGravity <- apply(as.matrix(dta), 2, mean)
filterVal <- sapply(1:nrow(dta), function(row) sum(centerOfGravity-dta[row, ]))
# ------------------------------------------------------
# ------------------------------------------------------ 
}


# ------------------------------------------------------
# DIVIDE DATA BY SUBINTERVALS  (código Fer)  -----------
# ------------------------------------------------------
#this section will create a data frame in which we will construct overlapping intervals
var_o <- filterVar
intervals_centers <- seq(min(var_o),max(var_o),length=n_int)  #basic partition = centers
interval_length <- intervals_centers[2]-intervals_centers[1]  #to create the overlaps of p% of this length
intervals <- data.frame(centers=intervals_centers)            #create a data frame
#create the overlapping intervals  
intervals$min <- intervals_centers - (0.5+p)*interval_length                     
intervals$max <- intervals_centers + (0.5+p)*interval_length
#decent name for the intervals e.g    [5.34;6.53)     [6.19;7.39)
intervals$interval <- seq(1,n_int)
intervals$name <- with(intervals, sprintf("[%.2f;%.2f)",min,max))
#function that will split the variable according to the invervals
splitting <- lapply(split(intervals,intervals$interval), function(x) which(var_o> x$min & var_o <= x$max)) 
# ------------------------------------------------------
# ------------------------------------------------------

# ------------------------------------------------------
# CLUSTERING PER GROUP   -------------------------------
# ------------------------------------------------------

# ------------------------------------------------------
# ------------------------------------------------------



# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------