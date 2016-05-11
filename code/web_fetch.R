
# matriz <- matrix(0, n, n)
  
flow_mat <- matrix(c(0, 2, 1, 3,
              -2, 0, 1, -1,
              -1, -1, 0, -2,
              -3, -1, -2, 0), ncol=4, nrow=4)
n <- nrow(flow_mat) # numero de vertices

X <- matrix(0, ncol=n, nrow=n*(n-1)/2+1)
y <- numeric(n*(n-1)/2+1)

X[1, ] <- 1 # llena de unos la primera fila de X
y[1] <- 0 # pone cero en la primera entrada de y
row <- 2
for(i in 1:(n-1)) for(j in (i+1):n) {
    X[row,i] <- -1
    X[row,j] <- 1
    y[row] <- mat[i, j]
    row <- row + 1
}

t(X) %*% X # vean que bonita queda la matriz

# y la solución por lo tanto es:
qr.solve(t(X)%*%X, t(X)%*%y)
# o más fácil
t(X)%*%y/n

mod <- lm(y~X-1, data=data.frame(X,y))
summary(mod)


# pagerank

transmat <- read.csv("MarkovChainExamples/data.csv")[ ,-1]
n <- nrow(transmat)
alpha <- .85
newtransmat <- alpha*transmat + (1-alpha)*matrix(1/n, n, n)

# 1 con page rank
library(expm)

v0 <- rep(1/n, n)
pi_vec <- t(v0) %*% (as.matrix(newtransmat) %^% 1000)
dd <- data.frame(university=names(newtransmat), probability=as.numeric(pi_vec))
View(dd[order(dd$probability, decreasing=TRUE), ])
names(transmat)[order(pi_vec, decreasing=TRUE)]

# HODGE RANK
newtransmat <- as.matrix(newtransmat) %^% 2

flow_mat <- matrix(0, n, n)
for(i in 1:(n-1)) for(j in (i+1):n) flow_mat[i, j] <- log(newtransmat[i, j]/newtransmat[j, i])
flow_mat <- flow_mat - t(flow_mat)

X <- matrix(0, ncol=n, nrow=n*(n-1)/2+1)
y <- numeric(n*(n-1)/2+1)

X[1, ] <- 1
y[1] <- 0
row <- 2
for(i in 1:(n-1)) for(j in (i+1):n) {
  X[row,i] <- -1
  X[row,j] <- 1
  y[row] <- flow_mat[i, j]
  row <- row + 1
}

# y la solución por lo tanto es:
qr.solve(t(X)%*%X, t(X)%*%y)
# o más fácil
sol <- t(X)%*%y/n

names(transmat)[order(sol, decreasing=TRUE)]





# RANKING FUTBOL
  
  data <- read.csv("torneo_clausura.csv", stringsAsFactors = FALSE)
  data[ ,1] <- gsub("[^ a-zA-Z']","", data[ ,1]) # quitemos espacios vacios
  data[ ,2]<- gsub("[^ a-zA-Z']","", data[ ,2])
  scores <- do.call("rbind", strsplit(data$Marcador, ":"))
  data$score_local <- as.numeric(scores[ ,1])
  data$score_visit <- as.numeric(scores[ ,2])
  
  nodos <- unique(data$Local)
  n <- length(nodos)
  flow_mat <- matrix(0, n, n, dimnames = list(nodos, nodos))
  for(i in 1:nrow(data)) 
    flow_mat[data$Local[i], data$Visitante[i]] <- 2*data$score_visit[i]-data$score_local[i]
  
  X <- matrix(0, ncol=n, nrow=n*(n-1)/2+1)
  y <- numeric(n*(n-1)/2+1)
  wts <- rep(1, n*(n-1)/2+1)
  X[1, ] <- 1
  
  y[1] <- 0
  wts[3] <- 0
  row <- 2
  for(i in 1:(n-1)) for(j in (i+1):n) {
    X[row,i] <- -1
    X[row,j] <- 1
    y[row] <- flow_mat[i, j]
    row <- row + 1
  }
  
  # y la solución por lo tanto es:
  sol <- qr.solve(t(X)%*% diag(wts) %*% X, t(X)%*% diag(wts) %*% y)
  dd <- data.frame(nodos, sol)
  View(dd[order(dd$sol, decreasing=TRUE), ])
  
  mod <- lm(y~X-1, data=data.frame(X,y), weights = wts)
  summary(mod)
  
  # con pesos
  