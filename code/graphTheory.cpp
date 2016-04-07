#include <Rcpp.h>
#include <queue>
using namespace Rcpp; // rcpp library connection
using namespace std; // standard library

// -------------------------------------------------------------------------
// --------------- BREADTH FIRST SEARCH ALGORITHM -------------------------
// -------------------------------------------------------------------------

// OBS: LAS CLASES NUMERICVECTOR NUMERICMATRIX NO LAS DA RCPP Y SON MUY COMODAS

// [[Rcpp::export]]
LogicalVector breadth_first_search(vector< vector<int> > adjlst, int root) {
  int n = adjlst.size(); // tamaño de la lista
  LogicalVector visited(n); // indicara si un nodo esta conectado o no
  // IntegerVector level(n); // la distancia al nodov 
  // IntegerVector parent(n); // el padre a traves del cual fue el nodo fue visitado
  queue<int> aux_queue; // la cola que vamos a usar
  int i; // iterador principal
  // STEP 0: BEFORE START
  for(i=0; i<n; i++){
    visited[i] = false;
    // level[i] = -1; // -1 significa infinito aquí
    // parent[i] = -1; // -1 significa no tiene padre
  }
  // STEP 1: START WITH ROOT
  visited[root] = true; //
  // level[root] = 0;
  // parent[root] = root;
  aux_queue.push(root);
  // STEP 2: BREADTH-SEARCH
  int top; // el elemento hasta arriba de la cola
  int nnbhs; // numero de vecinos en ese elemento de la lista
  while(!aux_queue.empty()){
    top = aux_queue.front(); // extraemos el nodo en el tope de la cola para explorar sus hijos
    nnbhs = adjlst[top].size(); // cantidad de vecinos de top
    aux_queue.pop(); // quitamos al elemento en el tope de la cola
    for(i=0; i<nnbhs; i++){
      if(visited[adjlst[top][i]]==false){
        visited[adjlst[top][i]]=true;
        // level[adjlst[top][i]] = level[top] + 1;
        // parent[adjlst[top][i]] = top;
        aux_queue.push(adjlst[top][i]);
      }
    }
  }
  return visited;
}




#include <Rcpp.h>
#include <queue>
using namespace Rcpp; // rcpp library connection
using namespace std; // standard library

// -------------------------------------------------------------------------
// -------- BREADTH FIRST SEARCH DENSITY CLUSTERING ALGORITHM ---------------
// -------------------------------------------------------------------------

// FUNCIÓN AYUDA, REGRESA EL ÍNDICE DEL J-ÉSIMO ELEMENTO MÁS PEQUEÑO DE UN VECTOR

// [[Rcpp::export]]

inline NumericVector smallest_k(NumericVector x, int k=1){
  NumericVector y=clone(x);
  NumericVector result(k);
  y.sort();
  int pos = 0;
  int it = 0;
  bool stop = false;
  while(pos < x.size() && !stop){
      if(x[pos] <= y[k-1]){
        result[it] = pos;
        it++;
      }
      if(it==k){
        stop=true;
      }
      pos++;
  }
  return result;
}


// [[Rcpp::export]]
LogicalVector bfs_density_component(NumericMatrix distmat, int root=0, double jump_factor=3, int nbs=10) {
  int n = distmat.nrow(); // cantidad de datos
  LogicalVector visited(n); // indicara si un nodo esta conectado o no
  double infinity=999999999;
  double density, sum_local_dists; // empezamos con un valor muy grande para la primera iter
  queue<int> aux_queue; // la cola que vamos a usar
  
  int i; // iterador principal

  // STEP 0: BEFORE START
  for(i=0; i<n; i++){
    visited[i] = false;
  }

  // STEP 1: INITIAL VALUES FOR BREADTH SEARCH
  NumericVector distances = distmat.row(root);
  NumericVector neighbours(nbs);
  distances[root] = infinity;
  visited[root] = true;
  int n_visited = 0; // para actualizar la densidad
  int nbs_queued;
  int top=root; // el frente de la cola
  aux_queue.push(root);
  density = infinity;
  
  // STEP 2: BREADTH-SEARCH
  do{
    // tomamos el frente de la cola
    // cout << "Padre: " << aux_queue.front();
    top = aux_queue.front();
    aux_queue.pop(); // quitamos el elemento tope
    
    // queremos calcular la distancias para la densidad y los vecinos
    distances = distmat.row(top);
    distances[top] = infinity;
    
    // vamos a ver que vecinos vamos a visitar
    neighbours = smallest_k(distances, nbs);

    nbs_queued=0;
    sum_local_dists=0;
    for(i=0; i < nbs; i++){
      // cout << "; vecino " << neighbours[i] << ": " << distances[neighbours[i]] << ", ";
      if(!visited[neighbours[i]] && (distances[neighbours[i]] < jump_factor*density) ){
        aux_queue.push(neighbours[i]);
        visited[neighbours[i]] = true;
        sum_local_dists += distances[neighbours[i]];
        nbs_queued++;
      }
    }
    density = (density*n_visited + sum_local_dists)/(n_visited + nbs_queued);
    n_visited += nbs_queued;
    
    // cout << "densidad: " << density;
    // cout << endl;
    
  } while (!aux_queue.empty());

  return visited;
}


/*** R
distmat <- as.matrix(dist(iris[ ,1:4]))

par(mfrow=c(2,2))

group <- bfs_density_component(distmat, nbs=5, jump_factor=2)
plot(iris[ ,1], iris[ ,2], bg=factor(group), pch=21, cex=2, main="View 1: factor brinco 2 y 5 vecinos")
group <- bfs_density_component(distmat, nbs=5, jump_factor=2)
plot(iris[ ,3], iris[ ,4], bg=factor(group), pch=21, cex=2, main="View 2: factor brinco 2 y 5 vecinos")
group <- bfs_density_component(distmat, nbs=10, jump_factor=8)
plot(iris[ ,1], iris[ ,2], bg=factor(group), pch=21, cex=2, main="View 1: factor brinco 8 y 10 vecinos")
group <- bfs_density_component(distmat, nbs=10, jump_factor=8)
plot(iris[ ,3], iris[ ,4], bg=factor(group), pch=21, cex=2, main="View 2: factor brinco 8 y 10 vecinos")


# USANDO SOLO LA VISTA 3,4

distmat <- as.matrix(dist(iris[ ,3:4]))

par(mfrow=c(1,2))

group <- bfs_density_component(distmat, nbs=5, jump_factor=2)
plot(iris[ ,3], iris[ ,4], bg=factor(group), pch=21, cex=2, main="factor brinco 2 y 5 vecinos")
group <- bfs_density_component(distmat, nbs=15, jump_factor=10)
plot(iris[ ,3], iris[ ,4], bg=factor(group), pch=21, cex=2, main="Vfactor brinco 10 y 10 vecinos")


# CLUSTERING ALGORITH!!!!

distmat <- as.matrix(dist(iris[ ,1:4]))

clustITAM <- function(distmat){
  mat <- distmat
  clustList <- list()
  while(nrow(mat)==0){
    indices <- which(bfs_density_component(mat))
    mat <- mat[-indices, -indices]
    clustlist <- c(clustList, list(indices))
  }
  return(clustList)
}


# adjacency2list <- function(adjmat){
#   lst <- vector("list", nrow(adjmat))
#   for(i in 1:nrow(adjmat)) lst[[i]] <- which(adjmat[i, ]!=0) - 1
#   return(lst)
# }
# adjmat <- matrix(c(0, 0, 1, 0, 0, 0, 1, 0, 0), nrow=3, ncol=3, byrow = TRUE)
#   adjlst <- adjacency2list(adjmat)
#   rBFS <- function(adjlst=NULL, root=1, adjmat=NULL){
#     if(is.null(adjlst)) adjlst = adjacency2list(adjmat)
#       component <- which(breadth_first_search(adjlst=adjlst, root=root-1))
#       return(component)
#   }
# 
# # 
# set.seed(110104)
#   nrows <- 1000
# nnbhs <- 5
# distmat <- matrix(runif(nrows^2), nrow=nrows, ncol=nrows)
#   distmat <- (distmat + t(distmat))/2
# diag(distmat) <- 0
# 
# # # rápidamente sacamos als ditancias de los k-vecinos más cercanos
# # densityKnbhs <- function(distmat, k=5){
# which_nbhs <- do.call("rbind", lapply(1:nrow(distmat), function(i) order(distmat[i, ])[2:(2+nnbhs)]))
#   nbhs_vals <- do.call("rbind", lapply(1:nrow(distmat), function(i) distmat[i, which_nbhs[i, ]]))
#   densities <- apply(nbhs_vals, 1, mean)
# # }
# # 
# # 
# # 
# # 
#   adjlst <- lapply(1:nrow(distmat), function(i) which_nbhs[i, ] - 1)
#   rBFS(adjlst, root=10)
#   rBFS(adjlst, root=12)
#   rBFS(adjlst, root=1)
#   
  */



