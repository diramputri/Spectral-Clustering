# generating data with mlbench package

library(mlbench)
set.seed(1)
data<-mlbench.spirals(100,1,0.025) #argument is (# of points, # of cycles, std. dev.)
data=(data$x)*4
plot(data)

# ----------------------------------------------------------------------------------------------- #

# creating a function in R

## compute similarity matrix using Gaussian/RBF kernel ##

s <- function(x1, x2, alpha=1) {     # the 's' function will define the Gaussian kernel
  exp(- alpha * norm(as.matrix(x1-x2), type="F")) 
}

similarity <- function(data, similarity) {    # compute weighted similarities
  N <- nrow(data)
  S <- matrix(rep(NA,N^2), ncol=N)
  for(i in 1:N) {
    for(j in 1:N) {
      S[i,j] <- similarity(data[i,], data[j,])
    }
  }
  S
}

S <- similarity(data, s)
S[1:8,1:8]

## adjacency matrix ##

affinity <- function(S, n.neighbors=2) {
  N <- length(S[,1])
  
  if (n.neighbors >= N) {  # fully connected
    A <- S
  } else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity 
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighbors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, i.e., the matrix becomes symmetric
      }
    }
  }
  A  
}

A <- affinity(S, 3)  # use 3 neighboors (includes self)
A[1:8,1:8]

## degree matrix ##

D <- diag(apply(A, 1, sum)) # sum rows
D[1:8,1:8]

## graph Laplacian ##

L <- D - A
round(L[1:12,1:12],1)

## eigenvalues ##

"%^%" <- function(M, power)
  with(eigen(M), vectors %*% (values^power * solve(vectors)))
L <- diag(nrow(my.data)) - solve(D) %*% A  # simple Laplacian
round(L[1:12,1:12],1)

## clustering ##

k   <- 2
evL <- eigen(L, symmetric=TRUE)
Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
plot(Z, pch=20) # all 50 points of each cluster are on top of each other

# ----------------------------------------------------------------------------------------------- #
