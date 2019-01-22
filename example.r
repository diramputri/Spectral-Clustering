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
# generating data with mlbench package

library(mlbench)
set.seed(1)
data<-mlbench.spirals(100,1,0.025) #argument is (# of points, # of cycles, std. dev.)
data=(data$x)*4
plot(data)

# ----------------------------------------------------------------------------------------------- #

# creating a function in R

## compute similarity matrix using Gaussian/RBF kernel ##

# The Gaussian kernel provides a measure of similarity. 
# When points x1 and x2 are close, similarity S estimates to 1. 
# If x1 and x2 are far, S is estimated to be 0. You can see why through the equation in line 21.
# We set alpha = 1 for simplicity.

s <- function(x1,x2,alpha=1) {     # the 's' function will define the Gaussian kernel, our similarity measure
  exp(-alpha*norm(as.matrix(x1-x2),type="F")) #type "F" = Frobenius/Euclidean norm
}

similarity <- function(data,similarity) {    # compute weighted similarities
  N <- nrow(data)
  S <- matrix(rep(NA,N^2),ncol=N)
  for(i in 1:N) {
    for(j in 1:N) {
      S[i,j] <- similarity(data[i,], data[j,])
    }
  }
  S
}

S <- similarity(data, s)
S[1:8,1:8] # here is our matrix showing similarities between points
# S is a 100x100 matrix; reduced to [1:8,1:8]

## adjacency matrix ##

# now using our similarity measure, we need to make an adjacency matrix
# conditons: symmetric and all entries are greater than or equal to 0
# below, we use a kNN filter to connect points that are close to each other

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
        A[j,i] <- S[i,j] # making our matrix symmetric i.e. an undirected graph
      }
    }
  }
  A  
}

A <- affinity(S, 3)  # use 3 neighbors (includes self)
A[1:8,1:8]

# with the adjacency matrix in place, we now have a graph partitioning problem
# partitioned such that edges in the same cluster have high weights
# edges connecting different clusters have low weight
# need degree matrix D, a diagonal matrix with entries relating to # of edges/weights of edges connected to a vertex

## degree matrix ##

D <- diag(apply(A, 1, sum)) # sum rows
D[1:8,1:8]

## graph Laplacian ##
# we just compute an unnormalized Laplacian, though we could do L_sym or L_rw as well

L <- D - A
#L is a 100x100 matrix
#round values to 1 decimal place, reduce L to 8x8
round(L[1:8,1:8],1)

## eigenvalues ##

"%^%" <- function(M, power)
  with(eigen(M), vectors %*% (values^power * solve(vectors))) 
# with() evaluates an expression in an environment constructed by data
# can possibly modify the data
# sovle() solves the equation a %*% x = b for x, where b can be either a vector or a matrix


## clustering ##

k=2 #number of clusters we want
# compute the k smallest eigenvectors, ignoring the constant eigenvector
evL <- eigen(L, symmetric=TRUE)
# cluster through k-means
Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
plot(Z, pch=20) # all 50 points of each cluster are on top of each other

# here is our original data visualized
# red points belong to one cluster, black points belong in another cluster
library(stats)
km <- kmeans(Z, centers=k, nstart=5)
plot(data, col=km$cluster)

# ----------------------------------------------------------------------------------------------- #

# there is already a package in R that does all of this
# located in the "kernlab" package
# we can use the same simulated data in the first block

library(kernlab)
spec <- specc(data,centers=2) # spectral clustering function
plot(data,col=spec,pch=4)         
points(data, pch=5) 

# we get the same result as above
