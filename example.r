# simulated data with mlbench package

library(mlbench)
set.seed(1)
data<-mlbench.spirals(100,1,0.025) #argument is (# of points, # of cycles, std. dev.)
data=(data$x)*4
plot(data)

# ----------------------------------------------------------------------------------------------- #

# regular function in R

## compute similarity matrix ##

s <- function(x1, x2, alpha=1) {
  exp(- alpha * norm(as.matrix(x1-x2), type="F"))
}

similarity <- function(data, similarity) {
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

## affinity matrix ##

affinity <- function(S, n.neighboors=2) {
  N <- length(S[,1])
  
  if (n.neighboors >= N) {  # fully connected
    A <- S
  } else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity 
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighboors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
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

# ----------------------------------------------------------------------------------------------- #

# using kernlab package
