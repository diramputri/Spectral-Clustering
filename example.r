library(mlbench)
set.seed(1)
obj=mlbench.spirals(100,1,0.025)
my.data=4 * obj$x
plot(my.data)
library(kernlab)
sc=specc(my.data, centers=2)
plot(my.data,col=sc,pch=4)        
points(my.data,col=obj$classes,pch=5)
