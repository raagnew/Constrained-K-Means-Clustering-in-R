rm(list=ls())
constrained.kmeans2 <- function(p,a,std,k,b,imax,seed,details){
# Bob Agnew (raagnew1@gmail.com, www.raagnew.com)
# Constrained k-means clustering using collapsed dual optimization
# Approximates Bradley-Bennett-Demiriz algorithm for large problems
# p = population vector associated with observations
# Population for each observation must be assigned to same cluster
# For a normal clustering problem, p is a vector of ones
# a = numeric attribute matrix with a row for each observation
# std = TRUE if attributes are to be standardized across observations,
# else FALSE
# k = specified number of k-means clusters
# b = vector of integer lower bounds for population clusters
# These are target lower bounds for cluster populations
# imax = maximum number of iterations in kmeans loop
# seed = random number seed for initial cluster assignments
# details = TRUE for printing of iteration details, else FALSE
m <- dim(a)[1] # Number of observations
n <- dim(a)[2] # Number of attributes
if (std) {a <- apply(a,2,function(u) (u-mean(u))/sd(u))}
set.seed(seed,kind=NULL,normal.kind=NULL) # Set random number seed
c <- sample.int(k,m,replace=TRUE,prob=NULL) # Initial cluster vector
num <- iter <- 1
while (iter <= imax & num > 0){
cost <- matrix(0,nrow=m,ncol=k)
for (j in 1:n){
cost <- cost + outer(a[,j],tapply(p*a[,j],c,sum)/tapply(p,c,sum),"-")^2}
cost <- matrix(p,nrow=m,ncol=k)*cost
# Collapsed dual function
dual <- function(v){-sum(apply(cost-outer(p,v,"*"),1,min)) - sum(b*v)}
# Approximate nonlinear solution of collapsed dual linear program
opt <- optim(rep(0,k),dual,method="L-BFGS-B",lower=rep(0,k),
control=list(maxit=1000))
# Generate new clusters
# Assign to minimum adjusted cost cluster
c1 <- apply(cost - outer(p,opt$par,"*"),1,which.min)
num <- sum(c1!=c)
if (details){
print(noquote("Iteration Number"))
print(iter)
print(noquote("Dual Objective"))
print(-opt$value)
print(noquote("Dual Variables"))
print(opt$par)
print(noquote("Number of revised cluster elements"))
print(num)}
c <- c1 #Clusters revised
iter <- iter + 1
}
means <- NULL
for (j in 1:n){
means <- cbind(means,tapply(p*a[,j],c,sum)/tapply(p,c,sum))}
cost <- matrix(0,nrow=m,ncol=k)
for (j in 1:n){
cost <- cost + outer(a[,j],tapply(p*a[,j],c,sum)/tapply(p,c,sum),"-")^2}
cost <- matrix(p,nrow=m,ncol=k)*cost
ss <- sum(cost*outer(c,1:k,"=="))
result <- list(iter-1,num==0,c,tapply(p,c,sum),means,ss)
names(result) <- c("iter","converge","cluster","population","centers",
"tot.withinss")
return(result)}

date() # Time
# Read in pulsar candidate emission data
# Downloaded from UCI Machine Learning Repository 
# (https://archive.ics.uci.edu/ml/datasets/HTRU2#) 
data <- read.csv("c:/Clustering/HTRU_2.csv",header=FALSE)
a <- as.matrix(data[,1:8]) # Eight numeric attributes for detected emissions
class <- as.vector(data[,9]) # Binary indicator of confirmed pulsar
m <- dim(a)[1] # Number of observations (emissions)
n <- dim(a)[2] # Number of cluster variables
# Standardize attributes for regular kmeans
as <- apply(a,2,function(u) (u-mean(u))/sd(u))
k <- 10 # Number of clusters
entropy <- function(t){
p <- t/rowSums(t)
l2 <- log2((p!=0)*p + (p==0))
q <- rowSums(t)/sum(t)
-sum(q*rowSums(p*l2))}

# Unconstrained clustering benchmark using R kmeans function
# Ignores cluster sizes
set.seed(23,kind=NULL,normal.kind=NULL) #Set random number seed
km <- kmeans(as,k,iter.max=100,nstart=1)
km$iter
km$size
km$tot.withinss
table(km$cluster,class)
entropy(table(km$cluster,class))

# Constrained clustering with lower bound of 500 on cluster size
ckm <- constrained.kmeans2(rep(1,m),a,TRUE,k,rep(500,k),100,23,FALSE)
ckm$iter
ckm$converge
ckm$population
ckm$tot.withinss
table(ckm$cluster,class)
entropy(table(ckm$cluster,class))
date() # Time

