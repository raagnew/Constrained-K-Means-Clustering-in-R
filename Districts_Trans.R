rm(list=ls())
constrained.kmeans <- function(p,a,std,k,b,imax,seed,details){
# Bob Agnew (raagnew1@gmail.com, www.raagnew.com)
# Constrained k-means clustering using lp.transport function
# in R package lpSolve
# Implements Bradley-Bennett-Demiriz algorithm
# p = population vector associated with observations
# Population for each observation must be assigned to same cluster
# For a normal clustering problem, p is a vector of ones
# a = numeric attribute matrix with a row for each observation
# std = TRUE if attributes are to be standardized across observations,
# else FALSE
# k = specified number of k-means clusters
# b = vector of integer lower bounds for population clusters
# If all populations are one, then these are exact lower bounds
# for cluster sizes
# If populations vary, then these are target lower bounds for
# cluster populations
# imax = maximum number of iterations in kmeans loop
# seed = random number seed for initial cluster assignments
# details = TRUE for printing of iteration details, else FALSE
m <- dim(a)[1] # Number of observations
n <- dim(a)[2] # Number of attributes
if (std) {a <- apply(a,2,function(u) (u-mean(u))/sd(u))}
set.seed(seed,kind=NULL,normal.kind=NULL) # Set random number seed
c <- sample.int(k,m,replace=TRUE,prob=NULL) # Initial cluster vector
require(lpSolve) # Load R linear programming package
num <- iter <- 1
while (iter <= imax & num > 0){
cost <- matrix(0,nrow=m,ncol=k)
for (j in 1:n){
cost <- cost + outer(a[,j],tapply(p*a[,j],c,sum)/tapply(p,c,sum),"-")^2}
# Solution of transportation linear program
trans <- lp.transport(cost,"min",rep("=",m),p,rep(">=",k),b)
# Generate new clusters
# For split allocation, assign to maximal cluster
c1 <- apply(trans$solution,1,which.max)
num <- sum(c1!=c)
if (details){
print(noquote("Iteration Number"))
print(iter)
print(trans)
print(noquote("Number of split allocations"))
print(sum(apply(trans$solution,1,max) < p))
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
# Read in Illinois census block group data
# Downloaded from https://www.census.gov/geo/reference/centersofpop.html
data <- read.table("c:/Clustering/CenPop2010_Mean_BG17.txt",sep=",",header=TRUE)
# Remove dummy block groups with no population 
data <- data[as.vector(data[,5]) > 0,]
p <- as.vector(data[,5]) #Census block group populations
# Number of census block groups with positive population
length(p)
# Block group population range
range(p)
# Block group population mean
mean(p)
# Block group population median
median(p)
a <- as.matrix(data[,c(7,6)]) # Block group longitudes and latitudes
k <- 18 # Number of clusters = number of congressional districts

# Constrained k-means clustering applied to IL congressional districting
# Target population lower bounds set tightly at average rounded down
ckm <- constrained.kmeans(p,a,FALSE,k,rep(floor(sum(p)/k),k),100,23,TRUE)
ckm$iter
ckm$converge
ckm$tot.withinss
ckm$population # Approximately equal
MEAN_LONGITUDE <- ckm$centers[,1]
MEAN_LATITUDE <- ckm$centers[,2]
c <- ckm$cluster

# Geoplot of IL census block groups and cluster centroids
pdf("c:/Clustering/IL_census_block_group_centroids.pdf")
plot(a[,1],a[,2],col="black",pch=19,cex=.25,asp=1.5,xlim=c(-94,-84),ylim=c(36,43),xlab="LONGITUDE",ylab="LATITUDE",
main="ILLINOIS CENSUS BLOCK GROUPS & CLUSTER CENTROIDS")
points(MEAN_LONGITUDE,MEAN_LATITUDE,col="red",pch=19,cex=.75)
dev.off()

# Geoplot of IL census block group clusters
pdf("c:/Clustering/IL_census_block_group_clusters.pdf")
plot(a[c==1,1],a[c==1,2],asp=1.5,cex=.25,col="black",pch=19,xlim=c(-94,-84),ylim=c(36,43),xlab="LONGITUDE",ylab="LATITUDE",
main="ILLINOIS CENSUS BLOCK GROUP CLUSTERS")
points(a[c==2,1],a[c==2,2],col="dimgray",pch=19,cex=.25)
points(a[c==3,1],a[c==3,2],col="gray70",pch=19,cex=.25)
points(a[c==4,1],a[c==4,2],col="red",pch=19,cex=.25)
points(a[c==5,1],a[c==5,2],col="darkred",pch=19,cex=.25)
points(a[c==6,1],a[c==6,2],col="darkorange1",pch=19,cex=.25)
points(a[c==7,1],a[c==7,2],col="purple1",pch=19,cex=.25)
points(a[c==8,1],a[c==8,2],col="blue",pch=19,cex=.25)
points(a[c==9,1],a[c==9,2],col="aquamarine1",pch=19,cex=.25)
points(a[c==10,1],a[c==10,2],col="dodgerblue1",pch=19,cex=.25)
points(a[c==11,1],a[c==11,2],col="cyan1",pch=19,cex=.25)
points(a[c==12,1],a[c==12,2],col="deeppink1",pch=19,cex=.25)
points(a[c==13,1],a[c==13,2],col="darkgoldenrod1",pch=19,cex=.25)
points(a[c==14,1],a[c==14,2],col="green1",pch=19,cex=.25)
points(a[c==15,1],a[c==15,2],col="darkgreen",pch=19,cex=.25)
points(a[c==16,1],a[c==16,2],col="maroon1",pch=19,cex=.25)
points(a[c==17,1],a[c==17,2],col="olivedrab1",pch=19,cex=.25)
points(a[c==18,1],a[c==18,2],col="darksalmon",pch=19,cex=.25)
dev.off()
date() # Time

