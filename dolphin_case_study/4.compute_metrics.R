# load package and data
library(rjags)
library(sna)
library(igraph)

load('dolphin_model.RData')

# pool chains
res <- rbind(jsample[[1]],jsample[[2]]) # join two lists of MCMC values
dim(res)
head(res)
colnames(res)

# get posterior summaries for model parameters
# values given in Table 2
apply(res[,1:8],2,mean)
apply(res[,1:8],2,quantile,probs=c(2.5,97.5)/100)

# extract states
grid <- seq(9,ncol(res)+1,by=2485)
x1 <- res[,grid[1]:(grid[2]-1)] # occ 1
x2 <- res[,grid[2]:(grid[3]-1)] # occ 2
x3 <- res[,grid[3]:(grid[4]-1)] # occ 3
x4 <- res[,grid[4]:(grid[5]-1)] # occ 4
x5 <- res[,grid[5]:(grid[6]-1)] # occ 5

# build all possible N(N-1)/2 dyads from N individuals
N <- 71
nb.dyads <- N*(N-1)/2
# matrix of dyads' labels
label.dyads <- matrix(0,ncol=3,nrow=nb.dyads)
counter <- 0
for (i in 1:N){ 
	for (j in (i+1):N){ 
		if ((i+1) > N) next # to overcome r-trap in the loop on previous line...
		counter <- counter + 1
		label.dyads[counter,1] <- counter
		label.dyads[counter,2] <- i
		label.dyads[counter,3] <- j
	}
}

#--- indiv metrics

nbMCMC <- 1000
sselect <- sample(nrow(res),nbMCMC)

degreeMCMC <- array(0,c(nbMCMC,5,71))
betweenMCMC <- array(0,c(nbMCMC,5,71)) 


# network at occasion 1
xx <- x1[sselect,]
for (i in 1:nrow(xx)){ 
network1 <- label.dyads
mask <- as.logical(2 - xx[i,]) # in x1: 1 = activated, 2 = non-activated 
network1 <- network1[,-1]
network1 <- network1[mask,]
nn <- as.network(network1)
# compute degree and betweeness
degreeMCMC[i,1,] <- degree(nn,gmod='graph') # node degree
betweenMCMC[i,1,] <- betweenness(nn,gmod='graph') # betweeness
	}

# network at occasion 2
xx <- x2[sselect,]
for (i in 1:nrow(xx)){ 
network1 <- label.dyads
mask <- as.logical(2 - xx[i,]) # in x1: 1 = activated, 2 = non-activated 
network1 <- network1[,-1]
network1 <- network1[mask,]
nn <- as.network(network1)
# compute degree and betweeness
degreeMCMC[i,2,] <- degree(nn,gmod='graph') # node degree
betweenMCMC[i,2,] <- betweenness(nn,gmod='graph') # betweeness
	}

# network at occasion 3
xx <- x3[sselect,]
for (i in 1:nrow(xx)){ 
network1 <- label.dyads
mask <- as.logical(2 - xx[i,]) # in x1: 1 = activated, 2 = non-activated 
network1 <- network1[,-1]
network1 <- network1[mask,]
nn <- as.network(network1)
# compute degree and betweeness
degreeMCMC[i,3,] <- degree(nn,gmod='graph') # node degree
betweenMCMC[i,3,] <- betweenness(nn,gmod='graph') # betweeness
	}

# network at occasion 4
xx <- x4[sselect,]
for (i in 1:nrow(xx)){ 
network1 <- label.dyads
mask <- as.logical(2 - xx[i,]) # in x1: 1 = activated, 2 = non-activated 
network1 <- network1[,-1]
network1 <- network1[mask,]
nn <- as.network(network1)
# compute degree and betweeness
degreeMCMC[i,4,] <- degree(nn,gmod='graph') # node degree
betweenMCMC[i,4,] <- betweenness(nn,gmod='graph') # betweeness
	}

# network at occasion 5
xx <- x5[sselect,]
for (i in 1:nrow(xx)){ 
network1 <- label.dyads
mask <- as.logical(2 - xx[i,]) # in x1: 1 = activated, 2 = non-activated 
network1 <- network1[,-1]
network1 <- network1[mask,]
nn <- as.network(network1)
# compute degree and betweeness
degreeMCMC[i,5,] <- degree(nn,gmod='graph') # node degree
betweenMCMC[i,5,] <- betweenness(nn,gmod='graph') # betweeness
	}


save(degreeMCMC,betweenMCMC,file='dolphin_indivmetrics.Rdata')



#------------- plot measures per indiv 

qdegree1=matrix(NA,71,5)
qbetweeness1=matrix(NA,71,5)
for(i in 1:71){
	qdegree1[i,]=quantile(degreeMCMC[,1,i],c(.025,.25,.5,.75,.975))
	qbetweeness1[i,]=quantile(betweenMCMC[,1,i],c(.025,.25,.5,.75,.975))
	}
	
qdegree2=matrix(NA,71,5)
qbetweeness2=matrix(NA,71,5)
for(i in 1:71){
	qdegree2[i,]=quantile(degreeMCMC[,2,i],c(.025,.25,.5,.75,.975))
	qbetweeness2[i,]=quantile(betweenMCMC[,2,i],c(.025,.25,.5,.75,.975))
	}

qdegree3=matrix(NA,71,5)
qbetweeness3=matrix(NA,71,5)
for(i in 1:71){
	qdegree3[i,]=quantile(degreeMCMC[,3,i],c(.025,.25,.5,.75,.975))
	qbetweeness3[i,]=quantile(betweenMCMC[,3,i],c(.025,.25,.5,.75,.975))
	}

qdegree4=matrix(NA,71,5)
qbetweeness4=matrix(NA,71,5)
for(i in 1:71){
	qdegree4[i,]=quantile(degreeMCMC[,4,i],c(.025,.25,.5,.75,.975))
	qbetweeness4[i,]=quantile(betweenMCMC[,4,i],c(.025,.25,.5,.75,.975))
	}

qdegree5=matrix(NA,71,5)
qbetweeness5=matrix(NA,71,5)
for(i in 1:71){
	qdegree5[i,]=quantile(degreeMCMC[,5,i],c(.025,.25,.5,.75,.975))
	qbetweeness5[i,]=quantile(betweenMCMC[,5,i],c(.025,.25,.5,.75,.975))
	}

ind <- 1:71

png("fig2.png",res=300,width=20,height=12,unit='in') 

par(mfrow=c(2,5))

for (j in 1:5){

degree <- paste0('qdegree',j)
xl <- order(eval(parse(text=degree))[,3])
qdegree_ord <- eval(parse(text=degree))[xl,]
plot(ind, qdegree_ord[,3],type="n",xlab="individual",ylab="degree",main=paste0("occasion ",j),ylim=c(0,max(qdegree_ord[,3])+10),cex.main=2,cex.lab=1.5)
points(ind, qdegree_ord[,3],pch=19,cex=0.5)
for(i in ind)
{   lines(rep(ind[i],2), qdegree_ord[i,c(1,5)])
    lines(rep(ind[i],2), qdegree_ord[i,c(2,4)],lwd=2)
}
}

for (j in 1:5){
between <- paste0('qbetweeness',j)
xl <- order(eval(parse(text=between))[,3])
qbetweeness_ord <- eval(parse(text=between))[xl,]
plot(ind, qbetweeness_ord[,3],type="n",xlab="individual",ylab="betweeness",main=paste0("occasion ",j),ylim=c(0,max(qbetweeness_ord[,3])+10),cex.main=2,cex.lab=1.5)
points(ind, qbetweeness_ord[,3],pch=19,cex=0.5)
for(i in ind){   
	lines(rep(i,2), qbetweeness_ord[i,c(1,5)])
    lines(rep(i,2), qbetweeness_ord[i,c(2,4)],lwd=2)
}

}

dev.off()

#--- group metrics

aplMCMC <- matrix(0,nbMCMC,5)
clusMCMC <- matrix(0,nbMCMC,5)

# network at occasion 1
xx <- x1[sselect,]
for (i in 1:nrow(xx)){ 
network1 <- label.dyads
mask <- as.logical(2 - xx[i,]) # in x1: 1 = activated, 2 = non-activated 
network1 <- network1[,-1]
network1 <- network1[mask,]
nn <- as.network(network1)
write(t(network1),file='klause.txt',ncol=2)
g <- read.graph('klause.txt')
# compute average path length and clustering
aplMCMC[i,1] <- average.path.length(g)
clusMCMC[i,1] <- transitivity(g)
	}

# network at occasion 2
xx <- x2[sselect,]
for (i in 1:nrow(xx)){ 
network1 <- label.dyads
mask <- as.logical(2 - xx[i,]) # in x1: 1 = activated, 2 = non-activated 
network1 <- network1[,-1]
network1 <- network1[mask,]
nn <- as.network(network1)
write(t(network1),file='klause.txt',ncol=2)
g <- read.graph('klause.txt')
# compute average path length and clustering
aplMCMC[i,2] <- average.path.length(g)
clusMCMC[i,2] <- transitivity(g)
	}

# network at occasion 3
xx <- x3[sselect,]
for (i in 1:nrow(xx)){ 
network1 <- label.dyads
mask <- as.logical(2 - xx[i,]) # in x1: 1 = activated, 2 = non-activated 
network1 <- network1[,-1]
network1 <- network1[mask,]
nn <- as.network(network1)
write(t(network1),file='klause.txt',ncol=2)
g <- read.graph('klause.txt')
# compute average path length and clustering
aplMCMC[i,3] <- average.path.length(g)
clusMCMC[i,3] <- transitivity(g)
	}

# network at occasion 4
xx <- x4[sselect,]
for (i in 1:nrow(xx)){ 
network1 <- label.dyads
mask <- as.logical(2 - xx[i,]) # in x1: 1 = activated, 2 = non-activated 
network1 <- network1[,-1]
network1 <- network1[mask,]
nn <- as.network(network1)
write(t(network1),file='klause.txt',ncol=2)
g <- read.graph('klause.txt')
# compute average path length and clustering
aplMCMC[i,4] <- average.path.length(g)
clusMCMC[i,4] <- transitivity(g)
	}

# network at occasion 5
xx <- x5[sselect,]
for (i in 1:nrow(xx)){ 
network1 <- label.dyads
mask <- as.logical(2 - xx[i,]) # in x1: 1 = activated, 2 = non-activated 
network1 <- network1[,-1]
network1 <- network1[mask,]
nn <- as.network(network1)
write(t(network1),file='klause.txt',ncol=2)
g <- read.graph('klause.txt')
# compute average path length and clustering
aplMCMC[i,5] <- average.path.length(g)
clusMCMC[i,5] <- transitivity(g)
	}


save(aplMCMC,clusMCMC,file='dolphin_groupmetrics.Rdata')

# values given in Table 2

apply(aplMCMC,2,mean)
apply(aplMCMC,2,quantile,probs=c(2.5/100,50/100,97.5/100))


#> apply(aplMCMC,2,mean)
#[1] 1.308181 1.649835 1.613416 1.604478 1.609010
#> apply(aplMCMC,2,quantile,probs=c(2.5/100,50/100,97.5/100))
#          [,1]     [,2]     [,3]     [,4]     [,5]
#2.5%  1.245970 1.541453 1.569874 1.554238 1.558922
#50%   1.306445 1.647114 1.613120 1.603133 1.609237
#97.5% 1.380640 1.785085 1.658459 1.654106 1.662782


# values given in Table 2

apply(clusMCMC,2,mean)
apply(clusMCMC,2,quantile,probs=c(2.5/100,50/100,97.5/100))

#> apply(clusMCMC,2,mean)
#[1] 0.6804506 0.3607442 0.4199761 0.3894948 0.3959426
#> apply(clusMCMC,2,quantile,probs=c(2.5/100,50/100,97.5/100))
#           [,1]      [,2]      [,3]      [,4]      [,5]
#2.5%  0.6103371 0.2689322 0.3891264 0.3485224 0.3551654
#50%   0.6818992 0.3606645 0.4200122 0.3903437 0.3958728
#97.5% 0.7429786 0.4467331 0.4529634 0.4297261 0.4343316


