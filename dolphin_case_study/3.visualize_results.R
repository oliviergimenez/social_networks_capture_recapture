# load packages
library(rjags)
library(network)
library(tidygraph)
library(ggraph)
library(png)
library(igraph)

# load parameter estimates
load('dolphin_model.RData')

# pool 2 chains
res <- rbind(jsample[[1]],jsample[[2]]) # join two lists of MCMC values
dim(res)
head(res)
colnames(res)

# get posterior summaries for model parameters
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

#---- plot network

## above which value we should keep the edge
#threshold <- 2400
#threshold/4000 # corresponding probability

probs <- NA

for (i in 1:5){

# get latent var estimates
mask <- grep(paste0("x",i),colnames(res))
x <- res[,mask]
x_associated <- apply(x==1,2,sum) # sum of 1's for each dyad
unique(x_associated)

# assign 1 to observed links
color_edge <- x_associated
color_edge <- ifelse(color_edge == 4000,1, color_edge)
unique(color_edge)

# assign 2 to links with non-null prob of existing
threshold <- quantile(x_associated,probs=90/100)
probs <- c(probs,threshold/4000)
color_edge <- ifelse(color_edge >= threshold,2,color_edge)
unique(color_edge)

# assign NA to all other links, ie missing links and links w/ prob < threshold
color_edge <- ifelse(color_edge==0,NA,color_edge)
unique(color_edge)
color_edge <- ifelse(color_edge>2,NA,color_edge)
unique(color_edge)

# get network
network <- label.dyads
network <- network[,-1]
mask <- !is.na(color_edge)
mon_graph <- data.frame(from = network[mask,1], to = network[mask,2], col = color_edge[mask])
dim(mon_graph)

# change default specs
nw <- graph_from_data_frame(mon_graph)
V(nw)$size <- 4
V(nw)$frame.color <- "white"
V(nw)$color <- "orange"
V(nw)$label <- "" 
E(nw)$arrow.mode <- 0
l <- layout_with_kk(nw)
#plot(nw,layout=l)

# observed edges are in grey
# I would like to use green 
# for edges that are estimated 
# to be existing but not observed

# compute how many such edges there are 
mask2 <- color_edge[mask] != 1
sum(mask2)

# generate edge color variable to plot the path
ecol <- rep("black", ecount(nw))
ecol[mask2] <- "green"
png(paste0("network",i,".png"),res=300,width=8,height=8,unit='in') 
plot(nw, edge.color=ecol,layout=l)
title(paste0('occasion ',i),cex.main=3)
box()
dev.off()

}


