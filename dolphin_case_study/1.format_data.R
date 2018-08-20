# 1. format data

# load package to handle dates and time
library(chron)

# data frame with 3 columns: 1-individual ids, 2-date of capture and 3-time of capture
datos.2007 <- read.table('individual_capture_times_Cephalorhynchus_commersonii_2007.txt',sep=",")

# get capture history matrix
# more precisely, nb of captures per individual (row; 71) per sampling occasion (column; 5)
# then, replace nb of occurrence by 1's
HC.mab <- as.matrix(table(datos.2007[,1:2]))
HC.mab[HC.mab>0] <- 1

# nb of indentified individuals 
n.ind <- length(unique(datos.2007[,1]))
ID <- 1:n.ind

# vector containg all ids
ids.name <- as.character(unique(datos.2007[,1]))

# nb of sampling ocassions
oc <- length(unique(datos.2007[,2]))

#Number of dyads given the total number of indentified individuals
N.diad <- choose(n.ind,2)

# set the time interval (in minutes) for which a dyad is considered as associated
# following Klaich et 2011, time interval is is 20 minutes
tiempo.asoc <- 20

# empty matrix where to fill the association data 
MAB <- as.data.frame(matrix(rep(0,(2+oc)*N.diad),N.diad)) 

# columns names
names(MAB) <- c("ind.a","ind.b",paste("ta",1:oc,"",sep=""))
    
# fill the two columns containg the ids combinations for all possible dyads
ind.1 <- 0
ind.2 <- 0
lab.1 <- 0
lab.2 <- 0
for(i in 1:(n.ind-1)){
	ind.1 <- c(ind.1,ids.name[rep(ID[i],length(seq(ID[i+1],n.ind)))])
	ind.2 <- c(ind.2,ids.name[seq(ID[i+1],n.ind)])
	lab.1 <- c(lab.1,ID[rep(ID[i],length(seq(ID[i+1],n.ind)))])
	lab.2 <- c(lab.2,ID[seq(ID[i+1],n.ind)])
}
ind.1 <- ind.1[-1]
ind.2 <- ind.2[-1]
lab.1 <- lab.1[-1]
lab.2 <- lab.2[-1]
MAB[,1] <- ind.1
MAB[,2] <- ind.2
    
# number code for associated and non-associated states. 2:associated, 3:non-associated
estado <- c(3,2)

# auxiliar matrix for build the association data
captures.h <- as.data.frame(matrix(rep(0,sum(HC.mab)*max(table(datos.2007[,1]))),max(table(datos.2007[,1]))))
                
count <- 0
dates <- unique(datos.2007[,2])
ids.cap <- unique(datos.2007[,1])
for(i in 1:length(dates)){
  for(ii in 1:length(ids.cap)){
    filt <- datos.2007[,1]==ids.cap[ii]&datos.2007[,2]==dates[i]
    if(sum(filt)>0){  
      count <- count+1
      h.cap.id <- times(datos.2007[filt,3])
      captures.h[1:length(h.cap.id),count] <- h.cap.id
    }
  }
} 
    
hora.capt <- captures.h[1:max(apply(captures.h>0,2,sum)),]

# list containing the ids for each sampling ocassion
tt <- as.list(tapply(datos.2007[,1],datos.2007[,2],unique))

# number of identified individual for each sampling ocassion
n.cap <- as.numeric(apply(HC.mab,2,sum))

# code for association states assigning given the time interval  
for(ii in 1:length(n.cap)){
	aux <- hora.capt[,rep(c(1:oc),n.cap)==ii]
	ids.aux <- tt[[ii]]	
	for(iii in 1:dim(aux)[2]){
		vec1 <- aux[aux[,iii]>0,iii]
		for(iv in 1:dim(aux)[2]){
			vec2 <- aux[aux[,iv]>0,iv]
			filt.ids <- MAB[,1]==ids.aux[iii]&MAB[,2]==ids.aux[iv]
			if(sum(filt.ids)>0){
				v.rest.h <- rep(vec1,length(vec2))-rep(vec2,each=length(vec1))
				v.est <- sum(sum(abs(v.rest.h)<=tiempo.asoc/(24*60))>0)
				MAB[filt.ids,ii+2] <- estado[v.est+1]
				apply(rbind(rep(vec1,length(vec2)),rep(vec2,each=length(vec1))),2,mean)[abs(v.rest.h)<=tiempo.asoc/(24*60)]
			}
		}
	}
}

# assigning the numeric code 1 when only one individual was sighted
HC.ids.I <- as.data.frame(cbind(row.names(HC.mab),HC.mab))  
for(v in 3:(oc+2)){
	for(vv in 1:length(row.names(HC.mab))){
		filt.hc.a <- MAB[MAB$ind.a==row.names(HC.mab)[vv],v]==0
		filt.hc.b <- MAB[MAB$ind.b==row.names(HC.mab)[vv],v]==0
		MAB[MAB$ind.a==row.names(HC.mab)[vv],v][filt.hc.a]=as.numeric(as.character(HC.ids.I[as.character(HC.ids.I[,1])==row.names(HC.mab)[vv],v-1]))
		MAB[MAB$ind.b==row.names(HC.mab)[vv],v][filt.hc.b]=as.numeric(as.character(HC.ids.I[as.character(HC.ids.I[,1])==row.names(HC.mab)[vv],v-1]))
	}
}
MAB <- cbind(MAB,lab.1,lab.2)

ah.dyads <- MAB[,3:7]
mask0 <- (ah.dyads == 0) # dyad non-observed
mask1 <- (ah.dyads == 1) # one of the two individuals observed
mask2 <- (ah.dyads == 2) # associated
mask3 <- (ah.dyads == 3) # non-associated

ah.dyads[mask0] <- 1 # dyad non-observed
ah.dyads[mask1] <- 2 # one of the two individuals non-observed
ah.dyads[mask2] <- 3 # dyad seen and activated
ah.dyads[mask3] <- 4 # dyad seen and non-activated

ah.dyads.final <- cbind(1:nrow(ah.dyads),MAB[,8:9],ah.dyads,MAB[,1:2])

# save data
save(ah.dyads,ah.dyads.final,file='dolphin_dyads.RData')
