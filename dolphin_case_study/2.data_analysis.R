# ah.dyads is a matrix with dyads in rows and sampling occasions in columns
# 1 = dyad non-observed, 
# 2 = one of the two individuals non-observed, 
# 3 = dyad seen and activated
# 4 = dyad seen and non-activated 

load('dolphin_dyads.RData')
ah.dyads
	
# model
sink("mnetwork.txt")
cat("
model{

# Pr(dyad's state)
px[1,1] <- psiAA # probability of staying associated
px[1,2] <- 1 - psiAA # probability of associated -> non-associated
px[2,1] <- 1 - psiBB # probability of non-associated -> associated
px[2,2] <- psiBB # probability of staying non-associated

# Pr(dyad's obs given dyad's state)
## pp is the individual detection probability
for (t in 1:J){
	po[1,1,t] <- (1-pp[t]) * (1-pp[t])
	po[1,2,t] <- 2 * pp[t] * (1-pp[t])
	po[1,3,t] <- pp[t] * pp[t]
	po[1,4,t] <- 0
	po[2,1,t] <- (1-pp[t]) * (1-pp[t])
	po[2,2,t] <- 2 * pp[t] * (1-pp[t])
	po[2,3,t] <- 0
	po[2,4,t] <- pp[t] * pp[t]
	}
			
# Pr(initial states)
px0[1] <- pi # prob. of being in initial state A
px0[2] <- 1-pi # prob. of being in initial state B

# Model likelihood
	for (i in 1:n){
		
		# record states for every sampling occasion
		x1[i] <- x[i,1]
		x2[i] <- x[i,2]
		x3[i] <- x[i,3]
		x4[i] <- x[i,4]
		x5[i] <- x[i,5]
		
		# for t = 1
		x[i,1] ~ dcat(px0[1:2])
		obs[i,1] ~ dcat(po[x[i,1],1:4,1])

		# for t > 1
		for (t in 2:J){
			
			#-- state equation
			# 1 = associated
			# 2 = non-associated
			x[i,t] ~ dcat(px[x[i,t-1],1:2]) 
			
			#-- observation equation
			# 1 = dyad non-observed, 
			# 2 = one of the two individuals non-observed, 
			# 3 = dyad seen and associated
			# 4 = dyad seen and non-associated 
			obs[i,t] ~ dcat(po[x[i,t],1:4,t])
							}
						}

# Priors
for (t in 1:J){
pp[t] ~ dunif(0,1) # detection pr
}
psiAA ~ dunif(0,1) # pr of staying associated
psiBB ~ dunif(0,1) # pr of staying non-associated
pi ~ dunif(0,1) # initial state pr
}
",fill=TRUE)
sink()

# load rjags package
library(rjags)

# nb iterations
ni=2000
# nb burn-in
nb=500
# nb thin
nt=1
# nb chains to be run in //
nc=2

# build list of data
obs <- ah.dyads
J <- ncol(obs) # nb capture occasions
n <- nrow(obs) # nb dyads
mydatax = list(J=J,n=n,obs=obs)

# initial values
x <- as.matrix(ah.dyads)
mask1 <- which(x==3)
x[mask1] <- 1
mask2 <- which(x==4)
x[mask2] <- 2
init1 = list(psiAA=0.3,pp=runif(J,0,1),x=x)
init2 = list(psiAA=0.8,pp=runif(J,0,1),x=x)
inits = list(init1,init2)

# parameters to be monitored
parameters <- c("psiAA","psiBB","pi","pp","x1","x2","x3","x4","x5")

# run jags
start<-as.POSIXlt(Sys.time())
jmodel <- jags.model("mnetwork.txt", mydatax, inits, n.chains = nc,n.adapt = nb)
jsample <- coda.samples(jmodel, parameters, n.iter=ni, thin = nt)
#m.dic <- dic.samples(jmodel,n.iter=ni)
end <-as.POSIXlt(Sys.time())
duration = end-start

save(jmodel,jsample,duration,file="dolphin_model.RData")

