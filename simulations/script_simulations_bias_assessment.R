# simulations to assess bias in parameters estimates

# required packages
library(R2jags)

# function to simulate data
sim_CRnetwork <- function(J = 5, n = 105, pp = 0.7, psiAA = 0.3, psiBB = 0.8, pi = 0.7){
library(runjags)
# code to simulate with jags, note the use of the data block
# parameters for simulations 
# J = nb occasions
# n = nb of dyads (= N(N-1)/2 where N is the number of individuals); by default, we consider N = 15 individuals, hence N(N-1)/2 = 105 possible dyads
# pp = detection
# psiAA = pr of staying associated
# psiBB = pr of staying non-associated
# pi = initial state pr
txtstring <- '					
data{
	
# A = associated
# B = non-associated
	
#Pr(dyads state)
px[1,1] <- psiAA # probability of staying associated
px[1,2] <- 1 - psiAA # probability of associated -> non-associated
px[2,1] <- 1 - psiBB # probability of non-associated -> associated
px[2,2] <- psiBB # probability of staying non-associated

# Pr(dyads obs given dyads state)
## pp is the individual detection probability
	po[1,1] <- (1-pp) * (1-pp)
	po[1,2] <- 2 * pp * (1-pp)
	po[1,3] <- pp * pp
	po[1,4] <- 0
	po[2,1] <- (1-pp) * (1-pp)
	po[2,2] <- 2 * pp * (1-pp)
	po[2,3] <- 0
	po[2,4] <- pp * pp
			
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
		obs[i,1] ~ dcat(po[x[i,1],1:4])

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
			obs[i,t] ~ dcat(po[x[i,t],1:4])
							}
						}
}
model{
fake <- 0
}
'

# parameters are treated as data for the simulation step
data<-list(n=n, J=J, pp=pp, psiAA=psiAA, psiBB=psiBB, pi=pi)

# run jags
out <- run.jags(txtstring, data = data, monitor=c("obs","x"), sample=1, n.chains=1, summarise=FALSE)

# reformat the outputs
Simulated <- coda::as.mcmc(out)
#Simulated
#dim(Simulated)
dat <- matrix(Simulated[1:(n*J)],ncol=J)
#dat
states <- matrix(Simulated[-(1:(n*J))],ncol=J)
#states
list(dat=dat,states=states) # outputs: dat = detections/non-detections; states = underlying states
}

# specify model that will be used to estimate network parameters
sink("sim_network.txt")
cat("
model{

# Pr(dyads state)
px[1,1] <- psiAA 			# probability of staying associated
px[1,2] <- 1 - psiAA 		# probability of associated -> non-associated
px[2,1] <- 1 - psiBB 		# probability of non-associated -> associated
px[2,2] <- psiBB 			# probability of staying non-associated

# Pr(dyads obs given dyads state)
## pp is the individual detection probability
po[1,1] <- (1-pp) * (1-pp)
po[1,2] <- 2 * pp * (1-pp)
po[1,3] <- pp * pp
po[1,4] <- 0
po[2,1] <- (1-pp) * (1-pp)
po[2,2] <- 2 * pp * (1-pp)
po[2,3] <- 0
po[2,4] <- pp * pp
			
# Pr(initial states)
px0[1] <- pi 				# prob. of being in initial state A
px0[2] <- 1-pi 				# prob. of being in initial state B

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
	obs[i,1] ~ dcat(po[x[i,1],1:4])

	# for t > 1
	for (t in 2:J){
		#-- state equation # 1 = associated # 2 = non-associated
		x[i,t] ~ dcat(px[x[i,t-1],1:2]) 
		#-- observation equation (1,2,3 ó 4)
		obs[i,t] ~ dcat(po[x[i,t],1:4])
			}
		}

# Priors
pp ~ dunif(0,1) # detection pr
psiAA ~ dunif(0,1) # pr of staying associated
psiBB ~ dunif(0,1) # pr of staying non-associated
pi ~ dunif(0,1) # initial state pr
}
",fill=TRUE)
sink()


# now consider the following 36 scenarios
#scenarios on pp = c(0.3, 0.8) # (2)
#scenarios on pi = c(0.2, 0.7) # (2)
#scenarios on psiAA = c(0.1, 0.4, 0.9) # (3)
#scenarios on psiBB = c(0.1, 0.4, 0.9) # (3)
grid <- expand.grid(pp=c(0.3,0.8),pi=c(0.2,0.7),psiAA=c(0.1,0.4,0.9),psiBB=c(0.1,0.4,0.9))

# nb of monte carlo iterations
nb_simulations <- 100

# matrix to store results with estimated values for pp, psiAA, psiBB, pi
res <- array(NA,dim=c(nrow(grid),nb_simulations,4))

# run simulation
for (index in 1:nrow(grid)){ # go through grid of scenarios
		for (i in 1:nb_simulations){
 
# 1. simulate
 
pp <- grid[index,1]
pi <- grid[index,2]
psiAA <- grid[index,3]
psiBB <- grid[index,4]

sim_data <- sim_CRnetwork(J=5,n=105, pp = pp, psiAA = psiAA, psiBB = psiBB, pi = pi)
dat <- sim_data[[1]]
states <- sim_data[[2]]

# 2. estimation

# initial values
init1 <- list(psiAA=grid[index,3],pp=grid[index,1],x=states)
inits <- list(init1)
 
# data
jags.data <- list(obs = dat, n = nrow(dat), J = ncol(dat)) 
 
# nb iterations
ni <- 2000
# nb burn-in
nb <- 1000
# nb thin
nt <- 1
# nb chains
nc <- 1

# parameters to be monitored
parameters_sim <- c("psiAA","psiBB","pi","pp","x1","x2","x3","x4","x5")

# call JAGS from R
mod <- jags(jags.data, inits, parameters_sim, 'sim_network.txt', n.chains = nc, n.thin = nt, 
n.iter = ni, n.burnin = nb, working.directory = getwd())

res[index,i,1] <- mean(mod$BUGSoutput$sims.matrix[,'pp']) # detection
res[index,i,2] <- mean(mod$BUGSoutput$sims.matrix[,'psiAA']) # associated
res[index,i,3] <- mean(mod$BUGSoutput$sims.matrix[,'psiBB']) # non-associated
res[index,i,4] <- mean(mod$BUGSoutput$sims.matrix[,'pi']) # prop of associated
   
}
 
}

save(res,file='simul_network_index36_sim100.RData')
# load("simul_network_index36_sim100.RData")

# compute relative bias in percent in pp, psi11, psiBB and pi
bias_param <- matrix(NA,nrow(grid),4)
for(i in 1:nrow(grid)){
	for (j in 1:4){ 
		bias_param[i,j] <- (mean(res[i,,c(1,4,2,3)[j]]) - grid[i,j])/grid[i,c(1,3,4,2)[j]]*100
	}
}
res_bias <- round(cbind(1:nrow(bias_param),grid,bias_param),2)
colnames(res_bias) <- c('scenario',names(grid),'bias_pp','bias_pi','bias_psiAA','bias_psiBB')
res_bias

#   scenario  pp  pi psiAA psiBB bias_pp bias_pi bias_psiAA bias_psiBB
#1         1 0.3 0.2   0.1   0.1    0.50   26.49     120.98      58.77
#2         2 0.8 0.2   0.1   0.1    0.08    9.46       4.37       1.36
#3         3 0.3 0.7   0.1   0.1   -0.33   -1.23     142.04      22.83
#4         4 0.8 0.7   0.1   0.1   -0.21   -3.23       8.91      -0.04
#5         5 0.3 0.2   0.4   0.1    0.07   14.73      27.30      53.03
#6         6 0.8 0.2   0.4   0.1   -0.04    1.02      -1.65       4.74
#7         7 0.3 0.7   0.4   0.1    0.60  -10.96      65.19      26.16
#8         8 0.8 0.7   0.4   0.1   -0.04   -0.37      -8.88       0.79
#9         9 0.3 0.2   0.9   0.1    0.29    4.46     -23.10      37.30
#10       10 0.8 0.2   0.9   0.1    0.11    2.29      -5.57       7.26
#11       11 0.3 0.7   0.9   0.1   -0.25    0.30     -14.44      28.57
#12       12 0.8 0.7   0.9   0.1    0.07   -0.55      -7.99       3.91
#13       13 0.3 0.2   0.1   0.4   -0.74   54.58      45.20      24.95
#14       14 0.8 0.2   0.1   0.4   -0.08    6.23       2.19       4.71
#15       15 0.3 0.7   0.1   0.4    0.27  -25.83      29.36       7.66
#16       16 0.8 0.7   0.1   0.4   -0.11  -11.71       3.05       1.72
#17       17 0.3 0.2   0.4   0.4    0.45   14.96      10.59      21.80
#18       18 0.8 0.2   0.4   0.4   -0.09    3.22      -1.45      -0.27
#19       19 0.3 0.7   0.4   0.4    0.64  -13.37       5.67       5.96
#20       20 0.8 0.7   0.4   0.4    0.02    0.24      -1.44      -0.71
#21       21 0.3 0.2   0.9   0.4   -0.26    8.35     -17.84     -28.74
#22       22 0.8 0.2   0.9   0.4    0.01    1.28      -1.62      -1.72
#23       23 0.3 0.7   0.9   0.4    0.45   -1.59     -10.12      -5.75
#24       24 0.8 0.7   0.9   0.4   -0.08   -0.52      -2.47      -0.54
#25       25 0.3 0.2   0.1   0.9    0.94   38.86      21.21      -1.08
#26       26 0.8 0.2   0.1   0.9    0.08    8.48       2.90       0.87
#27       27 0.3 0.7   0.1   0.9    0.11  -47.67      10.35      -2.45
#28       28 0.8 0.7   0.1   0.9   -0.34    2.48       1.29      -0.87
#29       29 0.3 0.2   0.4   0.9   -0.46   11.66      -4.68     -16.83
#30       30 0.8 0.2   0.4   0.9   -0.22    2.55      -0.36      -1.68
#31       31 0.3 0.7   0.4   0.9   -0.27   -6.82      -7.96      -7.75
#32       32 0.8 0.7   0.4   0.9    0.04   -1.00      -0.86      -1.29
#33       33 0.3 0.2   0.9   0.9    1.18    3.33     -30.90     -55.94
#34       34 0.8 0.2   0.9   0.9    0.12    0.47      -1.45      -2.53
#35       35 0.3 0.7   0.9   0.9   -0.74   -3.30     -20.09     -38.39
#36       36 0.8 0.7   0.9   0.9   -0.16   -1.19      -0.85      -1.26


# from parameter estimates in Table 2 of the paper, the dolphin case study 
# is in between scenario 19 and scenario 31
# bias is negligible on detection, around +5 percent on the transition 
# probabilities and around -13 percent on pi in scenario 19 with psiBB = 0.4
# when psiBB = 0.9 in scenario 31, the bias in pi decreases by a factor 2
