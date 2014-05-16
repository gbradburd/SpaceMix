#UPDATE FUNCTIONS ETC.
if(TRUE){
simulate.spacemix.dataset <- function(k,loci,admix.target,admix.source,admix.proportion,sim.a0,sim.aD,sim.a2,sample.sizes,nugget.option,generate.counts,boundary.option,filename){
	# recover()
	if(nugget.option == "sample.sizes" && generate.counts){
		stop("you are double-dipping your nugget!")
	}
	#set random seed
		random.seed <- sample(100:999,1)
			set.seed(random.seed)
	#simulate population locations
		sim.locations <- cbind(runif(2*k,-1,1),runif(2*k,-1,1))
		sim.locations[k+admix.target,] <- sim.locations[admix.source,]
		observed.X.coordinates <- sim.locations[(1:k),1]
		observed.Y.coordinates <- sim.locations[(1:k),2]
	#generate spatial covariance matrix
		sim.D <- fields::rdist(sim.locations)
		mean.sample.sizes <- rowMeans(sample.sizes)
		sim.covariance <- Covariance(sim.a0,sim.aD,sim.a2,sim.D)
	#generate admixed covariance matrix
		if(nugget.option == "none"){
			sim.nugget <- numeric(k)
		} else if(nugget.option == "sample.sizes"){
			sim.nugget <- 1/mean.sample.sizes
		}
		sim.admix.proportions <- numeric(k)
			sim.admix.proportions[admix.target] <- admix.proportion
		sim.admixed.covariance <- admixed.Covariance(sim.covariance,sim.admix.proportions,sim.nugget)
	#simulate MVN allele frequencies + counts
		sim.MVN.draws <- t(MASS::mvrnorm(n = loci, mu = numeric(k), Sigma = sim.admixed.covariance))
		sim.allele.freqs <- sim.MVN.draws + 0.5
		reduced.sample.sizes <- sample.sizes
		if(boundary.option == "drop"){
			loci.2.drop <- unique(c(unique(which(sim.allele.freqs>1,arr.ind=TRUE)[,2]), unique(which(sim.allele.freqs < 0,arr.ind=TRUE)[,2])))		
			if(sum(loci.2.drop > 0)){
				sim.allele.freqs <- sim.allele.freqs[,-loci.2.drop]
				reduced.sample.sizes <- sample.sizes[,-loci.2.drop]
			}
			cat("reduced loci = ",ncol(reduced.sample.sizes),"\n")
		} else if(boundary.option == "absorb"){
			sim.allele.freqs[which(sim.allele.freqs < 0)] <- 0
			sim.allele.freqs[which(sim.allele.freqs > 1)] <- 1
		}
		reduced.loci <- ncol(sim.allele.freqs)	
	if(generate.counts){
		sim.counts <- matrix(rbinom(n = (k*reduced.loci), size = sample.sizes, prob = sim.allele.freqs),nrow=k,ncol=reduced.loci)
	}
	save(list=objects(),file=paste(filename,".Robj",sep=""))
	return(0)
}

admix_target_location_and_nugget_gibbs_sampler <- function(last.params){
	new.params <- last.params
	pop.to.update <- sample(1:last.params$k,1)
	x.min <- min(last.params$population.coordinates[1:last.params$k,1]) - diff(range(last.params$population.coordinates[1:last.params$k,1]))*0.15
	x.max <- max(last.params$population.coordinates[1:last.params$k,1]) + diff(range(last.params$population.coordinates[1:last.params$k,1]))*0.15
	y.min <- min(last.params$population.coordinates[1:last.params$k,2]) - diff(range(last.params$population.coordinates[1:last.params$k,2]))*0.15
	y.max <- max(last.params$population.coordinates[1:last.params$k,2]) + diff(range(last.params$population.coordinates[1:last.params$k,2]))*0.15
	X.gridpoints <- seq(x.min,x.max,length.out=last.params$X.grid.fineness)
	Y.gridpoints <- seq(y.min,y.max,length=last.params$Y.grid.fineness)
	nugget.gridpoints <- -log(seq(1e-10,1,length.out=last.params$nugget.grid.fineness))
		nugget.gridpoints[which(nugget.gridpoints == 0)] <- 1e-5	
	lnL.array <- array(0,dim=c(last.params$X.grid.fineness,last.params$Y.grid.fineness,last.params$nugget.grid.fineness))
	prior.prob.nugget.array <- array(0,dim=c(last.params$X.grid.fineness,last.params$Y.grid.fineness,last.params$nugget.grid.fineness))
	prior.prob.admix.target.locations.array <- array(0,dim=c(last.params$X.grid.fineness,last.params$Y.grid.fineness,last.params$nugget.grid.fineness))
	post.prob.array <- array(0,dim=c(last.params$X.grid.fineness,last.params$Y.grid.fineness,last.params$nugget.grid.fineness))
	coords.prime <- matrix(last.params$population.coordinates,nrow=2*last.params$k,ncol=2)
	nugget.prime <- last.params$nugget
	covariance.prime <- last.params$covariance
	for(x in 1:length(X.gridpoints)){
		for(y in 1:length(Y.gridpoints)){
				coords.prime[pop.to.update,] <- c(X.gridpoints[x],Y.gridpoints[y])
				prior.prob.admix.target.locations.array[x,y,] <- Prior_prob_admix_target_locations(coords.prime[1:last.params$k,],last.params$observed.X.coordinates,last.params$observed.Y.coordinates,last.params$target.spatial.prior.scale)
				tmp.cov <- Covariance(last.params$a0,
									  last.params$aD,
								  	  last.params$a2,
									  fields::rdist(coords.prime[pop.to.update,1:2,drop=FALSE],coords.prime))
				covariance.prime[pop.to.update,] <- tmp.cov
				covariance.prime[,pop.to.update] <- tmp.cov
			for(z in 1:length(nugget.gridpoints)){
				nugget.prime[pop.to.update] <- nugget.gridpoints[z]
				admixed.covariance.prime <- admixed.Covariance(covariance.prime,last.params$admix.proportions,nugget.prime)
				transformed.covariance.prime <- transformed.Covariance(admixed.covariance.prime,last.params$projection.matrix)
				lnL.array[x,y,z] <- wishart.lnL(last.params$sample.cov,transformed.covariance.prime/last.params$loci,last.params$loci)
				prior.prob.nugget.array[x,y,z] <- Prior_prob_nugget(nugget.prime,last.params$mean.sample.sizes)
			}
		}
	}
	post.prob.array <- lnL.array + prior.prob.nugget.array + prior.prob.admix.target.locations.array
	current.lnl.pr <- last.params$LnL_freqs + last.params$prior_prob_nugget + last.params$prior_prob_admix_proportions
	tmp.prob <- c(post.prob.array,current.lnl.pr)
	sampling.probs <- exp(tmp.prob-max(tmp.prob))/sum(exp(tmp.prob-max(tmp.prob)))
	tmp.sampled.index <- sample(c(1:length(sampling.probs)),1,prob=sampling.probs)
	if(tmp.sampled.index != length(sampling.probs)){
		sampled.parameters.index <- which(post.prob.array == post.prob.array[tmp.sampled.index],arr.ind=TRUE)
		new.params$population.coordinates[pop.to.update,] <- c(X.gridpoints[sampled.parameters.index[1]],Y.gridpoints[sampled.parameters.index[2]])
		new.params$prior_prob_admix_target_locations <- prior.prob.admix.target.locations.array[sampled.parameters.index]
		new.params$nugget[pop.to.update] <- nugget.gridpoints[sampled.parameters.index[3]]
		new.params$prior_prob_nugget <- prior.prob.nugget.array[sampled.parameters.index]
		new.params$D <- fields::rdist(new.params$population.coordinates)
		new.params$covariance <- Covariance(last.params$a0,last.params$aD,last.params$a2,new.params$D)
		new.params$admixed.covariance <- admixed.Covariance(new.params$covariance,last.params$admix.proportions,new.params$nugget)
		new.params$transformed_covariance <- transformed.Covariance(new.params$admixed.covariance,last.params$projection.matrix)
		new.params$LnL_freqs <- lnL.array[sampled.parameters.index]
		new.params$admix_target_location_moves[pop.to.update] <- new.params$admix_target_location_moves[pop.to.update] + 1
		new.params$admix_target_location_accept[pop.to.update] <- new.params$admix_target_location_accept[pop.to.update] + 1
		new.params$admix_target_location_accept_rate[pop.to.update] <- new.params$admix_target_location_accept[pop.to.update]/new.params$admix_target_location_moves[pop.to.update]
		new.params$nugget_moves[pop.to.update] <- new.params$nugget_moves[pop.to.update] + 1
		new.params$nugget_accept[pop.to.update] <- new.params$nugget_accept[pop.to.update] + 1
		new.params$nugget_accept_rate[pop.to.update] <- new.params$nugget_accept[pop.to.update]/new.params$nugget_moves[pop.to.update]
	}
		return(new.params)
}

Update_admixture_target_location <- function(last.params){
	# recover()
	new.params <- last.params
	pop.to.update <- sample(1:last.params$k,1)
		population.coordinates_prime <- last.params$population.coordinates
		population.coordinates_prime[pop.to.update,] <- propose.new.location(population.coordinates_prime[pop.to.update,1],
																			population.coordinates_prime[pop.to.update,2],
																			exp(last.params$admix.target.location.lstp[pop.to.update]))
		prior_prob_admix_target_locations_prime <- Prior_prob_admix_target_locations(population.coordinates_prime[1:last.params$k,],
																						last.params$observed.X.coordinates,
																						last.params$observed.Y.coordinates,
																						last.params$target.spatial.prior.scale)
		D_prime <- fields::rdist(population.coordinates_prime)
		covariance_prime <- Covariance(last.params$a0,last.params$aD,last.params$a2,D_prime)
		admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions,last.params$nugget)
				transformed_covariance_prime <- transformed.Covariance(admixed.covariance_prime,
																		last.params$projection.matrix)
				LnL_freqs_prime <- wishart.lnL(last.params$sample.covariance,transformed_covariance_prime/last.params$loci,last.params$loci)
					if(metropolis_ratio(LnL_freqs_prime,prior_prob_admix_target_locations_prime,last.params$LnL_freqs,last.params$prior_prob_admix_target_locations)){
						new.params$population.coordinates <- population.coordinates_prime
						new.params$prior_prob_admix_target_locations <- prior_prob_admix_target_locations_prime
						new.params$D <- D_prime
						new.params$covariance <- covariance_prime
						new.params$admixed.covariance <- admixed.covariance_prime
						new.params$transformed_covariance <- transformed_covariance_prime
						new.params$LnL_freqs <- LnL_freqs_prime
						new.params$admix_target_location_accept[pop.to.update] <- new.params$admix_target_location_accept[pop.to.update] + 1
					}
	new.params$admix_target_location_moves[pop.to.update] <- new.params$admix_target_location_moves[pop.to.update] + 1
	new.params$admix_target_location_accept_rate[pop.to.update] <- new.params$admix_target_location_accept[pop.to.update]/new.params$admix_target_location_moves[pop.to.update]
	return(new.params)	
}

	
Update_admixture_source_location <- function(last.params){
	 # recover()
	new.params <- last.params
	pop.to.update <- sample(1:last.params$k,1)
		population.coordinates_prime <- last.params$population.coordinates
		population.coordinates_prime[pop.to.update + last.params$k,] <- propose.new.location(population.coordinates_prime[pop.to.update + last.params$k,1],
																							population.coordinates_prime[pop.to.update + last.params$k,2],
																							exp(last.params$admix.source.location.lstp[pop.to.update]))
	prior_prob_admix_source_locations_prime <- Prior_prob_admix_source_locations(population.coordinates_prime[(last.params$k+1):(2*last.params$k),],
																					last.params$centroid,
																					last.params$source.spatial.prior.scale)
	D_prime <- fields::rdist(population.coordinates_prime)
	covariance_prime <- Covariance(last.params$a0,last.params$aD,last.params$a2,D_prime)
	admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions,last.params$nugget)
			transformed_covariance_prime <- transformed.Covariance(admixed.covariance_prime,
																	last.params$projection.matrix)
				LnL_freqs_prime <- wishart.lnL(last.params$sample.covariance,transformed_covariance_prime/last.params$loci,last.params$loci)
				if(metropolis_ratio(LnL_freqs_prime,prior_prob_admix_source_locations_prime,last.params$LnL_freqs,last.params$prior_prob_admix_source_locations)){
					new.params$population.coordinates <- population.coordinates_prime
					new.params$prior_prob_admix_source_locations <- prior_prob_admix_source_locations_prime
					new.params$D <- D_prime
					new.params$covariance <- covariance_prime
					new.params$admixed.covariance <- admixed.covariance_prime
					new.params$transformed_covariance <- transformed_covariance_prime
					new.params$LnL_freqs <- LnL_freqs_prime
					new.params$admix_source_location_accept[pop.to.update] <- new.params$admix_source_location_accept[pop.to.update] + 1
				}
	new.params$admix_source_location_moves[pop.to.update] <- new.params$admix_source_location_moves[pop.to.update] + 1
	new.params$admix_source_location_accept_rate[pop.to.update] <- new.params$admix_source_location_accept[pop.to.update]/new.params$admix_source_location_moves[pop.to.update]
	return(new.params)	
}

Update_admixture_proportions <- function(last.params){
	new.params <- last.params
	pop.to.update <- sample(last.params$k,1)
	admix.proportions_prime <- last.params$admix.proportions
	admix.proportions_prime[pop.to.update] <- admix.proportions_prime[pop.to.update] + rnorm(1,0,exp(last.params$admix.proportions.lstp[pop.to.update]))
	prior_prob_admix_proportions_prime <- Prior_prob_admix_proportions(admix.proportions_prime)
	if(is.finite(prior_prob_admix_proportions_prime)){
		admixed.covariance_prime <- admixed.Covariance(last.params$covariance,admix.proportions_prime,last.params$nugget)
				transformed_covariance_prime <- transformed.Covariance(admixed.covariance_prime,
																		last.params$projection.matrix)
				LnL_freqs_prime <- wishart.lnL(last.params$sample.covariance,transformed_covariance_prime/last.params$loci,last.params$loci)
					if(metropolis_ratio(LnL_freqs_prime,prior_prob_admix_proportions_prime,last.params$LnL_freqs,last.params$prior_prob_admix_proportions)){
						new.params$admix.proportions <- admix.proportions_prime
						new.params$prior_prob_admix_proportions <- prior_prob_admix_proportions_prime
						new.params$admixed.covariance <- admixed.covariance_prime
						new.params$transformed_covariance <- transformed_covariance_prime
						new.params$LnL_freqs <- LnL_freqs_prime
						new.params$admix_proportions_accept[pop.to.update] <- new.params$admix_proportions_accept[pop.to.update] + 1
					}
			}
	new.params$admix_proportions_moves[pop.to.update] <- new.params$admix_proportions_moves[pop.to.update] + 1
	new.params$admix_proportions_accept_rate[pop.to.update] <- new.params$admix_proportions_accept[pop.to.update]/new.params$admix_proportions_moves[pop.to.update]
	return(new.params)
}

Update_a0 <- function(last.params){
	new.params <- last.params
	a0_prime <- last.params$a0 + rnorm(1,0,exp(last.params$a0.lstp))
	prior_prob_alpha0_prime <- Prior_prob_alpha0(a0_prime)
	if(prior_prob_alpha0_prime != -Inf){
		covariance_prime <- Covariance(a0_prime,last.params$aD,last.params$a2,last.params$D)
			admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions,last.params$nugget)
			transformed_covariance_prime <- transformed.Covariance(admixed.covariance_prime,last.params$projection.matrix)
				LnL_freqs_prime <- wishart.lnL(last.params$sample.covariance,transformed_covariance_prime/last.params$loci,last.params$loci)
				if(metropolis_ratio(LnL_freqs_prime,prior_prob_alpha0_prime,last.params$LnL_freqs,last.params$prior_prob_alpha0)){
					new.params$a0 <- a0_prime
					new.params$covariance <- covariance_prime
					new.params$admixed.covariance <- admixed.covariance_prime
					new.params$transformed_covariance <- transformed_covariance_prime					
					new.params$prior_prob_alpha0 <- prior_prob_alpha0_prime
					new.params$LnL_freqs <- LnL_freqs_prime
					new.params$a0_accept <- new.params$a0_accept + 1
				}
	}
	new.params$a0_moves <- new.params$a0_moves + 1
	new.params$a0_accept_rate <- new.params$a0_accept/new.params$a0_moves
	return(new.params)
}

Update_aD <- function(last.params){
	new.params <- last.params
	aD_prime <- last.params$aD + rnorm(1,0,exp(last.params$aD.lstp))
	prior_prob_alphaD_prime <- Prior_prob_alphaD(aD_prime)
	if(prior_prob_alphaD_prime != -Inf){
		covariance_prime <- Covariance(last.params$a0,aD_prime,last.params$a2,last.params$D)
			admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions,last.params$nugget)
			transformed_covariance_prime <- transformed.Covariance(admixed.covariance_prime,last.params$projection.matrix)
				LnL_freqs_prime <- wishart.lnL(last.params$sample.covariance,transformed_covariance_prime/last.params$loci,last.params$loci)
				if(metropolis_ratio(LnL_freqs_prime,prior_prob_alphaD_prime,last.params$LnL_freqs,last.params$prior_prob_alphaD)){
					new.params$aD <- aD_prime
					new.params$covariance <- covariance_prime
					new.params$admixed.covariance <- admixed.covariance_prime
					new.params$transformed_covariance <- transformed_covariance_prime					
					new.params$prior_prob_alphaD <- prior_prob_alphaD_prime
					new.params$LnL_freqs <- LnL_freqs_prime
					new.params$aD_accept <- new.params$aD_accept + 1
				}
	}
	new.params$aD_moves <- new.params$aD_moves + 1
	new.params$aD_accept_rate <- new.params$aD_accept/new.params$aD_moves
	return(new.params)
}

Update_a2 <- function(last.params){
		new.params <- last.params
		a2_prime <- last.params$a2 + rnorm(1,0,exp(last.params$a2.lstp))
		prior_prob_alpha2_prime <- Prior_prob_alpha2(a2_prime) 
		if(prior_prob_alpha2_prime != -Inf) {
			covariance_prime <- Covariance(last.params$a0,last.params$aD,a2_prime,last.params$D)
				admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions,last.params$nugget)
				transformed_covariance_prime <- transformed.Covariance(admixed.covariance_prime,last.params$projection.matrix)
				LnL_freqs_prime <- wishart.lnL(last.params$sample.covariance,transformed_covariance_prime/last.params$loci,last.params$loci)
					if(metropolis_ratio(LnL_freqs_prime,prior_prob_alpha2_prime,last.params$LnL_freqs,last.params$prior_prob_alpha2)){
						new.params$a2 <- a2_prime
						new.params$covariance <- covariance_prime
						new.params$admixed.covariance <- admixed.covariance_prime						
						new.params$transformed_covariance <- transformed_covariance_prime							
						new.params$prior_prob_alpha2 <- prior_prob_alpha2_prime
						new.params$LnL_freqs <- LnL_freqs_prime
						new.params$a2_accept <- new.params$a2_accept + 1
					}
		}
	new.params$a2_moves <- new.params$a2_moves + 1
	new.params$a2_accept_rate <- new.params$a2_accept/new.params$a2_moves
	return(new.params)
}

Update_nugget <- function(last.params){
	new.params <- last.params
		pop.to.update <- sample(1:last.params$k,1)
	nugget_prime <- last.params$nugget + c(rep(0,pop.to.update-1),rnorm(1,0,exp(last.params$nugget.lstp[pop.to.update])),rep(0,last.params$k-pop.to.update))
	prior_prob_nugget_prime <- Prior_prob_nugget(nugget_prime,last.params$mean.sample.sizes)
	if(prior_prob_nugget_prime != -Inf){
		covariance_prime <- Covariance(last.params$a0,last.params$aD,last.params$a2,last.params$D)
			admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions,nugget_prime)
			transformed_covariance_prime <- transformed.Covariance(admixed.covariance_prime,last.params$projection.matrix)
				LnL_freqs_prime <- wishart.lnL(last.params$sample.covariance,transformed_covariance_prime/last.params$loci,last.params$loci)
				if(metropolis_ratio(LnL_freqs_prime, prior_prob_nugget_prime,last.params$LnL_freqs,last.params$prior_prob_nugget)){
					new.params$nugget <- nugget_prime
					new.params$covariance <- covariance_prime
					new.params$admixed.covariance <- admixed.covariance_prime
					new.params$transformed_covariance <- transformed_covariance_prime
					new.params$prior_prob_nugget <- prior_prob_nugget_prime
					new.params$LnL_freqs <- LnL_freqs_prime
					new.params$nugget_accept[pop.to.update] <- new.params$nugget_accept[pop.to.update] + 1
				}
	}
	new.params$nugget_moves[pop.to.update] <- new.params$nugget_moves[pop.to.update] + 1
	new.params$nugget_accept_rate[pop.to.update] <- new.params$nugget_accept[pop.to.update]/new.params$nugget_moves[pop.to.update]
	return(new.params)
}

metropolis_ratio <- function(LnL_prime,prior_prime,LnL,prior){
	accept <- FALSE
	if( exp((LnL_prime + prior_prime) - (LnL+prior)) >= runif(1) ){
		accept <- TRUE
	}
	return(accept)
}

Prior_prob_admix_target_locations <- function(admix_target_locations,observed.X.coordinates,observed.Y.coordinates,target.spatial.prior.scale){
	prior_prob_admix_target_locations <- sum(-log(2*pi*target.spatial.prior.scale) + (-1/2) * 
			((admix_target_locations[,1] - observed.X.coordinates)^2 / target.spatial.prior.scale + 
				(admix_target_locations[,2] - observed.Y.coordinates)^2 / target.spatial.prior.scale))
	return(prior_prob_admix_target_locations)
}

Prior_prob_admix_source_locations <- function(admix_source_locations,centroid,source.spatial.prior.scale){
	prior_prob_admix_source_locations <- sum(-log(2*pi*source.spatial.prior.scale) + (-1/2) * 
			((admix_source_locations[,1] - centroid[1])^2 / source.spatial.prior.scale + 
				(admix_source_locations[,2] - centroid[2])^2 / source.spatial.prior.scale))
	return(prior_prob_admix_source_locations)
}

Prior_prob_alpha0 <- function(a0){
	log(dexp(a0,rate=1/100))
	# log(dgamma(a0,1,1))
}

Prior_prob_alphaD <- function(aD){
	log(dexp(aD,rate=1))
}


Prior_prob_alpha2 <- function(a2){
	log(dunif(a2,0.1,2))
}

Prior_prob_nugget <- function(nugget,mean.sample.sizes){
	sum(dexp(nugget,rate=mean.sample.sizes,log=TRUE))
}

Prior_prob_admix_proportions <- function(admix_proportions){
	sum(dbeta(admix_proportions,shape1=1,shape2=50,log=TRUE))
}

wishart.lnL <- function(sample.cov,par.cov,n){
	A <- solve(par.cov)
	lnL <- -0.5 * sum( A * sample.cov ) - (n/2)*determinant(par.cov,logarithm=TRUE)$modulus
	return(lnL)
}

total_likelihood_freqs <- function(freqs,covmat,loci) {
	cholcov <- chol(covmat)
	logsqrtdet <- sum(log(diag(cholcov)))
	x <- backsolve(cholcov,freqs,transpose=TRUE)
	return(-(1/2)*crossprod(as.vector(x))-loci*logsqrtdet)
}
	
Covariance <- function(a0,aD,a2,GeoDist) {
	covariance <- (1/a0)*exp(-(aD*GeoDist)^a2)
	return(covariance)
}

admixed.Covariance <- function(covariance,admix.proportions,nugget){
	# recover()
	admix.proportions <- admix.proportions/2
	k <- nrow(covariance)/2
	w_k <- matrix(admix.proportions,nrow=k,ncol=1)
	admixed.Covariance <- 	tcrossprod((1-w_k),(1-w_k)) * 	covariance[1:k,1:k] + 
							tcrossprod((1-w_k),(w_k)) 	* 	covariance[1:k,(k+1):(2*k)] +
							tcrossprod(w_k,(1-w_k)) 	*	covariance[(k+1):(2*k),1:k] +
							tcrossprod(w_k,w_k)			*	covariance[(k+1):(2*k),(k+1):(2*k)]
	diag(admixed.Covariance) <- diag(admixed.Covariance) + nugget
	return(admixed.Covariance)
}

transformed.Covariance <- function(covariance,projection.matrix){
	transformed.covariance <- 
		crossprod(projection.matrix,covariance) %*% projection.matrix
	return(transformed.covariance)		
}

get.mean.sample.size <- function(sample.sizes){
	mean.sample.size <- mean(sample.sizes[which(sample.sizes!=0)])
	return(mean.sample.size)
}

get.weighted.mean.frequency <- function(sample.frequencies,mean.sample.sizes){
	na.pops <- which(is.na(sample.frequencies))
	if(sum(na.pops) > 0){
		sample.frequencies <- sample.frequencies[-na.pops]
		mean.sample.sizes <- mean.sample.sizes[-na.pops]
	}
	weighted.sample.frequencies <- mean.sample.sizes*sample.frequencies
	sample.frequency.weighted.mean <- sum(weighted.sample.frequencies)/sum(mean.sample.sizes)
	return(sample.frequency.weighted.mean)
}

curate.count.data <- function(counts,sample.sizes,prefix){
	k <- nrow(counts)
	original.loci <- ncol(counts)
	sample.frequencies <- counts/sample.sizes
		no.samples <- which(colSums(sample.sizes)==0)
		fixed.alleles <- c(which(colSums(sample.frequencies)==0),which(colSums(sample.frequencies)==k))
		na.loci <- which(is.na(sample.frequencies),arr.ind=TRUE)[,2]
		loci.to.drop <- unique(c(no.samples,fixed.alleles,na.loci))
	if(sum(loci.to.drop) > 0){	
		counts <- counts[,-loci.to.drop]
		sample.sizes <- sample.sizes[,-loci.to.drop]
		sample.frequencies <- sample.frequencies[,-loci.to.drop]
	}	
	#get mean sample sizes, round 1
		mean.sample.sizes <- apply(sample.sizes,1,get.mean.sample.size)
	#get mean frequencies across loci
		mean.frequencies <- apply(sample.frequencies,2,get.weighted.mean.frequency,mean.sample.sizes=mean.sample.sizes)
	reduced.loci <- length(mean.frequencies)
		if(reduced.loci==0){
			stop("you have no loci remaining in your curated dataset.")
		}
	cat("\t",reduced.loci,"loci out of the original",original.loci,"left in curated dataset. \n")
	curated.count.data <- list(	"reduced.loci" = reduced.loci,
					 		"counts" = counts,
							"sample.sizes" = sample.sizes,
							"sample.frequencies" = sample.frequencies,
							"mean.sample.sizes" = mean.sample.sizes,
							"mean.frequencies" = mean.frequencies)
	save(curated.count.data,file=paste(prefix,"curated.count.data.Robj",sep=''))
	return(curated.count.data)
}	

curate.frequency.data <- function(sample.frequencies){
	k <- nrow(sample.frequencies)
	fixed.alleles <- c(which(colSums(sample.frequencies)==0),which(colSums(sample.frequencies)==k))
	na.loci <- which(is.na(sample.frequencies),arr.ind=TRUE)[,2]
	loci.to.drop <- unique(c(fixed.alleles,na.loci))
	if(sum(loci.to.drop) > 0){	
		sample.frequencies <- sample.frequencies[,-loci.to.drop]
	}
	if(ncol(sample.frequencies) == 0){
		stop("you have no loci remaining in your curated dataset.")
	}
	return(sample.frequencies)
}

mean.center.normalize.frequencies <- function(sample.frequencies,mean.sample.sizes,prefix){
	# recover()
	k <- nrow(sample.frequencies)
	#mean-center, normalize, and project sample frequencies
		mean.frequencies <- apply(sample.frequencies,2,get.weighted.mean.frequency,mean.sample.sizes=mean.sample.sizes)
		mean.freq.mat <- matrix(mean.frequencies,nrow=k,ncol=length(mean.frequencies),byrow=TRUE)
		normalized.sample.frequencies <- sample.frequencies / sqrt(mean.freq.mat * (1 - mean.freq.mat))
		mean.centered.sample.frequencies <- sample.frequencies - mean.freq.mat
		mean.centered.normalized.sample.frequencies <- mean.centered.sample.frequencies / sqrt(mean.freq.mat * (1 - mean.freq.mat))
	MCN.frequencies.list <- list(	"mean.sample.sizes" = mean.sample.sizes,
											"sample.frequencies" = sample.frequencies,
											"normalized.sample.frequencies" = normalized.sample.frequencies,
											"mean.centered.sample.frequencies" = mean.centered.sample.frequencies,
											"mean.centered.normalized.sample.frequencies" = mean.centered.normalized.sample.frequencies)
			save(MCN.frequencies.list,file=paste(prefix,"MCN.frequencies.list.Robj",sep=''))
	return(mean.centered.normalized.sample.frequencies)
}

get.projection.matrix <- function(mean.sample.sizes,proj.mat.option=NULL){
	if(is.null(proj.mat.option)){
		k <- length(mean.sample.sizes)
		transformation.matrix <- diag(k) - matrix(mean.sample.sizes/(sum(mean.sample.sizes)),nrow=k,ncol=k,byrow=TRUE)
			qr.transformation.matrix <- qr(t(transformation.matrix))
			projection.matrix <- qr.Q(qr.transformation.matrix)[,1:qr.transformation.matrix$rank]
			stopifnot(qr.transformation.matrix$rank == sum(abs(eigen(transformation.matrix)$values - 1) < 1e-2) )
	} else {
		projection.matrix <- diag(length(mean.sample.sizes))
	}
	return(projection.matrix)
}

spacemix.data <- function(data.type,proj.mat.option=NULL,sample.frequencies=NULL,loci,mean.sample.sizes=NULL,counts=NULL,sample.sizes=NULL,sample.covariance=NULL,prefix){
	# recover()
	if(data.type == "sample.covariance"){
		if(!matrixcalc::is.positive.definite(sample.covariance)){
			stop("the input sample.covariance matrix is not a valid covariance matrix")
		}
		if(is.null(mean.sample.sizes)){
			stop("you must specify the mean population sample sizes")
		}
		sample.frequencies <- NULL
		mean.centered.normalized.frequencies <- NULL
	} else if(data.type == "sample.frequencies"){
		if(is.null(mean.sample.sizes)){
			stop("you must specify the mean population sample sizes")
		}
		sample.frequencies <- curate.frequency.data(sample.frequencies)
		mean.centered.normalized.frequencies <- mean.center.normalize.frequencies(sample.frequencies,mean.sample.sizes,prefix)
		sample.covariance <- cov(t(mean.centered.normalized.frequencies))
		loci <- ncol(mean.centered.normalized.frequencies)
	} else if(data.type == "counts"){
		curated.count.data.list <- curate.count.data(counts,sample.sizes,prefix)
		sample.frequencies <- curate.frequency.data(curated.count.data.list$sample.frequencies)
		mean.sample.sizes <- curated.count.data.list$mean.sample.sizes
		mean.centered.normalized.frequencies <- mean.center.normalize.frequencies(curated.count.data.list$sample.frequencies,curated.count.data.list$mean.sample.sizes,prefix)
		sample.covariance <- cov(t(mean.centered.normalized.frequencies))
		loci <- ncol(mean.centered.normalized.frequencies)
	}	
	projection.matrix <- get.projection.matrix(mean.sample.sizes,proj.mat.option)
	sample.covariance <- t(projection.matrix) %*% sample.covariance %*% projection.matrix
	spacemix.data <- list(	"sample.frequencies" = sample.frequencies,
							"mean.centered.normalized.frequencies" = mean.centered.normalized.frequencies,
							"mean.sample.sizes" = mean.sample.sizes,
							"projection.matrix" = projection.matrix,
							"sample.covariance" = sample.covariance,
							"loci" = loci)
		save(spacemix.data,file=paste(prefix,"spacemix.data.Robj",sep=''))
	return(spacemix.data)
}

propose.new.location <- function(lat,long,dist.std){
	coords_prime <- c(lat,long) + rnorm(n = 2, mean = 0, sd = dist.std)
	return(coords_prime)
}

initiate.admix.proportions <- function(k,model.option){
	if(model.option == "no_movement"){
		admix.proportions <- numeric(k)
		prior_prob_admix_proportions <- 0
	} else if(model.option == "target"){
		admix.proportions <- numeric(k)
		prior_prob_admix_proportions <- 0
	} else if(model.option == "source"){
		admix.proportions <- rbeta(k,shape1=0.1,shape2=1)
		prior_prob_admix_proportions <- Prior_prob_admix_proportions(admix.proportions)
	} else if(model.option == "source_and_target"){
		admix.proportions <- rbeta(k,shape1=0.1,shape2=1)
		prior_prob_admix_proportions <- Prior_prob_admix_proportions(admix.proportions)
	}
	initaite.admix.proportions.list <- list("admix.proportions" = admix.proportions, 
											"prior_prob_admix_proportions" = prior_prob_admix_proportions)
	return(initaite.admix.proportions.list)
}

initiate.population.coordinates <- function(observed.X.coordinates,observed.Y.coordinates,k){
	population.coordinates <- 	rbind(
		cbind(observed.X.coordinates,
			observed.Y.coordinates),
		cbind(	runif(k, 
			min = min(observed.X.coordinates), 
			max = max(observed.X.coordinates)),
		runif(k, 
			min = min(observed.Y.coordinates), 
			max = max(observed.Y.coordinates)))
	)
	return(population.coordinates)
}

save.initial.parameters <- function(a0,aD,a2,nugget,admix.proportions,covariance,admixed.covariance,transformed_covariance,population.coordinates,D,projection.matrix,prefix){
	initial.parameters <- list("a0" = a0,"aD" = aD,"a2" = a2,"nugget" = nugget,"admix.proportions" = admix.proportions,
								"covariance" = covariance,"admixed.covariance" = admixed.covariance,"transformed_covariance" = transformed_covariance,
								"population.coordinates" = population.coordinates,"D" = D,"projection.matrix" = projection.matrix)
	save(initial.parameters,file=paste(prefix,"Initial.parameters.Robj",sep=''))
	return(0)
}

initiate.update.function.list <- function(model.option){
	if(model.option == "no_movement"){
		Updates <- list(Update_a0,Update_aD,Update_a2,Update_nugget)
	} else if(model.option == "target"){
		Updates <- list(Update_a0,Update_aD,Update_a2,Update_nugget,Update_admixture_target_location)
	} else if(model.option == "source"){
		Updates <- list(Update_a0,Update_aD,Update_a2,Update_nugget,Update_admixture_source_location,Update_admixture_proportions)
	} else if(model.option == "source_and_target"){
		Updates <- list(Update_a0,Update_aD,Update_a2,Update_nugget,Update_admixture_target_location,Update_admixture_source_location,Update_admixture_proportions)
	}
	return(Updates)
}

print.mcmc.update <- function(LnL_freqs,prior_prob_admix_proportions,prior_prob_nugget,prior_prob_alpha0,prior_prob_alphaD,prior_prob_alpha2,prior_prob_admix_target_locations,prior_prob_admix_source_locations){
	P <- LnL_freqs + prior_prob_admix_proportions + prior_prob_nugget + prior_prob_alpha0 + prior_prob_alphaD + prior_prob_alpha2 + prior_prob_admix_target_locations + prior_prob_admix_source_locations
	return(P)
}

update.lstp <- function(n,lstp,acceptance.fraction){
	if(length(lstp) > 1){
		for(i in 1:length(lstp)){
			if(acceptance.fraction[i] > 0.44){
				lstp[i] <- lstp[i] + min(0.01,n^(-0.5))
			} else if(acceptance.fraction[i] < 0.44){
				lstp[i] <- lstp[i] - min(0.01,n^(-0.5))
			}
		}
	} else {
		if(acceptance.fraction > 0.44){
			lstp <- lstp + min(0.01,n^(-0.5))
		} else if(acceptance.fraction < 0.44){
			lstp <- lstp - min(0.01,n^(-0.5))
		}
	}
	return(lstp)
}

get.diagn.step <- function(generation,mixing.diagn.freq){
	diagn.step <- generation%%mixing.diagn.freq
	if(diagn.step == 0){
		diagn.step <- mixing.diagn.freq
	}
	return(diagn.step)
}

initiate.last.params <- function(spacemix.data,population.coordinates,admix.proportions,a0,aD,a2,nugget,covariance,admixed.covariance,transformed_covariance,
						admix.proportions.lstp, admix.target.location.lstp,admix.source.location.lstp,nugget.lstp,a0.lstp,aD.lstp,a2.lstp,k,LnL_freqs,
						prior_prob_alpha0,prior_prob_alphaD,prior_prob_alpha2,prior_prob_nugget,prior_prob_admix_proportions, prior_prob_admix_target_locations,prior_prob_admix_source_locations,
						a0_accept_rate,aD_accept_rate,a2_accept_rate,nugget_accept_rate,admix_source_location_accept_rate,admix_proportions_accept_rate,admix_target_location_accept_rate,
						a0_moves,aD_moves,a2_moves,nugget_moves,admix_source_location_moves,admix_target_location_moves,admix_proportions_moves,
						a0_accept,aD_accept,a2_accept,nugget_accept,admix_source_location_accept,admix_target_location_accept,admix_proportions_accept,
						D,observed.X.coordinates,observed.Y.coordinates,target.spatial.prior.scale,source.spatial.prior.scale,centroid,gibbs.spatial.fineness,gibbs.nugget.fineness){
	last.params <- list("sample.covariance" = spacemix.data$sample.covariance,
						"projection.matrix" = spacemix.data$projection.matrix,						
						"population.coordinates" = population.coordinates,
						"admix.proportions" = admix.proportions,
						"a0" = a0,"aD" = aD,"a2" = a2,"nugget" = nugget,
						"covariance" = covariance,
						"admixed.covariance" = admixed.covariance,
						"transformed_covariance" = transformed_covariance,
						"admix.proportions.lstp" = admix.proportions.lstp, 
						"admix.target.location.lstp" = admix.target.location.lstp,
						"admix.source.location.lstp" = admix.source.location.lstp,
						"nugget.lstp" = nugget.lstp,
						"a0.lstp" = a0.lstp,"aD.lstp" = aD.lstp,"a2.lstp" = a2.lstp,
						"k" = k,
						"LnL_freqs" = LnL_freqs,
						"prior_prob_alpha0" = prior_prob_alpha0,
						"prior_prob_alphaD" = prior_prob_alphaD,
						"prior_prob_alpha2" = prior_prob_alpha2,
						"prior_prob_nugget" = prior_prob_nugget,
						"prior_prob_admix_proportions" = prior_prob_admix_proportions, 
						"prior_prob_admix_target_locations" = prior_prob_admix_target_locations,
						"prior_prob_admix_source_locations" = prior_prob_admix_source_locations,
						"a0_accept_rate" = a0_accept_rate,
						"aD_accept_rate" = aD_accept_rate,
						"a2_accept_rate" = a2_accept_rate,
						"nugget_accept_rate" = nugget_accept_rate,
						"admix_source_location_accept_rate" = admix_source_location_accept_rate,
						"admix_proportions_accept_rate" = admix_proportions_accept_rate,
						"admix_target_location_accept_rate" = admix_target_location_accept_rate,
						"a0_moves" = a0_moves,
						"aD_moves" = aD_moves, 
						"a2_moves" = a2_moves,
						"nugget_moves" = nugget_moves,
						"admix_source_location_moves" = admix_source_location_moves, 
						"admix_target_location_moves" = admix_target_location_moves,
						"admix_proportions_moves" = admix_proportions_moves,
						"a0_accept" = a0_accept,
						"aD_accept" = aD_accept,
						"a2_accept" = a2_accept,
						"nugget_accept" = nugget_accept,
						"admix_source_location_accept" = admix_source_location_accept,
						"admix_target_location_accept" = admix_target_location_accept,
						"admix_proportions_accept" = admix_proportions_accept,
						"loci" = spacemix.data$loci,
						"D" = D,
						"mean.sample.sizes" = spacemix.data$mean.sample.sizes,
						"observed.X.coordinates" = observed.X.coordinates,
						"observed.Y.coordinates" = observed.Y.coordinates,
						"target.spatial.prior.scale" = target.spatial.prior.scale,
						"source.spatial.prior.scale" = source.spatial.prior.scale,
						"centroid" = centroid,
						"X.grid.fineness" = gibbs.spatial.fineness, 
						"Y.grid.fineness" = gibbs.spatial.fineness, 
						"nugget.grid.fineness" = gibbs.nugget.fineness)
	return(last.params)
}


}	
	
	
MCMC <-function(model.option,				#no_movement, target, source, source_and_target
				data.type,					#sample.covariance, sample.frequencies, counts
				proj.mat.option = NULL,
				sample.frequencies = NULL,
				mean.sample.sizes = NULL,
				counts = NULL,
				sample.sizes = NULL,
				sample.covariance = NULL,
				target.spatial.prior.scale = NULL,
				source.spatial.prior.scale = NULL,
				observed.X.coordinates,
				observed.Y.coordinates,
				k,
				loci,
				ngen,
				printfreq,
				samplefreq,
				mixing.diagn.freq,
				gibbs.nugget.fineness=50,
				gibbs.spatial.fineness=50,
				gibbs.step.frequency = 10000,
				savefreq,
				directory=NULL,
                prefix="",
				continue=FALSE,
				continuing.params=NULL){
	# recover()
	if(!is.null(directory)){
		setwd(directory)
	}
	
	spacemix.data <- spacemix.data(data.type = data.type,
									proj.mat.option = proj.mat.option,
									sample.frequencies = sample.frequencies, 
									loci = loci,
									mean.sample.sizes = mean.sample.sizes,
									counts = counts,
									sample.sizes = sample.sizes,
									sample.covariance = sample.covariance,
									prefix = prefix)
		
	#Declare variables
	if(TRUE){
		LnL_freqs <- numeric(ngen/samplefreq)
		Prob <- numeric(ngen/samplefreq)
		population.coordinates <- vector("list",ngen/samplefreq)
		distances <- vector("list",ngen/samplefreq)
		transformed.covariance.list <- vector("list",ngen/samplefreq)
		admix.proportions <- vector("list",ngen/samplefreq)
		nugget <- matrix(0,nrow=k,ncol=ngen/samplefreq)
		a0 <- numeric(ngen/samplefreq)
		aD <- numeric(ngen/samplefreq)
		a2 <- numeric(ngen/samplefreq)
		a0_accept_rate <- numeric(ngen/samplefreq)
		aD_accept_rate <- numeric(ngen/samplefreq)
		a2_accept_rate <- numeric(ngen/samplefreq)
		nugget_accept_rate <- matrix(0,nrow=k,ncol=(ngen/samplefreq))
		admix_target_location_accept_rate <- matrix(0,nrow=k,ncol=(ngen/samplefreq))
		admix_source_location_accept_rate <- matrix(0,nrow=k,ncol=(ngen/samplefreq))
		admix_proportions_accept_rate <- matrix(0,nrow=k,ncol=(ngen/samplefreq))
		a0_lstp <- numeric(ngen/samplefreq)
		aD_lstp <- numeric(ngen/samplefreq)
		a2_lstp <- numeric(ngen/samplefreq)
		nugget_lstp <- matrix(0,nrow=k,ncol=(ngen/samplefreq))
		admix_target_location_lstp <- matrix(0,nrow=k,ncol=(ngen/samplefreq))
		admix_source_location_lstp <- matrix(0,nrow=k,ncol=(ngen/samplefreq))
		admix_proportions_lstp <- matrix(0,nrow=k,ncol=(ngen/samplefreq))
		a0_diagn <- numeric(mixing.diagn.freq)
		aD_diagn <- numeric(mixing.diagn.freq)
		a2_diagn <- numeric(mixing.diagn.freq)
		nugget_diagn <- matrix(0,nrow=k,ncol=mixing.diagn.freq)
		admix_target_location_diagn <- matrix(0,nrow=k,ncol=mixing.diagn.freq)
		admix_source_location_diagn <- matrix(0,nrow=k,ncol=mixing.diagn.freq)
		admix_proportions_diagn <- matrix(0,nrow=k,ncol=mixing.diagn.freq)
		#
		seed <- sample(111:999,1)
			save(seed,file=paste(prefix,"_seed.Robj",sep=''))
		set.seed(seed)
	}
		
	if(!continue) {
		#INITIALIZE MCMC
				Prob[1] <- -Inf
				covariance <- matrix(0,nrow=k*2,ncol=k*2)
				badness.counter <- 0

			while(Prob[1] == -Inf | !matrixcalc::is.positive.definite(covariance) && badness.counter < 100){
						nugget[,1] <- rexp(k,rate = spacemix.data$mean.sample.sizes)
						a0[1] <- rexp(1,1/100)
						aD[1] <- rexp(1,1)
						a2[1] <- runif(1,0.1,2)
						population.coordinates[[1]] <- 	initiate.population.coordinates(observed.X.coordinates,observed.Y.coordinates,k)
						initiate.admix.proportions.list <- initiate.admix.proportions(k,model.option)
						admix.proportions[[1]] <- numeric(k)	#initiate.admix.proportions.list$admix.proportions
						prior_prob_admix_proportions <- 0		#initiate.admix.proportions.list$prior_prob_admix_proportions
					distances[[1]] <- fields::rdist(population.coordinates[[1]])
						centroid <- c(mean(observed.X.coordinates),mean(observed.Y.coordinates))
						if(is.null(target.spatial.prior.scale)){
							target.spatial.prior.scale <- mean(distances[[1]][1:k,1:k]) / 2
						}
						if(is.null(source.spatial.prior.scale)){
							source.spatial.prior.scale <- mean(distances[[1]][1:k,1:k]) * 2
						}
					covariance <- Covariance(a0[1],aD[1],a2[1],distances[[1]])
					admixed.covariance <- admixed.Covariance(covariance,admix.proportions[[1]],nugget[,1])
					transformed_covariance <- transformed.Covariance(admixed.covariance,spacemix.data$projection.matrix)
					tmp <- save.initial.parameters(a0[1],aD[1],a2[1],nugget[,1],admix.proportions[[1]],
													covariance,admixed.covariance,transformed_covariance,
													population.coordinates[[1]],distances[[1]],spacemix.data$projection.matrix,prefix)
				LnL_freqs[1] <- wishart.lnL(spacemix.data$sample.covariance,transformed_covariance/spacemix.data$loci,spacemix.data$loci)
					cat("LnL: ",LnL_freqs[1],"\n")
				prior_prob_alpha0 <- Prior_prob_alpha0(a0[1])
					cat("Pr(a0): ",prior_prob_alpha0,"\n")
				prior_prob_alphaD <- Prior_prob_alphaD(aD[1])
					cat("Pr(aD): ",prior_prob_alphaD,"\n")
				prior_prob_alpha2 <- Prior_prob_alpha2(a2[1])
					cat("Pr(a2): ",prior_prob_alpha2,"\n")
				prior_prob_nugget <- Prior_prob_nugget(nugget[,1],spacemix.data$mean.sample.sizes)
					cat("Pr(nugget): ",prior_prob_nugget,"\n")
				prior_prob_admix_target_locations <- Prior_prob_admix_target_locations(population.coordinates[[1]][1:k,],observed.X.coordinates,observed.Y.coordinates,target.spatial.prior.scale)
					cat("Pr(admix_target_locations): ",prior_prob_admix_target_locations,"\n")
				prior_prob_admix_source_locations <- Prior_prob_admix_source_locations(population.coordinates[[1]][(k+1):(2*k),],centroid,source.spatial.prior.scale)
					cat("Pr(admix_source_locations): ",prior_prob_admix_source_locations,"\n")
					cat("Pr(admix_proportions): ",prior_prob_admix_proportions,"\n")
				Prob[1] <- LnL_freqs[1] + prior_prob_admix_proportions + prior_prob_nugget + prior_prob_alpha0 + prior_prob_alphaD + prior_prob_alpha2 + prior_prob_admix_target_locations + prior_prob_admix_source_locations
				cat("Prob: ",Prob[1],"\n")
			badness.counter <- badness.counter + 1
			if(badness.counter > 99){
				if(!is.finite(Prob[1])){													
					stop("Initial probability of model is NEGATIVE INFINITY! Please attempt to initiate chain again.")
				} else {
					stop("the initial transformed covariance matrix is not positive definite! Attempt to re-initialize MCMC!")
				}
			}
		}	
		last.params <- initiate.last.params(spacemix.data = spacemix.data,population.coordinates = population.coordinates[[1]],admix.proportions = admix.proportions[[1]],
											a0[1],aD[1],a2[1],nugget[,1],covariance,admixed.covariance,transformed_covariance,
											admix.proportions.lstp = numeric(k),admix.target.location.lstp = numeric(k),admix.source.location.lstp = numeric(k),nugget.lstp = numeric(k),
											a0.lstp = 0,aD.lstp = 0,a2.lstp = 0,k = k,LnL_freqs = LnL_freqs[1],
											prior_prob_alpha0 = prior_prob_alpha0,prior_prob_alphaD = prior_prob_alphaD,prior_prob_alpha2 = prior_prob_alpha2,
											prior_prob_nugget = prior_prob_nugget,prior_prob_admix_proportions = prior_prob_admix_proportions,
											prior_prob_admix_target_locations = prior_prob_admix_target_locations,prior_prob_admix_source_locations = prior_prob_admix_source_locations,
											a0_accept_rate = 0,aD_accept_rate = 0,a2_accept_rate = 0,nugget_accept_rate = numeric(k),admix_source_location_accept_rate = numeric(k),
											admix_proportions_accept_rate = numeric(k),admix_target_location_accept_rate = numeric(k),
											a0_moves = 0,aD_moves = 0,a2_moves = 0,nugget_moves = numeric(k),admix_source_location_moves = numeric(k),admix_target_location_moves = numeric(k),
											admix_proportions_moves = numeric(k),a0_accept = 0,aD_accept = 0,a2_accept = 0,nugget_accept = numeric(k),admix_source_location_accept = numeric(k),
											admix_target_location_accept = numeric(k),admix_proportions_accept = numeric(k),D = distances[[1]],
											observed.X.coordinates = observed.X.coordinates,observed.Y.coordinates = observed.Y.coordinates,
											target.spatial.prior.scale = target.spatial.prior.scale,source.spatial.prior.scale = source.spatial.prior.scale,
											centroid = centroid,gibbs.spatial.fineness = gibbs.spatial.fineness,gibbs.nugget.fineness = gibbs.nugget.fineness)
	} else {
		load(continuing.params)
		a0_diagn <- continuing.params$a0_diagn
		aD_diagn <- continuing.params$aD_diagn
		a2_diagn <- continuing.params$a2_diagn
		nugget_diagn <- continuing.params$nugget_diagn
		admix_target_location_diagn <- continuing.params$admix_target_location_diagn
		admix_source_location_diagn <- continuing.params$admix_source_location_diagn
		admix_proportions_diagn <- continuing.params$admix_proportions_diagn
		nugget[,1] <- continuing.params$nugget
		a0[1] <- continuing.params$a0
		aD[1] <- continuing.params$aD
		a2[1] <- continuing.params$a2
		population.coordinates[[1]] <- 	continuing.params$population.coordinates
		admix.proportions[[1]] <- continuing.params$admix.proportions
		prior_prob_admix_proportions <- continuing.params$prior_prob_admix_proportions
					distances[[1]] <- continuing.params$D
						centroid <- continuing.params$centroid
						target.spatial.prior.scale <- continuing.params$target.spatial.prior.scale
						source.spatial.prior.scale <- continuing.params$source.spatial.prior.scale
					covariance <- continuing.params$covariance
					admixed.covariance <- continuing.params$admixed.covariance
					transformed_covariance <- continuing.params$transformed_covariance
					tmp <- save.initial.parameters(a0[1],aD[1],a2[1],nugget[,1],admix.proportions[[1]],
													covariance,admixed.covariance,transformed_covariance,
													population.coordinates[[1]],distances[[1]],spacemix.data$projection.matrix,prefix)
				LnL_freqs[1] <- continuing.params$LnL_freqs
					cat("LnL: ",LnL_freqs[1],"\n")
				prior_prob_alpha0 <- continuing.params$prior_prob_alpha0
					cat("Pr(a0): ",prior_prob_alpha0,"\n")
				prior_prob_alphaD <- continuing.params$prior_prob_alphaD
					cat("Pr(aD): ",prior_prob_alphaD,"\n")
				prior_prob_alpha2 <- continuing.params$prior_prob_alpha2
					cat("Pr(a2): ",prior_prob_alpha2,"\n")
				prior_prob_nugget <- continuing.params$prior_prob_nugget
					cat("Pr(nugget): ",prior_prob_nugget,"\n")
				prior_prob_admix_target_locations <- continuing.params$prior_prob_admix_target_locations
					cat("Pr(admix_target_locations): ",prior_prob_admix_target_locations,"\n")
				prior_prob_admix_source_locations <- continuing.params$prior_prob_admix_source_locations
					cat("Pr(admix_source_locations): ",prior_prob_admix_source_locations,"\n")
					cat("Pr(admix_proportions): ",prior_prob_admix_proportions,"\n")
				Prob[1] <- LnL_freqs[1] + prior_prob_admix_proportions + prior_prob_nugget + prior_prob_alpha0 + prior_prob_alphaD + prior_prob_alpha2 + prior_prob_admix_target_locations + prior_prob_admix_source_locations
				cat("Prob: ",Prob[1],"\n")
		last.params <- initiate.last.params(spacemix.data = spacemix.data,
											population.coordinates = population.coordinates[[1]],
											admix.proportions = admix.proportions[[1]],
											a0[1],aD[1],a2[1],nugget[,1],covariance,admixed.covariance,transformed_covariance,
											admix.proportions.lstp = continuing.params$admix.proportions.lstp,
											admix.target.location.lstp = continuing.params$admix.target.location.lstp,
											admix.source.location.lstp = continuing.params$admix.source.location.lstp,
											nugget.lstp = continuing.params$nugget.lstp,
											a0.lstp = continuing.params$a0.lstp,
											aD.lstp = continuing.params$aD.lstp,
											a2.lstp = continuing.params$a2.lstp,
											k = k,LnL_freqs = LnL_freqs[1],
											prior_prob_alpha0 = prior_prob_alpha0,
											prior_prob_alphaD = prior_prob_alphaD,
											prior_prob_alpha2 = prior_prob_alpha2,
											prior_prob_nugget = prior_prob_nugget,
											prior_prob_admix_proportions = prior_prob_admix_proportions,
											prior_prob_admix_target_locations = prior_prob_admix_target_locations,
											prior_prob_admix_source_locations = prior_prob_admix_source_locations,
											a0_accept_rate = continuing.params$a0_accept_rate,
											aD_accept_rate = continuing.params$aD_accept_rate,
											a2_accept_rate = continuing.params$a2_accept_rate,
											nugget_accept_rate = continuing.params$nugget_accept_rate,
											admix_source_location_accept_rate = continuing.params$admix_source_location_accept_rate,
											admix_proportions_accept_rate = continuing.params$admix_proportions_accept_rate,
											admix_target_location_accept_rate = continuing.params$admix_target_location_accept_rate,
											a0_moves = continuing.params$a0_moves,
											aD_moves = continuing.params$aD_moves,
											a2_moves = continuing.params$a2_moves,
											nugget_moves = continuing.params$nugget_moves,
											admix_source_location_moves = continuing.params$admix_source_location_moves,
											admix_target_location_moves = continuing.params$admix_target_location_moves,
											admix_proportions_moves = continuing.params$admix_proportions_moves,
											a0_accept = continuing.params$a0_accept,
											aD_accept = continuing.params$aD_accept,
											a2_accept = continuing.params$a2_accept,
											nugget_accept = continuing.params$nugget_accept,
											admix_source_location_accept = continuing.params$admix_source_location_accept,
											admix_target_location_accept = continuing.params$admix_target_location_accept,
											admix_proportions_accept = continuing.params$admix_proportions_accept,
											D = distances[[1]],
											observed.X.coordinates = continuing.params$observed.X.coordinates,
											observed.Y.coordinates = continuing.params$observed.Y.coordinates,
											target.spatial.prior.scale = continuing.params$target.spatial.prior.scale,
											source.spatial.prior.scale = continuing.params$source.spatial.prior.scale,
											centroid = continuing.params$centroid,
											gibbs.spatial.fineness = continuing.params$gibbs.spatial.fineness,
											gibbs.nugget.fineness = continuing.params$gibbs.nugget.fineness)
	}
	
	#Run the MCMC
		Updates <- initiate.update.function.list(model.option)
	
	for(i in 2:ngen) {
		x <- sample(c(1:length(Updates)),1)
		new.params <- Updates[[x]](last.params)
		if(i%%gibbs.step.frequency == 0){
			new.params <- admix_target_location_and_nugget_gibbs_sampler(last.params)
		}

		if(i%%samplefreq == 0){
			j <- i/samplefreq
			population.coordinates[[j]] <- new.params$population.coordinates
			distances[[j]] <- new.params$D
			transformed.covariance.list[[j]] <- new.params$transformed_covariance
			admix.proportions[[j]] <- new.params$admix.proportions
			nugget[,j] <- new.params$nugget
			a0[j] <- new.params$a0
			aD[j] <- new.params$aD
			a2[j] <- new.params$a2
			LnL_freqs[j] <- new.params$LnL_freqs
			Prob[j] <- LnL_freqs[j] +
						new.params$prior_prob_admix_proportions +
						new.params$prior_prob_nugget +
						new.params$prior_prob_alpha0 +
						new.params$prior_prob_alphaD +
						new.params$prior_prob_alpha2 + 
						new.params$prior_prob_admix_target_locations + 
						new.params$prior_prob_admix_source_locations
			a0_accept_rate[j] <- new.params$a0_accept_rate
			aD_accept_rate[j] <- new.params$aD_accept_rate
			a2_accept_rate[j] <- new.params$a2_accept_rate
			nugget_accept_rate[,j] <- new.params$nugget_accept_rate
			admix_target_location_accept_rate[,j] <- new.params$admix_target_location_accept_rate
			admix_source_location_accept_rate[,j] <- new.params$admix_source_location_accept_rate
			admix_proportions_accept_rate[,j] <- new.params$admix_proportions_accept_rate
			a0_lstp[j] <- new.params$a0.lstp
			aD_lstp[j] <- new.params$aD.lstp
			a2_lstp[j] <- new.params$a2.lstp
			nugget_lstp[,j] <- new.params$nugget.lstp
			admix_target_location_lstp[,j] <- new.params$admix.target.location.lstp
			admix_source_location_lstp[,j] <- new.params$admix.source.location.lstp
			admix_proportions_lstp[,j] <- new.params$admix.proportions.lstp
		}
		
		diagn.step <- get.diagn.step(i,mixing.diagn.freq)
			a0_diagn[diagn.step] <- new.params$a0_accept_rate
			aD_diagn[diagn.step] <- new.params$aD_accept_rate
			a2_diagn[diagn.step] <- new.params$a2_accept_rate
			nugget_diagn[,diagn.step] <- new.params$nugget_accept_rate
			admix_target_location_diagn[,diagn.step] <- new.params$admix_target_location_accept_rate
			admix_source_location_diagn[,diagn.step] <- new.params$admix_source_location_accept_rate
			admix_proportions_diagn[,diagn.step] <- new.params$admix_proportions_accept_rate
		
		if(i%%mixing.diagn.freq == 0){
			n <- i %/% mixing.diagn.freq
			new.params$a0.lstp <- update.lstp(n,new.params$a0.lstp,mean(a0_diagn))
			new.params$aD.lstp <- update.lstp(n,new.params$aD.lstp,mean(aD_diagn))
			new.params$a2.lstp <- update.lstp(n,new.params$a2.lstp,mean(a2_diagn))
			new.params$nugget.lstp <- update.lstp(n,new.params$nugget.lstp,rowMeans(nugget_diagn))
			new.params$admix.target.location.lstp <- update.lstp(n,new.params$admix.target.location.lstp,rowMeans(admix_target_location_diagn))
			new.params$admix.source.location.lstp <- update.lstp(n,new.params$admix.source.location.lstp,rowMeans(admix_source_location_diagn))
			new.params$admix.proportions.lstp <- update.lstp(n,new.params$admix.proportions.lstp,rowMeans(admix_proportions_diagn))
		}

		last.params <- new.params
						
		if(i%%printfreq == 0){
			P <- print.mcmc.update(new.params$LnL_freqs,new.params$prior_prob_admix_proportions,
									new.params$prior_prob_nugget,new.params$prior_prob_alpha0,new.params$prior_prob_alphaD,
									new.params$prior_prob_alpha2,new.params$prior_prob_admix_target_locations,new.params$prior_prob_admix_source_locations)
			cat(i," ---- ",P,"\n")
		}
				
		if(i%%savefreq == 0){	
			save(last.params,
				LnL_freqs,Prob,covariance,admixed.covariance,transformed_covariance,distances,
				population.coordinates,transformed.covariance.list,admix.proportions,a0,aD,a2,nugget,samplefreq,ngen,
				admix_source_location_accept_rate,admix_target_location_accept_rate,admix_proportions_accept_rate,a0_accept_rate,aD_accept_rate,a2_accept_rate,nugget_accept_rate,
				admix_source_location_lstp,admix_proportions_lstp,a0_lstp,aD_lstp,a2_lstp,nugget_lstp,
				admix_target_location_lstp,a0_diagn,aD_diagn,a2_diagn,nugget_diagn,
				admix_target_location_diagn,admix_source_location_diagn,admix_proportions_diagn,target.spatial.prior.scale,source.spatial.prior.scale,
				file=paste(prefix,sprintf("space_MCMC_output%d.Robj",1),sep=''))
		}
	}
    return(paste("Output",i,"runs to",paste(prefix,"MCMC_output*.Robj",sep=''),"."))
}

make.continuing.params <- function(MCMC.output,file.name){
	load(MCMC.output)
		last.params <- c(last.params,
							"a0_diagn" = a0_diagn,
							"a2_diagn" = a2_diagn,
							"aD_diagn" = aD_diagn,
							"admix_proportions_diagn" = admix_proportions_diagn,
							"admix_source_location_diagn" = admix_source_location_diagn,
							"admix_target_location_diagn" = admix_target_location_diagn,
							"nugget_diagn" = nugget_diagn)
	with(last.params, {
		continuing.params <- list(	"a0" = a0,"aD" = aD,"a2" = a2,"nugget" = nugget,"D" = D,
									"population.coordinates" = population.coordinates,
									"admix.proportions" = admix.proportions,
									"prior_prob_alpha0" = prior_prob_alpha0,
									"prior_prob_alphaD" = prior_prob_alphaD,
									"prior_prob_alpha2" = prior_prob_alpha2,
									"prior_prob_nugget" = prior_prob_nugget,
									"prior_prob_admix_target_locations" = prior_prob_admix_target_locations,
									"prior_prob_admix_source_locations" = prior_prob_admix_source_locations,
									"prior_prob_admix_proportions" = prior_prob_admix_proportions,
									"centroid" = centroid,"target.spatial.prior.scale" = target.spatial.prior.scale,
									"source.spatial.prior.scale" = source.spatial.prior.scale,
									"covariance" = covariance,"admixed.covariance" = admixed.covariance,
									"transformed_covariance" = transformed_covariance,"LnL_freqs" = LnL_freqs,
									"a0_accept" = a0_accept,"aD_accept" = aD_accept,
									"a2_accept" = a2_accept,"nugget_accept" = nugget_accept,
									"admix_target_location_accept" = admix_target_location_accept,
									"admix_source_location_accept" = admix_source_location_accept,
									"admix_proportions_accept" = admix_proportions_accept,
									"a0_accept_rate" = a0_accept_rate,"aD_accept_rate" = aD_accept_rate,
									"a2_accept_rate" = a2_accept_rate,"nugget_accept_rate" = nugget_accept_rate,
									"admix_target_location_accept_rate" = admix_target_location_accept_rate,
									"admix_source_location_accept_rate" = admix_source_location_accept_rate,
									"admix_proportions_accept_rate" = admix_proportions_accept_rate,
									"a0_moves" = a0_moves,"aD_moves" = aD_moves,"a2_moves" = a2_moves,
									"nugget_moves" = nugget_moves,"admix_source_location_moves" = admix_source_location_moves, 
									"admix_target_location_moves" = admix_target_location_moves,"admix_proportions_moves" = admix_proportions_moves,
									"observed.X.coordinates" = observed.X.coordinates,"observed.Y.coordinates" = observed.Y.coordinates, 
									"X.grid.fineness" = X.grid.fineness,"Y.grid.fineness" = Y.grid.fineness,"nugget.grid.fineness" = nugget.grid.fineness,
									"a0.lstp" = a0.lstp,"aD.lstp" = aD.lstp,"a2.lstp" = a2.lstp,
									"nugget.lstp" = nugget.lstp,"admix.target.location.lstp" = admix.target.location.lstp,
									"admix.source.location.lstp" = admix.source.location.lstp,"admix.proportions.lstp" = admix.proportions.lstp,
									"a0_diagn" = a0_diagn,"a2_diagn" = a2_diagn,"aD_diagn" = aD_diagn,"admix_proportions_diagn" = admix_proportions_diagn,
									"admix_source_location_diagn" = admix_source_location_diagn,"admix_target_location_diagn" = admix_target_location_diagn,
									"nugget_diagn" = nugget_diagn)
	save(continuing.params,file=file.name)
	})
	return(0)
}

get.conditional.mean <- function(covariance,pop.to.sample,observations,drop.option=NULL){
 # recover()
	if(!is.null(drop.option)){
		observations <- observations[-pop.to.sample,]
	}
	omega_12 <- covariance[pop.to.sample,-pop.to.sample,drop=FALSE]
	omega_22 <- covariance[-pop.to.sample,-pop.to.sample,drop=FALSE]
	mu_bar <- omega_12%*%MASS::ginv(omega_22) %*% (observations) #0.5 + omega_12%*%MASS::ginv(omega_22) %*% (observations - 0.5)
	return(mu_bar)
}

get.spatial.fstat <- function(sampled.locations,proposed.location,observations,focal.population,covariance=NULL,a0,aD,a2,focal.population.mean=NULL,get.focal.population.mean=NULL){
	# recover()
	k <- nrow(sampled.locations)
	if(is.null(covariance)){
		covariance <- Covariance(a0,aD,a2,fields::rdist(sampled.locations))
	}
		covariance.no.focal <- matrix(0,nrow=k,ncol=k)
		covariance.no.focal[1:(k-1),1:(k-1)] <- covariance[-focal.population,-focal.population]
	if(is.null(focal.population.mean) && is.null(get.focal.population.mean)){
		stop("you must either specify the focal population conditional mean or you must estimate it")
	}
	if(!is.null(get.focal.population.mean)){
		focal.population.mean <- get.conditional.mean(covariance,focal.population,observations,drop.option=TRUE)
	}
	coords <- rbind(sampled.locations[-focal.population,],proposed.location)
	tmp.D <- fields::rdist(coords[k,1:2,drop=FALSE],coords)
	tmp.cov <- Covariance(a0,aD,a2,tmp.D)
	covariance.no.focal[k,] <- tmp.cov
	covariance.no.focal[,k] <- tmp.cov
	spatial.location.mean <- get.conditional.mean(covariance.no.focal,k,observations[-focal.population,])
	spatial.fstat <- mean((observations[focal.population,]-focal.population.mean)*spatial.location.mean)
	return(spatial.fstat)
}

procrusteez <- function(obs.locs,target.locs,k,source.locs = NULL,option){
	# recover()
	proc.loc <- procrustes(obs.locs,target.locs,scale=TRUE)
	if(option==1){
		proc.pop.loc <- proc.loc$scale * target.locs %*% proc.loc$rotation + matrix(proc.loc$translation,nrow=k,ncol=2,byrow=TRUE)
	} else if(option==2){
		proc.pop.loc <- proc.loc$scale * source.locs %*% proc.loc$rotation + matrix(proc.loc$translation,nrow=k,ncol=2,byrow=TRUE)
	}
	return(proc.pop.loc)	
}

plot.pop.coords <- function(observed.locations,target.coords,k,admix.proportions,source.coords=NULL,plot.range=NULL,segments=NULL){
	# recover()
	if(is.null(plot.range)){
		plot.range <- get.plot.range(target.coords,source.coords)
	}
	plot(observed.locations,col=rainbow(k),xlim=c(plot.range[[1]],plot.range[[2]]),ylim=c(plot.range[[3]],plot.range[[4]]),pch=13,cex=2)
		if(!is.null(segments)){
			if(is.null(source.coords)){
				stop("must specify admixture source coordinates")
			}
			for(i in 1:length(target.coords)){
				points(target.coords[[i]][,1],target.coords[[i]][,2],col=adjustcolor(rainbow(k),0.2),pch=19)
					points(source.coords[[i]][,1],source.coords[[i]][,2],col=adjustcolor(rainbow(k),0.2),pch=15)
					segments(x0 = source.coords[[i]][,1], y0 = source.coords[[i]][,2],
							x1 = target.coords[[i]][,1], y1 = target.coords[[i]][,2],
							col = adjustcolor(rainbow(k),0.2) , lwd = admix.proportions[[i]]*1)
			}
		} else if(is.null(segments)){
			for(i in 1:length(target.coords)){
				points(target.coords[[i]][,1],target.coords[[i]][,2],col=adjustcolor(rainbow(k),0.2),pch=19)
					points(source.coords[[i]][,1],source.coords[[i]][,2],col=adjustcolor(rainbow(k),0.2),pch=15)
			}
		}
	return(0)
}

visualize.admix.posterior <- function(MCMC.output.object,burnin,thinning,procrustes,source=NULL,plot.range=NULL,segments=NULL){
	# recover()
	load(MCMC.output.object)
	x <- seq(from = burnin, to = length(which(Prob!=0)), by = thinning)
		population.coordinates <- population.coordinates[x]
		admix.proportions <- admix.proportions[x]
		observed.locations <- cbind(last.params$observed.X.coordinates,last.params$observed.Y.coordinates)
	if(procrustes){
		proc.target.pop.coords <- lapply(1:length(population.coordinates),function(i){procrusteez(observed.locations,head(population.coordinates[[i]],n=last.params$k),last.params$k,option=1)})
		if(!is.null(source)){
			proc.source.pop.coords <- lapply(1:length(population.coordinates),function(i){procrusteez(observed.locations,head(population.coordinates[[i]],n=last.params$k),last.params$k,tail(population.coordinates[[i]],n=last.params$k),option=2)})
		} else if(is.null(source)){
			proc.source.pop.coords <- NULL
		}
	} else if(!procrustes){
		stop("uh oh")
	}
	plot.pop.coords(observed.locations = observed.locations,
					target.coords = proc.target.pop.coords,
					k = last.params$k,
					admix.proportions=admix.proportions,
					source.coords = proc.source.pop.coords,
					plot.range=plot.range,
					segments = segments)
	return(0)
}
	
if(FALSE){
#Prepping Globetrotter Data (Hellenthal et al)
a<-read.table("hellpopfreqs.txt.frq.strat",head=TRUE,nrow=150,as.is=TRUE)
populations<-unique(a$CLST)
npops<-length(populations)


num.snps<-130084 ##found by wc -l on file

read.this.many<-5e3
all.MAC<-matrix(NA,nrow=num.snps,ncol=npops)
all.sample.size<-matrix(NA,nrow=num.snps,ncol=npops)

snp.info<-character()
for(i in 1:floor((num.snps)/read.this.many)){
print(i)
my.snps<-read.table("hellpopfreqs.txt.frq.strat",skip=1+npops*(read.this.many*(i-1)), nrow=npops*read.this.many,as.is=TRUE)
colnames(my.snps)<-colnames(a)

MAC<-matrix(my.snps$MAC,ncol=npops,byrow=TRUE)
sample.size<-matrix(my.snps$NCHROBS,ncol=npops,byrow=TRUE)

these<-((read.this.many*(i-1)+1)):((read.this.many*i))
all.MAC[these,]<-MAC
all.sample.size[these,]<-sample.size
snp.info<-rbind(snp.info,my.snps[,c("CHR","SNP")])
}
i=27
my.snps<-read.table("hellpopfreqs.txt.frq.strat",skip=1+npops*(read.this.many*(i-1)),as.is=TRUE)
colnames(my.snps)<-colnames(a)
MAC<-matrix(my.snps$MAC,ncol=npops,byrow=TRUE)
sample.size<-matrix(my.snps$NCHROBS,ncol=npops,byrow=TRUE)
these<-((read.this.many*(i-1)+1)):((read.this.many*(i-1))+nrow(MAC))
all.MAC[these,]<-MAC
all.sample.size[these,]<-sample.size
colnames(all.MAC)<-populations
colnames(all.sample.size)<-populations
snp.info<-rbind(snp.info,my.snps[,c("CHR","SNP")])

save(file="hellpopfreqs.Robj",all.MAC,all.sample.size,snp.info)

		# a0_accept_rate[1] <- continuing.params$a0_accept_rate
		# aD_accept_rate[1] <- continuing.params$aD_accept_rate
		# a2_accept_rate[1] <- continuing.params$a2_accept_rate
		# nugget_accept_rate[,1] <- continuing.params$nugget_accept_rate
		# admix_target_location_accept_rate[,1] <- continuing.params$admix_target_location_accept_rate
		# admix_source_location_accept_rate[,1] <- continuing.params$admix_source_location_accept_rate
		# admix_proportions_accept_rate[,1] <- continuing.params$admix_proportions_accept_rate
		# a0_lstp[1] <- continuing.params$a0.lstp
		# aD_lstp[1] <- continuing.params$aD.lstp
		# a2_lstp[1] <- continuing.params$a2.lstp
		# nugget_lstp[,1] <- continuing.params$nugget.lstp
		# admix_target_location_lstp[,1] <- continuing.params$admix.target.location.lstp
		# admix_source_location_lstp[,1] <- continuing.params$admix.source.location.lstp
		# admix_proportions_lstp[,1] <- continuing.params$admix.proportions.lstp

}
	