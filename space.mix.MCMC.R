#UPDATE FUNCTIONS ETC.
if(TRUE){
simulate.spacemix.dataset <- function(k,loci,admix.target,admix.source,admix.proportion,sim.a0,sim.aD,sim.a2,sample.sizes,generate.counts,boundary.option,filename){
	#recover()
	#set random seed
		random.seed <- sample(100:999,1)
			set.seed(random.seed)
	#simulate population locations
		sim.locations <- cbind(runif(2*k,-1,1),runif(2*k,-1,1))
		sim.locations[k+admix.target,] <- sim.locations[admix.source,]
		spatial.prior.X.coordinates <- sim.locations[(1:k),1]
		spatial.prior.Y.coordinates <- sim.locations[(1:k),2]
	#generate spatial covariance matrix
		sim.D <- fields::rdist(sim.locations)
		mean.sample.sizes <- rowMeans(sample.sizes)
		sim.covariance <- Covariance(sim.a0,sim.aD,sim.a2,sim.D)
	#generate admixed covariance matrix
		sim.nugget <- numeric(k)
		sim.admix.proportions <- numeric(k)
			sim.admix.proportions[admix.target] <- admix.proportion
		sim.admixed.covariance <- admixed.Covariance(sim.covariance,sim.admix.proportions,sim.nugget,k,1/mean.sample.sizes)
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
				prior.prob.admix.target.locations.array[x,y,] <- Prior_prob_admix_target_locations(coords.prime[1:last.params$k,],last.params$spatial.prior.X.coordinates,last.params$spatial.prior.Y.coordinates,last.params$target.spatial.prior.scale)
				tmp.cov <- Covariance(last.params$a0,
									  last.params$aD,
								  	  last.params$a2,
									  spacemix.dist(coords.prime[pop.to.update,1:2,drop=FALSE],coords.prime))
				covariance.prime[pop.to.update,] <- tmp.cov
				covariance.prime[,pop.to.update] <- tmp.cov
			for(z in 1:length(nugget.gridpoints)){
				nugget.prime[pop.to.update] <- nugget.gridpoints[z]
				admixed.covariance.prime <- admixed.Covariance(covariance.prime,last.params$admix.proportions,nugget.prime,last.params$k,last.params$inv.mean.sample.sizes,last.params$identity.matrix)
				transformed.covariance.prime <- transformed.Covariance(admixed.covariance.prime,last.params$projection.matrix)
				lnL.array[x,y,z] <- likelihood(likelihood.option = last.params$likelihood.option,sample.cov = last.params$sample.cov,par.cov = transformed.covariance.prime,index.matrix = last.params$index.matrix,sd = last.params$sd,n = last.params$loci)
				prior.prob.nugget.array[x,y,z] <- Prior_prob_nugget(nugget.prime)
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
		new.params$D <- spacemix.dist(new.params$population.coordinates)
		new.params$covariance <- Covariance(last.params$a0,last.params$aD,last.params$a2,new.params$D)
		new.params$admixed.covariance <- admixed.Covariance(new.params$covariance,last.params$admix.proportions,new.params$nugget,last.params$k,last.params$inv.mean.sample.sizes,last.params$identity.matrix)
		new.params$transformed.covariance <- transformed.Covariance(new.params$admixed.covariance,last.params$projection.matrix)
		new.params$LnL_freqs <- lnL.array[sampled.parameters.index]
		new.params$moves$admix_target_location_moves[pop.to.update] <- new.params$moves$admix_target_location_moves[pop.to.update] + 1
		new.params$accpet$admix_target_location_accept[pop.to.update] <- new.params$accept$admix_target_location_accept[pop.to.update] + 1
		new.params$accept_rates$admix_target_location_accept_rate[pop.to.update] <- new.params$accept$admix_target_location_accept[pop.to.update]/new.params$moves$admix_target_location_moves[pop.to.update]
		new.params$moves$nugget_moves[pop.to.update] <- new.params$moves$nugget_moves[pop.to.update] + 1
		new.params$accept$nugget_accept[pop.to.update] <- new.params$accept$nugget_accept[pop.to.update] + 1
		new.params$accept_rates$nugget_accept_rate[pop.to.update] <- new.params$accept$nugget_accept[pop.to.update]/new.params$moves$nugget_moves[pop.to.update]
	}
		return(new.params)
}

update.matrix <- function(matrix,pop.to.update,update.vector){
	matrix[pop.to.update,] <- update.vector
	matrix[,pop.to.update] <- update.vector
	return(matrix)
}

Update_admixture_target_location <- function(last.params){
	# recover()
	new.params <- last.params
	pop.to.update <- sample(1:last.params$k,1)
		population.coordinates_prime <- last.params$population.coordinates
		population.coordinates_prime[pop.to.update,] <- propose.new.location(population.coordinates_prime[pop.to.update,1],
																			population.coordinates_prime[pop.to.update,2],
																			exp(last.params$lstps$admix_target_location_lstp[pop.to.update]))
		prior_prob_admix_target_locations_prime <- Prior_prob_admix_target_locations(population.coordinates_prime[1:last.params$k,],
																						last.params$spatial.prior.X.coordinates,
																						last.params$spatial.prior.Y.coordinates,
																						last.params$target.spatial.prior.scale)
		D_prime <- spacemix.dist(population.coordinates_prime[pop.to.update,1:2,drop=FALSE], population.coordinates_prime)
		tmp.covariance_prime <- Covariance(last.params$a0,last.params$aD,last.params$a2,D_prime)
		covariance_prime <- update.matrix(last.params$covariance,pop.to.update,tmp.covariance_prime)
		admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions,last.params$nugget,last.params$k,last.params$inv.mean.sample.sizes,last.params$identity.matrix)
				transformed.covariance.prime <- transformed.Covariance(admixed.covariance_prime,
																		last.params$projection.matrix)
				LnL_freqs_prime <- likelihood(likelihood.option = last.params$likelihood.option,sample.cov = last.params$sample.cov,par.cov = transformed.covariance.prime,index.matrix = last.params$index.matrix,sd = last.params$sd,n = last.params$loci)
					if(metropolis_ratio(LnL_freqs_prime,prior_prob_admix_target_locations_prime,last.params$LnL_freqs,last.params$prior_prob_admix_target_locations)){
						new.params$population.coordinates <- population.coordinates_prime
						new.params$prior_prob_admix_target_locations <- prior_prob_admix_target_locations_prime
						new.params$D <- update.matrix(new.params$D,pop.to.update,D_prime)
						new.params$covariance <- covariance_prime
						new.params$admixed.covariance <- admixed.covariance_prime
						new.params$transformed.covariance <- transformed.covariance.prime
						new.params$LnL_freqs <- LnL_freqs_prime
						new.params$accept$admix_target_location_accept[pop.to.update] <- new.params$accept$admix_target_location_accept[pop.to.update] + 1
					}
	new.params$moves$admix_target_location_moves[pop.to.update] <- new.params$moves$admix_target_location_moves[pop.to.update] + 1
	new.params$accept_rates$admix_target_location_accept_rate[pop.to.update] <- new.params$accept$admix_target_location_accept[pop.to.update]/new.params$moves$admix_target_location_moves[pop.to.update]
	return(new.params)	
}

Update_admixture_source_location <- function(last.params){
	 # recover()
	new.params <- last.params
	pop.to.update <- sample(1:last.params$k,1)
	ghost.to.update <- pop.to.update + last.params$k
		population.coordinates_prime <- last.params$population.coordinates
		population.coordinates_prime[ghost.to.update,] <- propose.new.location(population.coordinates_prime[ghost.to.update,1],
																							population.coordinates_prime[ghost.to.update,2],
																							exp(last.params$lstps$admix_source_location_lstp[pop.to.update]))
	prior_prob_admix_source_locations_prime <- Prior_prob_admix_source_locations(population.coordinates_prime[(last.params$k+1):(2*last.params$k),],
																					last.params$centroid,
																					last.params$source.spatial.prior.scale)
		D_prime <- spacemix.dist(population.coordinates_prime[ghost.to.update,1:2,drop=FALSE], population.coordinates_prime)
		tmp.covariance_prime <- Covariance(last.params$a0,last.params$aD,last.params$a2,D_prime)
		covariance_prime <- update.matrix(last.params$covariance,ghost.to.update,tmp.covariance_prime)
		admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions,last.params$nugget,last.params$k,last.params$inv.mean.sample.sizes,last.params$identity.matrix)
			transformed.covariance.prime <- transformed.Covariance(admixed.covariance_prime,
																	last.params$projection.matrix)
				LnL_freqs_prime <- likelihood(likelihood.option = last.params$likelihood.option,sample.cov = last.params$sample.cov,par.cov = transformed.covariance.prime,index.matrix = last.params$index.matrix,sd = last.params$sd,n = last.params$loci)
				if(metropolis_ratio(LnL_freqs_prime,prior_prob_admix_source_locations_prime,last.params$LnL_freqs,last.params$prior_prob_admix_source_locations)){
					new.params$population.coordinates <- population.coordinates_prime
					new.params$prior_prob_admix_source_locations <- prior_prob_admix_source_locations_prime
					new.params$D <- update.matrix(new.params$D,ghost.to.update,D_prime)
					new.params$covariance <- covariance_prime
					new.params$admixed.covariance <- admixed.covariance_prime
					new.params$transformed.covariance <- transformed.covariance.prime
					new.params$LnL_freqs <- LnL_freqs_prime
					new.params$accept$admix_source_location_accept[pop.to.update] <- new.params$accept$admix_source_location_accept[pop.to.update] + 1
				}
	new.params$moves$admix_source_location_moves[pop.to.update] <- new.params$moves$admix_source_location_moves[pop.to.update] + 1
	new.params$accept_rates$admix_source_location_accept_rate[pop.to.update] <- new.params$accept$admix_source_location_accept[pop.to.update]/new.params$moves$admix_source_location_moves[pop.to.update]
	return(new.params)	
}

Update_admixture_proportions <- function(last.params){
	new.params <- last.params
	pop.to.update <- sample(last.params$k,1)
	admix.proportions_prime <- last.params$admix.proportions
	admix.proportions_prime[pop.to.update] <- admix.proportions_prime[pop.to.update] + rnorm(1,0,exp(last.params$lstps$admix_proportions_lstp[pop.to.update]))
	prior_prob_admix_proportions_prime <- Prior_prob_admix_proportions(admix.proportions_prime)
	if(is.finite(prior_prob_admix_proportions_prime)){
		admixed.covariance_prime <- admixed.Covariance(last.params$covariance,admix.proportions_prime,last.params$nugget,last.params$k,last.params$inv.mean.sample.sizes,last.params$identity.matrix)
				transformed.covariance.prime <- transformed.Covariance(admixed.covariance_prime,
																		last.params$projection.matrix)
				LnL_freqs_prime <- likelihood(likelihood.option = last.params$likelihood.option,sample.cov = last.params$sample.cov,par.cov = transformed.covariance.prime,index.matrix = last.params$index.matrix,sd = last.params$sd,n = last.params$loci)
					if(metropolis_ratio(LnL_freqs_prime,prior_prob_admix_proportions_prime,last.params$LnL_freqs,last.params$prior_prob_admix_proportions)){
						new.params$admix.proportions <- admix.proportions_prime
						new.params$prior_prob_admix_proportions <- prior_prob_admix_proportions_prime
						new.params$admixed.covariance <- admixed.covariance_prime
						new.params$transformed.covariance <- transformed.covariance.prime
						new.params$LnL_freqs <- LnL_freqs_prime
						new.params$accept$admix_proportions_accept[pop.to.update] <- new.params$accept$admix_proportions_accept[pop.to.update] + 1
					}
			}
	new.params$moves$admix_proportions_moves[pop.to.update] <- new.params$moves$admix_proportions_moves[pop.to.update] + 1
	new.params$accept_rates$admix_proportions_accept_rate[pop.to.update] <- new.params$accept$admix_proportions_accept[pop.to.update]/new.params$moves$admix_proportions_moves[pop.to.update]
	return(new.params)
}

Update_a0 <- function(last.params){
	new.params <- last.params
	a0_prime <- last.params$a0 + rnorm(1,0,exp(last.params$lstps$a0_lstp))
	prior_prob_alpha0_prime <- Prior_prob_alpha0(a0_prime)
	if(prior_prob_alpha0_prime != -Inf){
		covariance_prime <- (last.params$covariance*last.params$a0)/a0_prime
			admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions,last.params$nugget,last.params$k,last.params$inv.mean.sample.sizes,last.params$identity.matrix)
			transformed.covariance.prime <- transformed.Covariance(admixed.covariance_prime,last.params$projection.matrix)
				LnL_freqs_prime <- likelihood(likelihood.option = last.params$likelihood.option,sample.cov = last.params$sample.cov,par.cov = transformed.covariance.prime,index.matrix = last.params$index.matrix,sd = last.params$sd,n = last.params$loci)
				if(metropolis_ratio(LnL_freqs_prime,prior_prob_alpha0_prime,last.params$LnL_freqs,last.params$prior_prob_alpha0)){
					new.params$a0 <- a0_prime
					new.params$covariance <- covariance_prime
					new.params$admixed.covariance <- admixed.covariance_prime
					new.params$transformed.covariance <- transformed.covariance.prime					
					new.params$prior_prob_alpha0 <- prior_prob_alpha0_prime
					new.params$LnL_freqs <- LnL_freqs_prime
					new.params$accept$a0_accept <- new.params$accept$a0_accept + 1
				}
	}
	new.params$moves$a0_moves <- new.params$moves$a0_moves + 1
	new.params$accept_rates$a0_accept_rate <- new.params$accept$a0_accept/new.params$moves$a0_moves
	return(new.params)
}

Update_aD <- function(last.params){
	new.params <- last.params
	aD_prime <- last.params$aD + rnorm(1,0,exp(last.params$lstps$aD_lstp))
	prior_prob_alphaD_prime <- Prior_prob_alphaD(aD_prime)
	if(prior_prob_alphaD_prime != -Inf){
		covariance_prime <- Covariance(last.params$a0,aD_prime,last.params$a2,last.params$D)
			admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions,last.params$nugget,last.params$k,last.params$inv.mean.sample.sizes,last.params$identity.matrix)
			transformed.covariance.prime <- transformed.Covariance(admixed.covariance_prime,last.params$projection.matrix)
				LnL_freqs_prime <- likelihood(likelihood.option = last.params$likelihood.option,sample.cov = last.params$sample.cov,par.cov = transformed.covariance.prime,index.matrix = last.params$index.matrix,sd = last.params$sd,n = last.params$loci)
				if(metropolis_ratio(LnL_freqs_prime,prior_prob_alphaD_prime,last.params$LnL_freqs,last.params$prior_prob_alphaD)){
					new.params$aD <- aD_prime
					new.params$covariance <- covariance_prime
					new.params$admixed.covariance <- admixed.covariance_prime
					new.params$transformed.covariance <- transformed.covariance.prime					
					new.params$prior_prob_alphaD <- prior_prob_alphaD_prime
					new.params$LnL_freqs <- LnL_freqs_prime
					new.params$accept$aD_accept <- new.params$accept$aD_accept + 1
				}
	}
	new.params$moves$aD_moves <- new.params$moves$aD_moves + 1
	new.params$accept_rates$aD_accept_rate <- new.params$accept$aD_accept/new.params$moves$aD_moves
	return(new.params)
}

Update_a2 <- function(last.params){
		new.params <- last.params
		a2_prime <- last.params$a2 + rnorm(1,0,exp(last.params$lstps$a2_lstp))
		prior_prob_alpha2_prime <- Prior_prob_alpha2(a2_prime) 
		if(prior_prob_alpha2_prime != -Inf) {
			covariance_prime <- Covariance(last.params$a0,last.params$aD,a2_prime,last.params$D)
				admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions,last.params$nugget,last.params$k,last.params$inv.mean.sample.sizes,last.params$identity.matrix)
				transformed.covariance.prime <- transformed.Covariance(admixed.covariance_prime,last.params$projection.matrix)
				LnL_freqs_prime <- likelihood(likelihood.option = last.params$likelihood.option,sample.cov = last.params$sample.cov,par.cov = transformed.covariance.prime,index.matrix = last.params$index.matrix,sd = last.params$sd,n = last.params$loci)
					if(metropolis_ratio(LnL_freqs_prime,prior_prob_alpha2_prime,last.params$LnL_freqs,last.params$prior_prob_alpha2)){
						new.params$a2 <- a2_prime
						new.params$covariance <- covariance_prime
						new.params$admixed.covariance <- admixed.covariance_prime						
						new.params$transformed.covariance <- transformed.covariance.prime							
						new.params$prior_prob_alpha2 <- prior_prob_alpha2_prime
						new.params$LnL_freqs <- LnL_freqs_prime
						new.params$accept$a2_accept <- new.params$accept$a2_accept + 1
					}
		}
	new.params$moves$a2_moves <- new.params$moves$a2_moves + 1
	new.params$accept_rates$a2_accept_rate <- new.params$accept$a2_accept/new.params$moves$a2_moves
	return(new.params)
}

Update_nugget <- function(last.params){
	new.params <- last.params
		pop.to.update <- sample(1:last.params$k,1)
	nugget_prime <- last.params$nugget + c(rep(0,pop.to.update-1),rnorm(1,0,exp(last.params$lstps$nugget_lstp[pop.to.update])),rep(0,last.params$k-pop.to.update))
	prior_prob_nugget_prime <- Prior_prob_nugget(nugget_prime)
	if(prior_prob_nugget_prime != -Inf){
		admixed.covariance_prime <- last.params$admixed.covariance
		diag(admixed.covariance_prime) <- diag(admixed.covariance_prime) - last.params$nugget + nugget_prime
		transformed.covariance.prime <- transformed.Covariance(admixed.covariance_prime,last.params$projection.matrix)
				LnL_freqs_prime <- likelihood(likelihood.option = last.params$likelihood.option,sample.cov = last.params$sample.cov,par.cov = transformed.covariance.prime,index.matrix = last.params$index.matrix,sd = last.params$sd,n = last.params$loci)
				if(metropolis_ratio(LnL_freqs_prime, prior_prob_nugget_prime,last.params$LnL_freqs,last.params$prior_prob_nugget)){
					new.params$nugget <- nugget_prime
					new.params$admixed.covariance <- admixed.covariance_prime
					new.params$transformed.covariance <- transformed.covariance.prime
					new.params$prior_prob_nugget <- prior_prob_nugget_prime
					new.params$LnL_freqs <- LnL_freqs_prime
					new.params$accept$nugget_accept[pop.to.update] <- new.params$accept$nugget_accept[pop.to.update] + 1
				}
	}
	new.params$moves$nugget_moves[pop.to.update] <- new.params$moves$nugget_moves[pop.to.update] + 1
	new.params$accept_rates$nugget_accept_rate[pop.to.update] <- new.params$accept$nugget_accept[pop.to.update]/new.params$moves$nugget_moves[pop.to.update]
	return(new.params)
}

metropolis_ratio <- function(LnL_prime,prior_prime,LnL,prior){
	accept <- FALSE
	if( exp((LnL_prime + prior_prime) - (LnL+prior)) >= runif(1) ){
		accept <- TRUE
	}
	return(accept)
}

Prior_prob_admix_target_locations <- function(admix_target_locations,spatial.prior.X.coordinates,spatial.prior.Y.coordinates,target.spatial.prior.scale){
	prior_prob_admix_target_locations <- sum(-log(2*pi*target.spatial.prior.scale) + (-1/2) * 
			((admix_target_locations[,1] - spatial.prior.X.coordinates)^2 / target.spatial.prior.scale + 
				(admix_target_locations[,2] - spatial.prior.Y.coordinates)^2 / target.spatial.prior.scale))
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

Prior_prob_nugget <- function(nugget){
	sum(dexp(nugget,log=TRUE))
}

Prior_prob_admix_proportions <- function(admix_proportions){
	sum(dbeta(admix_proportions,shape1=1,shape2=100,log=TRUE))
}

normal.lnL <- function(sample.cov,par.cov,sd,index.matrix){
	par.cov <- par.cov[index.matrix]
	lnL <- - (sample.cov-par.cov)^2/(2*sd^2)
	return(sum(lnL))
}

wishart.lnL <- function(sample.cov,par.cov,n){
	par.cov <- par.cov/n
	A <- solve(par.cov)
	lnL <- -0.5 * sum( A * sample.cov ) - (n/2)*determinant(par.cov,logarithm=TRUE)$modulus
	return(lnL)
}

likelihood <- function(likelihood.option,sample.cov,par.cov,index.matrix=NULL,sd=NULL,n=NULL){
	if(likelihood.option == "wishart"){
		lnL <- wishart.lnL(sample.cov,par.cov,n)
	} else if(likelihood.option == "normal_approx"){
		lnL <- normal.lnL(sample.cov,par.cov,sd,index.matrix)
	}
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

admixed.Covariance <- function(covariance,admix.proportions,nugget,k,inv.mean.sample.sizes,ident.mat){
	# recover()
	if(any(admix.proportions !=0)){
		w_k <- admix.proportions/2
		admixed.Covariance <- 	tcrossprod((1-w_k),(1-w_k)) * 	covariance[1:k,1:k] + 
								tcrossprod((1-w_k),(w_k)) 	* 	covariance[1:k,(k+1):(2*k)] +
								tcrossprod(w_k,(1-w_k)) 	*	covariance[(k+1):(2*k),1:k] +
								tcrossprod(w_k,w_k)			*	covariance[(k+1):(2*k),(k+1):(2*k)]
		admixed.Covariance <- admixed.Covariance + ident.mat * nugget + ident.mat * inv.mean.sample.sizes
	} else {
		admixed.Covariance <- covariance[1:k,1:k] + ident.mat * nugget + ident.mat * inv.mean.sample.sizes
	}
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

get.projection.matrix <- function(mean.sample.sizes,proj.mat.option=NULL,likelihood.option){
	if(likelihood.option == "wishart"){
		if(is.null(proj.mat.option)){
			k <- length(mean.sample.sizes)
			transformation.matrix <- get.transformation.matrix(mean.sample.sizes)
				qr.transformation.matrix <- qr(t(transformation.matrix))
				projection.matrix <- qr.Q(qr.transformation.matrix)[,1:qr.transformation.matrix$rank]
				stopifnot(qr.transformation.matrix$rank == sum(abs(eigen(transformation.matrix)$values - 1) < 1e-2) )
		} else if(proj.mat.option == "identity"){
			projection.matrix <- diag(length(mean.sample.sizes))
		}
	} else if(likelihood.option == "normal_approx"){
		projection.matrix <- get.transformation.matrix(mean.sample.sizes)
	}
	return(projection.matrix)
}

get.transformation.matrix <- function(mean.sample.sizes){
	k <- length(mean.sample.sizes)
	transformation.matrix <- diag(k) - matrix(mean.sample.sizes/(sum(mean.sample.sizes)),nrow=k,ncol=k,byrow=TRUE)
	return(transformation.matrix)
}

get.cov.variance <- function(k,mcn.freqs,mean.sample.cov,loci){
	var.sample.cov <- matrix(0,nrow=k,ncol=k)
		for(i in 1:k){
			for(j in 1:k){
				var.sample.cov[i,j] <- mean((mcn.freqs[i,]*mcn.freqs[j,] - mean.sample.cov[i,j])^2)
			}
		}
		var.sample.cov <- (1/(loci-1)) * var.sample.cov
	return(var.sample.cov)
}

get.cov.standard.error <- function(mcn.freqs){
	k <- nrow(mcn.freqs)
	loci <- ncol(mcn.freqs)
	mean.sample.cov <- cov(t(mcn.freqs))
	var.sample.cov <- get.cov.variance(k,mcn.freqs,mean.sample.cov,loci)
	std.error <- sqrt(var.sample.cov)/sqrt(loci-1)
	return(std.error)
}

project.sample.covariance <- function(sample.covariance,projection.matrix,likelihood.option){
	if(likelihood.option == "wishart"){
		sample.covariance <- t(projection.matrix) %*% sample.covariance %*% projection.matrix
	} else if(likelihood.option == "normal_approx"){
		sample.covariance <- sample.covariance
	}
	return(sample.covariance)
}

spacemix.data <- function(data.type,likelihood.option,proj.mat.option=NULL,sample.frequencies=NULL,loci,mean.sample.sizes=NULL,counts=NULL,sample.sizes=NULL,sample.covariance=NULL,cov.standard.error=NULL,prefix){
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
	projection.matrix <- get.projection.matrix(mean.sample.sizes,proj.mat.option,likelihood.option)
	sample.covariance <- project.sample.covariance(sample.covariance,projection.matrix,likelihood.option)
	sd <- NULL
	index.matrix <- NULL
	if(likelihood.option == "normal_approx"){
		if(is.null(mean.centered.normalized.frequencies) && is.null(cov.standard.error)){
			stop("you must either specify the standard error of the sample covariance\n provide SpaceMix with the allele counts or frequencies, so that it can calculate it for itself.")
		}
		if(is.null(cov.standard.error)){
			sd <- get.cov.standard.error(mean.centered.normalized.frequencies)
		}
		index.matrix <- upper.tri(sd,diag=TRUE)
		sd <- sd[index.matrix]
		sample.covariance <- sample.covariance[index.matrix]
	}
	spacemix.data <- list(	"sample.frequencies" = sample.frequencies,
							"mean.centered.normalized.frequencies" = mean.centered.normalized.frequencies,
							"inv.mean.sample.sizes" = 1/mean.sample.sizes,
							"projection.matrix" = projection.matrix,
							"sample.covariance" = sample.covariance,
							"loci" = loci,
							"sd" = sd,
							"index.matrix" = index.matrix)
		save(spacemix.data,file=paste(prefix,"spacemix.data.Robj",sep=''))
	return(spacemix.data)
}

propose.new.location.plane <- function(long,lat,dist.std){
	coords_prime <- c(long,lat) + rnorm(n = 2, mean = 0, sd = dist.std)
	return(coords_prime)
}

sphere.hop <- function(long,lat,distance,bearing){
	new.lat <- asin(sin(lat) * cos(distance) + cos(lat) * sin(distance) * cos(bearing))
	dlong <- atan2(sin(bearing)*sin(distance)*cos(lat),cos(distance)-sin(lat)*sin(new.lat))
	new.long <- (long - dlong + pi) %% (2*pi) - pi
	return(cbind(new.long,new.lat))
}

propose.new.location.sphere <- function(long,lat,dist.std){
	long <- degrees2radians(long)
	lat <- degrees2radians(lat)
	proposed.jump <- abs(rnorm(1,0,dist.std))
	jump.bearing <- runif(1,0,2*pi)
	coords_prime <- sphere.hop(long=long, lat=lat, distance=proposed.jump, bearing=jump.bearing)
	coords_prime <- radians2degrees(coords_prime)
	return(coords_prime)
}

degrees2radians <- function(degrees){
	radians <- (pi/180) * degrees
	return(radians)
}

radians2degrees <- function(radians){
	degrees <- (180/pi) * radians			
	return(degrees)
}

initiate.admix.proportions <- function(k,model.option){
	if(model.option == "no_movement"){
		admix.proportions <- numeric(k)
		prior_prob_admix_proportions <- 0
	} else if(model.option == "target"){
		admix.proportions <- numeric(k)
		prior_prob_admix_proportions <- 0
	} else if(model.option == "source"){
		admix.proportions <- matrix(rbeta(k,shape1=1,shape2=50),nrow=k,ncol=1)
		prior_prob_admix_proportions <- Prior_prob_admix_proportions(admix.proportions)
	} else if(model.option == "source_and_target"){
		admix.proportions <- matrix(rbeta(k,shape1=1,shape2=50),nrow=k,ncol=1)
		prior_prob_admix_proportions <- Prior_prob_admix_proportions(admix.proportions)
	}
	initiate.admix.proportions.list <- list("admix.proportions" = admix.proportions, 
											"prior_prob_admix_proportions" = prior_prob_admix_proportions)
	return(initiate.admix.proportions.list)
}

initiate.population.coordinates <- function(spatial.prior.X.coordinates,spatial.prior.Y.coordinates,k){
	population.coordinates <- 	rbind(
		cbind(spatial.prior.X.coordinates,
			spatial.prior.Y.coordinates),
		cbind(	runif(k, 
			min = min(spatial.prior.X.coordinates), 
			max = max(spatial.prior.X.coordinates)),
		runif(k, 
			min = min(spatial.prior.Y.coordinates), 
			max = max(spatial.prior.Y.coordinates)))
	)
	return(population.coordinates)
}

save.initial.parameters <- function(a0,aD,a2,nugget,admix.proportions,covariance,admixed.covariance,transformed.covariance,population.coordinates,D,projection.matrix,inv.mean.sample.sizes,prefix){
	initial.parameters <- list("a0" = a0,"aD" = aD,"a2" = a2,"nugget" = nugget,"admix.proportions" = admix.proportions,
								"covariance" = covariance,"admixed.covariance" = admixed.covariance,"transformed.covariance" = transformed.covariance,
								"population.coordinates" = population.coordinates,"D" = D,"projection.matrix" = projection.matrix,"inv.mean.sample.sizes"=inv.mean.sample.sizes)
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

initiate.lstps.list <- function(model.option,k,ngen,samplefreq,for.last.params){
	if(!for.last.params){
		if(model.option == "no_movement"){
			lstps <- list("a0_lstp" = numeric(ngen/samplefreq),
							"aD_lstp" = numeric(ngen/samplefreq),
							"a2_lstp" = numeric(ngen/samplefreq),
							"nugget_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		} else if(model.option == "target"){
			lstps <- list("a0_lstp" = numeric(ngen/samplefreq),
							"aD_lstp" = numeric(ngen/samplefreq),
							"a2_lstp" = numeric(ngen/samplefreq),
							"nugget_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_target_location_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		} else if(model.option == "source"){
			lstps <- list("a0_lstp" = numeric(ngen/samplefreq),
							"aD_lstp" = numeric(ngen/samplefreq),
							"a2_lstp" = numeric(ngen/samplefreq),
							"nugget_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_source_location_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_proportions_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		} else if(model.option == "source_and_target"){
			lstps <- list("a0_lstp" = numeric(ngen/samplefreq),
							"aD_lstp" = numeric(ngen/samplefreq),
							"a2_lstp" = numeric(ngen/samplefreq),
							"nugget_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_target_location_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_source_location_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_proportions_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		}
	} else {
		if(model.option == "no_movement"){
			lstps <- list("a0_lstp" = 0,
							"aD_lstp" = 0,
							"a2_lstp" = 0,
							"nugget_lstp" = numeric(k))
		} else if(model.option == "target"){
			lstps <- list("a0_lstp" = 0,
							"aD_lstp" = 0,
							"a2_lstp" = 0,
							"nugget_lstp" = numeric(k),
							"admix_target_location_lstp" = numeric(k))
		} else if(model.option == "source"){
			lstps <- list("a0_lstp" = 0,
							"aD_lstp" = 0,
							"a2_lstp" = 0,
							"nugget_lstp" = numeric(k),
							"admix_source_location_lstp" = numeric(k),
							"admix_proportions_lstp" = numeric(k))
		} else if(model.option == "source_and_target"){
			lstps <- list("a0_lstp" = 0,
							"aD_lstp" = 0,
							"a2_lstp" = 0,
							"nugget_lstp" = numeric(k),
							"admix_target_location_lstp" = numeric(k),
							"admix_source_location_lstp" = numeric(k),
							"admix_proportions_lstp" = numeric(k))
		}
	}
	return(lstps)
}

initiate.diagns.list <- function(model.option,k,mixing.diagn.freq){
	if(model.option == "no_movement"){
		diagns <- list("a0_diagn" = numeric(mixing.diagn.freq),
						"aD_diagn" = numeric(mixing.diagn.freq),
						"a2_diagn" = numeric(mixing.diagn.freq),
						"nugget_diagn" = matrix(0,nrow=k,ncol=mixing.diagn.freq))
	} else if(model.option == "target"){
		diagns <- list("a0_diagn" = numeric(mixing.diagn.freq),
						"aD_diagn" = numeric(mixing.diagn.freq),
						"a2_diagn" = numeric(mixing.diagn.freq),
						"nugget_diagn" = matrix(0,nrow=k,ncol=(mixing.diagn.freq)),
						"admix_target_location_diagn" = matrix(0,nrow=k,ncol=mixing.diagn.freq))
	} else if(model.option == "source"){
		diagns <- list("a0_diagn" = numeric(mixing.diagn.freq),
						"aD_diagn" = numeric(mixing.diagn.freq),
						"a2_diagn" = numeric(mixing.diagn.freq),
						"nugget_diagn" = matrix(0,nrow=k,ncol=mixing.diagn.freq),
						"admix_source_location_diagn" = matrix(0,nrow=k,ncol=mixing.diagn.freq),
						"admix_proportions_diagn" = matrix(0,nrow=k,ncol=mixing.diagn.freq))
	} else if(model.option == "source_and_target"){
		diagns <- list("a0_diagn" = numeric(mixing.diagn.freq),
						"aD_diagn" = numeric(mixing.diagn.freq),
						"a2_diagn" = numeric(mixing.diagn.freq),
						"nugget_diagn" = matrix(0,nrow=k,ncol=mixing.diagn.freq),
						"admix_target_location_diagn" = matrix(0,nrow=k,ncol=mixing.diagn.freq),
						"admix_source_location_diagn" = matrix(0,nrow=k,ncol=mixing.diagn.freq),
						"admix_proportions_diagn" = matrix(0,nrow=k,ncol=mixing.diagn.freq))
	}
	return(diagns)
}

initiate.moves.list <- function(model.option,k){
	if(model.option == "no_movement"){
		moves <- list("a0_moves" = 0,
						"aD_moves" = 0,
						"a2_moves" = 0,
						"nugget_moves" = numeric(k))
	} else if(model.option == "target"){
		moves <- list("a0_moves" = 0,
						"aD_moves" = 0,
						"a2_moves" = 0,
						"nugget_moves" = numeric(k),
						"admix_target_location_moves" = numeric(k))
	} else if(model.option == "source"){
		moves <- list("a0_moves" = 0,
						"aD_moves" = 0,
						"a2_moves" = 0,
						"nugget_moves" = numeric(k),
						"admix_source_location_moves" = numeric(k),
						"admix_proportions_moves" = numeric(k))
	} else if(model.option == "source_and_target"){
		moves <- list("a0_moves" = 0,
						"aD_moves" = 0,
						"a2_moves" = 0,
						"nugget_moves" = numeric(k),
						"admix_target_location_moves" = numeric(k),
						"admix_source_location_moves" = numeric(k),
						"admix_proportions_moves" = numeric(k))
	}
	return(moves)
}

initiate.accept.list <- function(model.option,k){
	if(model.option == "no_movement"){
		accept <- list("a0_accept" = 0,
						"aD_accept" = 0,
						"a2_accept" = 0,
						"nugget_accept" = numeric(k))
	} else if(model.option == "target"){
		accept <- list("a0_accept" = 0,
						"aD_accept" = 0,
						"a2_accept" = 0,
						"nugget_accept" = numeric(k),
						"admix_target_location_accept" = numeric(k))
	} else if(model.option == "source"){
		accept <- list("a0_accept" = 0,
						"aD_accept" = 0,
						"a2_accept" = 0,
						"nugget_accept" = numeric(k),
						"admix_source_location_accept" = numeric(k),
						"admix_proportions_accept" = numeric(k))
	} else if(model.option == "source_and_target"){
		accept <- list("a0_accept" = 0,
						"aD_accept" = 0,
						"a2_accept" = 0,
						"nugget_accept" = numeric(k),
						"admix_target_location_accept" = numeric(k),
						"admix_source_location_accept" = numeric(k),
						"admix_proportions_accept" = numeric(k))
	}
	return(accept)
}

initiate.accept.rates.list <- function(model.option,k,ngen,samplefreq,for.last.params){
	if(!for.last.params){
		if(model.option == "no_movement"){
			accept_rates <- list("a0_accept_rate" = numeric(ngen/samplefreq),
							"aD_accept_rate" = numeric(ngen/samplefreq),
							"a2_accept_rate" = numeric(ngen/samplefreq),
							"nugget_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		} else if(model.option == "target"){
			accept_rates <- list("a0_accept_rate" = numeric(ngen/samplefreq),
							"aD_accept_rate" = numeric(ngen/samplefreq),
							"a2_accept_rate" = numeric(ngen/samplefreq),
							"nugget_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_target_location_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		} else if(model.option == "source"){
			accept_rates <- list("a0_accept_rate" = numeric(ngen/samplefreq),
							"aD_accept_rate" = numeric(ngen/samplefreq),
							"a2_accept_rate" = numeric(ngen/samplefreq),
							"nugget_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_source_location_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_proportions_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		} else if(model.option == "source_and_target"){
			accept_rates <- list("a0_accept_rate" = numeric(ngen/samplefreq),
							"aD_accept_rate" = numeric(ngen/samplefreq),
							"a2_accept_rate" = numeric(ngen/samplefreq),
							"nugget_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_target_location_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_source_location_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_proportions_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		}
	} else {
		if(model.option == "no_movement"){
			accept_rates <- list("a0_accept_rate" = 0,
							"aD_accept_rate" = 0,
							"a2_accept_rate" = 0,
							"nugget_accept_rate" = numeric(k))
		} else if(model.option == "target"){
			accept_rates <- list("a0_accept_rate" = 0,
							"aD_accept_rate" = 0,
							"a2_accept_rate" = 0,
							"nugget_accept_rate" = numeric(k),
							"admix_target_location_accept_rate" = numeric(k))
		} else if(model.option == "source"){
			accept_rates <- list("a0_accept_rate" = 0,
							"aD_accept_rate" = 0,
							"a2_accept_rate" = 0,
							"nugget_accept_rate" = numeric(k),
							"admix_source_location_accept_rate" = numeric(k),
							"admix_proportions_accept_rate" = numeric(k))
		} else if(model.option == "source_and_target"){
			accept_rates <- list("a0_accept_rate" = 0,
							"aD_accept_rate" = 0,
							"a2_accept_rate" = 0,
							"nugget_accept_rate" = numeric(k),
							"admix_target_location_accept_rate" = numeric(k),
							"admix_source_location_accept_rate" = numeric(k),
							"admix_proportions_accept_rate" = numeric(k))
		}
	}
	return(accept_rates)
}

print.mcmc.update <- function(LnL_freqs,prior_prob_admix_proportions,prior_prob_nugget,prior_prob_alpha0,prior_prob_alphaD,prior_prob_alpha2,prior_prob_admix_target_locations,prior_prob_admix_source_locations){
	P <- LnL_freqs + prior_prob_admix_proportions + prior_prob_nugget + prior_prob_alpha0 + prior_prob_alphaD + prior_prob_alpha2 + prior_prob_admix_target_locations + prior_prob_admix_source_locations
	return(P)
}

make.update.sampled.accept.rates.function <- function(model.option){
	if(model.option == "no_movement"){
		update.sampled.accept.rates <<- function(accept_rates,j,new.params){
			accept_rates$a0_accept_rate[j] <- new.params$accept_rates$a0_accept_rate
			accept_rates$aD_accept_rate[j] <- new.params$accept_rates$aD_accept_rate
			accept_rates$a2_accept_rate[j] <- new.params$accept_rates$a2_accept_rate
			accept_rates$nugget_accept_rate[,j] <- new.params$accept_rates$nugget_accept_rate
			return(accept_rates)
		}
	} else if(model.option == "target"){
		update.sampled.accept.rates <<- function(accept_rates,j,new.params){
			accept_rates$a0_accept_rate[j] <- new.params$accept_rates$a0_accept_rate
			accept_rates$aD_accept_rate[j] <- new.params$accept_rates$aD_accept_rate
			accept_rates$a2_accept_rate[j] <- new.params$accept_rates$a2_accept_rate
			accept_rates$nugget_accept_rate[,j] <- new.params$accept_rates$nugget_accept_rate
			accept_rates$admix_target_location_accept_rate[,j] <- new.params$accept_rates$admix_target_location_accept_rate
			return(accept_rates)
		}
	} else if(model.option == "source"){
		update.sampled.accept.rates <<- function(accept_rates,j,new.params){
			accept_rates$a0_accept_rate[j] <- new.params$accept_rates$a0_accept_rate
			accept_rates$aD_accept_rate[j] <- new.params$accept_rates$aD_accept_rate
			accept_rates$a2_accept_rate[j] <- new.params$accept_rates$a2_accept_rate
			accept_rates$nugget_accept_rate[,j] <- new.params$accept_rates$nugget_accept_rate
			accept_rates$admix_source_location_accept_rate[,j] <- new.params$accept_rates$admix_source_location_accept_rate
			accept_rates$admix_proportions_accept_rate[,j] <- new.params$accept_rates$admix_proportions_accept_rate
			return(accept_rates)
		}
	} else if(model.option == "source_and_target"){
		update.sampled.accept.rates <<- function(accept_rates,j,new.params){
			accept_rates$a0_accept_rate[j] <- new.params$accept_rates$a0_accept_rate
			accept_rates$aD_accept_rate[j] <- new.params$accept_rates$aD_accept_rate
			accept_rates$a2_accept_rate[j] <- new.params$accept_rates$a2_accept_rate
			accept_rates$nugget_accept_rate[,j] <- new.params$accept_rates$nugget_accept_rate
			accept_rates$admix_target_location_accept_rate[,j] <- new.params$accept_rates$admix_target_location_accept_rate
			accept_rates$admix_source_location_accept_rate[,j] <- new.params$accept_rates$admix_source_location_accept_rate
			accept_rates$admix_proportions_accept_rate[,j] <- new.params$accept_rates$admix_proportions_accept_rate
			return(accept_rates)
		}
	}
	return(0)
}

make.update.sampled.lstps.function <- function(model.option){
	if(model.option == "no_movement"){
		update.sampled.lstps <<- function(lstps,j,new.params){
			lstps$a0_lstp[j] <- new.params$lstps$a0_lstp
			lstps$aD_lstp[j] <- new.params$lstps$aD_lstp
			lstps$a2_lstp[j] <- new.params$lstps$a2_lstp
			lstps$nugget_lstp[,j] <- new.params$lstps$nugget_lstp
			return(lstps)
		}
	} else if(model.option == "target"){
		update.sampled.lstps <<- function(lstps,j,new.params){
			lstps$a0_lstp[j] <- new.params$lstps$a0_lstp
			lstps$aD_lstp[j] <- new.params$lstps$aD_lstp
			lstps$a2_lstp[j] <- new.params$lstps$a2_lstp
			lstps$nugget_lstp[,j] <- new.params$lstps$nugget_lstp
			lstps$admix_target_location_lstp[,j] <- new.params$lstps$admix_target_location_lstp
			return(lstps)
		}	
	} else if(model.option == "source"){
		update.sampled.lstps <<- function(lstps,j,new.params){
			lstps$a0_lstp[j] <- new.params$lstps$a0_lstp
			lstps$aD_lstp[j] <- new.params$lstps$aD_lstp
			lstps$a2_lstp[j] <- new.params$lstps$a2_lstp
			lstps$nugget_lstp[,j] <- new.params$lstps$nugget_lstp
			lstps$admix_source_location_lstp[,j] <- new.params$lstps$admix_source_location_lstp
			lstps$admix_proportions_lstp[,j] <- new.params$lstps$admix_proportions_lstp
			return(lstps)
		}
	} else if(model.option == "source_and_target"){
		update.sampled.lstps <<- function(lstps,j,new.params){
			lstps$a0_lstp[j] <- new.params$lstps$a0_lstp
			lstps$aD_lstp[j] <- new.params$lstps$aD_lstp
			lstps$a2_lstp[j] <- new.params$lstps$a2_lstp
			lstps$nugget_lstp[,j] <- new.params$lstps$nugget_lstp
			lstps$admix_target_location_lstp[,j] <- new.params$lstps$admix_target_location_lstp
			lstps$admix_source_location_lstp[,j] <- new.params$lstps$admix_source_location_lstp
			lstps$admix_proportions_lstp[,j] <- new.params$lstps$admix_proportions_lstp
			return(lstps)
		}
	}
	return(0)
}

get.diagn.step <- function(generation,mixing.diagn.freq){
	diagn.step <- generation%%mixing.diagn.freq
	if(diagn.step == 0){
		diagn.step <- mixing.diagn.freq
	}
	return(diagn.step)
}

make.update.diagns.function <- function(model.option){
	if(model.option == "no_movement"){
		update.diagns <<- function(diagns,generation,mixing.diagn.freq,accept_rates){
			diagn.step <- get.diagn.step(generation,mixing.diagn.freq)
				diagns$a0_diagn[diagn.step] <- accept_rates$a0_accept_rate
				diagns$aD_diagn[diagn.step] <- accept_rates$aD_accept_rate
				diagns$a2_diagn[diagn.step] <- accept_rates$a2_accept_rate
				diagns$nugget_diagn[,diagn.step] <- accept_rates$nugget_accept_rate
			return(diagns)
		}
	} else if(model.option == "target"){
		update.diagns <<- function(diagns,generation,mixing.diagn.freq,accept_rates){
			diagn.step <- get.diagn.step(generation,mixing.diagn.freq)
				diagns$a0_diagn[diagn.step] <- accept_rates$a0_accept_rate
				diagns$aD_diagn[diagn.step] <- accept_rates$aD_accept_rate
				diagns$a2_diagn[diagn.step] <- accept_rates$a2_accept_rate
				diagns$nugget_diagn[,diagn.step] <- accept_rates$nugget_accept_rate
				diagns$admix_target_location_diagn[,diagn.step] <- accept_rates$admix_target_location_accept_rate
			return(diagns)
		}
	} else if(model.option == "source"){
		update.diagns <<- function(diagns,generation,mixing.diagn.freq,accept_rates){
			diagn.step <- get.diagn.step(generation,mixing.diagn.freq)
				diagns$a0_diagn[diagn.step] <- accept_rates$a0_accept_rate
				diagns$aD_diagn[diagn.step] <- accept_rates$aD_accept_rate
				diagns$a2_diagn[diagn.step] <- accept_rates$a2_accept_rate
				diagns$nugget_diagn[,diagn.step] <- accept_rates$nugget_accept_rate
				diagns$admix_source_location_diagn[,diagn.step] <- accept_rates$admix_source_location_accept_rate
				diagns$admix_proportions_diagn[,diagn.step] <- accept_rates$admix_proportions_accept_rate
			return(diagns)
		}
	} else if(model.option == "source_and_target"){
		update.diagns <<- function(diagns,generation,mixing.diagn.freq,accept_rates){
			diagn.step <- get.diagn.step(generation,mixing.diagn.freq)
				diagns$a0_diagn[diagn.step] <- accept_rates$a0_accept_rate
				diagns$aD_diagn[diagn.step] <- accept_rates$aD_accept_rate
				diagns$a2_diagn[diagn.step] <- accept_rates$a2_accept_rate
				diagns$nugget_diagn[,diagn.step] <- accept_rates$nugget_accept_rate
				diagns$admix_target_location_diagn[,diagn.step] <- accept_rates$admix_target_location_accept_rate
				diagns$admix_source_location_diagn[,diagn.step] <- accept_rates$admix_source_location_accept_rate
				diagns$admix_proportions_diagn[,diagn.step] <- accept_rates$admix_proportions_accept_rate
			return(diagns)
		}
	}
	return(0)
}

update.lstp <- function(n,lstp,acceptance.fraction){
	if(length(lstp) > 1){
		for(i in 1:length(lstp)){
			if(acceptance.fraction[i] > 0.44){
				lstp[i] <- min(lstp[i] + min(0.01,n^(-0.5)),20)
			} else if(acceptance.fraction[i] < 0.44){
				lstp[i] <- lstp[i] - min(0.01,n^(-0.5))
			}
		}
	} else {
		if(acceptance.fraction > 0.44){
			lstp <- min(lstp + min(0.01,n^(-0.5)),20)
		} else if(acceptance.fraction < 0.44){
			lstp <- lstp - min(0.01,n^(-0.5))
		}
	}
	return(lstp)
}

make.update.all.lstps.function <- function(model.option){
	if(model.option == "no_movement"){
		update.all.lstps <<- function(lstps,n,diagns){
			lstps$a0_lstp <- update.lstp(n,lstps$a0_lstp,mean(diagns$a0_diagn))
			lstps$aD_lstp <- update.lstp(n,lstps$aD_lstp,mean(diagns$aD_diagn))
			lstps$a2_lstp <- update.lstp(n,lstps$a2_lstp,mean(diagns$a2_diagn))
			lstps$nugget_lstp <- update.lstp(n,lstps$nugget_lstp,rowMeans(diagns$nugget_diagn))
			return(lstps)
		}
	} else if(model.option == "target"){
		update.all.lstps <<- function(lstps,n,diagns){
			lstps$a0_lstp <- update.lstp(n,lstps$a0_lstp,mean(diagns$a0_diagn))
			lstps$aD_lstp <- update.lstp(n,lstps$aD_lstp,mean(diagns$aD_diagn))
			lstps$a2_lstp <- update.lstp(n,lstps$a2_lstp,mean(diagns$a2_diagn))
			lstps$nugget_lstp <- update.lstp(n,lstps$nugget_lstp,rowMeans(diagns$nugget_diagn))
			lstps$admix_target_location_lstp <- update.lstp(n,lstps$admix_target_location_lstp,rowMeans(diagns$admix_target_location_diagn))
			return(lstps)
		}
	} else if(model.option == "source"){
		update.all.lstps <<- function(lstps,n,diagns){
			lstps$a0_lstp <- update.lstp(n,lstps$a0_lstp,mean(diagns$a0_diagn))
			lstps$aD_lstp <- update.lstp(n,lstps$aD_lstp,mean(diagns$aD_diagn))
			lstps$a2_lstp <- update.lstp(n,lstps$a2_lstp,mean(diagns$a2_diagn))
			lstps$nugget_lstp <- update.lstp(n,lstps$nugget_lstp,rowMeans(diagns$nugget_diagn))
			lstps$admix_source_location_lstp <- update.lstp(n,lstps$admix_source_location_lstp,rowMeans(diagns$admix_source_location_diagn))
			lstps$admix_proportions_lstp <- update.lstp(n,lstps$admix_proportions_lstp,rowMeans(diagns$admix_proportions_diagn))
			return(lstps)
		}
	} else if(model.option == "source_and_target"){
		update.all.lstps <<- function(lstps,n,diagns){
			lstps$a0_lstp <- update.lstp(n,lstps$a0_lstp,mean(diagns$a0_diagn))
			lstps$aD_lstp <- update.lstp(n,lstps$aD_lstp,mean(diagns$aD_diagn))
			lstps$a2_lstp <- update.lstp(n,lstps$a2_lstp,mean(diagns$a2_diagn))
			lstps$nugget_lstp <- update.lstp(n,lstps$nugget_lstp,rowMeans(diagns$nugget_diagn))
			lstps$admix_target_location_lstp <- update.lstp(n,lstps$admix_target_location_lstp,rowMeans(diagns$admix_target_location_diagn))
			lstps$admix_source_location_lstp <- update.lstp(n,lstps$admix_source_location_lstp,rowMeans(diagns$admix_source_location_diagn))
			lstps$admix_proportions_lstp <- update.lstp(n,lstps$admix_proportions_lstp,rowMeans(diagns$admix_proportions_diagn))
			return(lstps)
		}
	}
	return(0)
}

initiate.last.params <- function(model.option,likelihood.option,samplefreq,ngen,spacemix.data,population.coordinates,admix.proportions,a0,aD,a2,nugget,covariance,admixed.covariance,transformed.covariance,
						k,LnL_freqs,prior_prob_alpha0,prior_prob_alphaD,prior_prob_alpha2,prior_prob_nugget,prior_prob_admix_proportions,prior_prob_admix_target_locations,prior_prob_admix_source_locations,
						D,spatial.prior.X.coordinates,spatial.prior.Y.coordinates,target.spatial.prior.scale,source.spatial.prior.scale,centroid,gibbs.spatial.fineness,gibbs.nugget.fineness){
	last.params <- list("sample.covariance" = spacemix.data$sample.covariance,
						"projection.matrix" = spacemix.data$projection.matrix,						
						"population.coordinates" = population.coordinates,
						"admix.proportions" = admix.proportions,
						"a0" = a0,"aD" = aD,"a2" = a2,"nugget" = nugget,
						"covariance" = covariance,"admixed.covariance" = admixed.covariance,
						"transformed.covariance" = transformed.covariance,
						"lstps" = initiate.lstps.list(model.option,k,ngen,samplefreq,for.last.params=TRUE),
						"k" = k,"LnL_freqs" = LnL_freqs,
						"prior_prob_alpha0" = prior_prob_alpha0,
						"prior_prob_alphaD" = prior_prob_alphaD,
						"prior_prob_alpha2" = prior_prob_alpha2,
						"prior_prob_nugget" = prior_prob_nugget,
						"prior_prob_admix_proportions" = prior_prob_admix_proportions, 
						"prior_prob_admix_target_locations" = prior_prob_admix_target_locations,
						"prior_prob_admix_source_locations" = prior_prob_admix_source_locations,
						"accept_rates" = initiate.accept.rates.list(model.option,k,ngen,samplefreq,for.last.params=TRUE),
						"moves" = initiate.moves.list(model.option,k),"accept" = initiate.accept.list(model.option,k),
						"loci" = spacemix.data$loci,"D" = D,
						"inv.mean.sample.sizes" = spacemix.data$inv.mean.sample.sizes,
						"spatial.prior.X.coordinates" = spatial.prior.X.coordinates,
						"spatial.prior.Y.coordinates" = spatial.prior.Y.coordinates,
						"target.spatial.prior.scale" = target.spatial.prior.scale,
						"source.spatial.prior.scale" = source.spatial.prior.scale,
						"centroid" = centroid,
						"X.grid.fineness" = gibbs.spatial.fineness, 
						"Y.grid.fineness" = gibbs.spatial.fineness, 
						"nugget.grid.fineness" = gibbs.nugget.fineness,
						"sd" = spacemix.data$sd,
						"index.matrix" = spacemix.data$index.matrix,
						"likelihood.option" = likelihood.option,
						"identity.matrix" = diag(k))
	return(last.params)
}


}	
	
	
MCMC <-function(model.option,				#no_movement, target, source, source_and_target
				data.type,					#sample.covariance, sample.frequencies, counts
				likelihood.option,			#normal_approx,wishart
				proj.mat.option = NULL,
				sample.frequencies = NULL,
				mean.sample.sizes = NULL,
				counts = NULL,
				sample.sizes = NULL,
				sample.covariance = NULL,
				target.spatial.prior.scale = NULL,
				source.spatial.prior.scale = NULL,
				spatial.prior.X.coordinates,
				spatial.prior.Y.coordinates,
				initial.parameters = NULL,		#a0,aD,a2,population.coordinates,admix.proportions
				round.earth,
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
									likelihood.option = likelihood.option,
									sample.frequencies = sample.frequencies, 
									loci = loci,
									mean.sample.sizes = mean.sample.sizes,
									counts = counts,
									sample.sizes = sample.sizes,
									sample.covariance = sample.covariance,
									prefix = prefix)
	
	if(round.earth){
		propose.new.location <<- propose.new.location.sphere
		spacemix.dist <<- fields::rdist.earth
	} else {
		propose.new.location <<- propose.new.location.plane
		spacemix.dist <<- fields::rdist
	}
	
	#Declare variables
		LnL_freqs <- numeric(ngen/samplefreq)
		Prob <- numeric(ngen/samplefreq)
		population.coordinates <- vector("list",ngen/samplefreq)
		distances <- vector("list",ngen/samplefreq)
		transformed.covariance.list <- vector("list",ngen/samplefreq)
		admix.proportions <- matrix(0,nrow=k,ncol=ngen/samplefreq)
		nugget <- matrix(0,nrow=k,ncol=ngen/samplefreq)
		a0 <- numeric(ngen/samplefreq)
		aD <- numeric(ngen/samplefreq)
		a2 <- numeric(ngen/samplefreq)
		accept_rates <- initiate.accept.rates.list(model.option,k,ngen,samplefreq,for.last.params=FALSE)
		lstps <- initiate.lstps.list(model.option,k,ngen,samplefreq,for.last.params=FALSE)
		diagns <- initiate.diagns.list(model.option,k,mixing.diagn.freq)
		
		seed <- sample(0:1e6,1)
			save(seed,file=paste(prefix,"_seed.Robj",sep=''))
		set.seed(seed)
		
	if(!continue) {
		#INITIALIZE MCMC
				Prob[1] <- -Inf
				covariance <- matrix(0,nrow=k*2,ncol=k*2)
  				admixed.covariance <- matrix(0,nrow=k,ncol=k)
				transformed.covariance <- transformed.Covariance(admixed.covariance,spacemix.data$projection.matrix)
				badness.counter <- 0

			while(!is.finite(Prob[1]) && badness.counter < 100){
					if(!is.null(initial.parameters$nugget)){
						nugget[,1] <- initial.parameters$nugget
					} else { nugget[,1] <- rexp(k) }
					if(!is.null(initial.parameters$a0)){
						a0[1] <- initial.parameters$a0
					} else { a0[1] <- rexp(1,1/100) }
					if(!is.null(initial.parameters$aD)){
						aD[1] <- initial.parameters$aD
					} else { aD[1] <- rexp(1,1) }
					if(!is.null(initial.parameters$a2)){
						a2[1] <- initial.parameters$a2
					} else { a2[1] <- runif(1,0.1,2) }
					if(!is.null(initial.parameters$population.coordinates)){
						population.coordinates[[1]] <- 	initial.parameters$population.coordinates
					} else { population.coordinates[[1]] <-	initiate.population.coordinates(spatial.prior.X.coordinates,spatial.prior.Y.coordinates,k) }
					if(!is.null(initial.parameters$admix.proportions)){
						admix.proportions[,1] <- initial.parameters$admix.proportions
						prior_prob_admix_proportions <- Prior_prob_admix_proportions(admix.proportions[,1])
					} else { initiate.admix.proportions.list <- initiate.admix.proportions(k,model.option)
								admix.proportions[,1] <- initiate.admix.proportions.list$admix.proportions 
								prior_prob_admix_proportions <- initiate.admix.proportions.list$prior_prob_admix_proportions}
					distances[[1]] <- spacemix.dist(population.coordinates[[1]])
					centroid <- c(mean(spatial.prior.X.coordinates),mean(spatial.prior.Y.coordinates))
					if(is.null(target.spatial.prior.scale)){
						target.spatial.prior.scale <- mean(spacemix.dist(cbind(spatial.prior.X.coordinates,spatial.prior.Y.coordinates))) / 2
					}
					if(is.null(source.spatial.prior.scale)){
						source.spatial.prior.scale <- mean(spacemix.dist(cbind(spatial.prior.X.coordinates,spatial.prior.Y.coordinates))) * 2
					}
					covariance <- Covariance(a0[1],aD[1],a2[1],distances[[1]])
					admixed.covariance <- admixed.Covariance(covariance,admix.proportions[,1],nugget[,1],k,spacemix.data$inv.mean.sample.sizes,diag(k))
					transformed.covariance <- transformed.Covariance(admixed.covariance,spacemix.data$projection.matrix)
					tmp <- save.initial.parameters(a0[1],aD[1],a2[1],nugget[,1],admix.proportions[,1],covariance,admixed.covariance,transformed.covariance,population.coordinates[[1]],distances[[1]],spacemix.data$projection.matrix,spacemix.data$inv.mean.sample.sizes,prefix)
				LnL_freqs[1] <- likelihood(likelihood.option,spacemix.data$sample.covariance, transformed.covariance,spacemix.data$index.matrix,spacemix.data$sd,spacemix.data$loci)
					cat("LnL: ",LnL_freqs[1],"\n")
				prior_prob_alpha0 <- Prior_prob_alpha0(a0[1])
					cat("Pr(a0): ",prior_prob_alpha0,"\n")
				prior_prob_alphaD <- Prior_prob_alphaD(aD[1])
					cat("Pr(aD): ",prior_prob_alphaD,"\n")
				prior_prob_alpha2 <- Prior_prob_alpha2(a2[1])
					cat("Pr(a2): ",prior_prob_alpha2,"\n")
				prior_prob_nugget <- Prior_prob_nugget(nugget[,1])
					cat("Pr(nugget): ",prior_prob_nugget,"\n")
				prior_prob_admix_target_locations <- Prior_prob_admix_target_locations(population.coordinates[[1]][1:k,],spatial.prior.X.coordinates,spatial.prior.Y.coordinates,target.spatial.prior.scale)
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
					stop("Something has gone wrong with the attempt to initialize the MCMC! Please take a close look at your data and try again")
				}
			}
		}	
		last.params <- initiate.last.params(model.option = model.option,likelihood.option = likelihood.option,samplefreq = samplefreq,ngen = ngen,spacemix.data = spacemix.data,population.coordinates = population.coordinates[[1]],admix.proportions = admix.proportions[,1],
											a0[1],aD[1],a2[1],nugget[,1],covariance,admixed.covariance,transformed.covariance,k = k,LnL_freqs = LnL_freqs[1],
											prior_prob_alpha0 = prior_prob_alpha0,prior_prob_alphaD = prior_prob_alphaD,prior_prob_alpha2 = prior_prob_alpha2,
											prior_prob_nugget = prior_prob_nugget,prior_prob_admix_proportions = prior_prob_admix_proportions,
											prior_prob_admix_target_locations = prior_prob_admix_target_locations,prior_prob_admix_source_locations = prior_prob_admix_source_locations,
											D = distances[[1]],spatial.prior.X.coordinates = spatial.prior.X.coordinates,spatial.prior.Y.coordinates = spatial.prior.Y.coordinates,
											target.spatial.prior.scale = target.spatial.prior.scale,source.spatial.prior.scale = source.spatial.prior.scale,
											centroid = centroid,gibbs.spatial.fineness = gibbs.spatial.fineness,gibbs.nugget.fineness = gibbs.nugget.fineness)
		last.ngen <- 0
	} else {
		load(continuing.params)
		#FIX w/r/t diagn,lstp,accept_rate,moves,accept
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
		admix.proportions[,1] <- continuing.params$admix.proportions
		prior_prob_admix_proportions <- continuing.params$prior_prob_admix_proportions
					distances[[1]] <- continuing.params$D
						centroid <- continuing.params$centroid
						target.spatial.prior.scale <- continuing.params$target.spatial.prior.scale
						source.spatial.prior.scale <- continuing.params$source.spatial.prior.scale
					covariance <- continuing.params$covariance
					admixed.covariance <- continuing.params$admixed.covariance
					transformed.covariance <- continuing.params$transformed.covariance
					tmp <- save.initial.parameters(a0[1],aD[1],a2[1],nugget[,1],admix.proportions[,1],
													covariance,admixed.covariance,transformed.covariance,
													population.coordinates[[1]],distances[[1]],spacemix.data$projection.matrix,spacemix.data$inv.mean.sample.sizes,prefix)
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
											admix.proportions = admix.proportions[,1],
											a0[1],aD[1],a2[1],nugget[,1],covariance,admixed.covariance,transformed.covariance,
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
											spatial.prior.X.coordinates = continuing.params$spatial.prior.X.coordinates,
											spatial.prior.Y.coordinates = continuing.params$spatial.prior.Y.coordinates,
											target.spatial.prior.scale = continuing.params$target.spatial.prior.scale,
											source.spatial.prior.scale = continuing.params$source.spatial.prior.scale,
											centroid = continuing.params$centroid,
											gibbs.spatial.fineness = continuing.params$gibbs.spatial.fineness,
											gibbs.nugget.fineness = continuing.params$gibbs.nugget.fineness)
		last.ngen <- continuing.params$last.ngen
	}
	
	tmp <- make.update.sampled.accept.rates.function(model.option)
	tmp <- make.update.sampled.lstps.function(model.option)
	tmp <- make.update.diagns.function(model.option)
	tmp <- make.update.all.lstps.function(model.option)
	
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
			transformed.covariance.list[[j]] <- new.params$transformed.covariance
			admix.proportions[,j] <- new.params$admix.proportions
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
			accept_rates <- update.sampled.accept.rates(accept_rates,j,new.params)
			lstps <- update.sampled.lstps(lstps,j,new.params)
		}
		
		diagns <- update.diagns(diagns,i,mixing.diagn.freq,new.params$accept_rates)
		
		if(i%%mixing.diagn.freq == 0){
			n <- (i + last.ngen) %/% mixing.diagn.freq
			new.params$lstps <- update.all.lstps(new.params$lstps,n,diagns)
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
				LnL_freqs,Prob,covariance,admixed.covariance,transformed.covariance,distances,
				population.coordinates,transformed.covariance.list,admix.proportions,a0,aD,a2,nugget,samplefreq,ngen,
				accept_rates,lstps,diagns,
				target.spatial.prior.scale,source.spatial.prior.scale,last.ngen,
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
									"transformed.covariance" = transformed.covariance,"LnL_freqs" = LnL_freqs,
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
									"spatial.prior.X.coordinates" = spatial.prior.X.coordinates,"spatial.prior.Y.coordinates" = spatial.prior.Y.coordinates, 
									"X.grid.fineness" = X.grid.fineness,"Y.grid.fineness" = Y.grid.fineness,"nugget.grid.fineness" = nugget.grid.fineness,
									"a0.lstp" = a0.lstp,"aD.lstp" = aD.lstp,"a2.lstp" = a2.lstp,
									"nugget.lstp" = nugget.lstp,"admix.target.location.lstp" = admix.target.location.lstp,
									"admix.source.location.lstp" = admix.source.location.lstp,"admix.proportions.lstp" = admix.proportions.lstp,
									"a0_diagn" = a0_diagn,"a2_diagn" = a2_diagn,"aD_diagn" = aD_diagn,"admix_proportions_diagn" = admix_proportions_diagn,
									"admix_source_location_diagn" = admix_source_location_diagn,"admix_target_location_diagn" = admix_target_location_diagn,
									"nugget_diagn" = nugget_diagn, "last.ngen" = ngen)
	save(continuing.params,file=file.name)
	})
	return(0)
}

drop.objects.pre.link.up <- function(object,length){
	if(class(object)=="matrix"){
		object.dim <- dim(object)
	} else {
		object.dim <- length(object)
	}
	object.drop <- !any(grepl(length,object.dim))
	return(object.drop)
}

link.up.posteriors <- function(MCMC.output1,MCMC.output2,linked.up.output.file.name){
	recover()
	load(MCMC.output1)
		rm(list = setdiff(unique(c(	ls(pattern="diagn"),
							ls(pattern="last.params"),
							objects()[sapply(mget(objects()), drop.objects.pre.link.up,length=ngen/samplefreq)])),c("MCMC.output2","linked.up.output.file.name")))
    parameters <- objects()
    for(i in 1:length(parameters)){
        assign(sprintf("tmp.%s", parameters[i]), get(parameters[i]))
    }
    load(MCMC.output2)
    for (i in 1:length(parameters)) {
    	if(class(get(parameters[i])) == "numeric"){
    		assign(parameters[i],c(get(sprintf("tmp.%s", 
                    parameters[i])), get(parameters[i])))
    	} else if(class(get(parameters[i])) == "matrix"){
    		assign(parameters[i],cbind(get(sprintf("tmp.%s", 
                    parameters[i])), get(parameters[i])))
    	} else if(class(get(parameters[i])) == "list"){
    		assign(parameters[i],c(get(sprintf("tmp.%s", 
                    parameters[i])), get(parameters[i])))
    	}
    }
    rm(list = objects(pattern = "tmp."))
    save(list = setdiff(ls(all.names = TRUE), "linked.up.output.file.name"), 
        file = paste(linked.up.output.file.name, ".Robj", sep = ""))
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

get.mean.centering.matrix <- function(focal.pop,k,mean.sample.sizes=NULL){
	if(is.null(mean.sample.sizes)){
		mean.centering.matrix <- matrix(-1/(k-1),nrow=k,ncol=k,byrow=TRUE)
			diag(mean.centering.matrix) <- (k-2)/(k-1)
		mean.centering.matrix[,focal.pop] <- c(rep(0,focal.pop-1),1,rep(0,k-focal.pop))
	} else {
	mean.centering.matrix <- diag(k) - matrix(mean.sample.sizes/(sum(mean.sample.sizes[-focal.pop])),nrow=k,ncol=k,byrow=TRUE)
		mean.centering.matrix[,focal.pop] <- c(rep(0,focal.pop-1),1,rep(0,k-focal.pop))
	}
	return(mean.centering.matrix)		
}

get.uprank.mean.centering.matrix <- function(focal.pop,k,mean.sample.sizes=NULL){
	if(is.null(mean.sample.sizes)){
		mean.centering.matrix <- matrix(-1/(k-1),nrow=k,ncol=k,byrow=TRUE)
			diag(mean.centering.matrix) <- (k-2)/(k-1)
		mean.centering.matrix <- cbind(
									rbind(mean.centering.matrix,
											rep(-1/(k-1),k)),
									c(rep(0,k),1))
		mean.centering.matrix[,focal.pop] <- c(rep(0,focal.pop-1),1,rep(0,(k+1)-focal.pop))
	} else {
	mean.centering.matrix <- diag(k) - matrix(mean.sample.sizes/(sum(mean.sample.sizes[-focal.pop])),nrow=k,ncol=k,byrow=TRUE)
	mean.centering.matrix <- 	cbind(
									rbind(mean.centering.matrix,
											-mean.sample.sizes/(sum(mean.sample.sizes[-focal.pop]))),
									c(rep(0,k),1))
		mean.centering.matrix[,focal.pop] <- c(rep(0,focal.pop-1),1,rep(0,(k+1)-focal.pop))
	}
	return(mean.centering.matrix)		
}

calculate.fstat <- function(focal.population.observations,focal.pop.conditional.mean,admix.source.observations){
	fstat <- ((focal.population.observations - focal.pop.conditional.mean) %*% admix.source.observations)#/length(focal.population.observations)
	fstat <- fstat / sqrt(var(c(focal.population.observations - focal.pop.conditional.mean))*var(c(admix.source.observations)))
	fstat <- fstat/length(focal.population.observations)
	return(fstat)
}

get.fstat.pop <- function(focal.pop,source.pop,observations,covariance,mean.sample.sizes=NULL){
	k <- nrow(observations)
	mean.center.matrix <- get.mean.centering.matrix(focal.pop,k,mean.sample.sizes)
	mean.centered.covariance <- mean.center.matrix %*% covariance %*% t(mean.center.matrix)
	mean.centered.observations <- mean.center.matrix %*% observations
	pop.2.drop <- sample(c(1:k)[-c(focal.pop,source.pop)],1)
	if(pop.2.drop < focal.pop){
		focal.pop <- focal.pop - 1 
	}
	if(pop.2.drop < source.pop){
		source.pop <- source.pop - 1 
	}
	mean.centered.covariance <- mean.centered.covariance[-pop.2.drop,-pop.2.drop]
	mean.centered.observations <- mean.centered.observations[-pop.2.drop,]
	focal.pop.conditional.mean <- get.conditional.mean(mean.centered.covariance,focal.pop,mean.centered.observations,drop.option=1)
	fstat <- calculate.fstat(mean.centered.observations[focal.pop,],focal.pop.conditional.mean,mean.centered.observations[source.pop,])
	return(fstat)
}

get.focal.conditional.mean.fstat.location <- function(focal.pop,covariance,observations,mean.sample.sizes=NULL){
	k <- nrow(observations)
	mean.center.matrix <- get.mean.centering.matrix(focal.pop,k,mean.sample.sizes)
	mean.centered.covariance <- mean.center.matrix %*% covariance %*% t(mean.center.matrix)
	mean.centered.observations <- mean.center.matrix %*% observations
	dropped.focal.pop <- focal.pop
	pop.2.drop <- sample(c(1:k)[-c(focal.pop)],1)
	if(pop.2.drop < focal.pop){
		dropped.focal.pop <- focal.pop - 1 
	}
	dropped.mean.centered.covariance <- mean.centered.covariance[-pop.2.drop,-pop.2.drop]
	dropped.mean.centered.observations <- mean.centered.observations[-pop.2.drop,]
	focal.pop.conditional.mean <- get.conditional.mean(dropped.mean.centered.covariance,dropped.focal.pop,dropped.mean.centered.observations,drop.option=1)
	output <- list("mean.centered.observations" = mean.centered.observations,"focal.pop.conditional.mean" = focal.pop.conditional.mean)
	return(output)
}

get.fstat.location <- function(focal.pop,source.location,all.locations,covariance,observations,a0,aD,a2,mean.sample.sizes=NULL){
	# recover()
	prelims <- get.focal.conditional.mean.fstat.location(focal.pop,covariance,observations,mean.sample.sizes)
	focal.pop.conditional.mean <- prelims$focal.pop.conditional.mean
	mean.centered.observations <- prelims$mean.centered.observations
	spatial.covariance <- Covariance(a0,aD,a2,fields::rdist(rbind(all.locations,source.location)))
	mean.center.matrix <- get.uprank.mean.centering.matrix(focal.pop,k,mean.sample.sizes)
		mean.centered.spatial.covariance <- mean.center.matrix %*% spatial.covariance %*% t(mean.center.matrix)
	mean.centered.spatial.covariance <- mean.centered.spatial.covariance[-focal.pop,-focal.pop]
	source.location.conditional.mean <- get.conditional.mean(mean.centered.spatial.covariance,k,mean.centered.observations[-focal.pop,])
	fstat <- calculate.fstat(mean.centered.observations[focal.pop,1:ncol(mean.centered.observations),drop=FALSE],
								focal.pop.conditional.mean,
								t(source.location.conditional.mean))
	return(fstat)
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
							col = adjustcolor(rainbow(k),0.2) , lwd = admix.proportions[,i]*1)
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
		observed.locations <- cbind(last.params$spatial.prior.X.coordinates,last.params$spatial.prior.Y.coordinates)
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

if(FALSE){
get.spatial.fstat <- function(focal.pop,observations,source.pop,covariance=NULL,focal.pop.conditional.mean=NULL){
	if(is.null(focal.pop.conditional.mean)){
		if(is.null(covariance)){
			stop("you must specify a covariance matrix to calculate the conditional mean")
		}
		mean.centering.matrix <- get.mean.centering.matrix(focal.pop,nrow(observations))
		mean.centered.observations <- mean.centering.matrix %*% observations
		mean.centered.covariance <- mean.centering.matrix %*% covariance %*% t(mean.centering.matrix)
	focal.pop.conditional.mean <- get.conditional.mean(mean.centered.covariance,focal.pop,mean.centered.observations,drop.option=1)
	}
	fstat <- calculate.fstat(mean.centered.observations[focal.pop,],focal.pop.conditional.mean,mean.centered.observations[source.pop,])
	return(fstat)
}

get.all.spatial.fstat <- function(focal.pop,observations,covariance,focal.pop.conditional.mean){
	mean.centering.matrix <- get.mean.centering.matrix(focal.pop,nrow(observations))
	mean.centered.observations <- mean.centering.matrix %*% observations
	mean.centered.covariance <- mean.centering.matrix %*% covariance %*% t(mean.centering.matrix)
	focal.pop.conditional.mean <- get.conditional.mean(mean.centered.covariance,focal.pop,mean.centered.observations,drop.option=1)
	fstat <- numeric(nrow(observations)-1)
	for(i in 1:length(fstat)){
		fstat[i] <- calculate.fstat(mean.centered.observations[focal.pop,],focal.pop.conditional.mean,mean.centered.observations[i,])
	}
	return(fstat)
}

# # get.spatial.fstat <- function(sampled.locations,proposed.location,observations,focal.population,covariance=NULL,a0,aD,a2,focal.population.mean=NULL,get.focal.population.mean=NULL){
	# # recover()
	# k <- nrow(sampled.locations)
	# if(is.null(covariance)){
		# covariance <- Covariance(a0,aD,a2,fields::rdist(sampled.locations))
	# }
		# covariance.no.focal <- matrix(0,nrow=k,ncol=k)
		# covariance.no.focal[1:(k-1),1:(k-1)] <- covariance[-focal.population,-focal.population]
	# if(is.null(focal.population.mean) && is.null(get.focal.population.mean)){
		# stop("you must either specify the focal population conditional mean or you must estimate it")
	# }
	# if(!is.null(get.focal.population.mean)){
		# focal.population.mean <- get.conditional.mean(covariance,focal.population,observations,drop.option=TRUE)
	# }
	# coords <- rbind(sampled.locations[-focal.population,],proposed.location)
	# tmp.D <- fields::rdist(coords[k,1:2,drop=FALSE],coords)
	# tmp.cov <- Covariance(a0,aD,a2,tmp.D)
	# covariance.no.focal[k,] <- tmp.cov
	# covariance.no.focal[,k] <- tmp.cov
	# spatial.location.mean <- get.conditional.mean(covariance.no.focal,k,observations[-focal.population,])
	# spatial.fstat <- mean((observations[focal.population,]-focal.population.mean)*spatial.location.mean)
	# return(spatial.fstat)
# }
get.mean.centering.matrix2 <- function(focal.pop,k){
	mean.centering.matrix <- matrix(-1/(k-1),nrow=(k),ncol=k,byrow=TRUE)
		diag(mean.centering.matrix) <- (k-2)/(k-1)
		mean.centering.matrix[,focal.pop] <- c(rep(0,focal.pop-1),1,rep(0,k-focal.pop))
	return(mean.centering.matrix)		
}

get.fstat.pop2 <- function(focal.pop,source.pop,observations,covariance,mean.sample.sizes){
	k <- nrow(observations)
	mean.center.matrix <- get.mean.centering.matrix(focal.pop,k)
	mean.centered.covariance <- mean.center.matrix %*% covariance %*% t(mean.center.matrix)
	mean.centered.observations <- mean.center.matrix %*% observations
	pop.2.drop <- sample(c(1:k)[-c(focal.pop,source.pop)],1)
	if(pop.2.drop < focal.pop){
		focal.pop <- focal.pop - 1 
	}
	if(pop.2.drop < source.pop){
		source.pop <- source.pop - 1 
	}
	mean.centered.covariance <- mean.centered.covariance[-pop.2.drop,-pop.2.drop]
	mean.centered.observations <- mean.centered.observations[-pop.2.drop,]
	focal.pop.conditional.mean <- get.conditional.mean(mean.centered.covariance,focal.pop,mean.centered.observations,drop.option=1)
	fstat <- calculate.fstat(mean.centered.observations[focal.pop,],focal.pop.conditional.mean,mean.centered.observations[source.pop,])
	return(fstat)
}

}

query.MCMC.output <- function(MCMC.output,param.names=NULL,last.param.names=NULL){
	if(!is.null(param.names)){
		param.list <- vector("list",length=length(param.names))
		names(param.list) <- param.names
		load(MCMC.output)
			for(i in 1:length(param.names)){
				param.list[[param.names[i]]] <- get(param.names[i])
			}
	}
	if(!is.null(last.param.names)){
		last.param.list <- vector("list",length=length(last.param.names))
		names(last.param.list) <- last.param.names
		load(MCMC.output)
			for(i in 1:length(last.param.names)){
				last.param.list[[last.param.names[i]]] <- get(last.param.names[i],last.params)
			}
	}
	if(!is.null(param.names) & !is.null(last.param.names)){
		param.list <- c(param.list,last.param.list)
	}
	if(is.null(param.names) & !is.null(last.param.names)){
		param.list <- last.param.list
	}
	return(param.list)
}

run.spacemix.analysis <- function(n.fast.reps,
									fast.MCMC.ngen,
									fast.model.option,
									long.model.option,
									data.type,
									fast.likelihood.option,
									long.likelihood.option,
									proj.mat.option=NULL,
									sample.frequencies=NULL,
									mean.sample.sizes=NULL,
									counts=NULL,
									sample.sizes=NULL,
									sample.covariance=NULL,
									target.spatial.prior.scale=NULL,
									source.spatial.prior.scale=NULL,
									spatial.prior.X.coordinates,
									spatial.prior.Y.coordinates,
									round.earth,
									long.run.initial.parameters=NULL,
									k,
									loci,
									ngen,
									printfreq,
									samplefreq,
									mixing.diagn.freq,
									gibbs.nugget.fineness,
									gibbs.spatial.fineness,
									gibbs.step.frequency,
									savefreq,
									directory=NULL,
									prefix){
	if(n.fast.reps !=0){
		fast.run.dirs <- unlist(lapply(1:n.fast.reps,FUN=function(i){paste("fast_run_",i,sep="")}))
		for(i in 1:n.fast.reps){
			dir.create(fast.run.dirs[i])
			setwd(fast.run.dirs[i])
			random.initial.population.coordinates <- cbind(rep(0,2*k),rep(0,2*k))
#			random.initial.population.coordinates <- cbind( runif(2*k,min(spatial.prior.X.coordinates),
#																	max(spatial.prior.X.coordinates)),
#															runif(2*k,min(spatial.prior.Y.coordinates),
#																	max(spatial.prior.Y.coordinates)))
			initial.parameters <- list("population.coordinates" = random.initial.population.coordinates)
			tryCatch({
				MCMC(model.option = fast.model.option,
					data.type = data.type,
					likelihood.option = fast.likelihood.option,
					proj.mat.option = proj.mat.option,
					sample.frequencies = sample.frequencies,
					mean.sample.sizes = mean.sample.sizes,
					counts = counts,
					sample.sizes = sample.sizes,
					sample.covariance = sample.covariance,
					target.spatial.prior.scale = target.spatial.prior.scale,
					source.spatial.prior.scale = source.spatial.prior.scale,
					spatial.prior.X.coordinates = spatial.prior.X.coordinates,
					spatial.prior.Y.coordinates = spatial.prior.Y.coordinates,
					round.earth = round.earth,
					initial.parameters = initial.parameters,
					k = k,
					loci = loci,
					ngen = fast.MCMC.ngen,
					printfreq = 1e3,
					samplefreq = fast.MCMC.ngen/1e3,
					mixing.diagn.freq = mixing.diagn.freq,
					gibbs.nugget.fineness = gibbs.nugget.fineness,
					gibbs.spatial.fineness = gibbs.spatial.fineness,
					gibbs.step.frequency = gibbs.step.frequency,
					savefreq = fast.MCMC.ngen/2,
					directory = directory,
					prefix=paste(fast.run.dirs[i],"_",sep=""),
					continue=FALSE,
					continuing.params=NULL)
			},error=function(e){
							cat("fast run",i,"failed","\n")
							file.rename(
									paste("../",fast.run.dirs[i],sep=""),
									paste("../","failed_",fast.run.dirs[i],sep=""))
							})
			setwd("..")
		}
		last.probs <- c(rep(-Inf,n.fast.reps))
		for(i in 1:n.fast.reps){
			tryCatch({
				setwd(fast.run.dirs[i])
				prob.list <- query.MCMC.output(list.files()[grep("output",list.files())],param.names="Prob")
					last.probs[i] <- prob.list$Prob[length(which(prob.list$Prob!=0))]
				setwd("..")
				},error=function(e){cat("skipping",fast.run.dirs[i],"because it failed","\n")})
		}
		file.rename(fast.run.dirs[which.max(last.probs)],paste(fast.run.dirs[which.max(last.probs)],"_BestRun",sep=""))
		setwd(list.files()[grep("BestRun",list.files())])
			fast.run.initial.parameters <- query.MCMC.output(list.files()[grep("output",list.files())],last.param.names=c("a0","aD","a2","population.coordinates","nugget","admix.proportions"))
		setwd("..")
		dir.create(paste(prefix,"LongRun",sep="_"))
		setwd(list.files()[grep("LongRun",list.files())])
	} else {
		fast.run.initial.parameters <- long.run.initial.parameters
	}
			MCMC(model.option = long.model.option,
				data.type = data.type,
				likelihood.option = long.likelihood.option,
				proj.mat.option = proj.mat.option,
				sample.frequencies = sample.frequencies,
				mean.sample.sizes = mean.sample.sizes,
				counts = counts,
				sample.sizes = sample.sizes,
				sample.covariance = sample.covariance,
				target.spatial.prior.scale = target.spatial.prior.scale,
				source.spatial.prior.scale = source.spatial.prior.scale,
				spatial.prior.X.coordinates = spatial.prior.X.coordinates,
				spatial.prior.Y.coordinates = spatial.prior.Y.coordinates,
				initial.parameters = fast.run.initial.parameters,
				round.earth = round.earth,
				k = k,
				loci = loci,
				ngen = ngen,
				printfreq = printfreq,
				samplefreq = samplefreq,
				mixing.diagn.freq = mixing.diagn.freq,
				gibbs.nugget.fineness = gibbs.nugget.fineness, 
				gibbs.spatial.fineness = gibbs.spatial.fineness,
				gibbs.step.frequency = gibbs.step.frequency,
				savefreq = savefreq,
				directory = directory,
                prefix=prefix,
				continue=FALSE,
				continuing.params=NULL)
	setwd("..")
	return("analysis completed.")
}
