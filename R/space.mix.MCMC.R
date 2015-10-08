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
		tmp.covariance_prime <- Covariance(last.params$a0,last.params$a1,last.params$a2,D_prime)
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
		tmp.covariance_prime <- Covariance(last.params$a0,last.params$a1,last.params$a2,D_prime)
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

Update_a1 <- function(last.params){
	new.params <- last.params
	a1_prime <- last.params$a1 + rnorm(1,0,exp(last.params$lstps$a1_lstp))
	prior_prob_alpha1_prime <- Prior_prob_alpha1(a1_prime)
	if(prior_prob_alpha1_prime != -Inf){
		covariance_prime <- Covariance(last.params$a0,a1_prime,last.params$a2,last.params$D)
			admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions,last.params$nugget,last.params$k,last.params$inv.mean.sample.sizes,last.params$identity.matrix)
			transformed.covariance.prime <- transformed.Covariance(admixed.covariance_prime,last.params$projection.matrix)
				LnL_freqs_prime <- likelihood(likelihood.option = last.params$likelihood.option,sample.cov = last.params$sample.cov,par.cov = transformed.covariance.prime,index.matrix = last.params$index.matrix,sd = last.params$sd,n = last.params$loci)
				if(metropolis_ratio(LnL_freqs_prime,prior_prob_alpha1_prime,last.params$LnL_freqs,last.params$prior_prob_alpha1)){
					new.params$a1 <- a1_prime
					new.params$covariance <- covariance_prime
					new.params$admixed.covariance <- admixed.covariance_prime
					new.params$transformed.covariance <- transformed.covariance.prime					
					new.params$prior_prob_alpha1 <- prior_prob_alpha1_prime
					new.params$LnL_freqs <- LnL_freqs_prime
					new.params$accept$a1_accept <- new.params$accept$a1_accept + 1
				}
	}
	new.params$moves$a1_moves <- new.params$moves$a1_moves + 1
	new.params$accept_rates$a1_accept_rate <- new.params$accept$a1_accept/new.params$moves$a1_moves
	return(new.params)
}

Update_a2 <- function(last.params){
		new.params <- last.params
		a2_prime <- last.params$a2 + rnorm(1,0,exp(last.params$lstps$a2_lstp))
		prior_prob_alpha2_prime <- Prior_prob_alpha2(a2_prime) 
		if(prior_prob_alpha2_prime != -Inf) {
			covariance_prime <- Covariance(last.params$a0,last.params$a1,a2_prime,last.params$D)
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

Prior_prob_alpha1 <- function(a1){
	log(dexp(a1,rate=1))
}

Prior_prob_alpha2 <- function(a2){
	log(dunif(a2,0.1,2))
}

Prior_prob_nugget <- function(nugget){
	sum(dexp(nugget,log=TRUE))
}

Prior_prob_admix_proportions <- function(admix_proportions){
	admix_proportions <- admix_proportions * 2
	sum(dbeta(admix_proportions,shape1=1,shape2=100,log=TRUE))
}

normal.lnL <- function(sample.cov,par.cov,sd,index.matrix){
	par.cov <- par.cov[index.matrix]
	lnL <- - (sample.cov-par.cov)^2/(2*sd^2)
	return(sum(lnL))
}

wishart.lnL <- function(sample.cov,par.cov,n){
	A <- solve(par.cov)
	lnL <- -0.5 * n * sum( A * sample.cov ) - (n/2)*determinant(par.cov,logarithm=TRUE)$modulus
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
	
Covariance <- function(a0,a1,a2,GeoDist) {
	covariance <- (1/a0)*exp(-(a1*GeoDist)^a2)
	return(covariance)
}

admixed.Covariance <- function(covariance,admix.proportions,nugget,k,inv.mean.sample.sizes,ident.mat){
	# recover()
	if(any(admix.proportions !=0)){
		w_k <- admix.proportions
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
	save(curated.count.data,file=paste(prefix,"_curated.count.data.Robj",sep=''))
	return(curated.count.data)
}	

curate.frequency.data <- function(sample.frequencies){
	k <- nrow(sample.frequencies)
	fixed.alleles <- c(which(colSums(sample.frequencies)==0),which(colSums(sample.frequencies)==k))
	na.loci <- which(is.na(sample.frequencies),arr.ind=TRUE)[,2]
	missing.data <- which(!is.finite(sample.frequencies),arr.ind=TRUE)[,2]
	loci.to.drop <- unique(c(fixed.alleles,na.loci,missing.data))
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
			save(MCN.frequencies.list,file=paste(prefix,"_MCN.frequencies.list.Robj",sep=''))
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
			save(MCN.frequencies.list,file=paste(prefix,"_MCN.frequencies.list.Robj",sep=''))
	return(mean.centered.normalized.sample.frequencies)
}

get.sample.covariance <- function(counts,sample.sizes){
	sample.frequencies <- counts/sample.sizes
	mean.sample.sizes <- rowMeans(sample.sizes)
	mean.sample.frequencies <- matrix(apply(sample.frequencies,2,
											get.weighted.mean.frequency,
											mean.sample.sizes=mean.sample.sizes),
									nrow=length(mean.sample.sizes),ncol=ncol(sample.frequencies),byrow=TRUE)
	normalized.sample.frequencies <- sample.frequencies/sqrt(mean.sample.frequencies*(1-mean.sample.frequencies))
	sample.covariance <- cov(t(normalized.sample.frequencies),use="pairwise.complete.obs")
	loci <- ncol(sample.frequencies)
	return(list("sample.covariance" = sample.covariance,
				"mean.sample.sizes" = mean.sample.sizes,
				"loci" = loci))
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
		save(spacemix.data,file=paste(prefix,"_spacemix.data.Robj",sep=''))
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
		admix.proportions <- matrix(rbeta(k,shape1=1,shape2=50),nrow=k,ncol=1) / 2
		prior_prob_admix_proportions <- Prior_prob_admix_proportions(admix.proportions)
	} else if(model.option == "source_and_target"){
		admix.proportions <- matrix(rbeta(k,shape1=1,shape2=50),nrow=k,ncol=1) / 2
		prior_prob_admix_proportions <- Prior_prob_admix_proportions(admix.proportions)
	}
	initiate.admix.proportions.list <- list("admix.proportions" = admix.proportions, 
											"prior_prob_admix_proportions" = prior_prob_admix_proportions)
	return(initiate.admix.proportions.list)
}

initiate.population.coordinates <- function(spatial.prior.X.coordinates,spatial.prior.Y.coordinates,k){
	# recover()
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

save.initial.parameters <- function(a0,a1,a2,nugget,admix.proportions,covariance,admixed.covariance,transformed.covariance,population.coordinates,D,projection.matrix,inv.mean.sample.sizes,prefix){
	initial.parameters <- list("a0" = a0,"a1" = a1,"a2" = a2,"nugget" = nugget,"admix.proportions" = admix.proportions,
								"covariance" = covariance,"admixed.covariance" = admixed.covariance,"transformed.covariance" = transformed.covariance,
								"population.coordinates" = population.coordinates,"D" = D,"projection.matrix" = projection.matrix,"inv.mean.sample.sizes"=inv.mean.sample.sizes)
	save(initial.parameters,file=paste(prefix,"_Initial.parameters.Robj",sep=''))
	return(0)
}

initiate.update.function.list <- function(model.option){
	if(model.option == "no_movement"){
		Updates <- list(Update_a0,Update_a1,Update_a2,Update_nugget)
	} else if(model.option == "target"){
		Updates <- list(Update_a0,Update_a1,Update_a2,Update_nugget,Update_admixture_target_location)
	} else if(model.option == "source"){
		Updates <- list(Update_a0,Update_a1,Update_a2,Update_nugget,Update_admixture_source_location,Update_admixture_proportions)
	} else if(model.option == "source_and_target"){
		Updates <- list(Update_a0,Update_a1,Update_a2,Update_nugget,Update_admixture_target_location,Update_admixture_source_location,Update_admixture_proportions)
	}
	return(Updates)
}

initiate.lstps.list <- function(model.option,k,ngen,samplefreq,for.last.params){
	if(!for.last.params){
		if(model.option == "no_movement"){
			lstps <- list("a0_lstp" = numeric(ngen/samplefreq),
							"a1_lstp" = numeric(ngen/samplefreq),
							"a2_lstp" = numeric(ngen/samplefreq),
							"nugget_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		} else if(model.option == "target"){
			lstps <- list("a0_lstp" = numeric(ngen/samplefreq),
							"a1_lstp" = numeric(ngen/samplefreq),
							"a2_lstp" = numeric(ngen/samplefreq),
							"nugget_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_target_location_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		} else if(model.option == "source"){
			lstps <- list("a0_lstp" = numeric(ngen/samplefreq),
							"a1_lstp" = numeric(ngen/samplefreq),
							"a2_lstp" = numeric(ngen/samplefreq),
							"nugget_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_source_location_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_proportions_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		} else if(model.option == "source_and_target"){
			lstps <- list("a0_lstp" = numeric(ngen/samplefreq),
							"a1_lstp" = numeric(ngen/samplefreq),
							"a2_lstp" = numeric(ngen/samplefreq),
							"nugget_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_target_location_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_source_location_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_proportions_lstp" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		}
	} else {
		if(model.option == "no_movement"){
			lstps <- list("a0_lstp" = 0,
							"a1_lstp" = 0,
							"a2_lstp" = 0,
							"nugget_lstp" = numeric(k))
		} else if(model.option == "target"){
			lstps <- list("a0_lstp" = 0,
							"a1_lstp" = 0,
							"a2_lstp" = 0,
							"nugget_lstp" = numeric(k),
							"admix_target_location_lstp" = numeric(k))
		} else if(model.option == "source"){
			lstps <- list("a0_lstp" = 0,
							"a1_lstp" = 0,
							"a2_lstp" = 0,
							"nugget_lstp" = numeric(k),
							"admix_source_location_lstp" = numeric(k),
							"admix_proportions_lstp" = numeric(k))
		} else if(model.option == "source_and_target"){
			lstps <- list("a0_lstp" = 0,
							"a1_lstp" = 0,
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
						"a1_diagn" = numeric(mixing.diagn.freq),
						"a2_diagn" = numeric(mixing.diagn.freq),
						"nugget_diagn" = matrix(0,nrow=k,ncol=mixing.diagn.freq))
	} else if(model.option == "target"){
		diagns <- list("a0_diagn" = numeric(mixing.diagn.freq),
						"a1_diagn" = numeric(mixing.diagn.freq),
						"a2_diagn" = numeric(mixing.diagn.freq),
						"nugget_diagn" = matrix(0,nrow=k,ncol=(mixing.diagn.freq)),
						"admix_target_location_diagn" = matrix(0,nrow=k,ncol=mixing.diagn.freq))
	} else if(model.option == "source"){
		diagns <- list("a0_diagn" = numeric(mixing.diagn.freq),
						"a1_diagn" = numeric(mixing.diagn.freq),
						"a2_diagn" = numeric(mixing.diagn.freq),
						"nugget_diagn" = matrix(0,nrow=k,ncol=mixing.diagn.freq),
						"admix_source_location_diagn" = matrix(0,nrow=k,ncol=mixing.diagn.freq),
						"admix_proportions_diagn" = matrix(0,nrow=k,ncol=mixing.diagn.freq))
	} else if(model.option == "source_and_target"){
		diagns <- list("a0_diagn" = numeric(mixing.diagn.freq),
						"a1_diagn" = numeric(mixing.diagn.freq),
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
						"a1_moves" = 0,
						"a2_moves" = 0,
						"nugget_moves" = numeric(k))
	} else if(model.option == "target"){
		moves <- list("a0_moves" = 0,
						"a1_moves" = 0,
						"a2_moves" = 0,
						"nugget_moves" = numeric(k),
						"admix_target_location_moves" = numeric(k))
	} else if(model.option == "source"){
		moves <- list("a0_moves" = 0,
						"a1_moves" = 0,
						"a2_moves" = 0,
						"nugget_moves" = numeric(k),
						"admix_source_location_moves" = numeric(k),
						"admix_proportions_moves" = numeric(k))
	} else if(model.option == "source_and_target"){
		moves <- list("a0_moves" = 0,
						"a1_moves" = 0,
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
						"a1_accept" = 0,
						"a2_accept" = 0,
						"nugget_accept" = numeric(k))
	} else if(model.option == "target"){
		accept <- list("a0_accept" = 0,
						"a1_accept" = 0,
						"a2_accept" = 0,
						"nugget_accept" = numeric(k),
						"admix_target_location_accept" = numeric(k))
	} else if(model.option == "source"){
		accept <- list("a0_accept" = 0,
						"a1_accept" = 0,
						"a2_accept" = 0,
						"nugget_accept" = numeric(k),
						"admix_source_location_accept" = numeric(k),
						"admix_proportions_accept" = numeric(k))
	} else if(model.option == "source_and_target"){
		accept <- list("a0_accept" = 0,
						"a1_accept" = 0,
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
							"a1_accept_rate" = numeric(ngen/samplefreq),
							"a2_accept_rate" = numeric(ngen/samplefreq),
							"nugget_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		} else if(model.option == "target"){
			accept_rates <- list("a0_accept_rate" = numeric(ngen/samplefreq),
							"a1_accept_rate" = numeric(ngen/samplefreq),
							"a2_accept_rate" = numeric(ngen/samplefreq),
							"nugget_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_target_location_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		} else if(model.option == "source"){
			accept_rates <- list("a0_accept_rate" = numeric(ngen/samplefreq),
							"a1_accept_rate" = numeric(ngen/samplefreq),
							"a2_accept_rate" = numeric(ngen/samplefreq),
							"nugget_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_source_location_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_proportions_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		} else if(model.option == "source_and_target"){
			accept_rates <- list("a0_accept_rate" = numeric(ngen/samplefreq),
							"a1_accept_rate" = numeric(ngen/samplefreq),
							"a2_accept_rate" = numeric(ngen/samplefreq),
							"nugget_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_target_location_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_source_location_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)),
							"admix_proportions_accept_rate" = matrix(0,nrow=k,ncol=(ngen/samplefreq)))
		}
	} else {
		if(model.option == "no_movement"){
			accept_rates <- list("a0_accept_rate" = 0,
							"a1_accept_rate" = 0,
							"a2_accept_rate" = 0,
							"nugget_accept_rate" = numeric(k))
		} else if(model.option == "target"){
			accept_rates <- list("a0_accept_rate" = 0,
							"a1_accept_rate" = 0,
							"a2_accept_rate" = 0,
							"nugget_accept_rate" = numeric(k),
							"admix_target_location_accept_rate" = numeric(k))
		} else if(model.option == "source"){
			accept_rates <- list("a0_accept_rate" = 0,
							"a1_accept_rate" = 0,
							"a2_accept_rate" = 0,
							"nugget_accept_rate" = numeric(k),
							"admix_source_location_accept_rate" = numeric(k),
							"admix_proportions_accept_rate" = numeric(k))
		} else if(model.option == "source_and_target"){
			accept_rates <- list("a0_accept_rate" = 0,
							"a1_accept_rate" = 0,
							"a2_accept_rate" = 0,
							"nugget_accept_rate" = numeric(k),
							"admix_target_location_accept_rate" = numeric(k),
							"admix_source_location_accept_rate" = numeric(k),
							"admix_proportions_accept_rate" = numeric(k))
		}
	}
	return(accept_rates)
}

print.mcmc.update <- function(LnL_freqs,prior_prob_admix_proportions,prior_prob_nugget,prior_prob_alpha0,prior_prob_alpha1,prior_prob_alpha2,prior_prob_admix_target_locations,prior_prob_admix_source_locations){
	P <- LnL_freqs + prior_prob_admix_proportions + prior_prob_nugget + prior_prob_alpha0 + prior_prob_alpha1 + prior_prob_alpha2 + prior_prob_admix_target_locations + prior_prob_admix_source_locations
	return(P)
}

make.update.sampled.accept.rates.function <- function(model.option){
	if(model.option == "no_movement"){
		update.sampled.accept.rates <<- function(accept_rates,j,new.params){
			accept_rates$a0_accept_rate[j] <- new.params$accept_rates$a0_accept_rate
			accept_rates$a1_accept_rate[j] <- new.params$accept_rates$a1_accept_rate
			accept_rates$a2_accept_rate[j] <- new.params$accept_rates$a2_accept_rate
			accept_rates$nugget_accept_rate[,j] <- new.params$accept_rates$nugget_accept_rate
			return(accept_rates)
		}
	} else if(model.option == "target"){
		update.sampled.accept.rates <<- function(accept_rates,j,new.params){
			accept_rates$a0_accept_rate[j] <- new.params$accept_rates$a0_accept_rate
			accept_rates$a1_accept_rate[j] <- new.params$accept_rates$a1_accept_rate
			accept_rates$a2_accept_rate[j] <- new.params$accept_rates$a2_accept_rate
			accept_rates$nugget_accept_rate[,j] <- new.params$accept_rates$nugget_accept_rate
			accept_rates$admix_target_location_accept_rate[,j] <- new.params$accept_rates$admix_target_location_accept_rate
			return(accept_rates)
		}
	} else if(model.option == "source"){
		update.sampled.accept.rates <<- function(accept_rates,j,new.params){
			accept_rates$a0_accept_rate[j] <- new.params$accept_rates$a0_accept_rate
			accept_rates$a1_accept_rate[j] <- new.params$accept_rates$a1_accept_rate
			accept_rates$a2_accept_rate[j] <- new.params$accept_rates$a2_accept_rate
			accept_rates$nugget_accept_rate[,j] <- new.params$accept_rates$nugget_accept_rate
			accept_rates$admix_source_location_accept_rate[,j] <- new.params$accept_rates$admix_source_location_accept_rate
			accept_rates$admix_proportions_accept_rate[,j] <- new.params$accept_rates$admix_proportions_accept_rate
			return(accept_rates)
		}
	} else if(model.option == "source_and_target"){
		update.sampled.accept.rates <<- function(accept_rates,j,new.params){
			accept_rates$a0_accept_rate[j] <- new.params$accept_rates$a0_accept_rate
			accept_rates$a1_accept_rate[j] <- new.params$accept_rates$a1_accept_rate
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
			lstps$a1_lstp[j] <- new.params$lstps$a1_lstp
			lstps$a2_lstp[j] <- new.params$lstps$a2_lstp
			lstps$nugget_lstp[,j] <- new.params$lstps$nugget_lstp
			return(lstps)
		}
	} else if(model.option == "target"){
		update.sampled.lstps <<- function(lstps,j,new.params){
			lstps$a0_lstp[j] <- new.params$lstps$a0_lstp
			lstps$a1_lstp[j] <- new.params$lstps$a1_lstp
			lstps$a2_lstp[j] <- new.params$lstps$a2_lstp
			lstps$nugget_lstp[,j] <- new.params$lstps$nugget_lstp
			lstps$admix_target_location_lstp[,j] <- new.params$lstps$admix_target_location_lstp
			return(lstps)
		}	
	} else if(model.option == "source"){
		update.sampled.lstps <<- function(lstps,j,new.params){
			lstps$a0_lstp[j] <- new.params$lstps$a0_lstp
			lstps$a1_lstp[j] <- new.params$lstps$a1_lstp
			lstps$a2_lstp[j] <- new.params$lstps$a2_lstp
			lstps$nugget_lstp[,j] <- new.params$lstps$nugget_lstp
			lstps$admix_source_location_lstp[,j] <- new.params$lstps$admix_source_location_lstp
			lstps$admix_proportions_lstp[,j] <- new.params$lstps$admix_proportions_lstp
			return(lstps)
		}
	} else if(model.option == "source_and_target"){
		update.sampled.lstps <<- function(lstps,j,new.params){
			lstps$a0_lstp[j] <- new.params$lstps$a0_lstp
			lstps$a1_lstp[j] <- new.params$lstps$a1_lstp
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
				diagns$a1_diagn[diagn.step] <- accept_rates$a1_accept_rate
				diagns$a2_diagn[diagn.step] <- accept_rates$a2_accept_rate
				diagns$nugget_diagn[,diagn.step] <- accept_rates$nugget_accept_rate
			return(diagns)
		}
	} else if(model.option == "target"){
		update.diagns <<- function(diagns,generation,mixing.diagn.freq,accept_rates){
			diagn.step <- get.diagn.step(generation,mixing.diagn.freq)
				diagns$a0_diagn[diagn.step] <- accept_rates$a0_accept_rate
				diagns$a1_diagn[diagn.step] <- accept_rates$a1_accept_rate
				diagns$a2_diagn[diagn.step] <- accept_rates$a2_accept_rate
				diagns$nugget_diagn[,diagn.step] <- accept_rates$nugget_accept_rate
				diagns$admix_target_location_diagn[,diagn.step] <- accept_rates$admix_target_location_accept_rate
			return(diagns)
		}
	} else if(model.option == "source"){
		update.diagns <<- function(diagns,generation,mixing.diagn.freq,accept_rates){
			diagn.step <- get.diagn.step(generation,mixing.diagn.freq)
				diagns$a0_diagn[diagn.step] <- accept_rates$a0_accept_rate
				diagns$a1_diagn[diagn.step] <- accept_rates$a1_accept_rate
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
				diagns$a1_diagn[diagn.step] <- accept_rates$a1_accept_rate
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
			lstps$a1_lstp <- update.lstp(n,lstps$a1_lstp,mean(diagns$a1_diagn))
			lstps$a2_lstp <- update.lstp(n,lstps$a2_lstp,mean(diagns$a2_diagn))
			lstps$nugget_lstp <- update.lstp(n,lstps$nugget_lstp,rowMeans(diagns$nugget_diagn))
			return(lstps)
		}
	} else if(model.option == "target"){
		update.all.lstps <<- function(lstps,n,diagns){
			lstps$a0_lstp <- update.lstp(n,lstps$a0_lstp,mean(diagns$a0_diagn))
			lstps$a1_lstp <- update.lstp(n,lstps$a1_lstp,mean(diagns$a1_diagn))
			lstps$a2_lstp <- update.lstp(n,lstps$a2_lstp,mean(diagns$a2_diagn))
			lstps$nugget_lstp <- update.lstp(n,lstps$nugget_lstp,rowMeans(diagns$nugget_diagn))
			lstps$admix_target_location_lstp <- update.lstp(n,lstps$admix_target_location_lstp,rowMeans(diagns$admix_target_location_diagn))
			return(lstps)
		}
	} else if(model.option == "source"){
		update.all.lstps <<- function(lstps,n,diagns){
			lstps$a0_lstp <- update.lstp(n,lstps$a0_lstp,mean(diagns$a0_diagn))
			lstps$a1_lstp <- update.lstp(n,lstps$a1_lstp,mean(diagns$a1_diagn))
			lstps$a2_lstp <- update.lstp(n,lstps$a2_lstp,mean(diagns$a2_diagn))
			lstps$nugget_lstp <- update.lstp(n,lstps$nugget_lstp,rowMeans(diagns$nugget_diagn))
			lstps$admix_source_location_lstp <- update.lstp(n,lstps$admix_source_location_lstp,rowMeans(diagns$admix_source_location_diagn))
			lstps$admix_proportions_lstp <- update.lstp(n,lstps$admix_proportions_lstp,rowMeans(diagns$admix_proportions_diagn))
			return(lstps)
		}
	} else if(model.option == "source_and_target"){
		update.all.lstps <<- function(lstps,n,diagns){
			lstps$a0_lstp <- update.lstp(n,lstps$a0_lstp,mean(diagns$a0_diagn))
			lstps$a1_lstp <- update.lstp(n,lstps$a1_lstp,mean(diagns$a1_diagn))
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

initiate.last.params <- function(model.option,likelihood.option,samplefreq,ngen,spacemix.data,population.coordinates,admix.proportions,a0,a1,a2,nugget,covariance,admixed.covariance,transformed.covariance,
						k,LnL_freqs,prior_prob_alpha0,prior_prob_alpha1,prior_prob_alpha2,prior_prob_nugget,prior_prob_admix_proportions,prior_prob_admix_target_locations,prior_prob_admix_source_locations,
						D,spatial.prior.X.coordinates,spatial.prior.Y.coordinates,target.spatial.prior.scale,source.spatial.prior.scale,centroid){
	last.params <- list("sample.covariance" = spacemix.data$sample.covariance,
						"projection.matrix" = spacemix.data$projection.matrix,						
						"population.coordinates" = population.coordinates,
						"admix.proportions" = admix.proportions,
						"a0" = a0,"a1" = a1,"a2" = a2,"nugget" = nugget,
						"covariance" = covariance,"admixed.covariance" = admixed.covariance,
						"transformed.covariance" = transformed.covariance,
						"lstps" = initiate.lstps.list(model.option,k,ngen,samplefreq,for.last.params=TRUE),
						"k" = k,"LnL_freqs" = LnL_freqs,
						"prior_prob_alpha0" = prior_prob_alpha0,
						"prior_prob_alpha1" = prior_prob_alpha1,
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
						"sd" = spacemix.data$sd,
						"index.matrix" = spacemix.data$index.matrix,
						"likelihood.option" = likelihood.option,
						"identity.matrix" = diag(k))
	return(last.params)
}
	
	
MCMC <- function(model.option,
				data.type,
				sample.frequencies = NULL,
				mean.sample.sizes = NULL,
				counts = NULL,
				sample.sizes = NULL,
				sample.covariance = NULL,
				target.spatial.prior.scale = NULL,
				source.spatial.prior.scale = NULL,
				spatial.prior.X.coordinates,
				spatial.prior.Y.coordinates,
				initial.parameters = NULL,
				round.earth,
				k,
				loci,
				ngen,
				printfreq,
				samplefreq,
				mixing.diagn.freq = 50,
				savefreq,
				directory=NULL,
                prefix=""){
	# recover()
	if(!is.null(directory)){
		setwd(directory)
	}
	proj.mat.option <- NULL
	likelihood.option <- "wishart"
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
		a1 <- numeric(ngen/samplefreq)
		a2 <- numeric(ngen/samplefreq)
		accept_rates <- initiate.accept.rates.list(model.option,k,ngen,samplefreq,for.last.params=FALSE)
		lstps <- initiate.lstps.list(model.option,k,ngen,samplefreq,for.last.params=FALSE)
		diagns <- initiate.diagns.list(model.option,k,mixing.diagn.freq)
		
		seed <- sample(0:1e6,1)
			save(seed,file=paste(prefix,"_seed.Robj",sep=''))
		set.seed(seed)
		
	# if(!continue) {
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
					if(!is.null(initial.parameters$a1)){
						a1[1] <- initial.parameters$a1
					} else { a1[1] <- rexp(1,1) }
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
					covariance <- Covariance(a0[1],a1[1],a2[1],distances[[1]])
					admixed.covariance <- admixed.Covariance(covariance,admix.proportions[,1],nugget[,1],k,spacemix.data$inv.mean.sample.sizes,diag(k))
					transformed.covariance <- transformed.Covariance(admixed.covariance,spacemix.data$projection.matrix)
					tmp <- save.initial.parameters(a0[1],a1[1],a2[1],nugget[,1],admix.proportions[,1],covariance,admixed.covariance,transformed.covariance,population.coordinates[[1]],distances[[1]],spacemix.data$projection.matrix,spacemix.data$inv.mean.sample.sizes,prefix)
				LnL_freqs[1] <- likelihood(likelihood.option,spacemix.data$sample.covariance, transformed.covariance,spacemix.data$index.matrix,spacemix.data$sd,spacemix.data$loci)
					cat("LnL: ",LnL_freqs[1],"\n")
				prior_prob_alpha0 <- Prior_prob_alpha0(a0[1])
					cat("Pr(a0): ",prior_prob_alpha0,"\n")
				prior_prob_alpha1 <- Prior_prob_alpha1(a1[1])
					cat("Pr(a1): ",prior_prob_alpha1,"\n")
				prior_prob_alpha2 <- Prior_prob_alpha2(a2[1])
					cat("Pr(a2): ",prior_prob_alpha2,"\n")
				prior_prob_nugget <- Prior_prob_nugget(nugget[,1])
					cat("Pr(nugget): ",prior_prob_nugget,"\n")
				prior_prob_admix_target_locations <- Prior_prob_admix_target_locations(population.coordinates[[1]][1:k,],spatial.prior.X.coordinates,spatial.prior.Y.coordinates,target.spatial.prior.scale)
					cat("Pr(admix_target_locations): ",prior_prob_admix_target_locations,"\n")
				prior_prob_admix_source_locations <- Prior_prob_admix_source_locations(population.coordinates[[1]][(k+1):(2*k),],centroid,source.spatial.prior.scale)
					cat("Pr(admix_source_locations): ",prior_prob_admix_source_locations,"\n")
					cat("Pr(admix_proportions): ",prior_prob_admix_proportions,"\n")
				Prob[1] <- LnL_freqs[1] + prior_prob_admix_proportions + prior_prob_nugget + prior_prob_alpha0 + prior_prob_alpha1 + prior_prob_alpha2 + prior_prob_admix_target_locations + prior_prob_admix_source_locations
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
											a0[1],a1[1],a2[1],nugget[,1],covariance,admixed.covariance,transformed.covariance,k = k,LnL_freqs = LnL_freqs[1],
											prior_prob_alpha0 = prior_prob_alpha0,prior_prob_alpha1 = prior_prob_alpha1,prior_prob_alpha2 = prior_prob_alpha2,
											prior_prob_nugget = prior_prob_nugget,prior_prob_admix_proportions = prior_prob_admix_proportions,
											prior_prob_admix_target_locations = prior_prob_admix_target_locations,prior_prob_admix_source_locations = prior_prob_admix_source_locations,
											D = distances[[1]],spatial.prior.X.coordinates = spatial.prior.X.coordinates,spatial.prior.Y.coordinates = spatial.prior.Y.coordinates,
											target.spatial.prior.scale = target.spatial.prior.scale,source.spatial.prior.scale = source.spatial.prior.scale,
											centroid = centroid)
		last.ngen <- 0

	tmp <- make.update.sampled.accept.rates.function(model.option)
	tmp <- make.update.sampled.lstps.function(model.option)
	tmp <- make.update.diagns.function(model.option)
	tmp <- make.update.all.lstps.function(model.option)
	
	#Run the MCMC
		Updates <- initiate.update.function.list(model.option)
	
	for(i in 2:ngen) {
		x <- sample(c(1:length(Updates)),1)
		new.params <- Updates[[x]](last.params)

		if(i%%samplefreq == 0){
			j <- i/samplefreq
			population.coordinates[[j]] <- new.params$population.coordinates
			distances[[j]] <- new.params$D
			transformed.covariance.list[[j]] <- new.params$transformed.covariance
			admix.proportions[,j] <- new.params$admix.proportions
			nugget[,j] <- new.params$nugget
			a0[j] <- new.params$a0
			a1[j] <- new.params$a1
			a2[j] <- new.params$a2
			LnL_freqs[j] <- new.params$LnL_freqs
			Prob[j] <- LnL_freqs[j] +
						new.params$prior_prob_admix_proportions +
						new.params$prior_prob_nugget +
						new.params$prior_prob_alpha0 +
						new.params$prior_prob_alpha1 +
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
									new.params$prior_prob_nugget,new.params$prior_prob_alpha0,new.params$prior_prob_alpha1,
									new.params$prior_prob_alpha2,new.params$prior_prob_admix_target_locations,new.params$prior_prob_admix_source_locations)
			cat(i," ---- ",P,"\n")
		}
				
		if(i%%savefreq == 0){	
			save(last.params,
				LnL_freqs,Prob,distances,population.coordinates,
				transformed.covariance.list,admix.proportions,a0,a1,a2,nugget,samplefreq,ngen,
				accept_rates,lstps,diagns,
				target.spatial.prior.scale,source.spatial.prior.scale,
				file=paste(prefix,sprintf("_space_MCMC_output%d.Robj",1),sep=''))
		}
	}
    return(paste("Output",i,"runs to",paste(prefix,"_MCMC_output*.Robj",sep=''),"."))
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

#' Runs a SpaceMix analysis
#'
#' This function runs a Markov chain Monte Carlo to estimate a geogenetic map of your genotyped samples.
#' A SpaceMix analysis can be run with any of 4 separate models:
#' \enumerate{
#' \item "no_movement" - populations do not choose their own locations, nor can they draw admixture.
#' The only parameters to be estimated are: the alpha parameters of the spatial 
#' covariance function and the nugget parameters
#' \item "source" - populations do not choose their own locations, but they do draw admixture.
#' The parameters to be estimated are: the alpha parameters of the spatial 
#' covariance function, the nugget parameters, the locations of the sources of admixture,
#' the strength of that admixture.
#' \item "target" - populations choose their own locations, but no admixture.
#' The parameters to be estimated are: the alpha parameters of the spatial 
#' covariance function, the nugget parameters, the population locations.
#' \item "source_and_target" - populations choose their own locations AND they draw admixture.
#' The parameters to be estimated are: the alpha parameters of the spatial 
#' covariance function, the nugget parameters, the population locations,
#' the locations of the sources of admixture, and the strength of that admixture.
#' }
#'
#' The algorithm proceeds by running a user-specified number of short initial runs
#' from random locations in parameter space to find a generally good area, then one
#' long run from the final location in parameter space from the best short run,
#' the results of which are what the user cares about.  The user can also choose
#' to run only a single long analysis, for which the initial parameters may be specified.
#' 
#' @param n.fast.reps The number of short initial runs to perform.
#' @param fast.MCMC.ngen The number of generations to run each initial MCMC analysis.
#' @param fast.model.option The model to be used in the short runs:  may be "no_movement","source","target","source_and_target".
#' @param long.model.option The model to be used in the long run:  may be "no_movement","source","target","source_and_target".
#' @param data.type The data type to be used.  May be "sample.covariance","sample.frequencies","counts".
#'			Please see the vignette for a discussion of what these different data elements should look like.
#' @param sample.frequencies Data to be specified if "sample.frequencies" is chosen as data.type.
#' @param mean.sample.sizes Data to be specified if "sample.frequencies" or "sample.covariance" are chosen as data.type.
#' @param counts Data to be specified if "counts" is chosen as data.type.
#' @param sample.sizes Data to be specified if "counts" is chosen as data.type.
#' @param sample.covariance Data to be specified if "sample.covariance" is chosen as data.type.
#' @param target.spatial.prior.scale The variance on the spatial prior on population locations, default is half the pairwise observed distance.
#' @param source.spatial.prior.scale The variance on the spatial prior on sources of admixture, default is twice the pairwise observed distance.
#' @param spatial.prior.X.coordinates 'Observed' sample longitude, or, if you want to examine the influence of the prior, random values.
#' @param spatial.prior.Y.coordinates 'Observed' sample latitude, or, if you want to examine the influence of the prior, random values.
#' @param round.earth Option of whether you want to estimate locations on a plane (round.earth = FALSE) or a sphere (round.earth = TRUE).
#' @param long.run.initial.parameters List of parameter values that can be passed directly to the long run MCMC as initial parameter values.
#' 			The list should include values for \code{a0}, \code{a1}, \code{a2}, \code{population.coordinates}, \code{admix.proportions},
#'			and \code{nugget}, and each element of the list should named for the corresponding parameter 
#'			(e.g., \code{list("a0" = 1.07, "a1" = 0.5, ...)}).
#' @param k Number of samples.
#' @param loci Number of loci.
#' @param ngen Number of MCMC gnereations for the long MCMC.
#' @param printfreq Frequency with which updates are printed.  The updates consist of the current MCMC iteration followed by the posterior probability.
#' @param samplefreq Frequency with which samples are logged from the MCMC (basically the thinning).
#' @param mixing.diagn.freq Frequency of adaptive Metropolis-within-Gibbs updates do the tuning parameters of the proposal distributions.
#'			Default value is every 50 MCMC iterations.
#' @param savefreq Frequency with which MCMC_output object is saved.
#' @param directory Directory into which you want output to be saved.
#' @param prefix Prefix to be attached to all output files.
#'
#' @return This function saves an output R object (".Robj") which contains the results of the analysis.
#' The components of this R object are:
#' \itemize{
#' \item a0 - The posterior distribution on parameter \eqn{\alpha_0}.
#' \item a1 - The posterior distribution on parameter \eqn{\alpha_1}.
#' \item a2 - The posterior distribution on parameter \eqn{\alpha_2}.
#' \item accept_rates - The list of acceptance rates of different parameters over the course of the MCMC. 
#'		The total number of elements in each element of the list is equal to the number of sampled 
#'		MCMC iterations (i.e., the total number of generations divided by the sample frequency).
#' \item admix.proportions - The posterior distribution on admixture proportions.  This is a matrix
#'		in which the \eqn{i}th column is the vector of estimated admixture proportions from the 
#'		\eqn{i}th sampled generation of the MCMC.
#' \item diagns - The list of acceptance rates for each parameter over the last 50 MCMC iterations.
#' \item distances - The list of pairwise distances between all samples and their sources of admixture over the 
#'		course of the MCMC.  Each element of the list is a pairwise distance matrix of dimension \eqn{2*K} by
#'		\eqn{2*K}, where K is the number of samples.  The total number of elements in the list is equal to the
#'		number of sampled MCMC iterations (i.e., the total number of generations divided by the sample frequency).
#' \item last.params - The list of values passed between each iteration of the MCMC,
#' 		sampled at the last iteration of the MCMC (i.e., the location in parameter space from the
#'		very end of the analysis, along with other quantities passed between parameter update functions).
#' \item LnL_freqs - The vector of likelihood values sampled over the course of the MCMC.
#' \item lstps - A list giving the log of the scale of the tuning parameters, updated via an 
#'		adaptive MCMC procedure, for each model parameter. The total number of elements in each 
#'		element of the list is equal to the number of sampled MCMC iterations 
#'		(i.e., the total number of generations divided by the sample frequency).
#' \item ngen - The user-specified number of generations of the MCMC.
#' \item nugget - The posterior distribution on nugget parameters.  This is a matrix
#'	in which the \eqn{i}th column is the vector of estimated nuggets from the \eqn{i}th sampled 
#'	generation of the MCMC.
#' \item population.coordinates - The posterior distribution on sample coordinates in geogenetic space.  Each 
#'		element of the list is a matrix with 2 columns (Eastings and Northings, which correspond to Long and Lat 
#'		in the geogenetic space and \eqn{2*K} rows, where \eqn{K} is the number of samples in the dataset.  The first
#'		\eqn{K} rows correspond to the geogenetic coordinates of the samples themselves, and the \eqn{K+1}:\eqn{2*K}
#'		rows give the geogenetic coordinates of the source of admixture for each sample.
#' \item Prob - The vector of posterior probability values sampled over the course of the MCMC.
#' \item samplefreq - The number of iterations between each time the MCMC is sampled. A higher frequency (lower \code{samplefreq})
#'		result in more sampled iterations per analysis, with a higher autocorrelation between sampled parameter estimates.
#' \item source.spatial.prior.scale - The variance of the prior distribution on admixture source geogenetic locations.
#' \item target.spatial.prior.scale - The variance of the prior distribution on sample geogenetic locations.
#' \item transformed.covariance.list - The posterior distribution of the mean-centered and projected parametric covariance matrix.
#' 		This is of dimension \eqn{K-1} by \eqn{K-1}, where \eqn{K} is the number of samples.
#' }
#' 
#' @examples
#' # load example dataset
#' data(spacemix.example.dataset)
#' 
#' # run example analysis
#' run.spacemix.analysis(n.fast.reps = 2,
#' 			fast.MCMC.ngen = 100,
#' 			fast.model.option = "source_and_target",
#' 			long.model.option = "source_and_target",
#' 			data.type = "counts",
#' 			sample.frequencies = NULL,
#' 			mean.sample.sizes = NULL,
#' 			counts = spacemix.example.dataset$allele.counts,
#' 			sample.sizes = spacemix.example.dataset $sample.sizes,
#' 			sample.covariance = NULL,
#'			target.spatial.prior.scale = NULL,
#' 			source.spatial.prior.scale = NULL,
#' 			spatial.prior.X.coordinates = spacemix.example.dataset $population.coordinates[,1],
#' 			spatial.prior.Y.coordinates = spacemix.example.dataset $population.coordinates[,2],
#' 			round.earth = FALSE,
#' 			long.run.initial.parameters = NULL,
#' 			k = nrow(spacemix.example.dataset $allele.counts),
#' 			loci = ncol(spacemix.example.dataset $allele.counts),
#' 			ngen = 5000,
#'			printfreq = 50,
#' 			samplefreq = 5,
#' 			mixing.diagn.freq = 50,
#' 			savefreq = 5000,
#' 			directory = NULL,
#' 			prefix = "example_run")

		
run.spacemix.analysis <- function(n.fast.reps,
									fast.MCMC.ngen,
									fast.model.option,
									long.model.option,
									data.type,
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
									mixing.diagn.freq = 50,
									savefreq,
									directory=NULL,
									prefix){
	# recover()									
	if(n.fast.reps !=0){
		fast.run.dirs <- unlist(lapply(1:n.fast.reps,FUN=function(i){paste("fast_run_",i,sep="")}))
		for(i in 1:n.fast.reps){
			dir.create(fast.run.dirs[i])
			setwd(fast.run.dirs[i])
#			random.initial.population.coordinates <- cbind(rep(0,2*k),rep(0,2*k))
			random.initial.population.coordinates <- cbind( runif(2*k,min(spatial.prior.X.coordinates),
																	max(spatial.prior.X.coordinates)),
															runif(2*k,min(spatial.prior.Y.coordinates),
																	max(spatial.prior.Y.coordinates)))
			initial.parameters <- list("population.coordinates" = random.initial.population.coordinates)
			tryCatch({
				MCMC(model.option = fast.model.option,
					data.type = data.type,
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
					savefreq = fast.MCMC.ngen/2,
					directory = directory,
					prefix=paste(fast.run.dirs[i],"_",sep=""))
			},error=function(e){
							cat("fast run",i,"failed","\n") ;
							cat(paste(e,"\n",sep="")) ;
							file.rename(
									paste("../",fast.run.dirs[i],sep=""),
									paste("../","failed_",fast.run.dirs[i],sep="")) ;
							save.image(file="spacemix.dump.RData")
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
			fast.run.initial.parameters <- query.MCMC.output(list.files()[grep("output",list.files())],last.param.names=c("a0","a1","a2","population.coordinates","nugget","admix.proportions"))
		setwd("..")
		dir.create(paste(prefix,"LongRun",sep="_"))
		setwd(list.files()[grep("LongRun",list.files())])
	} else {
		fast.run.initial.parameters <- long.run.initial.parameters
	}
			MCMC(model.option = long.model.option,
				data.type = data.type,
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
				savefreq = savefreq,
				directory = directory,
                prefix=prefix)
	setwd("..")
	return("analysis completed.")
}

#' Performs a Procrustes superimposition
#'
#' This function performs a full Procrustes superimposition of the 
#' geogenetic locations estimated with SpaceMix onto a reference set 
#' of observed geographic sampling coordinates. This superimposition
#' consists of a rotation, reflection, scaling, and translation of the 
#' estimated locations so that they most closely match the reference
#' coordinates.  This function can also be used to transform the 
#' estimated sources of admixture along with the geogenetic locations,
#' so that both are in the same frame of reference.
#' 
#' @param X A k x 2 matrix of the reference coordinates.
#' @param Y A k x 2 matrix of the coordinates to be transformed.
#' @param k The number of samples.  This should be \code{nrow(X)}.
#' @param admix.source.locs The k x 2 matrix of estimated sources
#' 		of admixture. This may be left \code{NULL} if the geogenetic
#'		locations are being transformed.
#' @param option This argument, which takes the values 1 or 2, 
#'		specifies whether the function returns the transformed 
#' 		geogenetic locations (1) or transformed sources of admixture (2).
#'		
#' @return This function returns a k x 2 matrix of transformed coordinates,
#'		either those of the geogenetic locations (if option 1 is specified), 
#'		or those of the sources of admixture (if option 2 is specified).
#' 
#' @examples
#' 
#' # load example location data
#'	data(spacemix.location.data)
#' 
#' # perform Procrustes superimposition of 
#' #  geogenetic locations onto sampling coordinates
#' proc.geogen.coords <- spacemix.procrustes(X = spacemix.location.data$sample.coords,
#' 						Y = spacemix.location.data$geogen.coords,
#'  						k = nrow(spacemix.location.data$sample.coords),
#'  						admix.source.locs = NULL,
#'  						option = 1)
#' # perform Procrustes superimposition of 
#' #  admixture locations into same frame of
#' #  reference as the superimposition of the
#' #  geogenetic locations onto the sampling coordinates
#' proc.admix.coords <- spacemix.procrustes(X = spacemix.location.data$sample.coords,
#' 						Y = spacemix.location.data$geogen.coords,
#'  						k = nrow(spacemix.location.data$sample.coords),
#'  						admix.source.locs = spacemix.location.data$admix.coords,
#'  						option = 2)
#' par(mfrow=c(1,3))
#' 	plot(spacemix.location.data$sample.coords,
#' 			pch=20,xlab='',ylab='',main="Sample coordinates",
#'			col=spacemix.location.data$pop.colors)
#' 	plot(spacemix.location.data$geogen.coords,
#'			pch=20,xlab='',ylab='',main="Raw geogenetic and admixture coordinates",
#'			col=spacemix.location.data$pop.colors)
#'		points(spacemix.location.data$admix.coords,
#'			col=spacemix.location.data$admix.colors,pch=8)
#' 	plot(proc.geogen.coords,pch=20,xlab='',ylab='',
#'		main="Procrustes superimposition\n of geogenetic and admixture coordinates",
#'		col=spacemix.location.data$pop.colors)
#'		points(proc.admix.coords,col=spacemix.location.data$admix.colors,pch=8)

spacemix.procrustes <- function(X,Y,k,admix.source.locs = NULL,option){
	if(is.null(option)){
		stop("\nYou must specify an option of 1 or 2\n")
	}
	if(k != nrow(X) | k != nrow(Y) | nrow(X) != nrow(Y)){
		stop("\nYou have mis-specified the dimensions of the matrices\n")
	}
	proc.obj <- shapes::procOPA(X,Y,scale=TRUE,reflect=TRUE)
	rot.scal.locs <- proc.obj$s * Y %*% proc.obj$R
	translation.matrix <- matrix(colMeans(X) - colMeans(rot.scal.locs),nrow=k,ncol=2,byrow=TRUE)
	if(option==1){
		proc.locs <- rot.scal.locs + translation.matrix
	} else if(option==2){
		if(is.null(admix.source.locs)){
			stop("\nYou must specify admixture source locations for option 2.\n")
		}
		proc.locs <- proc.obj$s * admix.source.locs %*% proc.obj$R + translation.matrix
	}
	return(proc.locs)
}

get.procrustes.locations.posterior.list <- function(geographic.locations,population.coordinates.posterior){
	if(any(is.null(population.coordinates.posterior))){
		sampled.gen <- min(which(unlist(lapply(population.coordinates.posterior,is.null))))
	} else {
		sampled.gen <- length(population.coordinates.posterior)
	}
	geogen.coords.list <- vector(mode="list",length = sampled.gen-1)
	admix.source.coords.list <- vector(mode="list",length = sampled.gen-1)
	k <- nrow(geographic.locations)
	for(i in 1:length(geogen.coords.list)){
		geogen.coords.list[[i]] <- spacemix.procrustes(X = geographic.locations,
														Y = population.coordinates.posterior[[i]][1:k,],
														k = k,
														option = 1)
		admix.source.coords.list[[i]] <- spacemix.procrustes(X = geographic.locations,
														Y = population.coordinates.posterior[[i]][1:k,],
														k = k,
														admix.source.locs = population.coordinates.posterior[[i]][(k+1):(2*k),],
														option = 2)
	}
	return(list(geogen.coords.list = geogen.coords.list, admix.source.coords.list = admix.source.coords.list))
}

#' Fades a sample's color
#' 
#' This function fades a sample's color in proportion to its admixture proportion. 
#' A sample with a smaller admixture proportion is represented by a fainter color. 
#' When visualizing SpaceMix's output, it can be helpful to visualize the amount 
#' of admixture a sample is drawing by having samples that draw small amounts of 
#' admixture show up more faintly in figures.
#' 
#' @param pop.cols This is a vector of colors of length k in which each element
#' 		gives the color in which the corresponding sample should be plotted.
#' @param admix.proportions This is a vector of length k in which each element
#' 		is the estimated admixture proportion in the corresponding sample.
#' 
#' @return This function returns a vector of faded colors, for which the 
#' 		extent of fading is determined by the admixture proportion.
#' 
#' @examples
#' 
#' #generate example admixture proportions
#'  admix.values <- seq(0.5,0,length.out=100)
#' # make plotting colors
#'  admix.cols <- fade.admixture.source.points(rep("red",100),seq(1,0,length.out=100))
#' plot(admix.values,ylab="admix proportion",pch=19,col=admix.cols)

fade.admixture.source.points <- function(pop.cols,admix.proportions){
	faded.colors <- numeric(length(pop.cols))
	for(i in 1:length(pop.cols)){
		faded.colors[i] <- adjustcolor(pop.cols[i],admix.proportions[i])
	}
	return(faded.colors)
}

#' Plots admixture arrows
#' 
#' This function plots the admixture arrows, from the sources of admixture to their
#' targets, on the geogenetic map estimated as part of a SpaceMix analysis.
#' The width and opacity of these arrows can be proportional to their admixture
#' proportions, which provides an intuitive visualization of the amount of admixture
#' each sample is drawing.
#' 
#' @param admix.source.coords This is a k x 2 matrix of coordinates in which
#'		the two elements of the ith row give the x-coordinate and y-coordinate
#'		of the source admixture for sample i.
#' @param geogen.coords This is a k x 2 matrix of coordinates in which
#'		the two elements of the ith row give the x-coordinate and y-coordinate
#'		of the geogenetic location of the ith sample.
#' @param admix.proportions This is a vector of length k in which each element
#' 		is the estimated admixture proportion in the corresponding sample.
#' @param colors This is a vector of colors of length k in which each element
#' 		gives the color in which the corresponding sample should be plotted.
#' @param length This value determines the length of the edges of the arrow head
#'		(in inches).  The default value is 0.1.
#'
#' @return This function plots the admixture arrows on an existing plotting window.
#'		Its return value is invisible.
plot.admix.arrows <- function(admix.source.coords,geogen.coords,admix.proportions,colors=NULL,length=0.1){
	if(is.null(colors)){
		colors <- rep("black",nrow(admix.source.coords))
	}
	arrows(x0 = admix.source.coords[,1],
			y0 = admix.source.coords[,2],
			x1 = geogen.coords[,1],
			y1 = geogen.coords[,2],
			lwd = admix.proportions,
			col = colors,
			length = length)
	return(invisible("arrows!"))
}

load_MCMC_output <- function(MCMC.output.file){
    tmpenv <- environment()
	tmp <- load(MCMC.output.file,envir=tmpenv)
	mcmc.output <- lapply(tmp,get,envir=tmpenv)
	names(mcmc.output) <- tmp
	return(mcmc.output)
}

get.posterior.location.matrix.from.list <- function(posterior.list,population.index){
	post.location.matrix <- matrix(unlist(
								lapply(posterior.list,
									FUN=function(elem){elem[population.index,]})),
								nrow=length(posterior.list),ncol=2,byrow=TRUE)
	return(post.location.matrix)
}

get.credible.ellipse <- function(posterior.points,quantile){
	fit <- MASS::cov.mve(posterior.points, quantile.used = nrow(posterior.points) * quantile)
	points_in_ellipse <- posterior.points[fit$best,]
	ellipse_boundary <- stats::predict(cluster::ellipsoidhull(points_in_ellipse))
	return(ellipse_boundary)
}

plot.credible.ellipse <- function(ellipse_boundary,population.color,fading=0.3,lty=1){
	polygon(ellipse_boundary,col=adjustcolor(population.color,fading),border=1,lty=lty)
}

#' Post-processes SpaceMix output for plotting
#'
#' This function post-processes the output of a SpaceMix analysis to create
#' a list object that can be used or reference later for easy visualization 
#' of SpaceMix results.
#' 
#' @param MCMC.output.file This is the full path, given as a character argument 
#' 		(i.e., in quotes), of the MCMC_output R object generated by a SpaceMix 
#'		analysis.
#' @param geographic.locations This is a k x 2 matrix in which the ith row
#'		gives the geographic coordinates (i.e., longitude and latitude) of 
#'		the ith sample.
#' @param name.vector This is a character vector of length k in which each 
#'		element gives the name of the corresponding sample.
#' @param color.vector This is a vector of colors of length k in which each element
#' 		gives the color in which the corresponding sample should be plotted.
#' @param quantile This value determines the size of the credible interval
#' 		calculated for model parameters.  A value of 0.95 denotes the 95\% credible 
#'		interval.
#' @param burnin This value determines how many sampled iterations of the MCMC to discard
#' 		as `burn-in.' A value of 200 means that the first 200 sampled iterations of the MCMC 
#' 		will be discarded when calculating credible intervals.
#'
#' @return This function returns a list that can be used for plotting and visualization
#'		of SpaceMix output.  The components of this list are:
#' \itemize{
#' \item MCMC.output This is a list of the output of the SpaceMix analysis, 
#'			containing all the elements of the output .Robj file.
#' \item geographic.locations This is a k x 2 matrix in which the ith row
#'		gives the geographic coordinates (i.e., longitude and latitude) of 
#'		the ith sample.
#' \item name.vector This is a character vector of length k in which each 
#'		element gives the name of the corresponding sample.
#' \item color.vector This is a vector of colors of length k in which each element
#' 		gives the color in which the corresponding sample should be plotted.
#' \item quantile This value determines the size of the credible interval
#' 		calculated for model parameters.
#' \item best.iter This is the index of the sampled MCMC iteration with the largest
#'		posterior probability.  We refer to parameter estimates in that iteration as
#'		the maximum a posteriori (MAP) estimates.
#' \item admix.source.color.vector This is a vector of faded colors (the same as given
#'		in \code{color.vector}), for which the extent of fading is determined by the 
#' 		admixture proportion.  These colors, for which the opacity is proportional
#'		to the estimated admixture proportion, are used in plotting the admixture 
#'		sources and admixture arrows.
#' \item k This is the number of samples in the analysis.
#' \item MAPP.geogen.coords This is the Procrustes-transformed MAP geogenetic location 
#'		coordinates.
#' \item MAPP.admix.source.coords This is the Procrustes-transformed MAP admixture source 
#'		location coordinates.
#' \item procrustes.coord.posterior.lists This is a list of the Procrustes-transformed 
#'		location parameter coordinates.
#' 		\itemize{
#' 			\item geogen.coords.list A list of length N, where I is the number of sampled
#'					MCMC iterations.  The ith element of the list contains the Procrustes-
#'					transformed geogenetic location coordinates in the ith sampled iteration 
#'					of the MCMC.  As a whole, this list represents the posterior distribution 
#'					of geogenetic location parameters for all samples.
#' 			\item admix.source.coords.list A list of length N, where I is the number of sampled
#'					MCMC iterations.  The ith element of the list contains the Procrustes-
#'					transformed admixture source location coordinates in the ith sampled iteration 
#'					of the MCMC.  As a whole, this list represents the posterior distribution 
#'					of admixture source location parameters for all samples.
#' 		}
#' \item pp.geogen.location.matrices A list of length k in which the ith element is the Procrustes-
#' 		transformed posterior distribution of geogenetic location coordinates for the ith sample.
#' \item pp.admix.source.location.matrices A list of length k in which the ith element is the Procrustes-
#' 		transformed posterior distribution of admixture source location coordinates for the ith sample.
#' \item pp.geogen.ellipses A list of length k in which the ith element gives the boundaries of the 
#' 		95\% credible ellipse of the Procrustes-transformed posterior distribution of geogenetic 
#'		location coordinates of the ith sample.
#' \item pp.admix.source.ellipses A list of length k in which the ith element gives the boundaries of the 
#' 		95\% credible ellipse of the Procrustes-transformed posterior distribution of admixture source  
#'		location coordinates of the ith sample.
#' }

make.spacemix.map.list <- function(MCMC.output.file,geographic.locations,name.vector,color.vector,quantile=0.95,burnin=0){
	MCMC.output <- load_MCMC_output(MCMC.output.file)
	x <- seq((burnin+1),length(MCMC.output$Prob),by=1)
	best.iter <- which.max(MCMC.output$Prob[x])
	if(max(MCMC.output$Prob) > max(MCMC.output$Prob[x])){
		message("warning: the MCMC iteration with the\n 
				highest posterior probability is being\n 
				discarded as burn-in")
	}
	k <- MCMC.output$last.params$k
	if(is.null(color.vector)){
		color.vector <- rep(1,k)
	}
	admix.source.color.vector <- fade.admixture.source.points(color.vector,rowMeans(MCMC.output$admix.proportions[,x]))
	MAPP.geogen.coords <- spacemix.procrustes(X = geographic.locations,
											Y = MCMC.output$population.coordinates[[best.iter]][1:k,],
											k = k,option = 1)
	MAPP.admix.source.coords <- spacemix.procrustes(X = geographic.locations,
												Y = MCMC.output$population.coordinates[[best.iter]][1:k,],
												k = k,
												admix.source.locs = MCMC.output$population.coordinates[[best.iter]][(k+1):(2*k),],
												option=2)
	procrustes.coord.posterior.lists <- get.procrustes.locations.posterior.list(geographic.locations = geographic.locations,
																				population.coordinates.posterior = MCMC.output$population.coordinates[x])
	pp.geogen.location.matrices <- lapply(1:k,get.posterior.location.matrix.from.list,posterior.list=procrustes.coord.posterior.lists$geogen.coords.list)
	pp.admix.source.location.matrices <- lapply(1:k,get.posterior.location.matrix.from.list,posterior.list=procrustes.coord.posterior.lists$admix.source.coords.list)
	pp.geogen.ellipses <- lapply(pp.geogen.location.matrices,get.credible.ellipse,quantile)
	pp.admix.source.ellipses <- lapply(pp.admix.source.location.matrices,get.credible.ellipse,quantile)
	spacemix.map.list <- c(list(MCMC.output=MCMC.output),
							list(geographic.locations=geographic.locations),
								list(name.vector=name.vector),list(color.vector=color.vector),
								list(quantile=quantile),list(best.iter = best.iter),
								list(admix.source.color.vector = admix.source.color.vector),
								list(k = k),list(MAPP.geogen.coords = MAPP.geogen.coords),
								list(MAPP.admix.source.coords = MAPP.admix.source.coords),
								list(procrustes.coord.posterior.lists = procrustes.coord.posterior.lists),
								list(pp.geogen.location.matrices = pp.geogen.location.matrices),
								list(pp.admix.source.location.matrices = pp.admix.source.location.matrices),
								list(pp.geogen.ellipses = pp.geogen.ellipses),
								list(pp.admix.source.ellipses = pp.admix.source.ellipses))
	return(spacemix.map.list)
}

#' Plots the output of a SpaceMix analysis
#' 
#' This function plots the output of a SpaceMix analysis. Users can specify
#' whether they wish to also visualize the names of the samples, the 
#' credible ellipses of estimated sample geogenetic and admixture source 
#' location coordinates, as well as whether they wish to visualize the 
#' sources of admixture at all.
#' 
#' @param spacemix.map.list This is the name of the list produced using the function 
#'		\code{\link{make.spacemix.map.list}}.  This list should already be loaded into 
#'		R's working memory.
#' @param text This option (\code{TRUE} or \code{FALSE}) specifies whether sample 
#'		names should be printed on the plotted figure at the locations of the MAP 
#'		parameter estimates for each sample's geogenetic location and that of its 
#'		admixture source location.
#' @param ellipses This option (\code{TRUE} or \code{FALSE}) specifies whether 
#'		bivariate credible interval ellipses should be plotted for sample geogenetic 
#'		locations and admixture source locations.
#' @param source.option This option (\code{TRUE} or \code{FALSE}) specifies whether 
#'		the sources of admixture for each sample will be plotted, with accompanying 
#'		admixture arrows.
#' @param xlim This vector of length 2 gives the minimum and maximum x-value of the 
#'		plotting window.  Default is NULL, in which case the \code{xlim} used will 
#'		be the minimum and maximum credible ellipse x-coordinates.
#' @param ylim This vector of length 2 gives the minimum and maximum y-value of the 
#'		plotting window. Default is NULL, in which case the \code{xlim} used will 
#'		be the minimum and maximum credible ellipse y-coordinates.
#'
#' @return This function plots the output of a SpaceMix analysis.  If specified, 
#'		sample names are plotted at the MAP estimates of their geogenetic and admixture 
#'		source locations.  Ellipses with solid margins denote credible intervals on 
#'		the posterior distributions of sample geogenetic locations, while those with dashed 
#'		margins denote credible intervals on the posterior distributions of sample admixture 
#'		source locations. The return value of this function is invisible.
#' 
#' @examples
#' 
#' #load example spacemix.map.list object
#' data(example.spacemix.map.list)
#' 
#' #make SpaceMix output map
#' #	without plotting admixture sources, and with sample names plotted
#' make.spacemix.map(example.spacemix.map.list,source=FALSE,text=TRUE)
#' 
make.spacemix.map <- function(spacemix.map.list,text=FALSE,ellipses=TRUE,source.option=TRUE,xlim=NULL,ylim=NULL){
	with(spacemix.map.list,{ 
		plot(MAPP.geogen.coords,type='n',xlim=xlim,ylim=ylim,xlab="",ylab="")
			if(ellipses){
				lapply(1:k,FUN=function(i){plot.credible.ellipse(pp.geogen.ellipses[[i]],color.vector[i])})
			}
			if(text){
				text(MAPP.geogen.coords,col=color.vector,font=2,labels=name.vector,cex=0.7)
			}
			if(source.option){
				if(ellipses){
					lapply(1:k,FUN=function(i){plot.credible.ellipse(pp.admix.source.ellipses[[i]],admix.source.color.vector[i],fading=1,lty=2)})
				}
				text(MAPP.admix.source.coords,col= admix.source.color.vector,font=3,labels=name.vector,cex=0.7)
				plot.admix.arrows(MAPP.admix.source.coords, MAPP.geogen.coords,
									admix.proportions=MCMC.output$admix.proportions[,best.iter],
									colors=admix.source.color.vector,length=0.1)
			}
				box(lwd=2)
	})
	return(invisible("spacemix map!"))
}

#' Highlights specific samples in a geogenetic map
#' 
#' This function highlights the geogenetic locations and admixture 
#' source locations for a subset of the samples in a geogenetic map 
#' generated from a SpaceMix analysis. Users can specify the names of the 
#' samples to be highlighted, as well as whether they wish to visualize the 
#' sources of admixture for those samples.  This function should be called
#' after a a geogenetic map has already been plotted using 
#' \code{\link{make.spacemix.map}}
#' 
#' @param focal.pops This is a character vector specifying the names of the 
#'		samples to be highlighted on the existing SpaceMix map.
#' @param spacemix.map.list This is the name of the list produced using the function 
#'		\code{\link{make.spacemix.map.list}}.  This list should already be loaded into 
#'		R's working memory.
#' @param ellipses This option (\code{TRUE} or \code{FALSE}) specifies whether 
#'		bivariate credible interval ellipses should be plotted for sample geogenetic 
#'		locations and admixture source locations.
#' @param source.option This option (\code{TRUE} or \code{FALSE}) specifies whether 
#'		the sources of admixture for each sample will be plotted, with accompanying 
#'		admixture arrows.
#'
#' @return This function highlights the geogenetic locations and admixture source locations
#' 		for a specified subset of samples. The return value of this function is invisible.
#' 
#' @examples
#' 
#' #load example spacemix.map.list object
#' data(example.spacemix.map.list)
#' 
#' #make SpaceMix output map
#' #	without plotting admixture sources, and without sample names plotted
#' make.spacemix.map(example.spacemix.map.list,source=FALSE,text=FALSE)
#' 
#' #highlight two samples, without their admixture sources
#'  query.spacemix.map(focal.pops=c("Sample_9","Sample_23"),
#'						spacemix.map.list = example.spacemix.map.list,
#'						ellipses=TRUE,source.option=FALSE)
query.spacemix.map <- function(focal.pops,spacemix.map.list,ellipses=TRUE,source.option=TRUE){
	with(spacemix.map.list,{
		# browser()
		focal.indices <- match(focal.pops,name.vector)
		if(ellipses){
			for(i in 1:length(focal.indices)){
				plot.credible.ellipse(pp.geogen.ellipses[[focal.indices[i]]],color.vector[focal.indices[i]],fading=1)
			}
		}
			if(source.option){
				if(ellipses){
					for(i in 1:length(focal.indices)){
						plot.credible.ellipse(pp.admix.source.ellipses[[focal.indices[i]]], admix.source.color.vector[focal.indices[i]],fading=1,lty=2)
					}
				}
				text(MAPP.admix.source.coords[focal.indices,,drop=FALSE],col=1,font=3,labels=name.vector[focal.indices])
				arrows(	x0 = MAPP.admix.source.coords[focal.indices,1],
						y0 = MAPP.admix.source.coords[focal.indices,2],
						x1 = MAPP.geogen.coords[focal.indices,1],
						y1 = MAPP.geogen.coords[focal.indices,2],
						col= admix.source.color.vector[focal.indices],
						lwd=1,
						length=0.1)
			}
			text(MAPP.geogen.coords[focal.indices,,drop=FALSE],col=1,font=2,labels=name.vector[focal.indices],cex=1)
				box(lwd=2)
	})
	return(invisible("highlighted samples!"))
}