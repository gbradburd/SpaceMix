Geo_coord_correction <- function(population.coordinates){
	# recover()
	corrected.population.coordinates <- population.coordinates
		#correct longitude
			corrected.population.coordinates[,1][which(corrected.population.coordinates[,1] > 180)] <- -180 + abs(180 - corrected.population.coordinates[,1][which(corrected.population.coordinates[,1] > 180)])
			corrected.population.coordinates[,1][which(corrected.population.coordinates[,1] < -180)] <- 180 - abs(-180 - corrected.population.coordinates[,1][which(corrected.population.coordinates[,1] < -180)])
		#correct latitude
			corrected.population.coordinates[,2][which(corrected.population.coordinates[,2] > 90)]	<- 90 - abs(90 - corrected.population.coordinates[,2][which(corrected.population.coordinates[,2] > 90)])
			corrected.population.coordinates[,2][which(corrected.population.coordinates[,2] < -90)] <- 	-90 + abs(-90 - corrected.population.coordinates[,2][which(corrected.population.coordinates[,2] < -90)])
	return(corrected.population.coordinates)
}

Update_population_location <- function(last.params){
	# recover()
	new.params <- last.params
	pop.to.update <- sample(1:last.params$k,1)
		population.coordinates_prime <- last.params$population.coordinates + matrix(c(rep(0,(pop.to.update-1)*2),
																							rnorm(2,0,last.params$location_stp),
																								rep(0,(last.params$k*2-pop.to.update)*2)),
																							nrow=2*last.params$k,ncol=2,byrow=TRUE)
		population.coordinates_prime  <- Geo_coord_correction(population.coordinates_prime)
		D_prime <- fields::rdist.earth(population.coordinates_prime,R=1)
		covariance_prime <- Covariance(last.params$a0,last.params$aD,last.params$a2,D_prime,last.params$delta,last.params$nugget,last.params$mean.sample.sizes)		
		admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions)
			if(!any(eigen(admixed.covariance_prime)$values < 0)){
				transformed_covariance_prime <- transformed.Covariance(admixed.covariance_prime,
																		last.params$projection.matrix)
				LnL_freqs_prime <- total_likelihood_freqs(last.params$freqs,transformed_covariance_prime,last.params$loci)
					if(exp(LnL_freqs_prime - last.params$LnL_freqs) >= runif(1)){
						new.params$population.coordinates <- population.coordinates_prime
						new.params$D <- D_prime
						new.params$covariance <- covariance_prime
						new.params$admixed.covariance <- admixed.covariance_prime
						new.params$transformed_covariance <- transformed_covariance_prime
						new.params$LnL_freqs <- LnL_freqs_prime
						new.params$location_accept <- last.params$location_accept + 1
					}
			}
	new.params$location_moves <- new.params$location_moves+1
	return(new.params)	
}

	
Update_admixture_source_location <- function(last.params){
	 # recover()
	new.params <- last.params
	updated.admixture.source <- sample((last.params$k+1):(2*last.params$k),1)
	population.coordinates_prime <- last.params$population.coordinates + matrix(c(rep(0,(updated.admixture.source-1)*2),
																							rnorm(2,0,last.params$admix.source.location.stp),
																								rep(0,(last.params$k*2-updated.admixture.source)*2)),
																							nrow=2*last.params$k,ncol=2,byrow=TRUE)
	population.coordinates_prime  <- Geo_coord_correction(population.coordinates_prime)
	D_prime <- fields::rdist.earth(population.coordinates_prime,R=1)
	covariance_prime <- Covariance(last.params$a0,last.params$aD,last.params$a2,D_prime,last.params$delta,last.params$nugget,last.params$mean.sample.sizes)
	admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions)
		if(!any(eigen(admixed.covariance_prime)$values < 0)){
			transformed_covariance_prime <- transformed.Covariance(admixed.covariance_prime,
																	last.params$projection.matrix)
			LnL_freqs_prime <- total_likelihood_freqs(last.params$freqs,transformed_covariance_prime,last.params$loci)
				if(exp(LnL_freqs_prime - last.params$LnL_freqs) >= runif(1)){
					new.params$population.coordinates <- population.coordinates_prime
					new.params$D <- D_prime
					new.params$covariance <- covariance_prime
					new.params$admixed.covariance <- admixed.covariance_prime
					new.params$transformed_covariance <- transformed_covariance_prime
					new.params$LnL_freqs <- LnL_freqs_prime
					new.params$admix_source_location_accept <- last.params$admix_source_location_accept + 1
				}
		}
	new.params$admix_source_location_moves <- last.params$admix_source_location_moves + 1
	return(new.params)	
}

Update_admixture_proportion <- function(last.params){
	new.params <- last.params
	updated.population <- sample(last.params$k,1)
	admix.proportions_prime <- last.params$admix.proportions + c(rep(0,updated.population-1),
																	rnorm(1,0,last.params$admix.proportions.stp),
																		rep(0,last.params$k-updated.population))
	prior_prob_admix_proportions_prime <- Prior_prob_admix_proportions(admix.proportions_prime)
	if(!any(is.infinite(prior_prob_admix_proportions_prime))){
		admixed.covariance_prime <- admixed.Covariance(last.params$covariance,admix.proportions_prime)
			if(!any(eigen(admixed.covariance_prime)$values < 0)){
				transformed_covariance_prime <- transformed.Covariance(admixed.covariance_prime,
																		last.params$projection.matrix)
				LnL_freqs_prime <- total_likelihood_freqs(last.params$freqs,transformed_covariance_prime,last.params$loci)
					if(exp((sum(prior_prob_admix_proportions_prime) + LnL_freqs_prime) - 
								(sum(last.params$prior_prob_admix_proportions) + last.params$LnL_freqs)) >= runif(1)){
						new.params$admix.proportions <- admix.proportions_prime
						new.params$prior_prob_admix_proportions <- prior_prob_admix_proportions_prime
						new.params$admixed.covariance <- admixed.covariance_prime
						new.params$transformed_covariance <- transformed_covariance_prime
						new.params$LnL_freqs <- LnL_freqs_prime
						new.params$admix_proportion_accept <- last.params$admix_proportion_accept + 1
					}
			}
	}
	new.params$admix_proportion_moves <- last.params$admix_proportion_moves + 1
	return(new.params)
}
Update_a0 <- function(last.params){
	new.params <- last.params
	a0_prime <- last.params$a0 + rnorm(1,0,last.params$a0_stp)
	prior_prob_alpha0_prime <- Prior_prob_alpha0(a0_prime)
	if(prior_prob_alpha0_prime != -Inf){
		covariance_prime <- Covariance(a0_prime,last.params$aD,last.params$a2,last.params$D,last.params$delta,last.params$nugget,last.params$mean.sample.sizes)
			if(matrixcalc::is.positive.definite(covariance_prime)){
			admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions)
			transformed_covariance_prime <- transformed.Covariance(admixed.covariance_prime,last.params$projection.matrix)
			LnL_freqs_prime <- total_likelihood_freqs(last.params$freqs,transformed_covariance_prime,last.params$loci)
				if(exp((prior_prob_alpha0_prime + LnL_freqs_prime) - (last.params$prior_prob_alpha0 + last.params$LnL_freqs)) >= runif(1)){
					new.params$a0 <- a0_prime
					new.params$covariance <- covariance_prime
					new.params$admixed.covariance <- admixed.covariance_prime
					new.params$transformed_covariance <- transformed_covariance_prime					
					new.params$prior_prob_alpha0 <- prior_prob_alpha0_prime
					new.params$LnL_freqs <- LnL_freqs_prime
					new.params$a0_accept <- new.params$a0_accept + 1 				
				}
		}		
	}
	new.params$a0_moves <- last.params$a0_moves + 1
	return(new.params)
}

Update_aD <- function(last.params){
	new.params <- last.params
	aD_prime <- last.params$aD + rnorm(1,0,last.params$aD_stp)
	prior_prob_alphaD_prime <- Prior_prob_alphaD(aD_prime)
	if(prior_prob_alphaD_prime != -Inf){
		covariance_prime <- Covariance(last.params$a0,aD_prime,last.params$a2,last.params$D,last.params$delta,last.params$nugget,last.params$mean.sample.sizes)
			if(matrixcalc::is.positive.definite(covariance_prime)){
			admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions)
			transformed_covariance_prime <- transformed.Covariance(admixed.covariance_prime,last.params$projection.matrix)
			LnL_freqs_prime <- total_likelihood_freqs(last.params$freqs,transformed_covariance_prime,last.params$loci)
				if(exp((prior_prob_alphaD_prime + LnL_freqs_prime) - (last.params$prior_prob_alphaD + last.params$LnL_freqs)) >= runif(1)){
					new.params$aD <- aD_prime
					new.params$covariance <- covariance_prime
					new.params$admixed.covariance <- admixed.covariance_prime
					new.params$transformed_covariance <- transformed_covariance_prime					
					new.params$prior_prob_alphaD <- prior_prob_alphaD_prime
					new.params$LnL_freqs <- LnL_freqs_prime
					new.params$aD_accept <- new.params$aD_accept + 1 				
				}
		}		
	}
	new.params$aD_moves <- last.params$aD_moves + 1
	return(new.params)
}

Update_a2 <- function(last.params){
		new.params <- last.params
		a2_prime <- last.params$a2+rnorm(1,0,last.params$a2_stp)
		prior_prob_alpha2_prime <- Prior_prob_alpha2(a2_prime) 
		if(prior_prob_alpha2_prime != -Inf) {
			covariance_prime <- Covariance(last.params$a0,last.params$aD,a2_prime,last.params$D,last.params$delta,last.params$nugget,last.params$mean.sample.sizes)
			if(matrixcalc::is.positive.definite(covariance_prime)){
				admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions)
				transformed_covariance_prime <- transformed.Covariance(admixed.covariance_prime,last.params$projection.matrix)
				LnL_freqs_prime <- total_likelihood_freqs(last.params$freqs,transformed_covariance_prime,last.params$loci)
					if(exp((prior_prob_alpha2_prime + LnL_freqs_prime) - (last.params$prior_prob_alpha2 + last.params$LnL_freqs)) >= runif(1)){
						new.params$a2 <- a2_prime
						new.params$covariance <- covariance_prime
						new.params$admixed.covariance <- admixed.covariance_prime						
						new.params$transformed_covariance <- transformed_covariance_prime							
						new.params$prior_prob_alpha2 <- prior_prob_alpha2_prime
						new.params$LnL_freqs <- LnL_freqs_prime
						new.params$a2_accept <- new.params$a2_accept + 1 				
					}
			}		
		}
	new.params$a2_moves <- new.params$a2_moves + 1	
	return(new.params)
}

Update_nugget <- function(last.params){
	new.params <- last.params
		updated.pop <- sample(1:last.params$k,1)
	nugget_prime <- last.params$nugget + c(rep(0,updated.pop-1),rnorm(1,0,last.params$nugget_stp),rep(0,last.params$k-updated.pop))
	prior_prob_nugget_prime <- Prior_prob_nugget(nugget_prime)
	if(prior_prob_nugget_prime != -Inf){
		covariance_prime <- Covariance(last.params$a0,last.params$aD,last.params$a2,last.params$D,last.params$delta,nugget_prime,last.params$mean.sample.sizes)
			if(matrixcalc::is.positive.definite(covariance_prime)){
			admixed.covariance_prime <- admixed.Covariance(covariance_prime,last.params$admix.proportions)
			transformed_covariance_prime <- transformed.Covariance(admixed.covariance_prime,last.params$projection.matrix)
			LnL_freqs_prime <- total_likelihood_freqs(last.params$freqs,transformed_covariance_prime,last.params$loci)
				if(exp((prior_prob_nugget_prime + LnL_freqs_prime) - (last.params$prior_prob_nugget+last.params$LnL_freqs)) >= runif(1)){
					new.params$nugget <- nugget_prime
					new.params$covariance <- covariance_prime
					new.params$admixed.covariance <- admixed.covariance_prime
					new.params$transformed_covariance <- transformed_covariance_prime
					new.params$prior_prob_nugget <- prior_prob_nugget_prime
					new.params$LnL_freqs <- LnL_freqs_prime
					new.params$nugget_accept <- new.params$nugget_accept + 1 				
				}
		}		
	}
	new.params$nugget_moves <- last.params$nugget_moves + 1
	return(new.params)
}


Prior_prob_alpha0 <- function(a0){
		log(dgamma(a0,1,1))
	}

Prior_prob_alphaD <- function(aD){
		dexp(aD,log=TRUE)
	}

Prior_prob_alpha2 <- function(a2){
		log(dunif(a2,0.1,2))
	}

Prior_prob_nugget <- function(nugget){
	sum(dexp(nugget,log=TRUE))
}

Prior_prob_admix_proportions <- function(admix_proportions){
	dbeta(admix_proportions,shape1=0.05,shape2=1,log=TRUE)
}

total_likelihood_freqs <- function(freqs,covmat,loci) {
		cholcov <- chol(covmat)
		invcholcov <- MASS::ginv(cholcov)
		logsqrtdet <- sum(log(diag(cholcov)))
		return(-(1/2)*(sum((crossprod(invcholcov,freqs))^2))-loci*logsqrtdet)
	}
	
Covariance <- function(a0,aD,a2,GeoDist,delta,nugget,mean.sample.sizes) {
		nugget <- c(nugget,rep(0,nrow(GeoDist)/2))
		inv.mean.sample.sizes <- c(1/mean.sample.sizes,rep(0,nrow(GeoDist)/2))
		covariance <- (1/a0)*exp((-(aD*GeoDist)^a2)-delta)
			diag(covariance) <- 1/a0 + nugget + inv.mean.sample.sizes
		return(covariance)
	}

admixed.Covariance <- function(covariance,admix.proportions){
	# recover()
	admix.proportions <- admix.proportions/2
	k <- nrow(covariance)/2
	w_k <- matrix(admix.proportions,nrow=k,ncol=1)
	admixed.Covariance <- 	tcrossprod((1-w_k),(1-w_k)) * 	covariance[1:k,1:k] + 
							tcrossprod((1-w_k),(w_k)) 	* 	covariance[1:k,(k+1):(2*k)] +
							tcrossprod(w_k,(1-w_k)) 	*	covariance[(k+1):(2*k),1:k] +
							tcrossprod(w_k,w_k)			*	covariance[(k+1):(2*k),(k+1):(2*k)]
	return(admixed.Covariance)
}

transformed.Covariance <- function(covariance,projection.matrix){
	transformed.covariance <- 
		crossprod(projection.matrix,covariance) %*% projection.matrix
	return(transformed.covariance)		
}

transform.frequencies <- function(counts,sample_sizes){
	k <- nrow(counts)
	#data curation
		#remove alleles at sample frequency 0 or 1
			sample.freqs <- counts/sample_sizes
				no.samples <- which(colSums(sample_sizes)==0)		#should I drop a locus w/ 1 pop w/ sample.size=0?
				fixed.alleles <- c(which(colSums(sample.freqs)==0),which(colSums(sample.freqs)==k))
				loci.to.drop <- c(no.samples,fixed.alleles)
			counts <- counts[,-loci.to.drop]
			sample_sizes <- sample_sizes[,-loci.to.drop]
			sample.freqs <- sample.freqs[,-loci.to.drop]
	#mean center and normalize			
		mean.sample.sizes <- rowMeans(sample_sizes)
		#create mean-centering transformation matrix
			transformation.matrix <- diag(k) - matrix(mean.sample.sizes/(sum(mean.sample.sizes)),nrow=k,ncol=k,byrow=TRUE)
			mean.centered.sample.freqs <- transformation.matrix %*% sample.freqs
		#create projection matrix	
			qr.transformation.matrix <- qr(transformation.matrix)
			projection.matrix <- qr.Q(qr.transformation.matrix)[,1:qr.transformation.matrix$rank]
		#normalize sample frequencies by mean frequencies
			mean.sample.freqs <- sample.freqs - mean.centered.sample.freqs
			mean.centered.normalized.sample.freqs <- mean.centered.sample.freqs / sqrt(mean.sample.freqs * (1 - mean.sample.freqs))
	freqs <- t(projection.matrix) %*% mean.centered.normalized.sample.freqs
	transform.frequencies.list <- list("freqs" = freqs,
						"projection.matrix" = projection.matrix,
						"mean.sample.sizes" = mean.sample.sizes,
						"curated.counts" = counts,
						"curated.sample.sizes" = sample_sizes,
						"sample.frequencies" = sample.freqs,
						"mean.centered.normalized.sample.freqs" = mean.centered.normalized.sample.freqs)					
	return(transform.frequencies.list)
}


MCMC <-
function(		model.option,
				freqs = NULL,
				mean.sample.sizes = NULL,
				transform = NULL,
				counts = NULL,
				sample_sizes = NULL,
				projection.matrix = NULL,
				observed.X.coordinates,
				observed.Y.coordinates,
				k,
				loci,
				delta,
				nugget_stp,
				a0_stp,
				aD_stp,
				a2_stp,
				location_stp,
				admix.source.location.stp,
				admix.proportions.stp,
				ngen,
				printfreq,
				samplefreq,
				savefreq,
				directory=NULL,
                prefix="",
				continue=FALSE,
				continuing.params=NULL){
	
	if(!is.null(directory)){
		setwd(directory)
	}
	
	if(!is.null(transform)){
		if(is.null(counts) | is.null(sample_sizes)){
			stop("you have failed to specify allele count data and/or population sample sizes")
		}
		transformed.frequency.list <- transform.frequencies(counts,sample_sizes)
		mean.sample.sizes <- transformed.frequency.list$mean.sample.sizes
		freqs <- transformed.frequency.list$freqs
		projection.matrix <- transformed.frequency.list$projection.matrix
		save(transformed.frequency.list,file=paste(prefix,"transformed.frequency.list.Robj",sep=''))
	} else if(is.null(transform)){
		if(is.null(freqs) | is.null(projection.matrix) | is.null(mean.sample.sizes)){
			stop("you have failed to specify the transformed frequencies and/or a projection matrix and/or the vector of mean population sample sizes")
		}
	} 
		
	#Declare variables
		LnL_freqs <- numeric(ngen/samplefreq)
		Prob <- numeric(ngen/samplefreq)
		population.coordinates <- vector("list",ngen/samplefreq)
		transformed.covariance.list <- vector("list",ngen/samplefreq)
		admix.proportions <- vector("list",ngen/samplefreq)
		nugget <- matrix(0,nrow=k,ncol=ngen/samplefreq)
		a0 <- numeric(ngen/samplefreq)
		aD <- numeric(ngen/samplefreq)
		a2 <- numeric(ngen/samplefreq)
		location_moves <- numeric(ngen/samplefreq)
		nugget_moves <- numeric(ngen/samplefreq)
		admix_source_location_moves <- numeric(ngen/samplefreq)
		admix_proportion_moves <- numeric(ngen/samplefreq)
		a0_moves <- numeric(ngen/samplefreq)
		aD_moves <- numeric(ngen/samplefreq)		
		a2_moves <- numeric(ngen/samplefreq)
		a0_accept <- numeric(ngen/samplefreq)
		aD_accept <- numeric(ngen/samplefreq)		
		a2_accept <- numeric(ngen/samplefreq)
		location_accept <- numeric(ngen/samplefreq)
		nugget_accept <- numeric(ngen/samplefreq)
		admix_source_location_accept <- numeric(ngen/samplefreq)
		admix_proportion_accept <- numeric(ngen/samplefreq)


	if(!continue) {
		#INITIALIZE MCMC
				Prob[1] <- -Inf
				covariance <- matrix(0,nrow=k*2,ncol=k*2)
				pos.def.counter <- 0

			while(Prob[1] == -Inf){
				while(!matrixcalc::is.positive.definite(covariance) && pos.def.counter < 100){
						nugget[,1] <- rexp(k)
						a0[1] <- rgamma(1,1,2)
						aD[1] <- rexp(1)
						a2[1] <- runif(1,0,2)
						population.coordinates[[1]] <- rbind(cbind(observed.X.coordinates,
																observed.Y.coordinates),
														cbind(	runif(k, 
																	min = min(observed.X.coordinates), 
																	max = max(observed.X.coordinates)),
																runif(k, 
																	min = min(observed.Y.coordinates), 
																	max = max(observed.Y.coordinates))))
						population.coordinates[[1]] <- Geo_coord_correction(population.coordinates[[1]])
						if(model.option == "no_movement"){
							admix.proportions[[1]] <- numeric(k)
						} else if(model.option == "target"){
							admix.proportions[[1]] <- numeric(k)
						} else if(model.option == "source"){
							admix.proportions[[1]] <- rbeta(k,shape1=0.1,shape2=1)
						} else if(model.option == "source_and_target"){
							admix.proportions[[1]] <- rbeta(k,shape1=0.1,shape2=1)
						}
					D <- fields::rdist.earth(population.coordinates[[1]],R=1)
					covariance <- Covariance(a0[1],aD[1],a2[1],D,delta,nugget[,1],mean.sample.sizes)
					admixed.covariance <- admixed.Covariance(covariance,admix.proportions[[1]])
					transformed_covariance <- transformed.Covariance(admixed.covariance,projection.matrix)
					Initial.parameters <- list(a0 = a0[1],aD = aD[1], a2 = a2[1],
												nugget = nugget[,1],admix.proportions = admix.proportions,
												covariance = covariance, admixed.covariance = admixed.covariance,
												transformed_covariance = transformed_covariance,
												population.coordinates = population.coordinates,
												D = D, projection.matrix = projection.matrix)
					pos.def.counter <- pos.def.counter + 1			
				}
					if(pos.def.counter > 99){
						stop("the initial transformed covariance matrix is not positive definite! Either attempt to re-initialize MCMC, or increase the size of the delta shift")
					}			
                save(Initial.parameters,file=paste(prefix,"Initial.parameters.Robj",sep=''))
				LnL_freqs[1] <- total_likelihood_freqs(freqs,transformed_covariance,loci)
				prior_prob_alpha0 <- Prior_prob_alpha0(a0[1])
				prior_prob_alphaD <- Prior_prob_alphaD(aD[1])
				prior_prob_alpha2 <- Prior_prob_alpha2(a2[1])
				prior_prob_nugget <- Prior_prob_nugget(nugget[,1])
					if(model.option == "no_movement"){
						prior_prob_admix_proportions <- numeric(k)
					} else if(model.option == "target"){
						prior_prob_admix_proportions <- numeric(k)
					} else if(model.option == "source"){
						prior_prob_admix_proportions <- Prior_prob_admix_proportions(admix.proportions[[1]])
					} else if(model.option == "source_and_target"){
						prior_prob_admix_proportions <- Prior_prob_admix_proportions(admix.proportions[[1]])
					}
				Prob[1] <- LnL_freqs[1] + sum(prior_prob_admix_proportions) + prior_prob_nugget + prior_prob_alpha0 + prior_prob_alphaD + prior_prob_alpha2
			}
	} else {
		stop("uh oh")
	}
					
			if(!is.finite(Prob[1])){													
				stop("Initial probability of model is NEGATIVE INFINITY! Please attempt to initiate chain again.")
			}
			
	last.params <- list("population.coordinates" = population.coordinates[[1]],"admix.proportions" = admix.proportions[[1]],
						"a0" = a0[1],"aD" = aD[1],"a2" = a2[1],"nugget" = nugget[,1],"delta" = delta,
						"covariance" = covariance,"admixed.covariance" = admixed.covariance, "transformed_covariance" = transformed_covariance,
						"freqs" = freqs, "admix.proportions.stp" = admix.proportions.stp, "admix.source.location.stp" = admix.source.location.stp,
						"nugget_stp" = nugget_stp,"a0_stp" = a0_stp,"aD_stp" = aD_stp,"a2_stp" = a2_stp,"k" = k,"LnL_freqs" = LnL_freqs[1],
						"prior_prob_alpha0" = prior_prob_alpha0,"prior_prob_alphaD" = prior_prob_alphaD,
						"prior_prob_alpha2" = prior_prob_alpha2,"prior_prob_nugget" = prior_prob_nugget,
						"prior_prob_admix_proportions" = prior_prob_admix_proportions,						
						"a0_moves" = a0_moves[1],"aD_moves" = aD_moves[1],"a2_moves" = a2_moves[1],"nugget_moves" = nugget_moves[1],
						"admix_source_location_moves" = admix_source_location_moves[1],"admix_proportion_moves" = admix_proportion_moves[1],
						"a0_accept" = a0_accept[1],"aD_accept" = aD_accept[1],"a2_accept" = a2_accept[1],"nugget_accept" = nugget_accept[1],
						"admix_source_location_accept" = admix_source_location_accept[1],"admix_proportion_accept" = admix_proportion_accept[1],
						"loci" = loci,"D" = D,"projection.matrix" = projection.matrix,
						"location_stp" = location_stp,"location_moves" = location_moves[1],"location_accept" = location_accept[1],
						"observed.X.coordinates" = observed.X.coordinates,"observed.Y.coordinates" = observed.Y.coordinates,"mean.sample.sizes" = mean.sample.sizes)

	#Run the MCMC

	if(model.option == "no_movement"){
		Updates <- list(Update_a0,Update_aD,Update_a2,Update_nugget)
	} else if(model.option == "target"){
		Updates <- list(Update_a0,Update_aD,Update_a2,Update_nugget,Update_population_location)
	} else if(model.option == "source"){
		Updates <- list(Update_a0,Update_aD,Update_a2,Update_nugget,Update_admixture_source_location,Update_admixture_proportion)
	} else if(model.option == "source_and_target"){
		Updates <- list(Update_a0,Update_aD,Update_a2,Update_nugget,Update_population_location,Update_admixture_source_location,Update_admixture_proportion)
	}
	
	for(i in 2:ngen) {
		x <- sample(c(1:length(Updates)),1)
		new.params <- Updates[[x]](last.params)

		if(i%%samplefreq == 0){
			j <- i/samplefreq
			population.coordinates[[j]] <- new.params$population.coordinates
			transformed.covariance.list[[j]] <- new.params$transformed_covariance
			admix.proportions[[j]] <- new.params$admix.proportions
			nugget[,j] <- new.params$nugget
			a0[j] <- new.params$a0
			aD[j] <- new.params$aD
			a2[j] <- new.params$a2
			LnL_freqs[j] <- new.params$LnL_freqs
			Prob[j] <- LnL_freqs[j] +
						sum(new.params$prior_prob_admix_proportions) +
						new.params$prior_prob_nugget +
						new.params$prior_prob_alpha0 +
						new.params$prior_prob_alphaD +
						new.params$prior_prob_alpha2
			admix_source_location_moves[j] <- new.params$admix_source_location_moves
			admix_proportion_moves[j] <- new.params$admix_proportion_moves
			location_moves[j] <- new.params$location_moves		
			a0_moves[j] <- new.params$a0_moves
			aD_moves[j] <- new.params$aD_moves			
			a2_moves[j] <- new.params$a2_moves
			nugget_moves[j] <- new.params$nugget_moves			
			location_accept[j] <- new.params$location_accept
			admix_source_location_accept[j] <- new.params$admix_source_location_accept
			admix_proportion_accept[j] <- new.params$admix_proportion_accept
			a0_accept[j] <- new.params$a0_accept
			aD_accept[j] <- new.params$aD_accept			
			a2_accept[j] <- new.params$a2_accept
			nugget_accept[j] <- new.params$nugget_accept			
		}
		
		last.params <- new.params

		if(i%%printfreq == 0){
			P <- (	new.params$LnL_freqs
					+ sum(new.params$prior_prob_admix_proportions)
					+ new.params$prior_prob_nugget
					+ new.params$prior_prob_alpha0
					+ new.params$prior_prob_alphaD					
					+ new.params$prior_prob_alpha2)
			print(i)
			print(sprintf("Prob=%s___admix.source.acc.rate=%s",P,
				new.params$admix_source_location_accept/new.params$admix_source_location_moves))
		}
				
		if(i%%savefreq == 0){	
			save(last.params,
				LnL_freqs,Prob,covariance,admixed.covariance,transformed_covariance,
				population.coordinates,transformed.covariance.list,admix.proportions,
				a0,aD,a2,nugget,
				samplefreq,ngen,
				admix_source_location_moves,admix_proportion_moves,
				a0_moves,aD_moves,a2_moves,nugget_moves,
				admix_source_location_accept,admix_proportion_accept,a0_accept,aD_accept,a2_accept,nugget_accept,
				admix.source.location.stp,admix.proportions.stp,a0_stp,aD_stp,a2_stp,nugget_stp,
				location_moves, location_accept, location_stp,mean.sample.sizes,
				file=paste(prefix,sprintf("space_MCMC_output%d.Robj",1),sep=''))
		}
	}
    return(paste("Output",i,"runs to",paste(prefix,"MCMC_output*.Robj",sep=''),"."))
}