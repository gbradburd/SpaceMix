################################################################
################################################################
#	SpaceMix plotting functions
################################################################
################################################################


################################
#	MCMC diagnosis
################################

################
#	Trace Plots
################
#	plot posterior probability trace plot
plot(Prob)

#	plot trace plots for parameters of 
#		spatial covariance function
par(mfrow=c(1,3))
	plot(a0)
	plot(a1)
	plot(a2)

#	plot trace plot for nugget parameters
#		(showing population specific drift)
matplot(t(nugget),type='l')

#	plot trace plot for admixture parameters
#		(showing population admixture proportions)
matplot(t(admix.proportions),type='l')

################
#	Acceptance Rate Plots
################

#	plot acceptance rates for parameters of 
#		spatial covariance function
par(mfrow=c(1,3))
	plot(accept_rates$a0_accept_rate)
	plot(accept_rates$a1_accept_rate)
	plot(accept_rates$a2_accept_rate)

#	plot trace plot for nugget parameters
#		(showing population specific drift)
matplot(t(accept_rates$nugget_accept_rate),type='l')

#	plot trace plot for admixture parameters
#		(showing population admixture proportions)
matplot(t(accept_rates$admix_proportions_accept_rate),type='l')


################################
#	Visualizing model output
################################

################
#	Evaluating model fit
################

#	plot fit of estimated parametric covariance matrix to 
#		sample covariance over the course of the MCMC
plot(last.params$sample.covariance,
		last.params$transformed.covariance,type='n',
		ylim=c(range(unlist(transformed.covariance.list))),
		xlim=c(range(last.params$sample.covariance)))
	abline(0,1,col="red")
	lapply(1:length(transformed.covariance.list),
			FUN=function(i){points(last.params$sample.covariance,
									transformed.covariance.list[[i]],
									pch = 20,
									col=adjustcolor(1,0.05));return()})

#	plot sample covariance against distance, and plot
#		parametric covariance against distance around that, 
#		to observe model fit

# first, load the standardized (mean-centered and normalized)
#	allele frequency data object
load(list.files(pattern="MCN.frequencies.list"))
	
	transformation.matrix <- diag(last.params$k) - 
								matrix(1/last.params$inv.mean.sample.sizes / 
									(sum(1/last.params$inv.mean.sample.sizes)),
										nrow=last.params$k,ncol=last.params$k,
								byrow=TRUE)
	
	sample.covariance <- cov(t(MCN.frequencies.list$mean.centered.normalized.sample.frequencies))
	mean.centered.parametric.covariance <- (transformation.matrix) %*% last.params$admixed.covariance %*% t(transformation.matrix)
	
	plot(last.params$D[1:last.params$k,1:last.params$k], sample.covariance,pch=19,col="red")
		points(last.params$D[1:last.params$k,1:last.params$k], mean.centered.parametric.covariance,col=1,pch=20)

# note, this can also be applied over the posterior distribution of parametric covariances, 
#	to visualize fit over the whole MCMC
		
################
#	Inference maps
################

# with observed coordinates (i.e., the sampling locations)
	observed.coordinates <- cbind(c(1,1,1,1,1,3,3,3,3,3,5,5,5,5,5,7,7,7,7,7,9,9,9,9,9,11,11,11,11,11),
									c(1,3,5,7,9,1,3,5,7,9,1,3,5,7,9,1,3,5,7,9,1,3,5,7,9,1,3,5,7,9))

# A color scheme is always nice
	pop.cols <- rainbow(last.params$k,start=4/6,end=6/6)[as.numeric(cut(observed.coordinates[,1],last.params$k))]
	plot(observed.coordinates,col=pop.cols,pch=19,cex=2)

# Without Admixture:
	plot(last.params$population.coordinates[1:last.params$k,],type='n')
		text(last.params$population.coordinates[1:last.params$k,],col=pop.cols,labels=paste(1:last.params$k))

			
# With Admixture:
	x.min <- min(last.params$population.coordinates[,1])
	x.max <- max(last.params$population.coordinates[,1])
	y.min <- min(last.params$population.coordinates[,2])
	y.max <- max(last.params$population.coordinates[,2])
	
	admix.source.pop.cols <- unlist(lapply(1:last.params$k,FUN=function(i){adjustcolor(pop.cols[i],last.params$admix.proportions[i])}))
	
	plot(last.params$population.coordinates[1:last.params$k,],type='n',
			xlim=c(x.min,x.max),
			ylim=c(y.min,y.max))
		text(last.params$population.coordinates[1:last.params$k,],
				col=pop.cols,labels=paste(1:last.params$k))
		text(last.params$population.coordinates[(last.params$k+1):(2*last.params$k),],
				col=admix.source.pop.cols,labels=paste(1:last.params$k),font=3)
		arrows(x0 = last.params$population.coordinates[(last.params$k+1):(2*last.params$k),1],
				y0 = last.params$population.coordinates[(last.params$k+1):(2*last.params$k),2],
				x1 = last.params$population.coordinates[1:last.params$k,1],
				y1 = last.params$population.coordinates[1:last.params$k,2],
				lwd = last.params$admix.proportions,
				col = admix.source.pop.cols)
	
# again, this plotting procedure can also be applied over the posterior distribution of population locations, 
#	to visualize the map over the whole MCMC
#
# also, these locations can be transformed via a procrustes superimposition around the true coordinates, 
#	to improve clarity of presentation. this is definitely advised if you're plotting locations over the posterior
#	of the MCMC, as the absolute frame of reference may wander and rotate during the course of the MCMC (although the
#	relative positions of the populations will remain stable if it has converged on the stationary distribution).