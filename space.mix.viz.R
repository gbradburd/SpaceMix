require(maps)
require(vegan)
	
world.map <- function(xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL){
		if(all(is.null(c(xmin,xmax,ymin,ymax)))){
			tmp.map <- map(database="world",col="black",fill=TRUE,plot=FALSE)
				x.min <- min(tmp.map$x[which(!is.na(tmp.map$x))])
				x.max <- max(tmp.map$x[which(!is.na(tmp.map$x))])
				y.min <- min(tmp.map$y[which(!is.na(tmp.map$x))])
				y.max <- max(tmp.map$y[which(!is.na(tmp.map$x))])
			quartz(width=9,height=5,pointsize=9)
			map(database="world",col="black",fill=TRUE)
				polygon(x = c(x.min,x.max,x.max,x.min) , y = c(y.min,y.min,y.max,y.max), col="gray80")		
			map(database="world",col="black",fill=TRUE,add=TRUE)
				lines(tmp.map$x,tmp.map$y,col="white")
		} else {
			tmp.map <- map(database="world",ylim=c(ymin,ymax),xlim=c(xmin,xmax),col="black",fill=TRUE,plot=FALSE)
				x.min <- min(tmp.map$x[which(!is.na(tmp.map$x))])
				x.max <- max(tmp.map$x[which(!is.na(tmp.map$x))])
				y.min <- min(tmp.map$y[which(!is.na(tmp.map$x))])
				y.max <- max(tmp.map$y[which(!is.na(tmp.map$x))])
			quartz(width=9,height=5,pointsize=9)
			map(database="world",ylim=c(ymin,ymax),xlim=c(xmin,xmax),col="black",fill=TRUE)
				polygon(x = c(x.min,x.max,x.max,x.min) , y = c(y.min,y.min,y.max,y.max), col="gray80")		
			map(database="world",ylim=c(ymin,ymax),xlim=c(xmin,xmax),col="black",fill=TRUE,add=TRUE)
				lines(tmp.map$x,tmp.map$y,col="white")
			box(lwd=1)
		}
			return(0)
}

posterior.covariance.fit <- function(MCMC.output,delay=0.5,halt=NULL){
	load(MCMC.output)
		x <- length(which(a0!=0))
	if(!is.null(halt)){
		x <- halt
	}
		for(i in 1:x){
			plot(cov(t(last.params$freqs)),transformed.covariance.list[[i]],
				ylim=c(min(unlist(transformed.covariance.list)),max(unlist(transformed.covariance.list))))
			abline(0,1,col="red")
			Sys.sleep(delay)
		}
	return(0)
}

visualize.admix.posterior <- function(sim,hgdp,input.data.object,MCMC.output,world.map,populations,track.pops,thinning=1,arrows,target.arrows,source.arrows,show.inference,show.sims,procrustes){
	# recover()
	if(!is.numeric(world.map)){
		if(world.map == "big"){
				world.map()
		} else if(world.map == "little"){
			world.map(xmin = -30,xmax = 30,ymin = -20,ymax = 20)
		}
	} else {
			world.map(xmin = world.map[1],xmax = world.map[2],ymin = world.map[3],ymax = world.map[4])
		}
	load(input.data.object)
		if(sim){
			points(sim.locations[(1:k),],col="black",
				bg=rainbow(k,start=3/6,end=1/12),
				pch=21,cex=2)
			if(populations=="ALL"){
				pop.index <- c(1:k)
			} else {
				pop.index <- populations
			}
		} else if(hgdp){
			points(as.double(HGDP.bedassle.data$hgdp.metadata$Longitude),
				as.double(HGDP.bedassle.data$hgdp.metadata$Latitude),col="black",
				bg=rainbow(HGDP.bedassle.data$number.of.populations,start=3/6,end=1/12),
				pch=21,cex=2)
			if(populations=="ALL"){
				pop.index <- c(1:length(HGDP.bedassle.data$hgdp.metadata$Population))
			} else {
				pop.index <- c(match(populations,HGDP.bedassle.data$hgdp.metadata$Population))
			}
		}
		
	if(show.inference){
		load(MCMC.output)
				x <- seq(1,length(which(a0!=0)),by=thinning)
				obs.X <- last.params$observed.X.coordinates
				obs.Y <- last.params$observed.Y.coordinates
			if(!arrows){
						for(i in x){
							if(!track.pops){
								if(!procrustes){
									points(population.coordinates[[i]][(last.params$k+pop.index),1],
											population.coordinates[[i]][(last.params$k+pop.index),2],
												col=adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.3),pch=8)
									segments( 	x0 = obs.X[pop.index],
												y0 = obs.Y[pop.index],
												x1 = population.coordinates[[i]][(last.params$k+pop.index),1],
												y1 = population.coordinates[[i]][(last.params$k+pop.index),2],
												col = adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.3),
												lwd = 5*admix.proportions[[i]][pop.index])
								} else if(procrustes){
									proc.loc <-  procrustes(cbind(obs.X,obs.Y),
																population.coordinates[[i]][(1:last.params$k),],
																	scale = TRUE)
									proc.source.loc <- proc.loc$scale *	
															population.coordinates[[i]][(1:last.params$k),] %*% 
																proc.loc$rotation + 
														matrix(proc.loc$translation,nrow=last.params$k,ncol=2,byrow=TRUE)
									points(proc.source.loc[,1],proc.source.loc[,2],col=adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.3),pch=8,cex=1)
									segments(	x0 = obs.X[pop.index],
												y0 = obs.Y[pop.index],
												x1 = proc.source.loc[pop.index,1],
												y1 = proc.source.loc[pop.index,2],
												col = adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.3),
												lwd = 5*admix.proportions[[i]][pop.index])
								}
							} else if(track.pops){
								if(!procrustes){
									points(population.coordinates[[i]][(last.params$k+pop.index),1],
											population.coordinates[[i]][(last.params$k+pop.index),2],
											col=adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.3),pch=8)
									points(population.coordinates[[i]][(pop.index),1],
											population.coordinates[[i]][(pop.index),2],
											col=adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.3),pch=15)
									segments( 	population.coordinates[[i]][pop.index,1],
												population.coordinates[[i]][pop.index,2],
												x1 = population.coordinates[[i]][(last.params$k+pop.index),1],
												y1 = population.coordinates[[i]][(last.params$k+pop.index),2],
												col = adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.3),
												lwd = 5*admix.proportions[[i]][pop.index])
								} else if(procrustes){
									proc.loc <-  procrustes(cbind(obs.X,obs.Y),
																population.coordinates[[i]][(1:last.params$k),],
																	scale = TRUE)
									proc.target.loc <- proc.loc$scale *	
															population.coordinates[[i]][(1:last.params$k),] %*% 
																proc.loc$rotation + 
														matrix(proc.loc$translation,nrow=last.params$k,ncol=2,byrow=TRUE)
									proc.source.loc <- proc.loc$scale *	
															population.coordinates[[i]][((last.params$k+1):(2*last.params$k)),] %*% 
																proc.loc$rotation + 
														matrix(proc.loc$translation,nrow=last.params$k,ncol=2,byrow=TRUE)
									points(proc.source.loc[,1],proc.source.loc[,2],col=adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.3),pch=8,cex=1)
									points(proc.target.loc[,1],proc.target.loc[,2],col=adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.3),pch=15,cex=1)
									segments(	x0 = proc.target.loc[pop.index,1],
												y0 = proc.target.loc[pop.index,2],
												x1 = proc.source.loc[pop.index,1],
												y1 = proc.source.loc[pop.index,2],
												col = adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.3),
												lwd = 5*admix.proportions[[i]][pop.index])
								}
						}
					} 
			} else if(arrows){
				if(!track.pops){
					if(source.arrows){
						if(!procrustes){
							arrows(x0 = population.coordinates[[max(x)]][(last.params$k+pop.index),1],
									y0 = population.coordinates[[max(x)]][(last.params$k+pop.index),2],
									x1 = obs.X[pop.index],
									y1 = obs.Y[pop.index],
									col = adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.8),
									lwd = 5*admix.proportions[[max(x)]][pop.index],
									length=0.2)
						} else if(procrustes){
							proc.loc <- procrustes(cbind(obs.X,obs.Y),population.coordinates[[max(x)]][(1:last.params$k),],scale = TRUE)
								proc.source.loc <- proc.loc$scale *	
															population.coordinates[[max(x)]][((last.params$k+1):(2*last.params$k)),] %*% 
																proc.loc$rotation + 
														matrix(proc.loc$translation,nrow=last.params$k,ncol=2,byrow=TRUE)
							arrows(x0 = proc.source.loc[pop.index,1],
									y0 = proc.source.loc[pop.index,2],
									x1 = obs.X[pop.index],
									y1 = obs.Y[pop.index],
									col = adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.8),
									lwd = 5*admix.proportions[[max(x)]][pop.index],
									length=0.2)
						}
					}
					if(target.arrows){
						if(!procrustes){
							arrows(x0 = obs.X[pop.index],
									y0 = obs.Y[pop.index],
									x1 = population.coordinates[[max(x)]][pop.index,1],
									y1 = population.coordinates[[max(x)]][pop.index,2],
									col = adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.8),
									lwd = 2,
									lty = "dotted",
									length = 0.2)
						} else if(procrustes){
							proc.loc <- procrustes(cbind(obs.X,obs.Y),population.coordinates[[max(x)]][(1:last.params$k),],scale = TRUE)
								proc.target.loc <- proc.loc$scale *	
															population.coordinates[[max(x)]][(1:last.params$k),] %*% 
																proc.loc$rotation + 
														matrix(proc.loc$translation,nrow=last.params$k,ncol=2,byrow=TRUE)
							arrows(x0 = obs.X[pop.index],
									y0 = obs.Y[pop.index],
									x1 = proc.target.loc[pop.index,1],
									y1 = proc.target.loc[pop.index,2],
									col = adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.8),
									lwd = 5*admix.proportions[[max(x)]][pop.index],
									lty = "dotted",
									length=0.2)
						}
					}						
				} else if(track.pops){
					if(source.arrows){
						if(!procrustes){
							arrows(x0 = population.coordinates[[max(x)]][(last.params$k+pop.index),1],
									y0 = population.coordinates[[max(x)]][(last.params$k+pop.index),2],
									x1 = population.coordinates[[max(x)]][pop.index,1],
									y1 = population.coordinates[[max(x)]][pop.index,2],
									col = adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.8),
									lwd = 2,
									length=0.2)
						} else if(procrustes){
							proc.loc <- procrustes(cbind(obs.X,obs.Y),population.coordinates[[max(x)]][(1:last.params$k),],scale = TRUE)
								proc.source.loc <- proc.loc$scale *	
															population.coordinates[[max(x)]][((last.params$k+1):(2*last.params$k)),] %*% 
																proc.loc$rotation + 
														matrix(proc.loc$translation,nrow=last.params$k,ncol=2,byrow=TRUE)
								proc.target.loc <- proc.loc$scale *	
															population.coordinates[[max(x)]][(1:last.params$k),] %*% 
																proc.loc$rotation + 
														matrix(proc.loc$translation,nrow=last.params$k,ncol=2,byrow=TRUE)
							arrows(x0 = proc.source.loc[pop.index,1],
									y0 = proc.source.loc[pop.index,2],
									x1 = proc.target.loc[pop.index,1],
									y1 = proc.target.loc[pop.index,2],
									col = adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.8),
									lwd = 5*admix.proportions[[max(x)]][pop.index],
									length = 0.2)
						}
					}
					if(target.arrows){
						if(!procrustes){
							arrows(x0 = obs.X[pop.index],
									y0 = obs.Y[pop.index],
									x1 = population.coordinates[[max(x)]][pop.index,1],
									y1 = population.coordinates[[max(x)]][pop.index,2],
									col = adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.8),
									lwd = 2,
									lty = "dotted",
									length = 0.2)
						} else if(procrustes){
							proc.loc <- procrustes(cbind(obs.X,obs.Y),population.coordinates[[max(x)]][(1:last.params$k),],scale = TRUE)
								proc.target.loc <- proc.loc$scale *	
															population.coordinates[[max(x)]][(1:last.params$k),] %*% 
																proc.loc$rotation + 
														matrix(proc.loc$translation,nrow=last.params$k,ncol=2,byrow=TRUE)
							arrows(x0 = obs.X[pop.index],
									y0 = obs.Y[pop.index],
									x1 = proc.target.loc[pop.index,1],
									y1 = proc.target.loc[pop.index,2],
									col = adjustcolor(rainbow(last.params$k,start=3/6,end=1/12)[pop.index],0.8),
									lwd = 2,
									lty = "dotted",
									length = 0.2)
						}
					}
				}
			}
		}
		if(sim){
			if(show.sims){
				arrows(x0 = sim.locations[c((k+1):(2*k))[pop.index],1],
						y0 = sim.locations[c((k+1):(2*k))[pop.index],2],
						x1 = sim.locations[c(1:k)[pop.index],1],
						y1 = sim.locations[c(1:k)[pop.index],2],
						col = "white",
						lwd = 10*sim.admix.proportions[pop.index],
						length = 0.2)
				arrows(x0 = sim.locations[c((k+1):(2*k))[pop.index],1],
						y0 = sim.locations[c((k+1):(2*k))[pop.index],2],
						x1 = sim.locations[c(1:k)[pop.index],1],
						y1 = sim.locations[c(1:k)[pop.index],2],
						col = adjustcolor(rainbow(k,start=3/6,end=1/12)[pop.index],0.8),
						lwd = 5*sim.admix.proportions[pop.index],
						length = 0.2)
			}
		}
}

			
############################
#	admixture!
############################
																
																							
		
		visualize.admix.posterior(	sim = FALSE,
									hgdp = TRUE,
									input.data.object = "~/Desktop/sim16/hgdp16.4/hgdp.space.data.Robj",
									MCMC.output = "~/Desktop/sim16/hgdp16.4/hgdp_16.4space_MCMC_output1.Robj",
									world.map = "big", #c(-20,150,0,68)
									populations = "ALL",
									track.pops = TRUE,
									thinning = 20,
									arrows = FALSE,
									target.arrows = TRUE,
									source.arrows = TRUE,
									show.inference = TRUE,
									show.sims = FALSE,
									procrustes = FALSE)

		visualize.admix.posterior(	sim = FALSE,
									hgdp = TRUE,
									input.data.object = "~/Desktop/hgdp15.2/hgdp.space.data.Robj",
									MCMC.output = "~/Desktop/hgdp15.2/hgdp_15.2space_MCMC_output1.Robj",
									world.map = c(-20,150,8,68),
									populations = "ALL",
									track.pops = TRUE,
									thinning = 20,
									arrows = FALSE,
									target.arrows = TRUE,
									source.arrows = TRUE,
									show.inference = TRUE,
									show.sims = FALSE,
									procrustes = TRUE)
									
									
		load("~/Desktop/hgdp15.1/hgdp.space.data.Robj")
		load("~/Desktop/hgdp15.1/hgdp_15.1space_MCMC_output1.Robj")
			plot(last.params$admix.proportions,type="n")
			text(last.params$admix.proportions,
					labels=HGDP.bedassle.data$hgdp.metadata$Population)

		load("~/Desktop/hgdp15.2/hgdp_15.2space_MCMC_output1.Robj")
			plot(last.params$admix.proportions,type="n")
			text(last.params$admix.proportions,
					labels=HGDP.bedassle.data$hgdp.metadata$Population)

			plot(last.params$nugget,type="n")
			text(last.params$nugget,
					labels=HGDP.bedassle.data$hgdp.metadata$Population)

		posterior.covariance.fit("~/Desktop/hgdp15.1/hgdp_15.1space_MCMC_output1.Robj",delay=0.01,halt=446)
		posterior.covariance.fit("~/Desktop/Dropbox/space.mix/wu_tang/sims/sim14/sims14.5/sims_14.5space_MCMC_output1.Robj",halt=20)


########################################
#	sims
########################################

		visualize.admix.posterior(	sim = TRUE,
									hgdp = FALSE,
									input.data.object = "~/Desktop/sim16/sim16.1/spacemix_sim_data_14.5.Robj",
									MCMC.output = "~/Desktop/sim16/sim16.1/sims_16.1space_MCMC_output1.Robj",
									world.map = "little",
									populations = "ALL",
									track.pops = TRUE,
									thinning = 20,
									arrows = TRUE,
									target.arrows = TRUE,
									source.arrows = TRUE,
									show.inference = TRUE,
									show.sims = FALSE,
									procrustes = TRUE)

	visualize.admix.posterior(	sim = TRUE,
									hgdp = FALSE,
									input.data.object = "~/Desktop/sim16/sim16.1/spacemix_sim_data_14.5.Robj",
									MCMC.output = "~/Desktop/sim16/sim16.2/sims_16.1space_MCMC_output1.Robj",
									world.map = "little",
									populations = "ALL",
									track.pops = TRUE,
									thinning = 20,
									arrows = TRUE,
									target.arrows = TRUE,
									source.arrows = TRUE,
									show.inference = TRUE,
									show.sims = FALSE,
									procrustes = TRUE)

	visualize.admix.posterior(	sim = TRUE,
									hgdp = FALSE,
									input.data.object = "~/Desktop/sim16/sim16.1/spacemix_sim_data_14.5.Robj",
									MCMC.output = "~/Desktop/sim16/sim16.3/sims_16.1space_MCMC_output1.Robj",
									world.map = "little",
									populations = "ALL",
									track.pops = TRUE,
									thinning = 20,
									arrows = TRUE,
									target.arrows = TRUE,
									source.arrows = TRUE,
									show.inference = TRUE,
									show.sims = FALSE,
									procrustes = TRUE)

	visualize.admix.posterior(	sim = TRUE,
									hgdp = FALSE,
									input.data.object = "~/Desktop/sim16/sim16.1/spacemix_sim_data_14.5.Robj",
									MCMC.output = "~/Desktop/sim16/sim16.4/sims_16.1space_MCMC_output1.Robj",
									world.map = "little",
									populations = "ALL",
									track.pops = TRUE,
									thinning = 20,
									arrows = FALSE,
									target.arrows = TRUE,
									source.arrows = TRUE,
									show.inference = TRUE,
									show.sims = FALSE,
									procrustes = TRUE)									

	posterior.covariance.fit("~/Desktop/sim16/sim16.3/sims_16.1space_MCMC_output1.Robj",halt=10)


