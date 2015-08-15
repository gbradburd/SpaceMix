options(error=recover)
source(list.files()[grepl("space.mix.MCMC",list.files())])
load(list.files()[grepl("dataset",list.files())])

# A SpaceMix analysis can estimate any of 4 separate models:
#	1) "no_movement" - populations do not choose their own locations, nor can they draw admixture.
#						The only parameters to be estimated are: the alpha parameters of the spatial 
#						covariance function and the nugget parameters
#	2) "source" - populations do not choose their own locations, but they do draw admixture.
#						The parameters to be estimated are: the alpha parameters of the spatial 
#						covariance function, the nugget parameters, the locations of the sources of admixture,
#						the strength of that admixture.
#	3) "target" - populations choose their own locations, but no admixture.
#						The parameters to be estimated are: the alpha parameters of the spatial 
#						covariance function, the nugget parameters, the population locations.
#	4) "source_and_target" - populations choose their own locations AND they draw admixture.
#						The parameters to be estimated are: the alpha parameters of the spatial 
#						covariance function, the nugget parameters, the population locations,
#						the locations of the sources of admixture, and the strength of that admixture.


run.spacemix.analysis(	n.fast.reps = 0,	#the number of short, preliminary MCMC runs to do
						fast.MCMC.ngen = 1e5,	#the number of generations to run each short MCMC
						fast.model.option = "target",	#the model to be used in the short runs:  may be "no_movement","source","target","source_and_target"
						long.model.option = "source_and_target", #the model to be used in the long run:  may be "no_movement","source","target","source_and_target"
						data.type = "counts", #the data type to be used.  may be "sample.covariance","sample.frequencies","counts"
						fast.likelihood.option = "wishart", #the likelihood option to be used in the short run.  may be "wishart" or "normal_approx"
						long.likelihood.option = "wishart", #the likelihood option to be used in the long run.  may be "wishart" or "normal_approx"
						proj.mat.option=NULL, #the option for the projection matrix.  keep as NULL.
						sample.frequencies=NULL, #data to be specified if "sample.frequencies" is chosen as data.type
						mean.sample.sizes=NULL, #data to be specified if "sample.frequencies" or "sample.covariance" are chosen as data.type
						counts = spacemix.dataset$allele.counts, #data to be specified if "counts" is chosen as data.type
						sample.sizes = spacemix.dataset$sample.sizes, #data to be specified if "counts" is chosen as data.type
						sample.covariance=NULL,	#data to be specified if or "sample.covariance" is chosen as data.type
						target.spatial.prior.scale=NULL, #the variance on the spatial prior on population locations, default is half the pairwise observed distance
						source.spatial.prior.scale=NULL, #the variance on the spatial prior on sources of admixture, default is twice the pairwise observed distance
						spatial.prior.X.coordinates = spacemix.dataset$population.coordinates[,1], #'observed' sample longitude, or, if you want to examine the influence of the prior, random values
						spatial.prior.Y.coordinates = spacemix.dataset$population.coordinates[,2], #'observed' sample latitude, or, if you want to examine the influence of the prior, random values
						round.earth = FALSE, #option of whether you want to estimate locations on a plane (round.earth = FALSE) or a sphere (round.earth = TRUE)
						long.run.initial.parameters=NULL, #list of parameter values that can be passed directly to the long run MCMC as initial parameter values
						k = nrow(spacemix.dataset$population.coordinates), #number of samples
						loci = ncol(spacemix.dataset$allele.counts), #number of loci
						ngen = 1e7, #number of MCMC gnereations for the long MCMC
						printfreq = 1e3, #frequency with which updates are printed
						samplefreq = 1e4, #frequency with which samples are logged from the MCMC (basically the thinning)
						mixing.diagn.freq = 50, #frequency of adaptive Metropolis-within-Gibbs updates do the tuning parameters of the proposal distributions
						savefreq=5e5, #frequency with which MCMC_output object is saved
						directory=NULL, #directory into which you want output to be saved
						prefix="SpaceMix_example") #prefix to be attached to all output files

# optional:
#	the user can specify the sample covariance matrix in lieu of the 
#		matrices of allele counts and sample sizes.  This is recommended
#		if there is the potential for lots of missing data.  For example, 
#		if most samples are sequenced for all loci, but a few samples are 
#		sequenced at only a fraction of the loci, then the default SpaceMix
#		machinery will perform listwise deletion of missing loci, potentially
#		throwing out a huge fraction of the data.  Alternatively, the user
#		can calculate the sample covariance with pairwise deletion, which 
#		will preserve more of the data. If the user chooses this route, 
#		s/he must also specify the vector of mean sample sizes across all samples.

data.for.spacemix <- get.sample.covariance(spacemix.dataset$allele.counts,spacemix.dataset$sample.sizes)


run.spacemix.analysis(	n.fast.reps = 0,	#the number of short, preliminary MCMC runs to do
						fast.MCMC.ngen = 1e5,	#the number of generations to run each short MCMC
						fast.model.option = "target",	#the model to be used in the short runs:  may be "no_movement","source","target","source_and_target"
						long.model.option = "source_and_target", #the model to be used in the long run:  may be "no_movement","source","target","source_and_target"
						data.type = "sample.covariance", #the data type to be used.  may be "sample.covariance","sample.frequencies","counts"
						fast.likelihood.option = "wishart", #the likelihood option to be used in the short run.  may be "wishart" or "normal_approx"
						long.likelihood.option = "wishart", #the likelihood option to be used in the long run.  may be "wishart" or "normal_approx"
						proj.mat.option=NULL, #the option for the projection matrix.  keep as NULL.
						sample.frequencies=NULL, #data to be specified if "sample.frequencies" is chosen as data.type
						mean.sample.sizes=data.for.spacemix$mean.sample.sizes, #data to be specified if "sample.frequencies" or "sample.covariance" are chosen as data.type
						counts = NULL, #data to be specified if "counts" is chosen as data.type
						sample.sizes = NULL, #data to be specified if "counts" is chosen as data.type
						sample.covariance = data.for.spacemix$sample.covariance,	#data to be specified if or "sample.covariance" is chosen as data.type
						target.spatial.prior.scale=NULL, #the variance on the spatial prior on population locations, default is half the pairwise observed distance
						source.spatial.prior.scale=NULL, #the variance on the spatial prior on sources of admixture, default is twice the pairwise observed distance
						spatial.prior.X.coordinates = spacemix.dataset$population.coordinates[,1], #'observed' sample longitude, or, if you want to examine the influence of the prior, random values
						spatial.prior.Y.coordinates = spacemix.dataset$population.coordinates[,2], #'observed' sample latitude, or, if you want to examine the influence of the prior, random values
						round.earth = FALSE, #option of whether you want to estimate locations on a plane (round.earth = FALSE) or a sphere (round.earth = TRUE)
						long.run.initial.parameters=NULL, #list of parameter values that can be passed directly to the long run MCMC as initial parameter values
						k = nrow(spacemix.dataset$population.coordinates), #number of samples
						loci = data.for.spacemix$loci, #number of loci
						ngen = 1e7, #number of MCMC gnereations for the long MCMC
						printfreq = 1e3, #frequency with which updates are printed
						samplefreq = 1e4, #frequency with which samples are logged from the MCMC (basically the thinning)
						mixing.diagn.freq = 50, #frequency of adaptive Metropolis-within-Gibbs updates do the tuning parameters of the proposal distributions
						savefreq=5e5, #frequency with which MCMC_output object is saved
						directory=NULL, #directory into which you want output to be saved
						prefix="SpaceMix_example") #prefix to be attached to all output files