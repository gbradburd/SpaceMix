options(error=recover)
source(list.files()[grepl("space.mix.MCMC",list.files())])
load(list.files()[grepl("dataset",list.files())])

MCMC(
	model.option = "target",						#no_movement,source,target,source_and_target
	data.type = "counts",							#sample.covariance, sample.frequencies, counts
	likelihood.option = "wishart",								#wishart, normal_approx
	proj.mat.option = NULL,
	sample.frequencies = NULL,
	mean.sample.sizes = NULL,
	counts = spacemix.dataset$allele.counts,
	sample.sizes = spacemix.dataset$sample.sizes,
	sample.covariance = NULL,
	target.spatial.prior.scale = NULL,
	source.spatial.prior.scale = NULL,
	observed.X.coordinates = spacemix.dataset$population.coordinates[,1],
	observed.Y.coordinates = spacemix.dataset$population.coordinates[,2],
	round.earth = FALSE,
	k = nrow(spacemix.dataset$population.coordinates),
	loci = ncol(spacemix.dataset$allele.counts),
	ngen = 1e5,
	printfreq = 1e3,
	samplefreq = 1e4,
	mixing.diagn.freq = 100,
	gibbs.nugget.fineness=50,
	gibbs.spatial.fineness=50,
	gibbs.step.frequency = 1e10,
	savefreq = 1e5,
	directory = NULL,
	prefix = "",
	continue = FALSE,
	continuing.params=NULL)