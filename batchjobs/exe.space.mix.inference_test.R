options(error=recover)
source("space.mix.MCMC.R")
load("humangous_data.Robj")

#continue.analysis <- make.continuing.params("~/Desktop/Dropbox/space.mix/data/globetrotter/globe_spacemix/globe_2/globe_2_space_MCMC_output1.Robj","~/Desktop/Dropbox/space.mix/data/globetrotter/globe_spacemix/globe_3/globe_2_continuing.params.Robj")

MCMC(
	model.option = "source_and_target",
	data.type = "sample.covariance",
	proj.mat.option = NULL,
	sample.frequencies = NULL,
	mean.sample.sizes = humangous.mean.sample.sizes,
	counts = NULL,
	sample.sizes = NULL,
	sample.covariance = humangous.sample.cov,
	target.spatial.prior.scale = NULL,
	source.spatial.prior.scale = NULL,
	observed.X.coordinates = runif(nrow(humangous.sample.cov)),
	observed.Y.coordinates = rep(0,nrow(humangous.sample.cov)),
	round.earth = TRUE,
	k = nrow(humangous.sample.cov),
	loci = humangous.loci,
	ngen = 1e8,
	printfreq = 5e3,
	samplefreq = 1e5,
	mixing.diagn.freq = 50,
	gibbs.nugget.fineness=50,
	gibbs.spatial.fineness=50,
	gibbs.step.frequency = 2e8,
	savefreq = 1e7,
	directory = NULL,
	prefix = "humangous_1",
	continue = FALSE,
	continuing.params = NULL)

