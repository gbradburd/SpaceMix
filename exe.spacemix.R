setwd("~/Desktop/tmp")
options(error=recover)
source(list.files()[grepl("space.mix.MCMC",list.files())])
load(list.files()[grepl("dataset",list.files())])



run.spacemix.analysis(	n.fast.reps = 0,
						fast.MCMC.ngen = 1e3,
						model.option = "target",
						data.type = "counts",
						fast.likelihood.option = "normal_approx",
						long.likelihood.option = "wishart",
						proj.mat.option=NULL,
						sample.frequencies=NULL,
						mean.sample.sizes=NULL,
						counts = spacemix.dataset$allele.counts,
						sample.sizes = spacemix.dataset$sample.sizes,
						sample.covariance=NULL,
						target.spatial.prior.scale=NULL,
						source.spatial.prior.scale=NULL,
						spatial.prior.X.coordinates = spacemix.dataset$population.coordinates[,1],
						spatial.prior.Y.coordinates = spacemix.dataset$population.coordinates[,2],
						round.earth = FALSE,
						long.run.initial.parameters=NULL,
						k = nrow(spacemix.dataset$population.coordinates),
						loci = ncol(spacemix.dataset$allele.counts),
						ngen = 1e4,
						printfreq = 1e3,
						samplefreq = 1e1,
						mixing.diagn.freq = 50,
						gibbs.nugget.fineness=50,
						gibbs.spatial.fineness=50,
						gibbs.step.frequency=1e10,
						savefreq=1e4,
						directory=NULL,
						prefix="test")