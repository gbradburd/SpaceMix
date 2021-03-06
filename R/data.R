# Documenting example datasets

#' Example SpaceMix dataset
#'
#' Example genetic data and geographic metadata 
#' for running a SpaceMix analysis
#'
#' @format A list with 3 elements
#' \describe{
#'   \item{allele.counts}{
#' 		matrix of allele counts with 
#'		\code{nrow} = number of samples and 
#'		\code{ncol} = number of loci}
#'   \item{sample.sizes}{
#'		matrix of sample sizes with 
#'		\code{nrow} = number of samples and 
#'		\code{ncol} = number of loci}
#'   \item{population.coordinates}{
#'		matrix of longitude and latitude with
#'		\code{nrow} = number of samples and 
#'		\code{ncol} = 2}
#' }
"spacemix.example.dataset"


#' Example model adequacy dataset
#'
#' Example output of processed (mean-centered 
#' and normalized) allele frequency data produced 
#' by running a SpaceMix analysis
#'
#' @format A list with 5 elements
#' \describe{
#'   \item{mean.sample.sizes}{
#' 		vector of mean sample sizes for all samples
#'		with \code{length} = number of samples}
#'   \item{sample.frequencies}{
#'		matrix of sample frequencies 
#'		(counts/sample.sizes) with 
#'		\code{nrow} = number of samples and 
#'		\code{ncol} = number of loci.}
#'   \item{normalized.sample.frequencies}{
#'		matrix of sample frequencies 
#'		normalized by mean.f * (1-mean.f),
#'		where mean.f at a locus is the 
#'		mean alelle frequency at that locus}
#'   \item{mean.centered.sample.frequencies}{
#'		matrix of sample frequencies with the 
#'		mean allele frequency at each locus 
#'		subtracted from all entries at that locus}
#'   \item{mean.centered.normalized.sample.frequencies}{
#'		matrix of sample frequencies that is 
#'		both normalized and mean-centered as described above.}
#' }
"MCN.frequencies.list"


#' Example location data
#' 
#' Example location data for visualizing
#' how a Procrustes superimposition works
#' 
#' @format A list with 5 elements, 
#'		where \code{K} = the number of samples
#' \describe{
#'   \item{geogen.coords}{
#' 		a matrix (2 x \code{K}) of geogenetic location coordinates 
#'		generated for each sample in a SpaceMix run}
#'   \item{admix.coords}{
#' 		a matrix (2 x \code{K}) of admixture source location 
#'		coordinates generated for each sample in a SpaceMix run}
#'   \item{sample.coords}{
#' 		a matrix (\code{K} x 2) of geographic sampling  
#'		coordinates for the samples in the dataset}
#'   \item{pop.colors}{
#'		a vector (\code{length = K}) of colors for 
#'		pretty plotting of sample coordinates}
#'   \item{admix.colors}{
#'		a vector (\code{length = K}) of colors for 
#'		admixture sources, generated using 
#'		\code{fade.admixture.source.points},
#'		which makes the opacity of the \code{pop.colors} 
#'		proportional to the admixture proportion for each sample}
#' }
"spacemix.location.data"

#' Example spacemix.map.list object
#' 
#' Example list generated by \code{make.spacemix.map.list} 
#' to be used in visualizing the output of a SpaceMix analysis
#' 
#' @format A list with 15 elements, 
#'		using \code{K} as the number of samples
#' 
#' \describe{
#' \item{MCMC.output}{
#'		This is a list of the output of the SpaceMix analysis, 
#'		containing all the elements of the output .Robj file.}
#' \item{geographic.locations}{
#'		This is a \code{K} x 2 matrix in which the ith row
#'		gives the geographic coordinates (i.e., longitude and  
#'		latitude) of the ith sample.}
#' \item{name.vector}{
#'		This is a character vector of length \code{K} in which each 
#'		element gives the name of the corresponding sample.}
#' \item{color.vector}{
#'		This is a vector of colors of length \code{K} 
#'		in which each element gives the color in which 
#' 		the corresponding sample should be plotted.}
#' \item{quantile}{
#'		This value determines the size of the credible 
#' 		interval calculated for model parameters.}
#' \item{best.iter}{
#'		This is the index of the sampled MCMC iteration with the largest
#'		posterior probability.  We refer to parameter estimates in that iteration as
#'		the maximum a posteriori (MAP) estimates.}
#' \item{admix.source.color.vector}{
#'		This is a vector of faded colors (the same as given
#'		in \code{color.vector}), for which the extent of fading is determined by the 
#' 		admixture proportion.  These colors, for which the opacity is proportional
#'		to the estimated admixture proportion, are used in plotting the admixture 
#'		sources and admixture arrows.}
#' \item{k}{
#'		This is the number of samples in the analysis.}
#' \item{MAPP.geogen.coords}{
#'		This is the Procrustes-transformed MAP geogenetic location 
#'		coordinates.}
#' \item{MAPP.admix.source.coords}{
#'		This is the Procrustes-transformed MAP admixture source 
#'		location coordinates.}
#' \item{procrustes.coord.posterior.lists}{
#'	 	This is a list of the Procrustes-transformed 
#'		location parameter coordinates.}
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
#' \item{pp.geogen.location.matrices}{
#'		A list of length \code{K} in which the ith element is the Procrustes-
#' 		transformed posterior distribution of geogenetic location coordinates for the ith sample.}
#' \item{pp.admix.source.location.matrices}{
#'		A list of length \code{K} in which the ith element is the Procrustes-
#' 		transformed posterior distribution of admixture source location coordinates for the ith sample.}
#' \item{pp.geogen.ellipses}{
#'		A list of length \code{K} in which the ith element gives the boundaries of the 
#' 		95\% credible ellipse of the Procrustes-transformed posterior distribution of geogenetic 
#'		location coordinates of the ith sample.}
#' \item{pp.admix.source.ellipses}{
#'		A list of length \code{K} in which the ith element gives the boundaries of the 
#' 		95\% credible ellipse of the Procrustes-transformed posterior distribution of admixture source  
#'		location coordinates of the ith sample.}
#' }
"example.spacemix.map.list"


#' Example SpaceMix output object
#' 
#' List of the objects output by a SpaceMix analysis
#' 
#' @format A list with 22 elements, 
#'		using \code{K} as the number of samples
#' 
#' \describe{
#' \item{a0}{
#'		The posterior distribution on parameter \eqn{\alpha_0}.}
#' \item{a1}{
#'		The posterior distribution on parameter \eqn{\alpha_1}.}
#' \item{a2}{
#'		The posterior distribution on parameter \eqn{\alpha_2}.}
#' \item{accept_rates}{
#'		The list of acceptance rates of different parameters over the course of the MCMC. 
#'		The total number of elements in each element of the list is equal to the number of sampled 
#'		MCMC iterations (i.e., the total number of generations divided by the sample frequency).}
#' \item{admix.proportions}{
#'		The posterior distribution on admixture proportions.  This is a matrix
#'		in which the \eqn{i}th column is the vector of estimated admixture proportions from the 
#'		\eqn{i}th sampled generation of the MCMC.}
#' \item{diagns}{
#'		The list of acceptance rates for each parameter over the last 50 MCMC iterations.}
#' \item{distances}{
#'	The list of pairwise distances between all samples and their sources of admixture over the 
#'		course of the MCMC.  Each element of the list is a pairwise distance matrix of dimension \eqn{2*K} by
#'		\eqn{2*K}. The total number of elements in the list is equal to the number of sampled MCMC 
#'		iterations (i.e., the total number of generations divided by the sample frequency).}
#' \item{last.params}{
#'		The list of values passed between each iteration of the MCMC,
#' 		sampled at the last iteration of the MCMC (i.e., the location in parameter space from the
#'		very end of the analysis, along with other quantities passed between parameter update functions).}
#' \item{LnL_freqs}{
#'		The vector of likelihood values sampled over the course of the MCMC.}
#' \item{lstps}{
#'		A list giving the log of the scale of the tuning parameters, updated via an 
#'		adaptive MCMC procedure, for each model parameter. The total number of elements in each 
#'		element of the list is equal to the number of sampled MCMC iterations 
#'		(i.e., the total number of generations divided by the sample frequency).}
#' \item{ngen}{
#'		The user-specified number of generations of the MCMC.}
#' \item{nugget}{
#'		The posterior distribution on nugget parameters.  This is a matrix
#'		in which the \eqn{i}th column is the vector of estimated nuggets from the 
#'		\eqn{i}th sampled generation of the MCMC.}
#' \item{population.coordinates}{
#'		The posterior distribution on sample coordinates in geogenetic space.  Each 
#'		element of the list is a matrix with 2 columns (Eastings and Northings, which correspond to Long and Lat 
#'		in the geogenetic space and \eqn{2*K} rows.  The first \eqn{K} rows correspond to the geogenetic 
#'		coordinates of the samples themselves, and the \eqn{K+1}:\eqn{2*K}
#'		rows give the geogenetic coordinates of the source of admixture for each sample.}
#' \item{Prob}{
#'		The vector of posterior probability values sampled over the course of the MCMC.}
#' \item{samplefreq}{
#'		The number of iterations between each time the MCMC is sampled. A higher frequency (lower \code{samplefreq})
#'		result in more sampled iterations per analysis, with a higher autocorrelation between sampled parameter estimates.}
#' \item{source.spatial.prior.scale}{
#'		The variance of the prior distribution on admixture source geogenetic locations.}
#' \item{target.spatial.prior.scale}{
#'		The variance of the prior distribution on sample geogenetic locations.}
#' \item{transformed.covariance.list}{
#'		The posterior distribution of the mean-centered and projected parametric covariance matrix.
#' 		This is of dimension \eqn{K-1} by \eqn{K-1}.}