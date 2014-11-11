################################
################################
#	function for converting VCF to 
#	SpaceMix data format
################################
################################

################
#	Arguments
################
#	 counts.file.012 - the file with allele count numbers
#	 individual.population.assignment - vector of length K, 
#		where K is your number of individuals, where the kth 
#		entry is the population to which individual k belongs,
#		e.g. - c("Pop1","Pop2","Pop2","Pop3","Pop1") for K=5 from 3 populations
#	 remove.indvs - does not need to be specified, but, if desired, 
#		is the index of individuals to be removed from analysis
#	 geo.coords - K x 2 matrix of longitude and latitude for each population
#	 file.name - output file name to save Robj as.

vcf2spacemix <- function(counts.file.012,individual.population.assignment,remove.indvs=NULL,geo.coords,file.name){
	# read in genotype count data
	indv.geno.data <- as.matrix(read.table(counts.file.012))
		indv.geno.data <- indv.geno.data[,2:ncol(indv.geno.data)]
	# remove any individuals you want to drop from the analysis
		if(!is.null(remove.indvs)){
			indv.geno.data <- indv.geno.data[-remove.indvs,]
		}
	# create sample size matrix, assuming individuals are diploid	
	indv.sample.sizes <- matrix(2,nrow=nrow(indv.geno.data),ncol=ncol(indv.geno.data))
	indv.sample.sizes[which(indv.geno.data==-1,arr.ind=TRUE)] <- 0
	indv.geno.data[which(indv.geno.data==-1,arr.ind=TRUE)] <- 0
	
	# go through and lump individuals by population membership
	population.names <- unique(individual.population.assignment)
		pop.geno.data <- matrix(0,nrow=length(population.names),ncol=ncol(indv.geno.data))
		pop.sample.sizes <- matrix(0,nrow=length(population.names),ncol=ncol(indv.geno.data))
		pop.geo.coords <- numeric(length(population.names))
	for(i in 1:length(population.names)){
		pop.geno.data[i,] <- colSums(indv.geno.data[which(individual.population.assignment==population.names[i]),])
		pop.sample.sizes[i,] <- colSums(indv.sample.sizes[which(individual.population.assignment==population.names[i]),])
		pop.geo.coords[i] <- unique(geo.coords[which(individual.population.assignment==population.names[i])])
	}
	
	# return object for SpaceMix run
	save(pop.geno.data,pop.sample.sizes,geo.coords,file=paste(file.name,".Robj",sep=""))
	return(0)
}
