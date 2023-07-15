#' Code for simulating genealogies for a "Binary state speciation and extinction model" (BiSSE)
#' Trees are simulated using a structured coalescent framework implemented in the phydynR package
#' In this model, a wild type mutates at a constant rate to the mutant type which has lower transmission fitness. 
#'
#' @author Erik M Volz <erik.volz@gmail.com>
#' @date 14 January 2020 


#' Theoretical frequency of the wild type assuming mutation & selection balance 
#' with a mutant that has lower transmission fitness
#' 
#' @return Scalar frequency of the wild type 
#' @param mu Rate of de novo mutation wt -> mutant
#' @param omega Relative transmission fitness (<1) of the mutant type
#' @export 
pwt_mutsel_balance <- function( mu, omega ){
	min(1, max(0, (1 - omega*(1+mu)) / (1-omega*(1+mu) + mu ) ) )
}

#' Generates an epidmielogical history as input for the tree simulation functions
#'
#' See _phydynR_ package for details on output format 
#' 
#' @param mu Rate of de novo mutation wt -> mutant
#' @param omega Relative transmission fitness (<1) of the mutant type
#' @param mh Maximum time in the past to simulate
#' @param Net Effective population size through time. This is stored as a two column
#' matrix, such that the first column is time and the second column is population size
#' @param res Number of time points to simulate
.make.bisseco.culdesac.tfgy <- function( mu, omega, mh, Net, res = 200 )
{
	# Birth rate of wild type 
	fw <- 1 + mu 
	# Birht rate of mutant 
	fr <- omega * fw 
	# Proportion wild type 
	pw <- (1 - omega*(1+mu)) / (1-omega*(1+mu) + mu ) # todo 0-1?
	dx <- mh / res 
	
	times <- sort( Net[,1] )
	ne <- Net[,2] 
	
	# Population sizes 
	Ywt <- ne * 2 * fw
	if ( pw < 0 )
		Ywt <- rep(0, length(ne))
	Y <- pmax(0, Ywt / pw  )
	Ymutant <- Y * (1-pw )
	nt <- length( times )
	.Y <- lapply( 1:nt, function(k) c( wt = Ywt[k] , mutant = Ymutant[k]) )
	
	# Mutation rates from wild type to mutant type
	.G <- lapply( 1:nt,  function(k) matrix( nrow=2, byrow=TRUE
	 , c( 0  , mu*dx*Ywt[k],
	      0  , 0
	 )))
	 
	 # Transmission rates in both wild type and mutant type 
	.F <- lapply( 1:nt, function(k) matrix( nrow=2, byrow=TRUE
	 , c( fw*dx*Ywt[k]  , 0,
	      0  , fr*dx*Ymutant[k]
	 )))
	 
	list( times = times
	 , births = .F
	 , migrations = .G
	 , size = .Y 
	 )
}

#' Simulate a structured coalescent tree for the 'Binary State Speciation and Extinction Model' 
#'
#'
#' @param mu Rate of de novo mutation wt -> mutant
#' @param omega Relative transmission fitness (<1) of the mutant type
#' @param maxHeight Maximum time in the past to simulate
#' @param Net Effective population size through time. This is stored as a two column
#' matrix, such that the first column is time and the second column is population size
#' @param res Number of time points to simulate
#' @param sampleTimes A vector of sample times 
#' @param sampleStates A matrix describing the state (wt or mutant) for each sample
#' @return A simulated genealogy in ape::phylo format 
#' @export 
#' @examples 
#' set.seed(23) 
#' # Sample times: these will be distributed over 150 years to show role of heterochronous sampling
#' sts <- setNames( sample( seq( 300, 500, length=20 )), letters[1:20] )
#' 
#' # We will have the effective sample size increase linearly through time so that external branch lengths are longer: 
#' net <- cbind( seq(0, 501, length=500),  seq(1,200,length=500)  )
#' 
#' # Sample states: arranged in a matrix where each row gives the probability of being wildtype or mutant
#' ssts <- cbind( wt = c( rep(1,15), rep(0, 5))
#'  , mutant = c( rep( 0, 15), rep(1, 5 ) ) 
#'  )
#' rownames( ssts ) <- letters[1:20]
#' 
#' # Scenario A: Low denovo mutation rate and relatively high transmission fitness
#' tr1 = simulate_bisseco( .001, .995
#'  , maxHeight = 500
#'  , Net = net
#'  , res = 500
#'  , sampleTimes = sts
#'  , sampleStates = ssts 
#' ) 
#' 
#' # Scenario B: High de novo mutation rate and relatively low transission fitness
#' tr2 = simulate_bisseco( .05, .75
#'  , maxHeight = 500
#'  , Net = net 
#'  , res = 500
#'  , sampleTimes = sts
#'  , sampleStates = ssts 
#' ) 
#' 
#' # Plot the annotated trees for both scenarios
#' p1 <-bisseco_plot( tr1, legend = TRUE ) 
#' p2 <-bisseco_plot(tr2 , legend=TRUE)
simulate_bisseco <- function( mu, omega, maxHeight , Net, sampleTimes, sampleStates, res = 200,  ...)
{
	# library( phydynR )
	# library( ape )

	if (pwt_mutsel_balance( mu, omega ) <= 0 ) 
		return(NULL)
	tfgy <- .make.bisseco.culdesac.tfgy ( mu, omega, maxHeight, Net, res)
browser()
	phydynR::sim.co.tree.fgy( tfgy, sampleTimes, sampleStates, ...)
}

#' Likelihood of the BiSSeCo model given rates of variant mutation, selection, a dated tree, and effective population size over time
#'
#' @param mu_omega a 2-vector containing parameters mu (mutation rate) and omega (selection coef)
#' @param dtr A phydynR::DatedTree
#' @param maxHeight Maximum time in the past to simulate
#' @param Net Effective population size through time. This is stored as a two column
#' matrix, such that the first column is time and the second column is population size
likelihood_bisseco <- function( mu_omega, dtr, Net, res = 200, ... )
{
	mu <- mu_omega[1]
	omega <- mu_omega[2]
	maxHeigh <- diff( range( Net[,1] ))
	tfgy <- .make.bisseco.culdesac.tfgy ( mu, omega, maxHeight, Net, res)
	phydynR::colik.pik.fgy(dtr 
		  , tfgy
		  , ...
		)
}

set.seed( 23 ) 

# Sample times: these will be distributed over 150 years to show role of heterochronous sampling
sts <- setNames( sample( seq( 300, 500, length=20 )), letters[1:20] )

# We will have the effective sample size increase linearly through time so that external branch lengths are longer: 
net <- cbind( seq(0, 501, length=500),  seq(1,200,length=500)  )

# Sample states: arranged in a matrix where each row gives the probability of being wildtype or mutant
ssts <- cbind( wt = c( rep(1,15), rep(0, 5))
 , mutant = c( rep( 0, 15), rep(1, 5 ) ) 
 )
rownames( ssts ) <- letters[1:20]

# Scenario A: Low denovo mutation rate and relatively high transmission fitness
tr1 = simulate_bisseco( .001, .995
 , maxHeight = 500
 , Net = net
 , res = 500
 , sampleTimes = sts
 , sampleStates = ssts 
) 

# (l1 <- likelihood_bisseco( c(.001,.995), tr1, net ))


#' Make a plot of a bisseco tree that shows ancestral states at internal nodes estimated using a simple discrete trait analysis 
#'
#' @param tr A tree simulated using simulate_bisseco
#' @param ssts A matrix of sample states used in simulate_bisseco 
#' @export 
bisseco_plot <- function(tr, ssts, legend=FALSE)
{
	stopifnot(require( ggtree ))
	stopifnot(require( treeio ))

	cols <- c( wt = 'blue', mutant = 'red' )
	
	mtips <- rownames(ssts)[ sst[,'mutant']>.5 ]
	x <- setNames( rep('wt', ape::Ntip(tr)), tr$tip.label)
	x[ names(x) %in% mtips ] <- 'mutant' 
	tr2 <- treeio::full_join(tr, data.frame(label = names(x), stat = x ), by = 'label')
	p <- ggtree(tr2) +
		geom_tippoint(aes(color = stat)) + 
		ggplot2::scale_color_manual(values = cols)

	fitER <- ape::ace(x,tr,model="ER",type="discrete")
	ancstats <- as.data.frame(fitER$lik.anc)
	ancstats$node <- 1:tr$Nnode+Ntip(tr)

	pies <- nodepie(ancstats, cols = 1:2)
	pies <- lapply(pies, function(g) g+ ggplot2::scale_fill_manual(values = cols))
	p = p + geom_inset(pies, width = .1, height = .1) 
	if (!legend)
		p = p + theme(legend.position='none')
		
	p
}

