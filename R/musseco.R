#' Code for simulating genealogies for a "Binary state speciation and extinction model" (BiSSE)
#' Trees are simulated using a structured coalescent framework implemented in the phydynR package
#' In this model, a wild type mutates at a constant rate to the variant type which has lower transmission fitness. 
#'
#' @author Erik M Volz <erik.volz@gmail.com>
#' @date 14 January 2020 


#' Theoretical frequency of the wild type assuming mutation & selection balance 
#' with a variant that has lower transmission fitness
#' 
#' @return Scalar frequency of the wild type 
#' @param mu Rate of de novo mutation ancestral -> variant
#' @param omega Relative transmission fitness (<1) of the variant type
#' @export 
pancestral_mutsel_balance <- function( mu, omega ){
	min(1, max(0, (1 - omega*(1+mu)) / (1-omega*(1+mu) + mu ) ) )
}

#' Generates an epidemiological history as input for the tree simulation functions
#'
#' See _phydynR_ package for details on output format 
#' 
#' @param mu Rate of de novo mutation ancestral -> variant
#' @param omega Relative transmission fitness (<1) of the variant type
#' @param mh Maximum time in the past to simulate
#' @param Net Effective population size through time. This is stored as a two column
#' matrix, such that the first column is time and the second column is population size
#' @param res Number of time points to simulate
.make.bisseco.culdesac.tfgy <- function( mu, omega, mh, Net, res = 200 )
{
	stopifnot( omega < 1  )
	# Birth rate of wild type 
	fw <- 1 + mu 
	# Birth rate of variant 
	fr <- omega * fw 
	# Proportion wild type 
	pw <- min(1, max(0, (1 - omega*(1+mu)) / (1-omega*(1+mu) + mu ) ) )
	# print( pw )

	dx <- mh / res 
	
	times <- sort( Net[,1] )
	ne <- Net[,2] 
	
	# Population sizes 
	Yancestral <- ne * 2 * fw
	if ( pw < 0 )
		Yancestral <- rep(0, length(ne))
	Y <- pmax(0, Yancestral / pw  )
	Yvariant <- Y * (1-pw )
	nt <- length( times )
	.Y <- lapply( 1:nt, function(k) c( ancestral = Yancestral[k] , variant = Yvariant[k]) )
	
	# Mutation rates from wild type to variant type
	.G <- lapply( 1:nt,  function(k) matrix( nrow=2, byrow=TRUE
	 , c( 0  , mu*dx*Yancestral[k],
	      0  , 0
	 )))
	 
	 # Transmission rates in both wild type and variant type 
	.F <- lapply( 1:nt, function(k) matrix( nrow=2, byrow=TRUE
	 , c( fw*dx*Yancestral[k]  , 0,
	      0  , fr*dx*Yvariant[k]
	 )))
	 
	list( times = times
	 , births = .F
	 , migrations = .G
	 , sizes = .Y 
	 )
}

#' Simulate a structured coalescent tree for the 'Binary State Speciation and Extinction Model' 
#'
#'
#' @param mu Rate of de novo mutation ancestral -> variant
#' @param omega Relative transmission fitness (<1) of the variant type
#' @param maxHeight Maximum time in the past to simulate
#' @param Net Effective population size through time. This is stored as a two column
#' matrix, such that the first column is time and the second column is population size
#' @param res Number of time points to simulate
#' @param sampleTimes A vector of sample times 
#' @param sampleStates A matrix describing the state (ancestral or variant) for each sample
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
#' # Sample states: arranged in a matrix where each row gives the probability of being wildtype or variant
#' ssts <- cbind( ancestral = c( rep(1,15), rep(0, 5))
#'  , variant = c( rep( 0, 15), rep(1, 5 ) ) 
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
#' p1 <-bisseco_plot( tr1, legend=TRUE) 
#' p2 <-bisseco_plot( tr2, legend=TRUE)
simulate_bisseco <- function( mu, omega, maxHeight , Net, sampleTimes, sampleStates, res = 200,  ...)
{
	# library( phydynR )
	# library( ape )

	if (pancestral_mutsel_balance( mu, omega ) <= 0 ) 
		return(NULL)
	tfgy <- .make.bisseco.culdesac.tfgy ( mu, omega, maxHeight, Net, res)
	phydynR::sim.co.tree.fgy( tfgy, sampleTimes, sampleStates, ...)
}


DEFAULT_COLIK_PARMS <- list( 
			AgtY_penalty = 0 
			, likelihood='PL2'
	)

#' Likelihood of the BiSSeCo model given rates of variant mutation, selection, a dated tree, and effective population size over time
#'
#' @param mu_omega a 2-vector containing parameters mu (mutation rate) and omega (relative fitness = 1+s)
#' @param dtr A phydynR::DatedTree
#' @param Net Effective population size through time. This is stored as a two column
#' matrix, such that the first column is time and the second column is population size
#' @export 
loglikelihood_bisseco <- function( mu_omega, dtr, Net, res = 200, ... )
{
	mu <- mu_omega[1]
	omega <- mu_omega[2]
	if( omega >= 1 ) return(-Inf)
	maxHeight <- diff( range( Net[,1] ))
	tfgy <- .make.bisseco.culdesac.tfgy ( mu, omega, maxHeight, Net, res)
	coparms  <- modifyList( DEFAULT_COLIK_PARMS, list(...)  )
	coparms <- c( dtr, tfgy, coparms )
	do.call( phydynR::colik.pik.fgy, coparms ) 
}

#' Make a plot of a bisseco tree that shows ancestral states at internal nodes estimated using a simple discrete trait analysis 
#'
#' @param tr A tree (phydynR::DatedTree) simulated using simulate_bisseco
#' @export 
bisseco_plot <- function(tr, legend=FALSE)
{
	stopifnot(require( ggtree ))
	stopifnot(require( treeio ))

	ssts <- tr$sampleStates
	cols <- c( ancestral = 'blue', variant = 'red' )
	
	mtips <- rownames(ssts)[ ssts[,'variant']>.5 ]
	x <- setNames( rep('ancestral', ape::Ntip(tr)), tr$tip.label)
	x[ names(x) %in% mtips ] <- 'variant' 
	class( tr ) <- 'phylo'
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

