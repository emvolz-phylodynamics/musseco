#' Code for simulating genealogies for a "Binary state speciation and extinction model" (BiSSE)
#' Trees are simulated using a structured coalescent framework implemented in the phydynR package
#' In this model, a wild type mutates at a constant rate to the variant type which has lower transmission fitness. 
#'
#' @author Erik M Volz <erik.volz@gmail.com>
#' 14 January 2020 


DEFAULT_SGPARMS  <- list(  tau = NULL, tau_lower = 0.1, tau_upper = 1e6, ncpu = 1, model = 1  ) 
DEFAULT_OPTIMPARMS <- list( method = 'Nelder-Mead', control = list(fnscale=-1) )
DEFAULT_COLIK_PARMS <- list( 
			AgtY_penalty = 0 
			, PL2=TRUE 
	)

#' Theoretical frequency of the wild type assuming mutation & selection balance 
#' with a variant that has lower transmission fitness
#' 
#' @return Scalar frequency of the wild type 
#' @param mu Rate of de novo mutation ancestral -> variant
#' @param omega Relative transmission fitness (must be <1) of the variant type
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
#' @param maxst, maximum sample time ; time of last tip in the tree 
#' @param Net Effective population size through time. This is stored as a two column
#' matrix, such that the first column is time (forward, since some point in the past) and the second column is population size
#' @param res Number of time points to simulate
.make.bisseco.culdesac.tfgy <- function( mu, omega, yscale, mh, maxst, Net, res = 200 )
{
	stopifnot( omega < 1  )
	demes <- c('ancestral', 'variant' )
	
	# Proportion ancestral type 
	pa <- min(1, max(0, (1 - omega*(1+mu)) / (1-omega*(1+mu) + mu ) ) )
	
	times <- seq( maxst - mh, maxst, length = res ) 
	dx <- diff( times[1:2] )
	nt <- length( times )
	
	ix <- order( Net[,1] )
	neaf <- approxfun( Net[ix, 1] , Net[ix, 2] , rule = 2 )
	ne <- neaf( times ) 
	y = 2*ne /( pa^2*(1-omega) + pa*omega )
	
	# Population sizes 
	Yancestral <-  pmax(0, yscale * y * pa )
	Yvariant   <-  pmax(0, yscale * y * (1-pa) )
	if ( pa < 0 ){
		Yancestral <- rep(0, nt)
		Yvariant <- Yancestral 
	}
	
	.Y <- lapply( 1:nt, function(k) c( ancestral = Yancestral[k] , variant = Yvariant[k]) )
	
	# Mutation rates from ancestral type to variant type
	.G <- lapply( 1:nt,  function(k) matrix( nrow=2, byrow=TRUE
	 , c( 0			, mu*Yancestral[k],
	      mu*Yvariant[k]	, 0
	 )))
	for (k in seq_along(.G)){
		rownames(.G[[k]]) = colnames(.G[[k]] ) <- demes 
	}
	
	# Transmission rates in both ancestral type and variant type 
	.F <- lapply( 1:nt, function(k) matrix( nrow=2, byrow=TRUE
	 , c( Yancestral[k]  , 0,
	      0  , omega*Yvariant[k]
	)))
	for (k in seq_along(.F)){
		rownames(.F[[k]]) = colnames(.F[[k]] ) <- demes 
	}
	 
	list( times = rev( times )
	 , births = rev( .F )
	 , migrations = rev( .G )
	 , sizes = rev( .Y  )
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
	tfgy <- .make.bisseco.culdesac.tfgy( mu, omega, 1, maxHeight, max(sampleTimes), Net, res)
	phydynR::sim.co.tree.fgy( tfgy, sampleTimes, sampleStates, ...)
}


sample_loglikelihood <- function(mu, omega, dtr)
{ 
	isa <- dtr$sampleStates[, 'ancestral']==1 
	pa <- pancestral_mutsel_balance( mu, omega )
	# dbinom( isa, size = 1, prob = pa , log=TRUE ) |> sum()
	dbinom( sum(isa), size = length(isa), prob = pa , log=TRUE )
}

#' Likelihood of the BiSSeCo model given rates of variant mutation, selection, a dated tree, and effective population size over time
#'
#' @param mu_omega a 2-vector containing parameters mu (mutation rate) and omega (relative fitness = 1+s)
#' @param dtr A phydynR::DatedTree
#' @param Net Effective population size through time. This is stored as a two column
#' matrix, such that the first column is time and the second column is population size
#' @export 
loglikelihood_bisseco <- function( parms, dtr, Net, augment_likelihood=TRUE, res = 200, ... )
{
	mu <- parms[1]
	omega <- parms[2]
	yscale <- parms[3] 
	if( omega >= 1 ) return(-Inf)
	tfgy <- .make.bisseco.culdesac.tfgy( mu, omega, yscale, dtr$maxHeight, dtr$maxSampleTime,  Net, res)
	coparms  <- modifyList( DEFAULT_COLIK_PARMS, list(...)  )
	coparms <- c(list( dtr ), list( tfgy), coparms )
	sl <- ifelse( augment_likelihood, sample_loglikelihood( mu, omega, dtr) , 0 ) 
	do.call( phydynR::colik.pik.fgy, coparms ) + sl 
}

#' Fits the parameters of a BiSSeCo model by maximum likelihood 
#' 
#' This method will estimate a de novo mutation rate (mu) reflecting within-host fitness of a variant and between host relative fitness (omega) reflecting differences in transmissibility. 
#' 
#' @param tr an ape::phylo representing a time-scaled phylogeny. These can be computed with the `treedater` package 
#' @param isvariant vector of type boolean with length equal to ape::Ntip(tr). Each element is TRUE if the corresponding element in tr$tip.label is a variant type 
#' @param Net  NULL or a matrix with two columns giving the effective population size. If NULL, the effective population size will be computed with the `mlesky` package. The first column should be time since some point in the past, and the second column should be an estimate of the effective population size. Time zero should correspond to the root of the tree 
#' @param theta0 Initial guess of (log) parameter values: mu, omega, and Ne scale. 
#' @param augment_likelihood if TRUE (default), will combine the coalescent likelihood with a binomial likelihood of sampling variants or ancestral types under the assumption of random sampling and mutation-selection balance 
#' @param optim_parms optional list of parameters passed to optim when fitting the model 
#' @param mlesky_parms optional list of parameters passed to mlesky::mlskygrid if estimating Ne(t) 
#' @export 
fitbisseco <- function(tr, isvariant, Net = NULL, theta0 = log(c(.05, .95, 1)), augment_likelihood = TRUE, optim_parms = list(), mlesky_parms = list() )
{

	sts <- node.depth.edgelength( tr )[1:Ntip(tr)] |> setNames( tr$tip.label )
	maxsts <- max(sts) 
	ssts <- matrix( 0, nrow = Ntip(tr), ncol = 2 )
	colnames(ssts) <- c( 'ancestral', 'variant' )
	rownames(ssts) <- tr$tip.label 
	ssts[ !isvariant , 'ancestral'] <- 1.0 
	ssts[ isvariant , 'variant'] <- 1.0 
	
	bdt <- phydynR::DatedTree( tr, sts, sampleStates = ssts ) 
	
	fsg <- NULL 
	if ( is.null( Net ) )
	{
		sgparms  <- modifyList( DEFAULT_SGPARMS, mlesky_parms )
		sgparms <- c( list( tr, sampleTimes = sts ), sgparms )
		fsg = do.call( mlesky::mlskygrid, sgparms )
		Net <- cbind( fsg$time, fsg$ne )
	}
	
	lfun <- function( theta  )
	{

		l = tryCatch( loglikelihood_bisseco(  exp(theta) , bdt, Net, augment_likelihood = augment_likelihood  )
			, error = function(e) -Inf )

		print( Sys.time() )
		print( exp(theta) ) 
		print( l ) 
		l 

	}
	oparms <- modifyList( DEFAULT_OPTIMPARMS , optim_parms )
	oparms <- c( list( par = theta0, fn = lfun ), oparms )
	o = do.call( optim, oparms )
	# o = optim( par = theta0, fn = lfun, method = 'Nelder-Mead', control = list(trace=-6, fnscale=-1))

	etheta <- exp( o$par ) 
	structure( 
	list( 
	     mu = etheta[1]
	     , omega = etheta[2] 
	     , s = etheta[2] -1 
	     , optimoutput = o 
	     , Net = Net 
	     , mleskyfit = fsg 
	     )
	     , class = 'bissecofit' )
}


	     # # TODO 
	     # print and summary methods 
	     # profile CIs, incl. signif level of s < 0
	     # include reversions

#' @export 
print.bissecofit <- function(x) 
{
	stopifnot( inherits(x, 'bissecofit' ))
	print( as.data.frame( coef(x)))
	print( paste( 'Likelihood:', x$optimoutput$value ))
	invisible( x )
}

#' @export 
coef.bissecofit <- function(x)
{
	stopifnot( inherits(x, 'bissecofit' ))
	c( mu = x$mu, omega = x$omega, s = x$s )
}


#' Make a plot of a bisseco tree that shows ancestral states at internal nodes estimated using a simple discrete trait analysis 
#'
#' @param tr A tree (phydynR::DatedTree) simulated using simulate_bisseco
#' @export 
bissecoplot <- function(tr, legend=FALSE)
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

