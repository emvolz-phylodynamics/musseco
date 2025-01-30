#' Code for simulating genealogies for a "Binary state speciation and extinction model" (BiSSE)
#' Trees are simulated using a structured coalescent framework implemented in the phydynR package
#' In this model, a wild type mutates at a constant rate to the variant type which has lower transmission fitness. 
#'
#' @author Erik M Volz <erik.volz@gmail.com>
#' 14 January 2020 


DEFAULT_SGPARMS  <- list(  tau = NULL, tau_lower = 0.1, tau_upper = 1e7, ncpu = 1, model = 2  ) 
DEFAULT_OPTIMPARMS <- list( method = 'Nelder-Mead', control = list(fnscale=-1) , hessian=TRUE)
DEFAULT_COLIK_PARMS <- list( 
			AgtY_penalty = 0 
			, PL2=TRUE 
	)

#' Theoretical frequency of the wild type assuming mutation & selection balance 
#' with a variant that has lower transmission fitness
pancestral_mutsel_balance1 <- function( mu, gamma, alpha, omega, tol = 1e-3 ){
	proot <- function( p ){# 
		p*(gamma*(1-p)+mu*(1-p)-mu*alpha*p) - (omega*(1-p))*(gamma*p + mu*alpha*p - mu*(1-p))  
	}
	uniroot( proot, c(0,1), tol = tol )$root
}

# changing defaults to opt tau 
.mlesky <- function( tre, sampleTimes = NULL, res = 25, tau = NULL, tau_lower = .01, tau_upper = 1e7, tau_tol = 0.001, ncross = 5, ncpu = 1, quiet = TRUE, NeStartTimeBeforePresent = Inf, ne0 = NULL, adapt_time_axis = FALSE, model = 1, formula = NULL, data = NULL ) 
{
	mlesky::mlskygrid(  tre, sampleTimes, res, tau, tau_lower, tau_upper , tau_tol , ncross , ncpu , quiet , NeStartTimeBeforePresent , ne0 , adapt_time_axis , model , formula , data  ) 
}

#' Generates an epidemiological history as input for the tree simulation functions
#'
#' See _phydynR_ package for details on output format 
.make.bisseco.culdesac.tfgy <- function( mu, gamma, alpha, omega, yscale, mh, maxst, Net, pa, res = 200 )
{
	# stopifnot( omega < 1  )
	demes <- c('ancestral', 'variant' )
	
	# Proportion ancestral type 
	# pa <- pancestral_mutsel_balance1( mu, gamma, alpha, omega )
	
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
	 , c( 0			, alpha*mu*Yancestral[k],
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
	stop('Not implemented' ) # need to update this method 

	if (pancestral_mutsel_balance( mu, omega ) <= 0 ) 
		return(NULL)
	tfgy <- .make.bisseco.culdesac.tfgy( mu, omega, 1, maxHeight, max(sampleTimes), Net, res)
	phydynR::sim.co.tree.fgy( tfgy, sampleTimes, sampleStates, ...)
}


#' Likelihood of the BiSSeCo model given rates of variant mutation, selection, a dated tree, and effective population size over time
#'
#' @param parms a 3-vector containing parameters alpha (ratio substitution rate to variant relative to ancestral type), omega (relative fitness = 1+s) and yscale (population size adjustment)
#' @param mu Mean clock rate of evolution 
#' @param gamma 1 / generation time 
#' @param dtr A phydynR::DatedTree
#' @param Net Effective population size through time. This is stored as a two column
#' @param res integer number of time points used for coalescent approximation 
#' matrix, such that the first column is time and the second column is population size
#' @export 
loglikelihood_bisseco <- function( parms, mu, gamma, dtr, Net, augment_likelihood=TRUE, res=200, ... )
{
	alpha  <- parms[1]
	omega <- parms[2]
	yscale <- parms[3] 
	if( omega >= 1 ) return(-Inf)
	pa <- pancestral_mutsel_balance1( mu, gamma, alpha, omega )
	tfgy <- .make.bisseco.culdesac.tfgy( mu, gamma, alpha, omega, yscale, dtr$maxHeight, dtr$maxSampleTime,  Net, pa, res)
	coparms  <- modifyList( DEFAULT_COLIK_PARMS, list(...)  )
	coparms <- c(list( dtr ), list( tfgy), coparms )
	sl <- ifelse( augment_likelihood, {
			if ( is.null( dtr$nanc )){
				dtr$nanc <- sum( dtr$sampleStates[, 'ancestral'])
			}
			dbinom( dtr$nanc, size = ape::Ntip(dtr), prob = pa , log=TRUE )
		} , 0 ) # 
	l = do.call( phydynR::colik.pik.fgy, coparms ) 
	l + sl 
}

#' Fits the parameters of a BiSSeCo model by maximum likelihood 
#' 
#' This method estimates a pair of relative fitness coefficients from pathogen phylogenies representing differences in transmissibility and differences in within-host fitness. 
#' Relative fitness is paramaterised by 1) the parameter \alpha which represents the fold-change in within-host replicative fitness and 
#' 2) the parameter \omega which represents the fold-change in between-host replicative fitness 
#'
#' Pathogen phylogenies are assumed to be reconstructed from population-based random samples of pathogen genomes and at most one sequence per host. 
#' Phylogenies should be time-scaled, and an estimate of the molecular clock rate of evolution should be provided, such as estimated with the `treedater` package. 
#' If not provided, the effective population size over time is estimated with the `mlesky` package, and this estimate is also provided in the returned BiSSeCo fit. 
#' 
#' 
#' @param tr an ape::phylo representing a time-scaled phylogeny. These can be computed with the `treedater` package 
#' @param isvariant vector of type boolean with length equal to ape::Ntip(tr). Each element is TRUE if the corresponding element in tr$tip.label is a variant type 
#' @param Tg Generation time, i.e. the mean time elapsed between generations. Units of this variable should match those used in tr$edge.length and 1/mu, e.g. days or years  
#' @param mu Molecular clock rate of evolution 
#' @param Net  NULL or a matrix with two columns giving the effective population size. If NULL, the effective population size will be computed with the `mlesky` package. The first column should be time since some point in the past, and the second column should be an estimate of the effective population size. Time zero should correspond to the root of the tree 
#' @param theta0 Initial guess of (log) parameter values: alpha, omega, and Ne scale. Can be a named vector or in the aforementioned order. 
#' @param augment_likelihood if TRUE (default), will combine the coalescent likelihood with a binomial likelihood of sampling variants or ancestral types under the assumption of random sampling and mutation-selection balance 
#' @param optim_parms optional list of parameters passed to optim when fitting the model 
#' @param mlesky_parms optional list of parameters passed to mlesky::mlskygrid if estimating Ne(t) 
#' @param res Integer number time steps used in coalescent likelihood 
#' @param ... Additional parameters are passed to phydynR::colik 
#' @return A fitted BiSSeCo model with coef and summary methods 
#' @export 
fitbisseco <- function(tr, isvariant, Tg, mu, Net = NULL
		       , theta0 = log(c(alpha=15, omega=.95, yscale=1))
		       , augment_likelihood = TRUE, optim_parms = list(), mlesky_parms = list(), res=200, ... )
{
	stopifnot( length( theta0 ) == 3 )
	pnames <- c( 'alpha', 'omega', 'yscale' )
	if (!is.null( names(theta0)))
		theta0 <- theta0[pnames]
	names( theta0 ) <- pnames 

	sts <- node.depth.edgelength( tr )[1:Ntip(tr)] |> setNames( tr$tip.label )
	maxsts <- max(sts) 
	ssts <- matrix( 0, nrow = Ntip(tr), ncol = 2 )
	colnames(ssts) <- c( 'ancestral', 'variant' )
	rownames(ssts) <- tr$tip.label 
	ssts[ !isvariant , 'ancestral'] <- 1.0 
	ssts[ isvariant , 'variant'] <- 1.0 
	
	bdt <- phydynR::DatedTree( tr, sts, sampleStates = ssts ) 
	bdt$nanc <- sum( bdt$sampleStates[, 'ancestral'])
	
	gamma = 1/Tg 
	
	fsg <- NULL 
	if ( is.null( Net ) )
	{
		sgparms  <- modifyList( DEFAULT_SGPARMS, mlesky_parms )
		sgparms <- c( list( tr, sampleTimes = sts ), sgparms )
		# if ( 'tau' %in% names(mlesky_parms) & is.null( mlesky_parms[['tau']] ) ) # need to modify any NULL parms manually 
			# sgparms$tau <- NULL 
		fsg = do.call( .mlesky, sgparms )
		Net <- cbind( fsg$time, fsg$ne )
	}
	
	lfun <- function( theta )
	{
		l = tryCatch( loglikelihood_bisseco(  exp(theta) , mu, gamma , bdt, Net, augment_likelihood = augment_likelihood, res = res , ... )
			, error = function(e) -Inf )
		print( Sys.time() )
		print( exp(theta) ) 
		print( l ) 
		l 
	}
	oparms <- modifyList( DEFAULT_OPTIMPARMS , optim_parms )
	oparms <- c( list( par = theta0, fn = lfun ), oparms )
	o = do.call( optim, oparms )

	fberr <- rep(NA, 3 )
	if (!is.null( o$hessian ))
		fberr <- (-o$hessian) |> solve() |> diag() |> abs() |> sqrt()
	
	etheta <- unname( exp( o$par )  )
	structure( 
		list( 
	     	     mu = mu
	     	     , alpha = etheta[1] 
	     	     , omega = etheta[2] 
	     	     , yscale = etheta[3] 
	     	     , s = etheta[2] -1 
	     	     , optimoutput = o 
	     	     , Net = Net 
	     	     , mleskyfit = fsg 
	     	     , err = fberr 
	     )
	     , class = 'bissecofit' )
}


	     # # TODO 
	     # profile CIs, incl. signif level of s < 0

# threshold for approximate CIs 
MAXCOEFVAR <- 1.5 

#' @export 
summary.bissecofit <- function(x)
{
	stopifnot( inherits(x, 'bissecofit' ))
	vnames <-  c('alpha', 'omega', 'yscale') 
	cx <- coef(x)[vnames]
	odf = as.data.frame( cx )
	colnames( odf ) <- ''
	odf$`2.5%` <- exp( log(coef(x)[vnames])-x$err*1.96 )
	odf$`97.5%` <- exp( log(coef(x)[vnames])+ x$err*1.96 )
	odf$`97.5%`[ x$err > MAXCOEFVAR ] <- Inf 
	odf <- round( odf, 3 )

	cat( 'Binary state speciation and extinction coalescent fit\n' )
	cat( '\n' )
	print( odf )
	cat( '\n' )
	cat( paste( 'Likelihood:', round(x$optimoutput$value, 4 ), '\n' ))
	invisible( list(
		coef = coef(x)
		, confint = odf 
		, fit= x$optimoutput 
		, phylodynamics = x$Net 
	))
}

#' @export 
print.bissecofit <- function(x) 
{
	summary(x)
	invisible( x )
}

#' @export 
coef.bissecofit <- function(x)
{
	stopifnot( inherits(x, 'bissecofit' ))
	c( alpha=x$alpha, omega = x$omega, yscale = x$yscale, s = x$s )
}


