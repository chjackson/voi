## From the ldr package, version 1.3.3
## By Kofi Placid Adragni and Andrew M. Raim
## https://www.jstatsoft.org/article/view/v061i03

## No modifications 

ldr <-
function(X, y=NULL, fy=NULL, Sigmas=NULL, ns=NULL, numdir=NULL, nslices=NULL, 
model=c("core", "lad", "pfc"), numdir.test=FALSE, ...)
{
	if (model=="pfc")
	{	
		if (is.null(fy)){stop("fy is not provided"); return()}
 
 		return(invisible(pfc(X=X, y=y, fy=fy, numdir=numdir, numdir.test=numdir.test, ...)))
	}
}
