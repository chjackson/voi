## DONE

* Basic package structure

* EVPPI methods copied in from BCEA and tidied up

* A few unit tests added for these 

* Nonparametric regression and importance sampling methods added for EVSI, and rough tests written


## TODO

EVPPI

* GP.  Compare with SAVI code, and use that if it does more 

* INLA: error handling, check non-default args, plots, 

* Test earth gam method more, compare against other methods to suggest advantages 

* Add single parameter methods for historical interest? 

* Standard errors for all nonparametric regression methods at least 

* Diagnostics for fit of nonparametric regression models 

* MLMC and QMC methods

* Standard 2 level Monte Carlo? [ not recommended for use, but as a convenience for methods development ] 

* Anything to facilitate algebraic methods?

* Any other plots for communication?

EVSI 

* Complete doc and unit tests for regression and IS methods 

* Add the other EVSI methods 

* Material from Anna's EVSI package 

GENERAL

* Package vignettes explaining merits of different methods - linking to our papers.  Integrate this with the work on the books. 

* Interface to heemod package.  Any other packages to interface with? 


## DISCUSSION 

* Does the baseline option matter for the purpose of regression-based EVPPI computation - can it always be treatment 1?   For consistency with what BCEA does, incremental for regressions is calculated as reference - active.  That seems backwards, but not sure it matters.

* Suggests or Imports for packages enabling specific methods?   Imports for default methods, suggests for others? 

* BCEA allows different regression methods for effects and costs, handle this? 
