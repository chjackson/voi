## DONE

* Basic package structure

* EVPPI methods copied in from BCEA/SAVI and tidied up

* A few unit tests added for these 

* Nonparametric regression and importance sampling methods added for EVSI, and rough tests written

* Study design for EVSI can be specified as an R function to simulate data from it, or as one of the built-in standard designs - currently only have 1-arm or 2 arm trials of a binary outcome. 

* Chemotherapy model code added to facilitate examples. Note exact format of this code may change as we implement different methods that need access to the decision model function. 

* Standard 2-level Monte Carlo 


## TODO


Doc for extra args in built in datagen fn eg sigma in normal_known.  Is it even handled? 

It's not in help(evsi)  neither is binary


EVPPI

* Different regression methods and options for effects and costs

* INLA: check plots work

* Test earth gam method more, compare against other methods to suggest advantages, while developing book material. 

* Standard errors for GP, earth.  INLA too?

* Diagnostics for fit of nonparametric regression models 

* MLMC and QMC methods

* Parallel processing for standard 2 level Monte Carlo, and partial results for convergence checking 

* Anything to facilitate algebraic methods?

* Any other plots for communication?




EVSI

* More built-in study designs - which ones? 

* Any need to handle designs controlled by more than one sample size?  If so, make sure multiple arguments to datagen_fn are handled nicely

* Guidance for users to specify appropriate GAM formulae for their applications when the default formula is too intensive.  Relatedly, appropriate number of PSA samples to use. 

* Add Heath EVSI method: Tricky bit is updating the parameters given new data and rerunning the model - code this more abstractly than is currently done.

* Add Jalal EVSI method

* 2-level Monte Carlo again

* Output analysis material from Anna's EVSI package.



GENERAL

* Interface to heemod package.  Any other packages to interface with?





## PRINCIPLES OF PACKAGE DEVELOPMENT


### Using

* Most important!

* Just use it any time you want to do a VoI calculation, and give feedback.
Describe what you had to do to make it work in your example.  If you need to do anything tedious, this may suggest how the package could be more helpful.

* Does other software do anything better?


### Design

* Decide what the package should do, and what is better done with other tools

* What functions should look like: argument formats, consistency between different parts of the code

* Identifying where we do the same task multiple times, therefore should have a function for that task


### Coding 

* Clean, modular and consistent style.  Each function does one thing that can be described concisely.

* Descriptive/concise variable and function names.

* Code should be understandable as much as possible by itself without the need for commenting.

* Comment in cases where it won't be instantly obvious what the code is doing. 

* Any time you add code, add a unit test - a concise example where that code is executed - to ensure that test gives the expected result every time the code is modified. 

* House style: 

	- underscores not dots to separate words in function or variable names:  my_variable_name

	- use spaces to clearly separate elements of code, e.g. 
       foo <- fn(a, b, c)  is much easier to read than  foo<-fn(a,b,c) 



### Documentation 

* Vignettes with worked examples

* Check can reproduce results from ConVoI 1 four case studies

* Any new ConVoI examples (book, papers) should come with code to do them with this package. 

* Note these may raise "further research is needed" questions!
