## PRINCIPLES OF PACKAGE DEVELOPMENT


### Using

* Most important!

* Just use it any time you want to do a VoI calculation, and give feedback.  Describe what you had to do to make it work in your example.  If you need to do anything tedious, this may suggest how the package could be more helpful.

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

	- underscores not dots to separate words in function or variable names:  `my_variable_name`

	- use spaces to clearly separate elements of code, e.g. 
       `x <- fn(a, b, c)`  is much easier to read than  `x<-fn(a,b,c)` 
