# Generate a string with all interactions of a certain degree, to be used in a GAM formula

Generate a string with all interactions of a certain degree, to be used
in a GAM formula

## Usage

``` r
all_interactions(x, degree = 2)
```

## Arguments

- x:

  Character vector of variable names

- degree:

  Maximum interaction degree

## Value

A string looking like the right hand side of a GAM formula with tensor
product interactions.

For example, if `x` is `c("x1","x2","x3")`, then
`all_interactions(x, degree=2)` should return

`"te(x1,x2) + te(x1,x3) + te(x1,x3)"`

## Examples

``` r
x <- c("x1","x2","x3")
all_interactions(x, 2)
#> [1] "te(x1,x2) + te(x1,x3) + te(x2,x3)"
```
