# musseco

**musseco** stands for [Mu]lti-[S]tate [S]peciation and [E]xtinction [Co]alescent and
is an [R](https://www.r-project.org/) package to estimate parameters for 
a speciation and extinction coalescent process.


We have implemented the **BiSSECo** ([Bi]nary-[S]tate [S]peciation and [E]xtinction [Co]alescent) model in which we can estimate the within-host and between-host 
replicative fitness using a fixed time-scaled tree.



## Installation

```r
# You will need to install the R package devtools 
# (https://github.com/r-lib/devtools)

install.packages("devtools")
devtools::install_github("emvolz-phylodynamics/musseco")
```



## Tutorials

* We recommend that you read the [Get started](http://emvolz-phylodynamics.github.io/musseco/articles/musseco.html) to 
understand the implementation of BiSSECo.



## Author

musseco has been developed by [Erik Volz](https://profiles.imperial.ac.uk/e.volz)



## Additional information

musseco works on a fixed phylogenetic tree and therefore is **not** a program that 
will estimate the tree for you.

