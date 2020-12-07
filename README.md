# OpenCLCalculus

These julia scripts are intended to help constructing, deriving and fitting mathematical formulas using OpenCL, mainly probability distribution functions.

The basic constituents are expressions:
  * sum reduction
  * product reduction
  * dot product
  * scalar division
  * scalar logarithm
  * scalar negation
  * parameteric Gauss, log Gauss, parameteric error function
  * logistic function
  * normalized double logistic pdf, its logarithm and its cumulative probability density
 
These can be arbitrarily combined and their derivatives calculated. Some higher level simplification is also applied.

For cases where known repetition occurs by a number unknown at compilation time, 
e.g. fitting a sum of Gaussians, two additional formulas are available:
  * parameteric weighted averege
  * derived paremeteric weighted average
  
One possible use case for these is the calculation of conditional probabilities and their derivatives.

The symbolic expressions and the parameteric formulas allow the automatic creation of julia and openCL code, 
and the code set contains regression and pdf estimation OpenCL kernels, wrapped into a NLOpt julia optimizer.
