nwSpGr: Matlab code for the calculation of nodes and weights 
for multivariate integration on sparse grids
(c) by Florian Heiss & Viktor Winschel, April 4, 2006

See "Likelihood Approximation by Numerical Integration on Sparse Grids" 
for details on the method.

Installation:
1) unpack "nwSpGr.zip" to some path
2) add this path with subfolders to the search path of Matlab
   * In Matlab, run "addpath(genpath('directory'))", where
     "directory" is the place where you unzipped the files
   * or add the path (with subdirectories) in the "File" => "Set
     path" dialog

The zip file contains:
* nwSpGr.m:       the code for the calculation of nodes and weights
* nwSpGr_demo.m:  code for demonstration and testing
* NodesWeights:   directory with precalculated univariate rules. Used
                  by nwSpGr.m

The syntax is
[n w] = nwSpGr(string scalar type, real scalar dim, real scalar k)

where 
dim  = dimension of the integration problem
k    = Accuracy level. The rule will be exact for polynomial up to
       total order 2k-1
type = String for type of 1D integration rule:
       The code can be used for any underlying univariate quadrature
       rule. Four rules are readily implemented:
       "KPU": Delayed Kronrod-Patterson for unweighted integral over
              [0,1]
       "KPN": Delayed Kronrod-Patterson for integral with Gaussian
              weight
       "GQU": Gaussian quadrature for unweighted integral over [0,1]
              (Gauss-Legendre)
       "GQN": Gaussian quadrature for integral with Gaussian weight
              (Gauss-Hermite)
       Any other rule can be defined. A function (call it func) has
       to be defined (and added to the Matlab search path). Its
       syntax has to be: "[n1d w1d] = func(level)"
       where n1d and w1d are univariate nodes and weights and level
       is the accuracy level. If this rule is exact for all univar.
       polynomials of order 2level-1, then the multivariate rule has
       exactness 2k-1 according to theorem 1 in the paper.
       Example: If you have installed the CompEcon Toolbox of Paul 
              Fackler and Mario Miranda defining e.g. "qnwlogn" for
              integration with a lognormal weight function 
              (http://www4.ncsu.edu/~pfackler/compecon/toolbox.html)
              you can use nwSpGr('qnwlogn', dim, k). 
Output:
n    = matrix of nodes with dim columns 
w    = row vector of corresponding weights
