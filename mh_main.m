clear all;

% metropolis-hastings basics ...
numOfDraws = 100000;
numOfPars = 10;
numOfSeqs = 10;
burnIn = 0.2 * numOfDraws;
% ... and parameters
gammaGE = 1.0;
gammaRW = 0.2;

% true parameters
mu   = 10.0 * ones(numOfPars,1);
vars = 5.0 * ones(numOfPars,1);
var  = diag(vars);																						
ivar = eye(numOfPars) / var;																							
dvar = sqrt( det(2 * pi * var) );								
% of the target density
evalTarget    = @(x)   exp( - 0.5 * (x-mu)' * ivar * (x-mu)) / dvar;
evalTargetOne = @(x,p) exp( - (x-mu(p)) .^ 2 / (2 * vars(p) )) / sqrt(2 * pi * vars(p) );

% allocate sequences
seqs.pars = zeros(numOfPars,numOfDraws,numOfSeqs);
seqs.targ = zeros(        1,numOfDraws,numOfSeqs);
for s=1:numOfSeqs
    seqs.pars(:,1,s) = randn(numOfPars,1) .* sqrt(vars) + mu;
    seqs.targ(1,1,s) = evalTarget(seqs.pars(:,1,s));
end

% estimate
[seq,accRatio] = mh(seqs,gammaGE,gammaRW,evalTarget,burnIn);

% plot histogram (kernel density estimator)
mh_plot(seq,evalTargetOne);

% print some descriptive statistics
mh_print(seq,mu,vars)

