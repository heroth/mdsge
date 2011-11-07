%*****************************************************************************
% testnwSpGr: test & demo for nwSpGr
% integrate function using sparse grids and simulation
%*****************************************************************************

% dimensions:
D = 10;

% max. accuracy level (pol. exactness wil be 2k-1):
maxk = 4;

% integrand: some function that evaluates g(x), where x is (R times D) and g is (R times 1)
func = 'prod( exp(-(x/2).^2/2)/2/sqrt(2*pi), 2)';

% calculate result of integration between 0 and 1:
trueval=(.5*(1+erf(1./sqrt(2)/2))-.5).^D;
% note: if unknown, replace it with simulated value:
%x=rand(1e+6,D);trueval=mean(eval(func));

for k=2:maxk
    % sparse grids integration:
    [x w] = nwSpGr('KPU', D, k);
	g = eval(func);
	SGappr = g'*w;
    SGerror = abs(SGappr - trueval)/trueval;
    % simulation
    numnodes = length(w);
    sim = zeros(1000,1);
    for r=1:1000
        x = rand(numnodes,D);
        g = eval(func);
        sim(r) = mean(g);
    end
    Simerror = sqrt(mean((sim-trueval).^2))/trueval;
    fprintf('\nD=%2.0f, k=%2.0f (nodes=%4.0f): SG error=%8.5f, Sim. error=%8.5f', ...
        [D,k,numnodes,SGerror,Simerror] )
end
