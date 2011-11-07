%*********************************************************************************************************************
% nwSpGr: Nodes and weights for numerical integration on sparse grids (Smolyak)
% (c) 2006 Florian Heiss, Viktor Winschel, Feb 24, 2006
%*********************************************************************************************************************
% nwSpGr(): main function for generating nodes & weights for sparse grids intergration 
% Syntax: 
%    [n w] = nwSpGr(string scalar type, real scalar dim, real scalar k)
% Input: 
%    type = String for type of 1D integration rule:
%           "KPU": Delayed Kronrod-Patterson for unweighted integral over [0,1]
%           "KPN": Delayed Kronrod-Patterson for integral with Gaussian weight
%           "GQU": Gaussian quadrature for unweighted integral over [0,1] (Gauss-Legendre)
%           "GQN": Gaussian quadrature for integral with Gaussian weight (Gauss-Hermite)
%           func:  any function name. Function must accept level l and
%                  return nodes n and weights w for univariate quadrature rule with
%                  polynomial exactness 2l-1 as [n w] = feval(func,level)
%                  Example: If you have installed the CompEcon Toolbox of Paul Fackler and Mario Miranda
%                           (http://www4.ncsu.edu/~pfackler/compecon/toolbox.html)
%                           and you want to integrate with a lognormal weight function, you can use
%                           nwSpGr('qnwlogn', dim, k). 
%    dim  = dimension of the integration problem
%    k    = Accuracy level. The rule will be exact for polynomial up to total order 2k-1
% Output:
%    n    = matrix of nodes with dim columns 
%    w    = row vector of corresponding weights
%*********************************************************************************************************************
function [nodes, weights] = nwSpGr(type, dim, k)

  switch type
      case {'KPU','KPN','GQU','GQN'}
          if k<=25
              nw1d = nwload(type, k);
          else    
              error([type 'implemented only up to accuracy level 25']);
          end
      otherwise 
          try
              nw1d = nwcalc(type, k);
          catch
              error(['Unknown 1D integration rule or error in function']);
          end
  end              

  % initialization
  minq = max(0,k-dim);
  maxq = k-1;
  nodes = [];
  weights = [];

  % outer loop over q
  for q = minq:maxq
    bq = (-1)^(maxq-q) * nchoosek(dim-1,dim+q-k);
    % matrix of all rowvectors in N^D_{q}
    is = SpGrGetSeq(dim,dim+q);
    % inner loop collecting product rules
    for j=1:size(is,1)
      newnw = SpGrKronProd(nw1d(is(j,:),:));
      nodes = [nodes ; newnw(:,1:dim)];
      weights = [weights ; bq .* newnw(:,dim+1)];
    end
  end
  
  % collect equal nodes: first sort
  [nodes sortvec] = sortrows(nodes);
  weights = weights(sortvec);
  keep = 1; 
  lastkeep = 1;
  % then make list of rows to keep and sum weights of equal nodes
  for j=2:size(nodes,1)
    if nodes(j,:)==nodes(j-1,:) 
        weights(lastkeep)=weights(lastkeep)+weights(j);
    else
      lastkeep = j;
      keep = [keep ; j ];
    end
  end
  % return matrix of nodes & weights (normalization to account for rounding errors)*/
  nodes = nodes(keep,:);
  weights = weights(keep) / sum(weights(keep));
end

%**************************************************************************************
%SpGrGetSeq(): function for generating matrix of all rowvectors in N^D_{norm} 
%Syntax: 
%    out = nwSpGr(d,norm)
%Input: 
%    d    = dimension, will be #columns in output
%    norm = row sum of elements each of the rows has to have
%Output:
%    out  = matrix with d columns. Each row represents one vector with all elements >=1
%           and the sum of elements == norm
%**************************************************************************************
function fs = SpGrGetSeq(d, norm)
  seq = zeros(1,d);
  a=norm-d;
  seq(1)=a;
  fs = seq;
  c=1;
  while seq(d)<a
    if (c==d) 
        for i=(c-1):-1:1
            c=i;
            if seq(i)~=0, break, end;
        end
    end
    seq(c) = seq(c)-1;
    c=c+1;
    seq(c) = a - sum(seq(1:(c-1)));
    if (c<d) 
        seq((c+1):d)=zeros(1,d-c);
    end
    fs = [fs;seq];
  end 
  fs = fs+1;
end

%************************************************************************************
%SpGrKronProd(): function for generating tensor product quadrature rule 
%Syntax: 
%    out  = SpGrKronProd(pointer matrix nw1d)
%Input: 
%    nw1d = Dx2 dimensionl matrix of pointers to matrices of nodes (column 1) 
%           and weights (column 2) for univariate integration
%Output:
%    out  = matrix with d+1 columns. First D columns: nodes, last column: weights
%*************************************************************************************

function nw = SpGrKronProd(nw1d)
  nodes = nw1d{1,1} ; weights = nw1d{1,2};
  for j=2:size(nw1d,1)
    newnodes = nw1d{j,1};
    nodes = [kron(nodes,ones(size(newnodes,1),1)) kron(ones(size(nodes,1),1),newnodes)];
    weights = kron(weights,nw1d{j,2});
  end
  nw = [nodes weights];
end


%************************************************************************************
%nwload(): Read 1D nodes & weights and store pointers to them
%Syntax: 
%    out  = nwGQ(string scalar filen, real scalar maxlevel)
%Input: 
%    filen    = Start of filename of precalculated Gaussian nodes & weights (mata file)
%    maxlevel = maximum accuracy level
%Output:
%    out      = maxlevelx2 matrix of pointers to nodes (column 1) and weights (col. 2)
%************************************************************************************
function nw = nwload(filen, maxlevel)
  nodes = cell(0,0);
  weights = cell(0,0);
  for level=1:maxlevel
    nw = load(['NodesWeights/nw' filen '_' num2str(level) '.asc']);
    nodes = [nodes ; {nw(:,1)}];
    weights = [weights ; {nw(:,2)}];
  end
  nw = [nodes weights];
end

%************************************************************************************
%calc(): Calculate 1D nodes & weights using supplied function and store pointers to them
%Syntax: 
%    out  = nwGQ(string scalar filen, real scalar maxlevel)
%Input: 
%    filen    = Start of filename of precalculated Gaussian nodes & weights (mata file)
%    maxlevel = maximum accuracy level
%Output:
%    out      = maxlevelx2 matrix of pointers to nodes (column 1) and weights (col. 2)
%************************************************************************************
function nw = nwcalc(func, maxlevel)
  nodes = cell(0,0);
  weights = cell(0,0);
  for level=1:maxlevel
    [n w] = feval(func,level);
    nodes = [nodes ;  {n}];
    weights = [weights ; {w}];
  end
  nw = [nodes weights];
end
