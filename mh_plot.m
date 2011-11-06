function mh_plot(seq,evalTargetOne)
    [numOfPars,numOfDraws] = size(seq.pars);
    numOfBins = ceil(log2(numOfDraws)+1);
    minpar = min(seq.pars,[],2);
    maxpar = max(seq.pars,[],2);
    for p=1:numOfPars
        figure;
        [b,d,x] = kde(seq.pars(p,:),numOfBins,minpar(p),maxpar(p));
        plot(x,d,'-b',x,evalTargetOne(x,p),'-.b');
    end
end