function [seq,accRatio] = mh(seqs,gammaGE,gammaRW,evalTarget,burnIn)
    [numOfPars,numOfDraws,numOfSeqs] = size(seqs.pars);
    accepted = 0;
    for d=2:numOfDraws
        print_progress_bar(d,numOfDraws);
        for s=1:numOfSeqs
            [s1,s2] = randSeqs(s,numOfSeqs);
            candidate =   seqs.pars(:,d-1,s) ...
                        + gammaGE * (seqs.pars(:,d-1,s1) - seqs.pars(:,d-1,s2)) ...
                        + gammaRW * randn(numOfPars,1);            
            targetVal = evalTarget(candidate);
            if  rand(1,1) < targetVal / seqs.targ(1,d-1,s); 
                seqs.pars(:,d,s) = candidate;
                seqs.targ(:,d,s) = targetVal;
                accepted = accepted + 1;
            else
                seqs.pars(:,d,s) = seqs.pars(:,d-1,s);
                seqs.targ(:,d,s) = seqs.targ(:,d-1,s);
            end
        end
    end
    accRatio = accepted / (numOfDraws * numOfSeqs);
    seq.pars = seqs.pars(:,burnIn:end,:);
    seq.targ = seqs.targ(1,burnIn:end,:);
    seq.pars = seq.pars(:,:);
    seq.pars = seq.pars(:,:);    
end

function print_progress_bar(d,numOfDraws)
    switch d
        case 2 
            fprintf('simulate');
        case numOfDraws       
            fprintf('end\n');
        otherwise
            if mod(mod(d,numOfDraws)/numOfDraws,0.1)==0
                fprintf('%s','.');
            end
                
    end
end

function [s1, s2] = randSeqs(s,numOfSeqs)
    s1 = s;
    s2 = s;
    while s1 == s
        s1 = floor(rand(1) * numOfSeqs) + 1;
    end
    while s2 == s || s2 == s1
        s2 = floor(rand(1) * numOfSeqs) + 1;
    end
end
