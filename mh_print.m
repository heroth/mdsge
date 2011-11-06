function mh_print(seq,mu,vars)
    estmu  = mean(seq.pars,2);
    estvar =  var(seq.pars,[],2);
    rermu  = 100 * abs((estmu - mu) ./ mu);
    rervar = 100 * abs((estvar - vars) ./ vars);
    
    fprintf('\n');
    fprintf('estimated parameters from %i draws\n',size(seq.pars,2));
    fprintf('--------------------------------------\n');
    fprintf('est. par | true par | rel.err || est. var | true var | rel.err | \n');
    fprintf('----------------------------------------------------------------\n');
    for p=1:length(estmu)
        fprintf('%9.4f|%10.4f|%7.2f %%||%10.4f|%10.4f|%7.2f %%|\n',...
        estmu(p),mu(p),rermu(p),estvar(p),vars(p),rervar(p))
    end 
    fprintf('\n');
end    
