function [X,Xname,Fc] = spm_Volterra_sim(U,bf,V)
    % Ali's note: This is writtent based on spm_Volterra()
    % V option is not supported in this simulation. 
    
    %-1st order terms
    %==========================================================================
    X     = [];
    for i = 1:numel(U)
        for k = 1:size(U(i).u,2)
        for p = 1:size(bf,2)
            x = U(i).u(:,k);
            d = 1:length(x);
            x = conv(full(x),bf(:,p));
            x = x(d);
            X = [X x];
        end
        end
        
        % Ali's note: These outputs are there in the orignal function plus a
        % few other things but we don't need them here in this simulation.
        Fc = [];
        Xname = [];
    end
end
