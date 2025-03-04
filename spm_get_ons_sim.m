function U = spm_get_ons_sim(SPM)
    % Ali's note: This is written based on spm_get_ons() function.
    TR = 1;

    k     = SPM.nscan;
    T     = SPM.xBF.T;
    dt    = SPM.xBF.dt;
    
    U = SPM.U;

    %-Onsets
    %----------------------------------------------------------------------
    ons   = U.ons(:);
    
    %-Durations
    %----------------------------------------------------------------------
    dur   = U.dur(:);
    
    u     = ons.^0;
    
    %-Create stimulus functions (32 bin offset)
    %======================================================================
    ton       = round(ons*TR/dt) + 33;               % onsets
    tof       = round(dur*TR/dt) + ton + 1;          % offset
    sf        = sparse((k*T + 128),size(u,2));
    ton       = max(ton,1);
    tof       = max(tof,1);
    for j = 1:length(ton)
        if size(sf,1) > ton(j)
            sf(ton(j),:) = sf(ton(j),:) + u(j,:);
        end
        if size(sf,1) > tof(j)
            sf(tof(j),:) = sf(tof(j),:) - u(j,:);
        end
    end
    sf        = cumsum(sf);                         % integrate
    sf        = sf(1:(k*T + 32),:);                 % stimulus
    
    %-Place in ouputs structure U
    %----------------------------------------------------------------------
    U.dt   = dt;         % - time bin {seconds}
    U.u    = sf;         % - stimulus function matrix
end
