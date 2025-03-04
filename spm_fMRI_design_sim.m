function X = spm_fMRI_design_sim(SPM)
    % Ali's note: this function tries to match the 
    % spm_fMRI_design() function.
    fMRI_T     = SPM.xBF.T;
    fMRI_T0    = SPM.xBF.T0;

    %-Time units, dt = time bin {secs}
    %--------------------------------------------------------------------------
    SPM.xBF.dt     = SPM.xY.RT/SPM.xBF.T;
    
    %-Get basis functions
    %--------------------------------------------------------------------------
    SPM.xBF        = spm_get_bf(SPM.xBF);
    
    %-Create convolved stimulus functions or inputs
    %======================================================================
    %-Get inputs, neuronal causes or stimulus functions U
    %----------------------------------------------------------------------
    U = spm_get_ons_sim(SPM);

    %-Convolve stimulus functions with basis functions
    %----------------------------------------------------------------------
    [X,~,~] = spm_Volterra_sim(U, SPM.xBF.bf, SPM.xBF.Volterra);

    %-Number of scans for this session
    %----------------------------------------------------------------------
    k = SPM.nscan;
    
    %-Resample regressors at acquisition times (32 bin offset)
    %----------------------------------------------------------------------
    if ~isempty(X)
        X = X((0:(k - 1))*fMRI_T + fMRI_T0 + 32,:);
    end
end