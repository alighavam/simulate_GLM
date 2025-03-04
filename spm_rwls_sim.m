function SPM = spm_rwls_sim(job)
    % Ali's note: this function tries to match the 
    % spm_rwls_run_fmri_specs() function.

    SPM = []; % empty SPM structure. This is supposed to be similar to 
              % the SPM.mat
  
    % Ali's note: The SPM structure has many many values. The values that 
    % I am adding here are needed for this simulation.
    SPM.nscan = job.nscan;

    %-Timing parameters & Basis functions
    %==========================================================================
    %-Repeat time
    %--------------------------------------------------------------------------
    SPM.xY.RT = job.timing.RT;
    
    %-Basis function parameters
    %--------------------------------------------------------------------------
    SPM.xBF.UNITS = job.timing.units;
    SPM.xBF.T     = job.timing.fmri_t;
    SPM.xBF.T0    = job.timing.fmri_t0;
    
    %-Basis functions
    %--------------------------------------------------------------------------
    bf = char(fieldnames(job.bases));
    if strcmp(bf,'hrf') % if hrf was params were defiend in job:
        if all(job.bases.hrf.derivs == [0 0])
            SPM.xBF.name = 'hrf';
        elseif all(job.bases.hrf.derivs == [1 0])
            SPM.xBF.name = 'hrf (with time derivative)';
        elseif all(job.bases.hrf.derivs == [1 1])
            SPM.xBF.name = 'hrf (with time and dispersion derivatives)';
        else
            error('Unknown HRF derivative choices.');
        end
        if (isfield(job.bases.hrf,'params'))            % Hand over params for hrf: JD 2/2013
            SPM.xBF.params = job.bases.hrf.params;
        end
    else
        error("this simulation will not work without defining hrf params\n")
    end
    
    %-Model interactions (Volterra)
    %--------------------------------------------------------------------------
    SPM.xBF.Volterra = job.volt;

    %-Events
    %----------------------------------------------------------------------
    % Ali's note: This is defined like this for the sake of this 
    % simulation. SPM's definition of events is much more complete
    % and detailed. Check spm_rwls_run_fmri_spec() for the details
    % of how U is made and how it looks like in real SPM. 
    SPM.U = job.U;
end