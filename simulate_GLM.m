function X = simulate_GLM(model, nscan, onsets, hrf_params, varargin)
% INPUTS:
%       model: the model to simulate
%       nscan: number of fMRI volumes scanned (i.e. ,numTR, nscan,
%       numScans, etc.)
%       onsets: onset of the events in seconds. e.g., [10, 25, 40, 55, ...]
%       hrf_params: parameters for hrf, i.e., [p1, p2 p3 p4 p5 p6]. Check
%       spm_hrf() to see what these parameters are. You can use
%       spm_hrf(TR,hrf_params) to plot the hrf with spm.
% Case: 
%       'data': Uses SPM-like and rwls-like functions to generate simulated
%       BOLD signal. It adds a noise to the data.
%       'GLM': Same as 'data' case but generates GLM design matrix.

% example simulation:
add_noise = 0; 
mu = 0; % noise mean
sigma = 1; % noise std
event_dur = 10; % 10(ms) default Dirac delta event
vararginoptions(varargin, {'add_noise','mu','sigma','event_dur'})

switch model
    case 'simulate_GLM'
        ons = onsets;
        dur = event_dur*ones(size(ons))/1000;
        
        % ====================================== 
        global defaults; 
        if (isempty(defaults)) 
            spm_defaults;
        end
        
        % Empty job structure:
        J = [];
        
        % onets of the events:
        % Ali's note: U is defined like this for the sake of this 
                    % simulation. SPM's definition of events is much more complete
                    % and detailed. Check spm_rwls_run_fmri_spec() for the details
                    % of how U is made and how it looks like in real SPM. 
        for i = 1:size(ons,2)
            J.U(i).ons = ons(i,:);  
            J.U(i).dur = dur(i,:);
        end
        J.nscan = nscan;
        
        J.timing.units = 'secs';
        J.timing.RT = 1;
        
        % number of temporal bins in which the TR is divided,
        % defines the discrtization of the HRF inside each TR
        J.timing.fmri_t = 16;
        
        % slice number that corresponds to that acquired halfway in
        % each TR
        J.timing.fmri_t0 = 1;
        
        % Specify hrf parameters for convolution with
        % regressors
        J.bases.hrf.derivs = [0,0];
        J.bases.hrf.params = hrf_params;  % positive and negative peak of HRF
        defaults.stats.fmri.hrf = J.bases.hrf.params; 
        
        % Specify the order of the Volterra series expansion 
        % for modeling nonlinear interactions in the BOLD response
        % *Example Usage*: Most analyses use 1, assuming a linear
        % relationship between neural activity and the BOLD
        % signal.
        J.volt = 1;

        % SIMULATING DATA:
        SPM = spm_rwls_sim(J);
        X = spm_fMRI_design_sim(SPM);
        
        % Add noise:
        if add_noise
            X = X + mu + sigma * randn(size(X));
        end
end

