%% Maximum-Likelihood estimation by fminsearch optimization.
function [EstimatedParameters, syntheticInput, parameter_strings] = ModelParameterOptimization(Input, SubjObj_flag);
    
    % Setting the initialized values of optimizer
    options = optimset('fminsearch');
    options.Display = 'iter'; % 'off'
    options.Iter = 1000000;
    options.TolFun = 1e-6;
    options.TolX = 1e-6;

    % parameters.pDrift_Initialization = 1;
    parameters.pWm_Initialization = 0.04;
    parameters.pWm_offset_Initialization = 0.04; % for evaluating how scalable the noisy sample interval measurement is, here and other parts were set to zero for final test of model as it was near zero.
    parameters.pBound_PrAn_Initialization = 0.5;
    parameters.pBound_AnPr_Initialization = 0.5;
    parameters.pAlphaErr_Initialization = 0.01;
    parameters.lapseRate_Initialization = 0.01;


    % initialization of parameters of psychometric function:
    StartPointInitializedValues = [ parameters.pWm_Initialization, parameters.pWm_offset_Initialization, parameters.pBound_PrAn_Initialization, parameters.pBound_AnPr_Initialization, parameters.pAlphaErr_Initialization, parameters.lapseRate_Initialization];
    parameter_strings = { 'local_parameters.pWm = parameters(1);', 'local_parameters.pWm_offset = parameters(2);', 'local_parameters.pBound_PrAn = parameters(3);', 'local_parameters.pBound_AnPr = parameters(4);', 'local_parameters.pAlphaErr = parameters(5);', 'local_parameters.lapseRate = parameters(6);'};
    % Scalable model (also change in the probalitity of antisaccade code if we uncomment this)
    % StartPointInitializedValues = [ parameters.pWm_Initialization, parameters.pBound_PrAn_Initialization, parameters.pBound_AnPr_Initialization, parameters.pAlphaErr_Initialization, parameters.lapseRate_Initialization];
    % parameter_strings = { 'local_parameters.pWm = parameters(1);', 'local_parameters.pBound_PrAn = parameters(3);', 'local_parameters.pBound_AnPr = parameters(4);', 'local_parameters.pAlphaErr = parameters(5);', 'local_parameters.lapseRate = parameters(6);'};
    
    
    
    % fit the psychometric model:
    [EstimatedParameters, FVAL] = fminsearch(@(parameters) logLikelihood_Of_BernouliDist_p_Anti_Given_td_C(Input, parameters, parameter_strings, SubjObj_flag), StartPointInitializedValues, options);
    
    % make synthetic data for ploting continuous psychometric function in
    % terms of continuous tDev and according to the fitted parameters "EstimatedParameters"
    for iCue = 0: 1
        flag = 0;
        for tDev = Input.DevValues(1): 1 : Input.DevValues(9)
            flag = flag+1;
            syntheticInput.RuleChoice = iCue;
            syntheticInput.Cue = iCue;
            syntheticInput.tDev = tDev;
            syntheticInput.tdMean = Input.tdMean;
            syntheticInput.DevValues = Input.DevValues;
            [syntheticInput.p_anti(iCue+1, flag)] = pr_anti_c_td(syntheticInput, EstimatedParameters, parameter_strings, SubjObj_flag);
            syntheticInput.x_axis(flag) = tDev;
        end
    end
    
end
