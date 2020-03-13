%% psychometric function Pr(Anti | rule, sample interval)
function [p_anti] = pr_anti_c_td(Input, parameters, parameter_strings, SubjObj_flag);
% Gaussian pdf
pdf_x = @(x, mu, sigma) ((1/sqrt(2*pi*(sigma^2))) .* exp(-((x - mu).^2) / (2*(sigma^2))) );
Pr_I_Given_C_td = @(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, x_min, x_max, x_resolution) ((C==0)*sum(x_resolution*pdf_x(pBound_PrAn:x_resolution:x_max, pDrift*td, pWm*td+pWm_offset)) + (C==1)*sum(x_resolution*pdf_x(x_min:x_resolution:pBound_AnPr, pDrift*td, pWm*td+pWm_offset)) );
Pr_Anti_Given_C_td = @(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate, x_min, x_max, x_resolution) ( (1-abs(pAlphaErr)) * ((1-lapseRate)*Pr_I_Given_C_td(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, x_min, x_max, x_resolution)   + 0.5*lapseRate)  * (lapseRate<=1 && lapseRate>=0) * (pAlphaErr<=1 && pAlphaErr>=0) * (pWm>0) * (pWm_offset>0) );


    for iParameters = 1: length(parameter_strings)
       eval(parameter_strings{iParameters}); 
    end
    local_parameters.pDrift = 1; % fixed

    for iTrial = 1: length(Input.tDev)
        if strcmp(SubjObj_flag, 'Subj')
        p_anti(iTrial) = Pr_Anti_Given_C_td(Input.RuleChoice(iTrial), 0.001*(Input.tDev(iTrial)+Input.tdMean), local_parameters.pDrift, local_parameters.pWm, local_parameters.pWm_offset, local_parameters.pBound_PrAn, local_parameters.pBound_AnPr, local_parameters.pAlphaErr, local_parameters.lapseRate, -1, 10, 0.01);
        % scalable model (uncommment for scalable model)
        % p_anti(iTrial) = Pr_Anti_Given_C_td(Input.RuleChoice(iTrial), 0.001*(Input.tDev(iTrial)+Input.tdMean), local_parameters.pDrift, local_parameters.pWm, 0, local_parameters.pBound_PrAn, local_parameters.pBound_AnPr, local_parameters.pAlphaErr, local_parameters.lapseRate, -1, 10, 0.01);
        elseif strcmp(SubjObj_flag, 'Obj')
        p_anti(iTrial) = Pr_Anti_Given_C_td(Input.Cue(iTrial), 0.001*(Input.tDev(iTrial)+Input.tdMean), local_parameters.pDrift, local_parameters.pWm, local_parameters.pWm_offset, local_parameters.pBound_PrAn, local_parameters.pBound_AnPr, local_parameters.pAlphaErr, local_parameters.lapseRate, -1, 10, 0.01);
        % scalable model (uncommment for scalable model)
        % p_anti(iTrial) = Pr_Anti_Given_C_td(Input.Cue(iTrial), 0.001*(Input.tDev(iTrial)+Input.tdMean), local_parameters.pDrift, local_parameters.pWm, 0, local_parameters.pBound_PrAn, local_parameters.pBound_AnPr, local_parameters.pAlphaErr, local_parameters.lapseRate, -1, 10, 0.01);
        end
    end
end
