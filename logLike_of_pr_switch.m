  
%% loglikelihood function
function [logLikeSwitch] = logLike_of_pr_switch(Input, alpha_transition, sigma_switch, pam3, psych_parameters, psych_parameter_strings, psych_SubjObj_flag, expectedAccuracy_Benchmark, mOfT);

[Output_pr_of_switch, Output_tDev_lastOne, Output_RuleChoice_lastOne, Output_T, Output_SW, mu_switch_estimated] = pr_switch_func(Input, alpha_transition, sigma_switch, pam3, psych_parameters, psych_parameter_strings, psych_SubjObj_flag, expectedAccuracy_Benchmark);

    % option #1:
        %index = find( (Output_T == 1) + (Output_T == 2)); %optimize for all 1-back and 2-back trials as there are majority of error trials which reported in paper (%97)
        %logLikeSwitch = - sum(  Output_SW(index) .* log(Output_pr_of_switch(index)+0.0001) + (1-Output_SW(index)) .* log(1-Output_pr_of_switch(index)+0.0001)  ) / length(Output_SW(index));
    
    if mOfT == 1        % the option for optimizing only for 1-back trials
        index = find( (Output_T == 1));
        logLikeSwitch = - sum(  Output_SW(index) .* log(Output_pr_of_switch(index)+0.0001) + (1-Output_SW(index)) .* log(1-Output_pr_of_switch(index)+0.0001)  ) / length(Output_SW(index));
    elseif mOfT == 2    % the option for optimizing only for 2-back trials
        index = find( (Output_T == 2));
        logLikeSwitch = - sum(  Output_SW(index) .* log(Output_pr_of_switch(index)+0.0001) + (1-Output_SW(index)) .* log(1-Output_pr_of_switch(index)+0.0001)  ) / length(Output_SW(index));
    elseif mOfT == 3    % the option for optimizing 1-back & 2-back trials as there are the majority of error trials
        index = find( (Output_T == 1) + (Output_T == 2));
        logLikeSwitch = - sum(  Output_SW(index) .* log(Output_pr_of_switch(index)+0.0001) + (1-Output_SW(index)) .* log(1-Output_pr_of_switch(index)+0.0001)  ) / length(Output_SW(index));
    end

end
