%% modeling the switching probability.
function [Output_pr_of_switch, Output_tDev_lastOne, Output_RuleChoice_lastOne, Output_T, Output_SW, mu_switch_estimated] = pr_switch_func(Input, alpha_transition, sigma_switch, pam3,  psych_parameters, psych_parameter_strings, psych_SubjObj_flag, expectedAccuracy_Benchmark);

    % Setting the initialized values of optimizer
    options = optimset('fminsearch');
    options.Display = 'off'; % 'off'
    options.Iter = 1000000;
    options.TolFun = 1e-10;
    options.TolX = 1e-10;


% some reminders:
    % The array "Input" includes only the error trials (sorted 1B, 2B, ...)
    % variables of switch_model.alpha_transition, switch_model.sigma_switch
    
    % How to contruct the Posterior odds
        % for 1 error: H(T) / [(1-H(T))(1-A(T))]   % arrayOfInput => (T):Array(1)
        % for 2 error: ( H(T-1) + H(T)[1-H(T-1)][1-A(T-1)] ) / ( [1-H(T-1)][1-A(T-1)][1-H(T)][1-A(T)] )   % arrayOfInput => (T-1):Array(1)  ,   (T):Array(2)
        % for 3 errors: ( H(T-2) + H(T-1)[1-H(T-2)][1-A(T-2)] + H(T)[1-H(T-1)][1-A(T-1)][1-H(T-2)][1-A(T-2)] )   /   ( [1-H(T-2)][1-A(T-2)][1-H(T-1)][1-A(T-1)][1-H(T)][1-A(T)] )   % arrayOfInput => (T-2):Array(1)  ,   (T-1):Array(2)   ,   (T):Array(3)
        % for n errors: Q_IOM in supplementary



    
    % (this function p_of_cr is not used in the code):
    % just for reminder:  [p_anti] = pr_anti_c_td(Input, psych_parameters, psych_parameter_strings, psych_SubjObj_flag);
    % Here, I wrote the probability of correct trial
    p_of_cr = @(Input, psych_model, psych_model_strings, psych_model_SubjObj_flag)(     (Input.RuleChoice==0)*(Input.tDev<0)*(1-pr_anti_c_td(Input, psych_model, psych_model_strings, psych_model_SubjObj_flag)) +  ...
                                                                                        (Input.RuleChoice==0)*(Input.tDev>0)*(pr_anti_c_td(Input, psych_model, psych_model_strings, psych_model_SubjObj_flag)) +  ...
                                                                                        (Input.RuleChoice==0)*(Input.tDev==0)*0.5 +  ...
                                                                                        (Input.RuleChoice==1)*(Input.tDev<0)*(pr_anti_c_td(Input, psych_model, psych_model_strings, psych_model_SubjObj_flag)) +  ...
                                                                                        (Input.RuleChoice==1)*(Input.tDev>0)*(1-pr_anti_c_td(Input, psych_model, psych_model_strings, psych_model_SubjObj_flag)) +  ...
                                                                                        (Input.RuleChoice==1)*(Input.tDev==0)*0.5);
    
    
    
    % finding the probability of switch (cosidering gaussian distribution)
    pr_of_switch = @(x, mu_switch, sigma_switch, T ) (  (0.5 * (1-erf((x-mu_switch)/(sqrt(2)*(sigma_switch*sqrt(T)) )))) );
    %pr_of_switch = @(x, mu_switch, sigma_switch, T ) (  (0.5 * (1-erf((x-mu_switch)/(sqrt(2)*(sigma_switch ) )))) );
    
    % cost function for finding the mean of distribution
    loss_function_pr_of_switch = @(switch_bound, mu_switch, sigma_switch, posteriorRation, T) ( (pr_of_switch(switch_bound, mu_switch, sigma_switch, T) - (posteriorRation/(posteriorRation+1)) )^2 );
    
    
  for iTrial = 1: length(Input)
    T = Input(iTrial).T;
    for iT = 1: T % saving the history for each trial in a new variable
        syntheticInput(iT).RuleChoice = Input(iTrial).RuleChoice(iT);
        syntheticInput(iT).Cue = Input(iTrial).Cue(iT);
        syntheticInput(iT).tDev = Input(iTrial).tDev(iT);
        syntheticInput(iT).Anti = (-1*Input(iTrial).PrAn(iT) +1)/2; % PrAn(-1: antisaccade, +1: prosaccade) => converted to ( 0: prosaccade, and 1: antisaccade )
        syntheticInput(iT).TF = Input(iTrial).TF(iT);
        syntheticInput(iT).tdMean = Input(iTrial).tdMean;
        syntheticInput(iT).DevValues =Input(iTrial).DevValues;
        
        iCT_index = find( (expectedAccuracy_Benchmark.RuleChoice == Input(iTrial).RuleChoice(iT) ) .* (expectedAccuracy_Benchmark.tDev == Input(iTrial).tDev(iT) )  );
        A(iT) = expectedAccuracy_Benchmark.expectedAccuracy(iCT_index);
        H(iT) = alpha_transition;
        
    end
    
    StartPointInitializedValues = [0.2];
    switch_bound = 1;
    %option #1:
     [mu_switch_estimated(iTrial), FVAL(iTrial)] = fminsearch(@(mu)  loss_function_pr_of_switch(switch_bound, mu, sigma_switch, PO(A, H, T), T), StartPointInitializedValues, options);
    %option #2:
     % mu_switch_estimated(iTrial) = abs(switch_bound - sqrt(2) * sigma_switch * sqrt(T) * erfinv( (1 - PO(A, H, T)) ./ (1 + PO(A, H, T))  ) );
    
    % implementing "alpha * mu" in the simple model
    mu_switch_estimated(iTrial) = mu_switch_estimated(iTrial) * pam3;
    
    % to have a bound on the value of parameters:
    if sigma_switch < 4 
    Output_pr_of_switch(iTrial) = pr_of_switch(switch_bound, mu_switch_estimated(iTrial), sigma_switch, T);
    else
    Output_pr_of_switch(iTrial) = 0.001; % This would help so algorithm tries to not exceed larger values (>=4)
    end
    Output_tDev_lastOne(iTrial) = Input(iTrial).tDev(T);
    Output_RuleChoice_lastOne(iTrial) = Input(iTrial).RuleChoice(T);
    Output_T(iTrial) = T;
    if sigma_switch < 4
        Output_SW(iTrial) = Input(iTrial).SW;
    else
        Output_SW(iTrial) = 0; % This would make log-likelihood so small so that it helps the algorithm tries to not exceed larger values (<=4)
    end

    clear A; clear H; clear syntheticInput;
  end

 
end

