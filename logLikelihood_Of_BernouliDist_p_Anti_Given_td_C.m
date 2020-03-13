%% to fit psychometric, here we construct the log-likelihood of Bernouli distribution of [p^a * (1-p)^(1-a)] which p=p(Anti|td, c)
function [logLikelihoodValue] = logLikelihood_Of_BernouliDist_p_Anti_Given_td_C(Input, parameters, parameter_strings, SubjObj_flag);

% Model parameters
[p_anti] = pr_anti_c_td(Input, parameters, parameter_strings, SubjObj_flag); % this function computes the Pr(Anti|c, td)
Anti = (-1*Input.PrAn +1)/2; % PrAn(-1: antisaccade, +1: prosaccade) => converted to ( 0: prosaccade, and 1: antisaccade ) for likelihood
TF = Input.TF;

% how to compute log-likelihood => LogLikelihood = log(p^a * (1-p)^(1-a)) = a * log(p) +
% (1-a) * log(1-a), where a is response and p is probability of response
% given observations
logLikelihoodValue = -sum(Anti.*log(p_anti+0.0001) + (1-Anti).*log(1 - p_anti+0.0001));
end
