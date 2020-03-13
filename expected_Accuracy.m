%% Expected accuracy function
function [A] = expected_Accuracy(Input, parameters, parameter_strings, SubjObj_flag);


% These equations are for measuring the expected accuracy. please refer to the supplementary and also Braden et al PNAS 2016 paper



% Gaussian pdf
pdf_x = @(x, mu, sigma) ((1/sqrt(2*pi*(sigma^2))) .* exp(-((x - mu).^2) ./ (2*(sigma^2))) );
Pr_I_Given_C_td = @(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, x_min, x_max, x_resolution) ((C==0)*sum(x_resolution*pdf_x(pBound_PrAn:x_resolution:x_max, pDrift*td, pWm*td+pWm_offset)) + (C==1)*sum(x_resolution*pdf_x(x_min:x_resolution:pBound_AnPr, pDrift*td, pWm*td+pWm_offset)) );
Pr_Anti_Given_C_td = @(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate, x_min, x_max, x_resolution) (  ((1-lapseRate)*((1-abs(pAlphaErr))*Pr_I_Given_C_td(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, x_min, x_max, x_resolution) )  + 0.5*lapseRate)  * (lapseRate<=1 && lapseRate>=0) * (pAlphaErr<=1 && pAlphaErr>=0) );


td_DevValues = 0.001*([-320, -160, -80, -40, 0, 40, 80, 160, 320] + 850); % range of sample intervals converted to unit of seconds.


Pr_tdm_given_I_C_td = @(tdm, reward, C, td, td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate) (...
                                     (C==0).* (reward==1 || reward==2).*(td>td_DevValues(5))  .*  ( (1-abs(pAlphaErr)).*(tdm>pBound_PrAn).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)  ) +...
                                     (C==0).* (reward==1 || reward==2).*(td<td_DevValues(5))  .*  ( (1-abs(pAlphaErr)).*(tdm<=pBound_PrAn).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)  +   abs(pAlphaErr).*(tdm>pBound_PrAn).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)  ) +...
                                     (C==0).* (reward==1 || reward==2).*(td==td_DevValues(5))  .* (0.5*( (1-abs(pAlphaErr)).*(tdm>pBound_PrAn).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)  )   + 0.5*( (1-abs(pAlphaErr)).*(tdm<=pBound_PrAn).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)  +   abs(pAlphaErr).*(tdm>pBound_PrAn).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)  )  ) +...
                                     (C==1).* (reward==1 || reward==2).*(td<td_DevValues(5))  .*  ( (1-abs(pAlphaErr)).*(tdm<=pBound_AnPr).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)  ) +...
                                     (C==1).* (reward==1 || reward==2).*(td>td_DevValues(5))  .*  ( (1-abs(pAlphaErr)).*(tdm>pBound_AnPr).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)  +   abs(pAlphaErr).*(tdm<=pBound_AnPr).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)  ) +...
                                     (C==1).* (reward==1 || reward==2).*(td==td_DevValues(5))  .* ( 0.5*( (1-abs(pAlphaErr)).*(tdm<=pBound_AnPr).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)  )   + 0.5*( (1-abs(pAlphaErr)).*(tdm>pBound_AnPr).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)  +   abs(pAlphaErr).*(tdm<=pBound_AnPr).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)  )   ) ...
                                     +...
                                     (C==0).* (reward==2).*(td>td_DevValues(5))  .*  ( (1-abs(pAlphaErr)).*(tdm<=pBound_PrAn).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)  +   abs(pAlphaErr) ) + ...
                                     (C==0).* (reward==2).*(td<td_DevValues(5))  .*  ( (tdm>pBound_PrAn).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate) ) + ...
                                     (C==0).* (reward==2).*(td==td_DevValues(5))  .*  ( 0.5*( (1-abs(pAlphaErr)).*(tdm<=pBound_PrAn).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)  +   abs(pAlphaErr) )  + 0.5*( (tdm>pBound_PrAn).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate) )      ) +...
                                     (C==1).* (reward==2).*(td<td_DevValues(5))  .*  ( (1-abs(pAlphaErr)).*(tdm>pBound_AnPr).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate) +   abs(pAlphaErr) ) + ...
                                     (C==1).* (reward==2).*(td>td_DevValues(5))  .*  ( (tdm<=pBound_AnPr).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)   ) +...
                                     (C==1).* (reward==2).*(td==td_DevValues(5))  .* (0.5*( (1-abs(pAlphaErr)).*(tdm>pBound_AnPr).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate) +   abs(pAlphaErr) )  +0.5*( (tdm<=pBound_AnPr).*((1-lapseRate).*pdf_x(tdm, pDrift*td, pWm*(td)+pWm_offset) + 0.5*lapseRate)   )) ...
                                 );

                          
Pr_F_given_C_tdm = @(tdm, C, td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate ) (...
                        ( (( Pr_tdm_given_I_C_td(tdm, 1, C, td_DevValues(1), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 1, C, td_DevValues(2), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 1, C, td_DevValues(3), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 1, C, td_DevValues(4), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 1, C, td_DevValues(5), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 1, C, td_DevValues(6), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 1, C, td_DevValues(7), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 1, C, td_DevValues(8), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 1, C, td_DevValues(9), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9 )) )  ...
                            ./ ...
                        ( (( Pr_tdm_given_I_C_td(tdm, 2, C, td_DevValues(1), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 2, C, td_DevValues(2), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 2, C, td_DevValues(3), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 2, C, td_DevValues(4), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 2, C, td_DevValues(5), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 2, C, td_DevValues(6), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 2, C, td_DevValues(7), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 2, C, td_DevValues(8), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9  +  Pr_tdm_given_I_C_td(tdm, 2, C, td_DevValues(9), td_DevValues, pWm, pWm_offset, pDrift, pBound_PrAn, pBound_AnPr, pAlphaErr, lapseRate)/9 )) ) );

    for iParameters = 1: length(parameter_strings)
       eval(parameter_strings{iParameters}); 
    end
    local_parameters.pDrift = 1; % The drift variable is no longer a free parameter in the psychometric function optimization and it's set to 1 everywhere.

    td_min = -4; td_resolution = 0.01; td_max = 10;
    for iTrial = 1: length(Input.tDev) % integrating over tm (measured sample interval). 
        tdm = td_min:td_resolution:td_max; 
        
        for itdm = 1: length(tdm)
           A_1(itdm) = Pr_F_given_C_tdm(tdm(itdm), Input.RuleChoice(iTrial), td_DevValues, local_parameters.pWm, local_parameters.pWm_offset, local_parameters.pDrift, local_parameters.pBound_PrAn, local_parameters.pBound_AnPr, local_parameters.pAlphaErr, local_parameters.lapseRate) ;
           A_2(itdm) = pdf_x(tdm(itdm),    local_parameters.pDrift * 0.001*(Input.tDev(iTrial)+Input.tdMean),   local_parameters.pWm * (0.001*(Input.tDev(iTrial)+Input.tdMean)) + local_parameters.pWm_offset );
        end
        A(iTrial) = sum(A_1 .* A_2)*td_resolution;
    end

   
end
