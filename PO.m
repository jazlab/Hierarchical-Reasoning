%% posterior odds in ideal observer
function [output] = PO(A, H, T);
% equation of posterior odds in ideal observer, which is function of expected accuracy, and the number of previous errors.
    if T == 1
        output = ( abs(H(T)) / (abs(1-H(T))*abs(1-A(T))) );
    elseif T == 2
        output = ( abs(H(T-1)) + abs(H(T))*abs(1-H(T-1))*abs(1-A(T-1)) ) / ( abs(1-H(T-1))*abs(1-A(T-1))*abs(1-H(T))*abs(1-A(T)) );
    elseif T == 3
        output = ( abs(H(T-2)) + abs(H(T-1))*abs(1-H(T-2))*abs(1-A(T-2)) + H(T)*abs(1-H(T-1))*abs(1-A(T-1))*abs(1-H(T-2))*abs(1-A(T-2)) ) / ( abs(1-H(T-2))*abs(1-A(T-2))*abs(1-H(T-1))*abs(1-A(T-1))*abs(1-H(T))*abs(1-A(T)) );
    elseif T>3
        Numinator_output = abs(H(1));
        for iT = 2: T
            weight_help = abs(H(iT));
            for j = 1: iT-1
                weight_help = weight_help * ( abs(1-H(j))*abs(1-A(j)) );
            end
            Numinator_output = Numinator_output + weight_help;
        end
        
        Denominator_output = 1;
        for j = 1: iT-1
            Denominator_output = Denominator_output * ( abs(1-H(j))*abs(1-A(j)) );
        end
        output = Numinator_output ./ Denominator_output;
        
    end
end
