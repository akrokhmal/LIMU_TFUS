function min_factor = find_min_factor(min_number, max_number)
% Return minumum number from the range with whe minimum prime factor
%
% DESCRIPTION:
%     find_odd_factor loops through the given range of numbers and finds the
%     odd number with the smallest maximum prime factors. This allows suitable
%     grid sizes to be selected to maximise the speed of the FFT. 
%     The prime factor is printed to the command line.
%    
%
% INPUTS:
%     min_number    - integer specifying the lower bound of values to test
%     max_number    - integer specifying the upper bound of values to test
%
% OUTPUT:
%     min_factor    - integer number from the given range with the minimum prime factor
%
%
% ABOUT:
%     author        - Alisa Krokhmal
%     last update   - 26th August 2024

%
% extract factors
facs = zeros(1, max_number - min_number);
fac_max = facs;
for index = min_number:max_number
    facs(index - min_number + 1) = length(factor(index));
    fac_max(index - min_number + 1) = max(factor(index));
end

% compute best factors in range
disp('Prime factor is');
ind = min(fac_max);
disp(num2str(ind));

min_factor =  min_number + find(fac_max == min(fac_max), 1) - 1;


