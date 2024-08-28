function min_factor = find_odd_factor(min_number, max_number)
% Return minumum odd number from the range with whe minimum prime factor
%
% DESCRIPTION:
%     find_odd_factor loops through the given range of numbers and finds the
%     odd number with the smallest maximum prime factors. This allows suitable
%     grid sizes to be selected to maximise the speed of the FFT. Odd number is
%     necessary for the field calculation of the axial symmetric source to compare 
%     axial profiles. The prime factor is printed to the command line.
%    
%
% INPUTS:
%     min_number    - integer specifying the lower bound of values to test
%     max_number    - integer specifying the upper bound of values to test
%
% OUTPUT:
%     min_factor    - integer odd number from the given range with the minimum prime factor
%
%
% ABOUT:
%     author        - Alisa Krokhmal
%     last update   - 26th August 2024

% extract factors
nums = min_number:max_number;
odd_nums = nums(mod(nums,2)==1);

facs = zeros(size(odd_nums));
fac_max = facs;
min_odd = min(odd_nums);

for index = 1:length(odd_nums)
    facs(index) = length(factor(odd_nums(index)));
    fac_max(index) = max(factor(odd_nums(index)));
end

% compute best factors in range
disp('Prime factor is');
ind = min(fac_max);
disp(num2str(ind));

min_factor =  odd_nums(find(fac_max == min(fac_max), 1));


