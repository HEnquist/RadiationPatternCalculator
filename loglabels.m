function [numbers, labels] = loglabels(start,stop)
%Radiation Pattern Calculator
%Copyright Henrik Enquist 2017
%This file is part of Radiation Pattern Calculator.
% 
%Radiation Pattern Calculator is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
% 
%Radiation Pattern Calculator is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
% 
%You should have received a copy of the GNU General Public License
%along with Radiation Pattern Calculator.  If not, see <http://www.gnu.org/licenses/>.

evenstart = floor(log10(start));
evenend = ceil(log10(stop));
nbrdecades = evenend-evenstart;

numbers = round(logspace(evenstart,evenend,3*nbrdecades+1),1,'significant');

if sum(numbers<=start)>1
    numbers = numbers(sum(numbers<=start):end);
end
if sum(numbers>=stop)>1
    numbers = numbers(1:(end-(sum(numbers>=stop)-1)));
end
for n = 1:length(numbers)
    labels{n} = num2sip(numbers(n));
end

