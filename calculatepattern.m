function [ world ] = calculatepattern( plane,source,world,h)
%calculatepattern(plane,source,world): Calculates the radiation pattern
%Copyright Henrik Enquist 2016
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

for nn=1:source.n
    if h~=0
        waitbar(nn/source.n,h)
    end
    r_x=plane.xg-source.x(nn);
    r_y=plane.yg-source.y(nn);
    r_z=plane.zg-source.z(nn);
    r=sqrt(r_x.^2+r_y.^2+r_z.^2);
    if source.dipole(nn)
        alpha = acos((r_x.*source.dir(1) + r_y.*source.dir(2) + r_z.*source.dir(3))./r) ;
        world=world+source.ints(nn).*1i.*cos(alpha).*exp(1i*2*pi*r/source.lambda)./r;
    else
        world=world+source.ints(nn).*exp(1i*2*pi*r/source.lambda)./r;
    end
end
end

