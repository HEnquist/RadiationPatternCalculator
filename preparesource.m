

function [ source ] = preparesource( source ,fignum)
%preparesource( source ): Calculates the parameters needed in the source structure
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

scrsz = get(groot,'ScreenSize');

source.nz=round(source.Lz/source.dx); %number of source points, z
tempz=linspace(-source.Lz/2,source.Lz/2,source.nz);
source.ny=round(source.Ly/source.dx); %number of source points, y
tempy=linspace(-source.Ly/2,source.Ly/2,source.ny);
[source.z, source.y] = meshgrid(tempz,tempy);
source.y=source.y(:);
source.z=source.z(:);


source.r=sqrt(source.y.^2 + source.z.^2); %radius
source.x=sqrt(source.radcurv^2-source.r.^2)-source.radcurv + (source.r-source.radius)*source.conedepth/source.radius; %calculate depth
source.ints = ones(size(source.y)); %source point intensities
mask=source.r<source.radius; %circular source

source.r=source.r(mask);
source.x=source.x(mask);
source.y=source.y(mask);
source.z=source.z(mask);
source.ints=source.ints(mask);
source.ints=1e5*source.ints/sum(source.ints);
source.n=length(source.ints);

%rotate and move the source
if source.dipole
    source.dir=[cosd(source.rotz)  sind(source.rotz)  0;
        -sind(source.rotz)  cosd(source.rotz) 0;
        0              0           1]*source.dir;
    source.dir=source.dir/norm(source.dir); %normalize source direction vector
end
xs_rot=cosd(source.rotz)*source.x + sind(source.rotz)*source.y;
zs_rot=source.z;
ys_rot=-sind(source.rotz)*source.x + cosd(source.rotz)*source.y;

source.x=xs_rot+source.xpos;
source.y=ys_rot+source.ypos;
source.z=zs_rot+source.zpos;


%plot source points
h=figure(fignum);
set(h,'OuterPosition',[1 2*scrsz(4)/3 scrsz(3)/3 scrsz(4)/3.2]);
plot3(source.x,source.y,source.z,'.')
axis equal
xlabel('X, m')
ylabel('Y, m')
zlabel('Z, m')
end

