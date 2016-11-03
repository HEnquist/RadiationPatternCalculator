

function [ sourceout ] = preparesource( sources ,fignum)
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

isoctave = (exist ('OCTAVE_VERSION', 'builtin') > 0); %running in matlab or octave?
if isoctave 
  scrsz = get(0,'screensize');
else
  scrsz = get(groot,'ScreenSize');
end


sourceout = sources{1};
c=340;
sourceout.lambda=c/sourceout.f;


stemp = sources{1};
sourceout.x = [];
sourceout.y = [];
sourceout.z = [];
%sourceout.r = [];
sourceout.ints = [];
sourceout.n= 0;
sourceout.dipole = [];
if ~isfield(sources{1},'stereo')
   	sourceout.stereo=false;
end

for m=1:length(sources)
    % complete missing parameters
    if ~(isfield(sources{m},'Ly')&& isfield(sources{m},'Lz'))
        if isfield(sources{m},'radius')
            sources{m}.Lz=2*sources{m}.radius;
            sources{m}.Ly=2*sources{m}.radius;
        else
            error('Incomplete source!')
        end
    end
    if ~isfield(sources{m},'radius')
        if (isfield(sources{m},'Ly')&& isfield(sources{m},'Lz'))
            sources{m}.radius=max([sources{m}.Ly sources{m}.Lz]);
        else
            error('Incomplete source!')
        end
    end
    if ~isfield(sources{m},'radcurv')
        sources{m}.radcurv=1000*sources{m}.radius;
    end
    if ~isfield(sources{m},'conedepth')
        sources{m}.conedepth=0;
    end
    if ~isfield(sources{m},'dipole')
        sources{m}.dipole=false;
    elseif sources{m}.dipole && ~isfield(sources{m},'dir')
        sources{m}.dir=[1;0;0];
    end
    if ~isfield(sources{m},'xpos')
        sources{m}.xpos=0;
    end
    if ~isfield(sources{m},'ypos')
        sources{m}.ypos=0;
    end
    if ~isfield(sources{m},'zpos')
        sources{m}.zpos=0;
    end
    if ~isfield(sources{m},'roty')
        sources{m}.roty=0;
    end
    if ~isfield(sources{m},'rotz')
        sources{m}.rotz=0;
    end
    if ~isfield(sources{m},'level')
        sources{m}.level=0;
    end
    if ~isfield(sources{m},'phase')
        sources{m}.phase=0;
    end
    if ~isfield(sources{m},'wavespeed')
        sources{m}.wavespeed=1e9;
    end
    if ~isfield(sources{m},'xo')
        sources{m}.xo=@(s) 1;
    end
    
    stemp.nz=round(sources{m}.Lz/sources{m}.dx); %number of source points, z
    tempz=linspace(-sources{m}.Lz/2,sources{m}.Lz/2,stemp.nz);
    stemp.ny=round(sources{m}.Ly/sources{m}.dx); %number of source points, y
    tempy=linspace(-sources{m}.Ly/2,sources{m}.Ly/2,stemp.ny);
    [stemp.z, stemp.y] = meshgrid(tempz,tempy);
    stemp.y=stemp.y(:);
    stemp.z=stemp.z(:);
    
    
    stemp.r=sqrt(stemp.y.^2 + stemp.z.^2); %radius
    stemp.x=sqrt(sources{m}.radcurv^2-stemp.r.^2)-sources{m}.radcurv + (stemp.r-sources{m}.radius)*sources{m}.conedepth/sources{m}.radius; %calculate depth
    stemp.ints = ones(size(stemp.y)); %source point intensities
    mask=stemp.r<sources{m}.radius; %circular source
    
    stemp.r=stemp.r(mask);
    stemp.x=stemp.x(mask);
    stemp.y=stemp.y(mask);
    stemp.z=stemp.z(mask);
    stemp.ints=stemp.ints(mask);
    Anorm = 1e5*stemp.ints/sum(stemp.ints);
    Alevel = 10^(sources{m}.level/20);
    Aphase = exp(1i*sources{m}.phase*pi/180);
    Awave = exp(-1i*2*pi*stemp.f*stemp.r/sources{m}.wavespeed);
    Axo = sources{m}.xo(1i*2*pi*stemp.f);
    stemp.ints = Anorm .* Aphase .* Alevel .* Awave .* Axo ;
    stemp.n=length(stemp.ints);
    stemp.dipole = logical(sources{m}.dipole*ones(size(stemp.ints)));
    
    %rotate and move the source
    if sources{m}.dipole
        stemp.dir=[cosd(sources{m}.roty) 0 sind(sources{m}.roty);
            0              1           0;
            -sind(sources{m}.roty) 0 cosd(sources{m}.roty)]*sources{m}.dir;
        
        stemp.dir=[cosd(sources{m}.rotz)  sind(sources{m}.rotz)  0;
            -sind(sources{m}.rotz)  cosd(sources{m}.rotz) 0;
            0              0           1]*sources{m}.dir;
        stemp.dir=stemp.dir/norm(stemp.dir); %normalize source direction vector
    end
    
    xs_rot=cosd(sources{m}.roty)*stemp.x + sind(sources{m}.roty)*stemp.z;
    ys_rot=stemp.y;
    zs_rot=-sind(sources{m}.roty)*stemp.x + cosd(sources{m}.roty)*stemp.z;
    
    xs2_rot=cosd(sources{m}.rotz)*xs_rot + sind(sources{m}.rotz)*ys_rot;
    zs2_rot=zs_rot;
    ys2_rot=-sind(sources{m}.rotz)*xs_rot + cosd(sources{m}.rotz)*ys_rot;
    
    stemp.x=xs2_rot+sources{m}.xpos;
    stemp.y=ys2_rot+sources{m}.ypos;
    stemp.z=zs2_rot+sources{m}.zpos;
    
    sourceout.x = [sourceout.x; stemp.x];
    sourceout.y = [sourceout.y; stemp.y];
    sourceout.z = [sourceout.z; stemp.z];
    %sourceout.r = [];
    sourceout.ints = [sourceout.ints; stemp.ints];
    sourceout.dipole = [sourceout.dipole; stemp.dipole];
    sourceout.n= sourceout.n + stemp.n;
    
end


%plot source points
if fignum>0
    h=figure(fignum);
    
    if isoctave
        set(h,'position',[1 2*scrsz(4)/3 scrsz(3)/3-30 scrsz(4)/3.2-90]);
    else
        set(h,'OuterPosition',[1 2*scrsz(4)/3 scrsz(3)/3 scrsz(4)/3.2]);
    end
    plot3(sourceout.x,sourceout.y,sourceout.z,'.')
    axis equal
    xlabel('X, m')
    ylabel('Y, m')
    zlabel('Z, m')
end
end

