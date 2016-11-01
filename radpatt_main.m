%Radiation Pattern Calculator
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

clear all
close all

isoctave = (exist ('OCTAVE_VERSION', 'builtin') > 0); %running in matlab or octave?
if isoctave 
  scrsz = get(0,'screensize');
else
  scrsz = get(groot,'ScreenSize');
end

%% define the room
room = struct;
sweetspot = struct;

room.dx = 10e-3; %pixelsize in meters
room.lx = 5.1; %room length, m
room.lz = 2.4; %room height, m
room.ly = 3.2; %room width, m

sweetspot.l = 200e-3; %sweetspot size
sweetspot.n = 10; %sweetspot number of points per side
sweetspot.xp = 3.6; %sweetspot distance from front wall
sweetspot.yp = 0; %sweetspot sideways position
sweetspot.zp = 0.0; %sweetspot vertical position



%% define the source

% --Definition of the source parameters--
% sources{} - array of sources

% Mandatory
% sources{}.f  : frequency, Hz
% sources{}.dx : source point spacing, m
% sources{}.Lz : source height, m, must be defined for rectangular sources
% sources{}.Ly : source width, m, must be defined for rectangular sources
% sources{}.radius : source radius, m, must be defined for round sources

% Optional, if left out these are disabled or put to zero
% sources{}.radcurv : source radius of curvature, m
% sources{}.conedepth : source cone depth
% sources{}.dir : [x;y;z] vector for lobe direction before rotation by rotz, only used when dipole=true
% sources{}.dipole : false for monopoles, true for dipoles
% sources{}.xpos : source distance from front wall
% sources{}.zpos : vertical position
% sources{}.ypos : sideways position
% sources{}.roty : rotation around the y axis (tilt up/down), deg
% sources{}.rotz : rotation around the vertical axis (toe-in), deg
% sources{}.level : source level adjust, dB
% sources{}.phase : source phase, only relevant when using more than one source
% sources{}.wavespeed : propagation speed in the cone, outwards from voice coil, m/s
% sources{}.stereo : true to mirror source sideways for stereo, needs ypos>0 to make sense
%          only used for horizontal and back wall patterns

%example: dipole line source 300x20mm
% sources{1}.f=5000;
% sources{1}.dx=5e-3;
% sources{1}.Lz=0.3;
% sources{1}.Ly=0.02;
% sources{1}.dir=[1;0;0];
% sources{1}.dipole = true;
% sources{1}.xpos = 1;
% sources{1}.ypos = 0;
% sources{1}.zpos = 0;
% sources{1}.roty = 0;
% sources{1}.rotz = 0;
% sources{1}.level = 0;
% sources{1}.phase = 0;


%example: plane rectangular source 100x60mm with symmetric stereo enabled
% sources{1}.f=5000; 
% sources{1}.dx=10e-3;
% sources{1}.Lz=100e-3;
% sources{1}.Ly=60e-3;
% sources{1}.dir=[1;0;0];
% sources{1}.dipole = true;
% sources{1}.xpos = 1;
% sources{1}.ypos = 1.5;
% sources{1}.roty = 0;
% sources{1}.zpos = 0;
% sources{1}.rotz = 25;
% sources{1}.level = 0;
% sources{1}.stereo=true;

%example: dome, 26mm diameter, 30mm radius of curvature
% sources{1}.f=15000;
% sources{1}.dx=2e-3;
% sources{1}.radius=13e-3;
% sources{1}.radcurv=30e-3;
% sources{1}.xpos = 0;
% sources{1}.ypos = 0;
% sources{1}.zpos = 0;
% sources{1}.roty = 0;
% sources{1}.rotz = 0;
% sources{1}.level = 0;


%example: cone, 200mm diameter, 60mm deep, minimum info needed
sources{1}.f=5000;
sources{1}.dx=10e-3;
sources{1}.radius=100e-3;
sources{1}.conedepth=60e-3;
%sources{1}.wavespeed = 500;



%example: MTM, dome tweeter plus double woofers
% sources{1}.f=2000;
% sources{1}.dx=10e-3;
% sources{1}.radius=80e-3;
% sources{1}.conedepth=40e-3;
% sources{1}.xpos = 0;
% sources{1}.zpos = -150e-3;
% sources{1}.ypos = 0;
% sources{1}.roty = -5;
% sources{1}.rotz = 0;
% sources{1}.level = 0;
% sources{1}.phase = 0;
% 
% sources{2}.dx=10e-3;
% sources{2}.radius=80e-3;
% sources{2}.conedepth=40e-3;
% sources{2}.xpos = 0;
% sources{2}.zpos = 150e-3;
% sources{2}.ypos = 0;
% sources{2}.roty = 5;
% sources{2}.rotz = 0;
% sources{2}.level = 0;
% sources{2}.phase = 0;
% 
% sources{3}.dx=3e-3;
% sources{3}.radius=12.5e-3;
% sources{3}.radcurv=20e-3;
% sources{3}.xpos = 0;
% sources{3}.zpos = 0;
% sources{3}.ypos = 0;
% sources{3}.roty = 0;
% sources{3}.rotz = 0;
% sources{3}.level = 0;
% sources{3}.phase = 0;

%% calculate source points
sourcestruct=preparesource(sources,10);
fprintf('The source consists of %d points. The area is %4.2f cm^2.\n',sourcestruct.n, sourcestruct.n*sourcestruct.dx^2*1e4)


%% create matrises och vectors
sp_per_wl=sourcestruct.lambda/sourcestruct.dx;
fprintf('Wavelenght is %4.2f mm. The source has %4.2f points per wavelength.\n',1000*sourcestruct.lambda, sp_per_wl)
if sp_per_wl<5 %5 seems like a reasonable limit
    fprintf('Warning: the source has too few points!\n')
end



room.nz=round(room.lz/room.dx); %vertical number of pixels (front wall height)
room.nx=round(room.lx/room.dx); %horizontal number of pixels (room length)
room.ny=round(room.ly/room.dx); %number of pixels for depth (room width)
room.z=linspace(-room.nz/2*room.dx,room.nz/2*room.dx,room.nz);
room.x=linspace(0,room.nx*room.dx,room.nx);
room.y=linspace(-room.ny/2*room.dx,room.ny/2*room.dx,room.ny);
fprintf('Room is %d pixels high, %d pixels wide and %d pixels long.\n',room.nz,room.ny,room.nx)

plotroom(sourcestruct,room,sweetspot,11);

% XZ plane (vertical) 
plane_xz = struct;
plane_xz.y0 = 0;
[plane_xz.xg, plane_xz.zg]=meshgrid(room.x,room.z); 
plane_xz.yg=plane_xz.y0*ones(size(plane_xz.xg));

% XY plane (horizontal)
plane_xy = struct;
plane_xy.z0=0; 
[plane_xy.xg, plane_xy.yg]=meshgrid(room.x,room.y); 
plane_xy.zg=plane_xy.z0*ones(size(plane_xy.xg));

% YZ plane (back wall)
plane_yz = struct;
plane_yz.x0=room.nx*room.dx; %x-coordinate of the plane
[plane_yz.yg, plane_yz.zg]=meshgrid(room.y,room.z); 
plane_yz.xg=plane_yz.x0*ones(size(plane_yz.zg));

% sweetspot plane
plane_ss = struct;
plane_ss.x0=sweetspot.xp; %x-coordinate of the plane
sweetspot.z = linspace(-sweetspot.l/2,sweetspot.l/2,sweetspot.n)+sweetspot.zp;
sweetspot.y = linspace(-sweetspot.l/2,sweetspot.l/2,sweetspot.n)+sweetspot.yp;
[plane_ss.zg, plane_ss.yg]=meshgrid(sweetspot.z,sweetspot.y); 
plane_ss.xg=plane_ss.x0*ones(size(plane_ss.zg));


%% calculate vertical pattern
world_v = zeros(room.nz,room.nx);
h = waitbar(0,'wait..');
world_v=calculatepattern(plane_xz,sourcestruct,world_v,h);
close(h)

logimage_v=20*log10(abs((world_v)));
h=figure(20);
if isoctave
    set(h,'position',[scrsz(3)/3 2*scrsz(4)/3 scrsz(3)/3-30 scrsz(4)/3.2-90]);
else
    set(h,'OuterPosition',[scrsz(3)/3 2*scrsz(4)/3 scrsz(3)/3 scrsz(4)/3.2]);
end

imagesc(room.x,room.z,logimage_v,[0 140]);
set(gca,'YDir','normal')
colormap(jet(256))
colorbar
axis image
title(['Vertical radiation pattern at ',num2str(sourcestruct.f), ' Hz'])
xlabel('X, m')
ylabel('Z, m')



%% calculate horizontal pattern
world_h = zeros(room.ny,room.nx);
h = waitbar(0,'wait..');
world_h=calculatepattern(plane_xy,sourcestruct,world_h,h);
close(h)

if sourcestruct.stereo
    image_lr=makestereoimage(world_h,'ud',1);
    h=figure(31);
    if isoctave
        imshow(image_lr);
    else
        imshow(image_lr,'InitialMagnification','fit')
    end
    set(gca,'YDir','normal')
    title(['Horizontal channel balance, ',num2str(sourcestruct.f), ' Hz'])
    if isoctave
        set(h,'position',[2*scrsz(3)/3 2*scrsz(4)/3 scrsz(3)/3-30 scrsz(4)/3.2-90]);
    else
        set(h,'OuterPosition',[2*scrsz(3)/3 2*scrsz(4)/3 scrsz(3)/3 scrsz(4)/3.2]);
    end
    world_h=world_h + flipud(world_h);
end


logimage_h=20*log10(abs((world_h)));
h=figure(30);

if isoctave
    set(h,'position',[scrsz(3)/3 scrsz(4)/3+20 scrsz(3)/3-30 scrsz(4)/3.2-90]);
else
    set(h,'OuterPosition',[scrsz(3)/3 scrsz(4)/3+20 scrsz(3)/3 scrsz(4)/3.2]);
end
imagesc(room.x,room.y,logimage_h,[0 140]);
set(gca,'YDir','normal')
colormap(jet(256))
colorbar
axis image
title(['Horizontal radiation pattern at ',num2str(sourcestruct.f), ' Hz'])
xlabel('X, m')
ylabel('Y, m')



%% calculate pattern on back wall
world_bak = zeros(room.nz,room.ny);
h = waitbar(0,'wait..');
world_bak=calculatepattern(plane_yz,sourcestruct,world_bak,h);
close(h)
if sourcestruct.stereo
    image_blr = makestereoimage(world_bak,'lr',1);
    h=figure(41);
    if isoctave
        imshow(image_blr);
    else
        imshow(image_blr,'InitialMagnification','fit');
    end
    set(gca,'YDir','normal')
    title(['Back wall channel balance, ',num2str(sourcestruct.f), ' Hz'])
    
    if isoctave
        set(h,'position',[2*scrsz(3)/3 scrsz(4)/3+20 scrsz(3)/3-30 scrsz(4)/3.2-90]);
    else
        set(h,'OuterPosition',[2*scrsz(3)/3 scrsz(4)/3+20 scrsz(3)/3 scrsz(4)/3.2]);
    end
    world_bak=world_bak + fliplr(world_bak);
end
logimage_bak=20*log10(abs((world_bak)));
h=figure(40);

if isoctave
    set(h,'position',[scrsz(3)/3 40 scrsz(3)/3-30 scrsz(4)/3.2-90]);
else
    set(h,'OuterPosition',[scrsz(3)/3 40 scrsz(3)/3 scrsz(4)/3.2]);
end
imagesc(room.y,room.z,logimage_bak,[0 140]);
set(gca,'YDir','normal')
colormap(jet(256))
colorbar
axis image
title(['Radiation pattern on back wall, ',num2str(sourcestruct.f), ' Hz'])
xlabel('Y, m')
ylabel('Z, m')


%% calculate frequency response in sweetspot
h = waitbar(0,'wait..');

nbrf=100;


fscan=logspace(log10(20),log10(20000),nbrf);
ss_average = zeros(size(fscan));

for nf=1:nbrf
    
    sources{1}.f=fscan(nf);
    sourcestruct=preparesource(sources,0);
    waitbar(nf/nbrf,h)
    world_ss = zeros(sweetspot.n,sweetspot.n);
    world_ss=calculatepattern(plane_ss,sourcestruct,world_ss,0);
    ss_average(nf)=20*log10(mean(abs(world_ss(:))));
end
close(h)

h=figure(50);

if isoctave
    set(h,'position',[1 40 scrsz(3)/3-30 scrsz(4)/3.2-90]);
else
    set(h,'OuterPosition',[1 40 scrsz(3)/3 scrsz(4)/3.2]);
end
semilogx(fscan,ss_average)
grid on
title(['Averaged level in sweetspot. Diff is: ',num2str(ss_average(nbrf)-ss_average(1),'%2.1f'),' dB']) 
xlabel('Frequency, Hz')
ylabel('SPL, dB')