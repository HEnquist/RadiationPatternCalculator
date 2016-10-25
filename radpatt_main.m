%Radiation Pattern Calculator v5
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
scrsz = get(groot,'ScreenSize');
%isoctave = (exist ('OCTAVE_VERSION', 'builtin') > 0); %running in matlab or octave?

%% define the room
room = struct;
sweetspot = struct;

room.dx = 50e-3; %pixelsize in meters
room.lx = 5.1; %room length, m
room.lz = 2.4; %room height, m
room.ly = 3.2; %room width, m

sweetspot.l = 200e-3; %sweetspot size
sweetspot.n = 10; %sweetspot number of points per side
sweetspot.xp = 3.6; %sweetspot distance from front wall
sweetspot.yp = 0; %sweetspot sideways position
sweetspot.zp = 0.0; %sweetspot vertical position



%% define the source
source = struct;

source.f=5000; %frequency, Hz
source.dx=10e-3; %source point spacing

% --Definition of the source parameters--
% source.Lz : source height, m
% source.Ly : source width, m
% source.radius : source radius, m, make large to disable
% source.radcurv : source radius of curvature, m, make large to disable
% source.conedepth : source cone depth, put to 0 to disable
% source.dir : [x;y;z] vector for lobe direction before rotation by rotz, only used when dipole=true
% source.dipole : false fot monopoles, true for dipoles
% source.xpos : source distance from front wall
% source.zpos : vertical position
% source.ypos : sideways position
% source.rotz : rotation around the vertical axis (toe-in), deg
% source.stereo : true to mirror source sideways for stereo, needs Syp>0 to make sense
%          only used for horizontal and back wall petterns

%example: dipole line source 300x20mm
source.Lz=1.32;
source.Ly=0.025;
source.radius=1000;
source.radcurv=1000;
source.conedepth=0;
source.dir=[1;0;0];
source.dipole = true;
source.xpos = 1;
source.ypos = 1;
source.zpos = 0;
source.rotz = 0;
source.stereo=true;

%example: plane rectangular source 100x60mm with stereo enabled
% source.Lz=100e-3;
% source.Ly=60e-3;
% source.radius=1000;
% source.radcurv=1000;
% source.conedepth=0;
% source.dir=[1;0;0];
% source.dipole = true;
% source.xpos = 1;
% source.ypos = 1.5;
% source.zpos = 0;
% source.rotz = 25;
% source.stereo=true;

%example: dome, 26mm diameter, 30mm radius of curvature
% source.Lz=26e-3;
% source.Ly=26e-3;
% source.radius=13e-3;
% source.radcurv=30e-3;
% source.conedepth=0;
% source.dipole = false;
% source.xpos = 0;
% source.ypos = 0;
% source.zpos = 0;
% source.rotz = 0;
% source.stereo=false;

%example: cone, 200mm diameter, 60mm deep
% source.Lz=200e-3;
% source.Ly=200e-3;
% source.radius=100e-3;
% source.radcurv=100;
% source.conedepth=60e-3;
% source.dir=[1;0;0];
% source.dipole = false;
% source.xpos = 0;
% source.zpos = 0;
% source.ypos = 0;
% source.rotz = 0;
% source.stereo=false;


%% calculate source points
source=preparesource(source,10);
fprintf('The source consists of %d points. The area is %4.2f cm^2.\n',source.n, source.n*source.dx^2*1e4)


%% create matrises och vectors
c=340;
source.lambda=c/source.f;
sp_per_wl=source.lambda/source.dx;
fprintf('Wavelenght is %4.2f mm. The source has %4.2f points per wavelength.\n',1000*source.lambda, sp_per_wl)
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

plotroom(source,room,sweetspot,11);

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
world_v=calculatepattern(plane_xz,source,world_v,h);
close(h)

logimage_v=20*log10(abs((world_v)));
h=figure(20);
set(h,'OuterPosition',[scrsz(3)/3 2*scrsz(4)/3 scrsz(3)/3 scrsz(4)/3.2]);
imagesc(room.x,room.z,logimage_v,[0 140]);
colormap(jet(256))
colorbar
axis image
title(['Vertical radiation pattern, ',num2str(source.Lz*100),' x ',num2str(source.Ly*100),' cm @',num2str(source.f), ' Hz'])
xlabel('X, m')
ylabel('Z, m')



%% calculate horizontal pattern
world_h = zeros(room.ny,room.nx);
h = waitbar(0,'wait..');
world_h=calculatepattern(plane_xy,source,world_h,h);
close(h)

if source.stereo
    image_lr=makestereoimage(world_h,'ud',1);
    h=figure(31);
    imshow(image_lr,'InitialMagnification','fit') %/max(image_lr(:))
    title(['Horizontal channel balance, ',num2str(source.f), ' Hz'])
    set(h,'OuterPosition',[2*scrsz(3)/3 2*scrsz(4)/3 scrsz(3)/3 scrsz(4)/3.2]);
    world_h=world_h + flipud(world_h);
end


logimage_h=20*log10(abs((world_h)));
h=figure(30);
set(h,'OuterPosition',[scrsz(3)/3 scrsz(4)/3+20 scrsz(3)/3 scrsz(4)/3.2]);
imagesc(room.x,room.y,logimage_h,[0 140]);
colormap(jet(256))
colorbar
axis image
title(['Horizontal radiation pattern, ',num2str(source.Lz*100),' x ',num2str(source.Ly*100),' cm @',num2str(source.f), ' Hz'])
xlabel('X, m')
ylabel('Y, m')



%% calculate pattern on back wall
world_bak = zeros(room.nz,room.ny);
h = waitbar(0,'wait..');
world_bak=calculatepattern(plane_yz,source,world_bak,h);
close(h)
if source.stereo
    image_blr = makestereoimage(world_bak,'lr',1);
    h=figure(41);
    imshow(image_blr,'InitialMagnification','fit') %/max(image_lr(:))
    title(['Back wall channel balance, ',num2str(source.f), ' Hz'])
    set(h,'OuterPosition',[2*scrsz(3)/3 scrsz(4)/3+20 scrsz(3)/3 scrsz(4)/3.2]);
    world_bak=world_bak + fliplr(world_bak);
end
logimage_bak=20*log10(abs((world_bak)));
h=figure(40);
set(h,'OuterPosition',[scrsz(3)/3 40 scrsz(3)/3 scrsz(4)/3.2]);
imagesc(room.y,room.z,logimage_bak,[0 140]);
colormap(jet(256))
colorbar
axis image
title(['Radiation pattern on back wall, ',num2str(source.f), ' Hz'])
xlabel('Y, m')
ylabel('Z, m')


%% calculate frequency response in sweetspot
h = waitbar(0,'wait..');

nbrf=100;
db20 = 1;
db20k = 1;

fscan=logspace(log10(20),log10(20000),nbrf);
ss_average = zeros(size(fscan));

for nf=1:nbrf
    source.lambda=c./fscan(nf);
    waitbar(nf/nbrf,h)
    world_ss = zeros(sweetspot.n,sweetspot.n);
    world_ss=calculatepattern(plane_ss,source,world_ss,0);
    ss_average(nf)=20*log10(mean(abs(world_ss(:))));
end
close(h)

h=figure(50);
set(h,'OuterPosition',[1 40 scrsz(3)/3 scrsz(4)/3.2]);
semilogx(fscan,ss_average)
grid on
title(['Averaged level in sweetspot. Diff is: ',num2str(ss_average(nbrf)-ss_average(1),'%2.1f'),' dB']) 
xlabel('Frequency, Hz')
ylabel('SPL, dB')