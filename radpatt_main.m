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

clear all
close all

isoctave = (exist ('OCTAVE_VERSION', 'builtin') > 0); %running in matlab or octave?
if isoctave 
  scrsz = get(0,'screensize');
else
  scrsz = get(groot,'ScreenSize');
end

%% define the room and source
%run('examples\talldipole.m');
%run('examples\dome.m');
%run('examples\cone.m');
%run('examples\array.m');
%run('examples\stereo.m');
run('examples\twoway.m');

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

% directivity pattern
%L_arc = 5; %m, distance to arcs
L_arc = sweetspot.xp - mean(sourcestruct.x); %%m, distance to arcs, default is distance between first source and sweetspot
np_arc = 100; %number of points
nf_arc = 200; %number of frequencies
arc_angles = linspace(-pi/2,pi/2,np_arc);
arc_hor.xg = L_arc*cos(arc_angles)+mean(sourcestruct.x);
arc_hor.yg = L_arc*sin(arc_angles)+mean(sourcestruct.y);
arc_hor.zg = zeros(size(arc_hor.xg))+mean(sourcestruct.z);
arc_vert.zg = L_arc*sin(arc_angles)+mean(sourcestruct.z);
arc_vert.xg = L_arc*cos(arc_angles)+mean(sourcestruct.x);
arc_vert.yg = zeros(size(arc_hor.xg))+mean(sourcestruct.y);



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
if ~isoctave
    [ticks,labels]=loglabels(fscan(1),fscan(end));
    set(gca,'XTick',ticks);
    set(gca,'XTickLabel',labels);
end
grid on
title(['Averaged level in sweetspot. Diff is: ',num2str(max(ss_average)-min(ss_average),'%2.1f'),' dB']) 
xlabel('Frequency, Hz')
ylabel('SPL, dB')


%% calculate directivity patterns (disabled when stereo is enabled)
if sourcestruct.stereo~=true;
    h = waitbar(0,'wait..');
    
    fscan=logspace(log10(100),log10(20000),nf_arc);
    world_dir_hor = zeros(np_arc,nf_arc);
    world_dir_vert = zeros(np_arc,nf_arc);
    world_temp = zeros(1,np_arc);
    for nf=1:nf_arc
        
        sources{1}.f=fscan(nf);
        sourcestruct=preparesource(sources,0);
        waitbar(nf/nf_arc,h)
        world_temp = zeros(1,np_arc);
        world_temp=calculatepattern(arc_hor,sourcestruct,world_temp,0);
        world_dir_hor(:,nf)=world_temp;
        world_temp = zeros(1,np_arc);
        world_temp=calculatepattern(arc_vert,sourcestruct,world_temp,0);
        world_dir_vert(:,nf)=world_temp;
    end
    close(h)
    
    h=figure(60);
    
    logimage_dir_hor=20*log10(abs((world_dir_hor)));
    
    
    if isoctave
        set(h,'position',[2*scrsz(3)/3 scrsz(4)/3+20 scrsz(3)/3-30 scrsz(4)/3.2-90]);
    else
        set(h,'OuterPosition',[2*scrsz(3)/3 scrsz(4)/3+20 scrsz(3)/3 scrsz(4)/3.2]);
    end
    %imagesc( [fscan(10), fscan(end)],arc_angles*180/pi,logimage_dir_hor,[0 140]);
    %imagesc(logimage_dir_hor,[max(logimage_dir_hor(:))-30 max(logimage_dir_hor(:))]);
    if isoctave
        imagesc(([log10(fscan(1)),log10(fscan(end))]),arc_angles*180/pi,logimage_dir_hor,[max(logimage_dir_hor(:))-20 max(logimage_dir_hor(:))]);
        set(gca,'YDir','normal')
        ticks = get(gca,'XTick');
        set(gca, 'XTickLabel',10.^ticks);
    else    
        imagesc([0,1],arc_angles*180/pi,logimage_dir_hor,[max(logimage_dir_hor(:))-20 max(logimage_dir_hor(:))]);
        set(gca,'YDir','normal')
        [ticks,labels]=loglabels(fscan(1),fscan(end));
        ticks = (log10(ticks)-(log10(fscan(1))))/(log10(fscan(end))-log10(fscan(1)));
        set(gca, 'XTick',ticks)
        set(gca, 'XTickLabel',labels);
    end
    colormap(jet(256))
    colorbar
    %axis image
    title(['Horizontal directivity'])
    xlabel('f, Hz')
    ylabel('angle, deg')
    
    logimage_dir_vert=20*log10(abs((world_dir_vert)));
    h=figure(61);
    if isoctave
        set(h,'position',[2*scrsz(3)/3 40 scrsz(3)/3-30 scrsz(4)/3.2-90]);
    else
        set(h,'OuterPosition',[2*scrsz(3)/3 40 scrsz(3)/3 scrsz(4)/3.2]);
    end
    

    if isoctave
        imagesc(([log10(fscan(1)),log10(fscan(end))]),arc_angles*180/pi,logimage_dir_vert,[max(logimage_dir_vert(:))-20 max(logimage_dir_vert(:))]);
        set(gca,'YDir','normal')
        ticks = get(gca,'XTick');
        set(gca, 'XTickLabel',10.^ticks);
    else
        imagesc(([0,1]),arc_angles*180/pi,logimage_dir_vert,[max(logimage_dir_vert(:))-20 max(logimage_dir_vert(:))]);
        set(gca,'YDir','normal')
        [ticks,labels]=loglabels(fscan(1),fscan(end));
        ticks = (log10(ticks)-(log10(fscan(1))))/(log10(fscan(end))-log10(fscan(1)));
        set(gca, 'XTick',ticks)
        set(gca, 'XTickLabel',labels);
    end
    colormap(jet(256))
    colorbar
    %axis image
    title(['Vertical directivity'])
    xlabel('f, Hz')
    ylabel('angle, deg')
end