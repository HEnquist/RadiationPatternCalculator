%% define the room
room = struct;
sweetspot = struct;

room.dx = 50e-3; %pixelsize in meters
room.lx = 10; %room length, m
room.lz = 10; %room height, m
room.ly = 3.2; %room width, m

sweetspot.l = 200e-3; %sweetspot size
sweetspot.n = 10; %sweetspot number of points per side
sweetspot.xp = 7; %sweetspot distance from front wall
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
% sources{}.xo : transfer function of a crossover filter,
%                example 1st order lowpass at 100 Hz: @(s) (2*pi*100)./(s + 2*pi*100)
% sources{}.stereo : true to mirror source sideways for stereo, needs ypos>0 to make sense
%          only used for horizontal and back wall patterns

Ndr = 200;
Ldr = 10e-3;
Rcurv  = 3.6;
angsep = 180/pi * Ldr/Rcurv;
%example: line array with Ndr drivers,

sources{1}.zpos = -Ndr/2 * Ldr + 0.5*Ldr;
for mm = 1:Ndr
    if mm==1
        sources{mm}.f=15000;
        sources{mm}.dx=5e-3;
        sources{mm}.Lz=Ldr;
        sources{mm}.Ly=5e-3;
        %sources{mm}.conedepth=10e-3;
    else
        sources{mm} = sources{mm-1};
    end
    temppos = -Ndr/2 - 0.5 + mm ;
    alpha = temppos*angsep;
    sources{mm}.zpos = Rcurv*sin(pi/180*alpha);
    sources{mm}.xpos = Rcurv*(cos(pi/180*alpha)-1);
    sources{mm}.roty = -alpha;
    sources{mm}.level = -abs(12*(2*temppos/Ndr)^2); %shading
    %sources{mm}.zpos = sources{mm-1}.zpos + Ldr;
end