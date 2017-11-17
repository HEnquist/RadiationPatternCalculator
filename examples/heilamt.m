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
% sources{}.xo : transfer function of a crossover filter,
%                example 1st order lowpass at 100 Hz: @(s) (2*pi*100)./(s + 2*pi*100)
% sources{}.stereo : true to mirror source sideways for stereo, needs ypos>0 to make sense
%          only used for horizontal and back wall patterns


%example: Heil AMT  + 8-inch woofer, ideal brickwall filter
sources{1}.f=15000;
sources{1}.dx=10e-3;
sources{1}.radius=85e-3;
sources{1}.conedepth=40e-3;
sources{1}.xpos = 0;
sources{1}.zpos = -150e-3;
sources{1}.ypos = 0;
sources{1}.roty = 0;
sources{1}.rotz = 0;
sources{1}.level = 0;
sources{1}.phase = 0;
w0 = 2*pi*2000;
%sources{1}.xo = @(s) w0^2./(s.^2 + w0/Q*s + w0^2);  
sources{1}.xo = @(s) imag(s)<w0; 

sources{2}.dx=3e-3;
sources{2}.Lz=130e-3;
sources{2}.Ly=30e-3;
sources{2}.dir=[1;0;0];
%sources{2}.dipole = true;
sources{2}.xpos = -20e-3;
sources{2}.zpos = 100e-3;
sources{2}.ypos = 0;
sources{2}.roty = 0;
sources{2}.rotz = 0;
sources{2}.level = -0.5;
sources{2}.phase = 0;
%sources{2}.xo = @(s) -s.^2./(s.^2 + w0/Q*s + w0^2); 
sources{2}.xo = @(s) imag(s)>w0; 