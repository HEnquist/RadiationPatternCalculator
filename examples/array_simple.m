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


%example: line array with 5 drivers,
sources{1}.f=5000;
sources{1}.dx=5e-3;
sources{1}.radius=40e-3;
sources{1}.conedepth=10e-3;
sources{1}.zpos = -200e-3;
for mm = 2:5
    sources{mm} = sources{mm-1};
    sources{mm}.zpos = sources{mm-1}.zpos + 100e-3;
end