function [] = plotroom( source,room,sweetspot,fignum )
% plotroom( source,room,sweetspot,fignum ) : show the source in the room
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

corners=[0 -room.ny*room.dx/2 -room.nz*room.dx/2;
        room.nx*room.dx -room.ny*room.dx/2 -room.nz*room.dx/2;
        room.nx*room.dx room.ny*room.dx/2 -room.nz*room.dx/2;
        0  room.ny*room.dx/2 -room.nz*room.dx/2;
        0 -room.ny*room.dx/2 -room.nz*room.dx/2;
        0 -room.ny*room.dx/2 room.nz*room.dx/2;
        room.nx*room.dx -room.ny*room.dx/2 room.nz*room.dx/2;
        room.nx*room.dx room.ny*room.dx/2 room.nz*room.dx/2;
        0  room.ny*room.dx/2 room.nz*room.dx/2;
        0 -room.ny*room.dx/2 room.nz*room.dx/2;
        nan nan nan;
        0 room.ny*room.dx/2 -room.nz*room.dx/2;
        0 room.ny*room.dx/2 room.nz*room.dx/2;
        nan nan nan;
        room.nx*room.dx room.ny*room.dx/2 -room.nz*room.dx/2;
        room.nx*room.dx room.ny*room.dx/2 room.nz*room.dx/2
        nan nan nan;
        room.nx*room.dx -room.ny*room.dx/2 -room.nz*room.dx/2;
        room.nx*room.dx -room.ny*room.dx/2 room.nz*room.dx/2 ];
    
corners_ss=[sweetspot.xp   -sweetspot.l/2+sweetspot.yp   -sweetspot.l/2+sweetspot.zp;
            sweetspot.xp   -sweetspot.l/2+sweetspot.yp    sweetspot.l/2+sweetspot.zp;
            sweetspot.xp    sweetspot.l/2+sweetspot.yp    sweetspot.l/2+sweetspot.zp;
            sweetspot.xp    sweetspot.l/2+sweetspot.yp   -sweetspot.l/2+sweetspot.zp;
            sweetspot.xp   -sweetspot.l/2+sweetspot.yp   -sweetspot.l/2+sweetspot.zp];

h=figure(fignum);

if isoctave
    set(h,'position',[1 scrsz(4)/3+20 scrsz(3)/3-30 scrsz(4)/3.2-90]);
else
    set(h,'OuterPosition',[1 scrsz(4)/3+20 scrsz(3)/3 scrsz(4)/3.2]);
end
plot3(source.x,source.y,source.z,'.',corners(:,1),corners(:,2),corners(:,3),corners_ss(:,1),corners_ss(:,2),corners_ss(:,3))
axis equal
xlabel('X, m')
ylabel('Y, m')
zlabel('Z, m')

end

