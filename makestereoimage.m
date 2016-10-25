function [ image_stereo ] = makestereoimage( world,flipdir,style )
%makestereoimage(world,flipdir,style): Makes an image showing stereo balance
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

image_stereo=zeros([size(world),3]);
image_stereo(:,:,1)=20*log10(abs((world)))/60-0.5;
if strcmp(flipdir,'lr')
    image_stereo(:,:,3)=20*log10(fliplr(abs((world))))/60-0.5;
elseif strcmp(flipdir,'ud')
    image_stereo(:,:,3)=20*log10(flipud(abs((world))))/60-0.5;
end
if style == 1
    image_stereo(:,:,2)=min(image_stereo(:,:,[1 3]),[],3);
elseif style == 2
    image_stereo(:,:,2)=0.6*(image_stereo(:,:,1)+image_stereo(:,:,3));
end
image_stereo(image_stereo<0)=0;
image_stereo(image_stereo>1)=1;

end

