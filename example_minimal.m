%% MINIMAL EXAMPLE FOR USING YASSES
%   Step 1: read the siemens star image
%   Step 2: provide information about magnification and pixel size of the camera
%   Step 3: run yasses
%   Step 4: enjoy.
%
%   For a more complex example on how to extract a PSF from the MTF, check
%   out example_extended.m

clear
close all
clc

%% STEP 1: get an example star image
star = imread('Star.tif');

%% STEP 2: provide information about magnification and pixel size of the camera
% Pixel size and magnification to convert px grid to mm
% If image-side resolution is desired, use magnification = 1
% For other magnifications, the object-side MTF is calculated which is
% probably what you want.
magnification=1;  
pxsize=6.5;         % pixel size in µm


%% STEP 3: run yasses
% get a handle on yasses' subfunctions
ya = yasses();
% show plots? 1=on, 0=off
showplots=1;
% get the MTF curve form the siemens star and show the results
[lpmm,mtf] = ya.getMTFfromStar(star,magnification,pxsize,showplots);