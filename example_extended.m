%% EXTENDED EXAMPLE FOR USING YASSES
%   Step 1: read the siemens star image
%   Step 2: provide information about magnification and pixel size of the camera
%   Step 3: run yasses
%   Step 4: extract a PSF from your MTF
%   Step 5: create an artifical additional MTF/PSF that is used to blur the
%   siemens star image. 
%   Step 6: evaluate the blurred siemens star image and evaluate it to see
%   if it matches the predicted total MTF/PSF 


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

%% STEP 4: EXTRACT PSF FROM MTF
% get the PSF from the MTF by inverse FFT
% as the MTF will probably have no values for low lpmm values, either a
% gaussian has to be fit to the data or extrapolation to lpmm -> 0 is
% necessary
% methods may either be 'gaussian' which fits a gaussian function to the
% data. this will force the extracted PSF to be also a Gaussian (because
% FFT or iFFT of a Gaussian is a Gaussian)
% other methods are passed to fillmissin. try 'linear',
% 'makima','pchip',...
mtf1 = ya.extrapolateLowLPMM(lpmm,mtf,'gaussian');
[x,psf] = ya.createPSFfromMTF(lpmm,mtf1);

% how broad is the psf?
psf_fwhm = ya.getFWHM(x,psf);
% how broad is the mtf? use MTF1 for this because prependen NaNs causes
% failure
mtf_fwhm = ya.getFWHM(lpmm,mtf1);

%% STEP 5: ADDITIONAL MTF
% Now, to verify the code, let's assume the siemens star image above is now
% recorded with a lens that has a more limited performance.
% As all components of the imaging system have their own MTF and degrade
% optical resolution by convolution with the corresponding PSF, the total
% MTF may be retrieved by multiplication of the individual MTFs (because
% convolution is multiplication after FFT)
% construct a gaussian MTF for the lens
mtf_lens_fwhm = 50; % play around with this value to change the optical quality of the lens
mtf_lens = ya.createMTFfromFHWM(lpmm,mtf_lens_fwhm);
% create the PSF for this lens
[~,psf_lens] = ya.createPSFfromMTF(lpmm,mtf_lens);
psf_lens_fwhm = ya.getFWHM(x,psf_lens);
% also create a 2D version of this in image space to blurr the siemens star
% by 2D convolution
% this needs the desired image size and also the pxsize and magnification
% to know how large the kernel is supposed to be
% attention: the FWHM is in mm, the pxsize in µm.
imsize = size(star); % provide the image size to the function
psf_kernel = ya.create2DPSFfromFWHM(imsize,psf_lens_fwhm,pxsize,magnification);

% predict the resulting PSF and MTF for the newly blurred siemens star
mtf_total = mtf1.*mtf_lens;
[~,psf_total] = ya.createPSFfromMTF(lpmm,mtf_total);

% now see what the siemens star evaluation looks like after convolution
star_lens = conv2(star,psf_kernel,'same');
% ... and get the MTF from the new star
[lpmm,mtf_star_total] = ya.getMTFfromStar(star_lens,magnification,pxsize,showplots);
% ... and make a gaussian extrapolation
mtf_star_total = ya.extrapolateLowLPMM(lpmm,mtf_star_total,'gaussian');
[~,psf_star_total] = ya.createPSFfromMTF(lpmm,mtf_star_total);


%% SHOW THE PLOTS

if showplots
    
    figure(1);subplot(2,2,4);hold all;
    plot(lpmm,mtf1,'r-.','LineWidth',2)
    
    figure(2);subplot(2,2,4);hold all;
    plot(lpmm,mtf_star_total,'r-.','LineWidth',2)
    
    
    figure(3);clf;
    subplot(2,2,1);hold all;
    plot(x*1000,psf,'r-.','LineWidth',2);
    plot(x*1000,psf_lens,'b-','LineWidth',2);
    plot(x*1000,psf_total,'g-','LineWidth',2);
    plot(x*1000,psf_star_total,'k-.','LineWidth',2,'DisplayName','Evaluated total');
    xlabel('x in \mum');
    ylabel('psf')
    xlim([-3*1000*psf_fwhm 3*1000*psf_fwhm])
    box on;
    grid on;
    
    subplot(2,2,2);hold all;
    plot(lpmm,mtf,'k-','LineWidth',2,'DisplayName','Original MTF');
    plot(lpmm,mtf1,'r-.','LineWidth',2,'DisplayName','Gaussian Fit MTF/PSF');
    plot(lpmm,mtf_lens,'b-','LineWidth',2,'DisplayName','Lens MTF/PSF');
    plot(lpmm,mtf_total,'g-','LineWidth',2,'DisplayName','Predicted total MTF/PSF');
    plot(lpmm,mtf_star_total,'k-.','LineWidth',2,'DisplayName','Evaluated total MTF/PSF');
    xlabel('lpmm');
    ylabel('mtf')
    xlim([0 2*mtf_fwhm])
    box on;
    grid on;
    legend('show')
    
    subplot(2,2,3);hold all;
    imagesc(psf_kernel)
    title('2D lens PSF')
    box on;
    grid on;
end

