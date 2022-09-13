function h = yasses()
%YASSES some functions for yasses.
% see description in subfunctions. to add more, just add functions and
% corresponding handles as below.
h.createMTFfromFHWM = @createMTFfromFHWM;
h.createPSFfromFHWM = @createPSFfromFHWM;
h.createPSFfromMTF = @createPSFfromMTF;
h.createMTFfromPSF = @createMTFfromPSF;
h.getFWHM = @getFWHM;
h.getMTFfromStar = @getMTFfromStar;
h.extrapolateLowLPMM = @extrapolateLowLPMM;
h.create2DPSFfromFWHM = @create2DPSFfromFWHM;
end

function [lpmm,mtf,ctf] = getMTFfromStar(star,magnification,pxsize,showplots,cx,cy)
% getMTF
% function to retrieve an MTF curve and object-side resolution
% from a Siemens Star Image. Magnification and pxsize are needed to
% retrieve object side resolution and nyquist frequency (aliasing).
% if parameter showplots is set to 1, results are displayed immediately
% the script assumes that the siemens star is of type square wave, not sine
% wave. thus, actually the CTF is measured, not the MTF. Conversion from
% CTF to MTF is done using Coltman's formula, see below.

% number of spokes (divide by two for linepairs) of the used star
nspokes=72;

% try to find center of star if not provided manually
if nargin == 4
    [~,cx] = max(sum(star,1));
    [~,cy] = max(sum(star,2));
end
% convert to polar coordinates to easily access the radius. this includes
% bilinear interpolation which should not tamper with the resolution
% result too much. oversampling may be specified to check for sensitivity.
[rows, cols] = size(star);
rmax = min([cx-1, cols-cx, cy-1, rows-cy]); % same definition as in polartrans
supsample_r = 5; % supersampling in radial direction
supsample_theta = 20; % number of pixels per spoke in theta direction
polstar = polartrans(star,supsample_r*rmax,supsample_theta*nspokes,cx,cy,'linear','valid');

% reshape the result (spokes are now vertical lines) to a 3D array with one
% linepair per 2D frame. it does not really matter, if it is split up at an
% edge, as the contrast later is defined on the maximum and minimum
% intensity in the given frame anyway.
% the resulting 36 (half of nspokes) frames are averaged (dimension 3)
linepairs = reshape(polstar,size(polstar,1),[],nspokes/2);
linepair = mean(linepairs,3,'omitnan');

% now derive the contrast (I_max-I_min)/(I_max+I_min) for every radius
I_max = max(linepair,[],2);
I_min = min(linepair,[],2);
ctf = (I_max-I_min)./(I_max+I_min);
% normalize to DC contrast
% ATTENTION: this assumes that DC contrast is actually visible in the
% image! This may not be the case (e.g. for badly focused systems)
ctf_raw = ctf./max(ctf);

% convert axis from r in super pixel to lp/mm
r_px = linspace(0,rmax-1,length(ctf));      % radial axis in original px dimensions
r_mm = r_px*(pxsize*1e-3)/magnification;    % radial axis in object-side(!) mm
lpmm_raw = nspokes/2./(2*pi*r_mm);          % 36 linepairs on a circle with circumference 2pi*r

% interpolate the result on a uniform grid for later (i)FFT
lpmm = 0:0.1:max(lpmm_raw(~isinf(lpmm_raw)),[],'omitnan');
idxvalid = (~isinf(lpmm_raw) & ~isnan(lpmm_raw)); % get valid indices without NaNs and Infs
ctf = interp1(lpmm_raw(idxvalid),ctf_raw(idxvalid),lpmm);

% get the MTF from coltmanns formula:
% Coltman, J. W., June 1954, The Specification of Imaging Properties by Response to a Sine WaveTarget, Journal of the Optical Society of America, Vol. 44, No. 6, pp. 468-471
% MTF(f) = pi/4*[CTF(f)+1/3CTF(3F)-1/5CTF(5f) + 1/7(CTF(7f) + ...]
% that is: + + - + + - - - + for order 0, 3, 5, 7, 11, 13, 15, 17, 19
sign = [1,-1,1,1,-1,-1,-1,1];
c = zeros(length(ctf),9);   % preallocate array for all 9 terms in Coltman equation
c(:,1)=ctf;                 % first order term
term = [3,5,7,11,13,15,17,19];  % higher order terms
for i = 1:length(term)
    % get all available correction terms.
    % due to finite signal length, higher terms can only be computed up to
    % a certain lpmm value. use 0 for those (handled by preallocation with
    % 0)
    termval = sign(i)*ctf(1:term(i):end)/term(i); 
    c(1:length(termval),i+1)=termval;        
end
% sum all terms
mtf = sum(c,2);
% normalize the MTF
mtf = mtf/max(mtf);


if showplots
    figure;
    subplot(2,2,1);hold on;
    imagesc(star);
    plot(cx,cy,'Marker','x','MarkerSize',10,'LineWidth',5,'Color','r')
    title('raw image with detected center')
    xlabel('x in px')
    ylabel('y in px')
    colormap('bone')
    axis image
    subplot(2,2,2);
    imagesc(polstar);
    title('polar transform of image')
    xlabel('\Theta')
    ylabel('Radius in px')
    subplot(2,2,3);
    plot(r_mm,ctf_raw,'LineWidth',3,'Color','k');
    xlabel('Object-side radius in mm')
    ylabel('CTF')
    grid on
    box on
    title('Radial contrast')
    set(gca,'YTick',0:0.25:1);
    subplot(2,2,4);hold all;
    plot(lpmm,ctf,'LineWidth',3,'Color','k');
    plot(lpmm,mtf,'LineWidth',3,'Color','r');
    xlabel('Object-side lp/mm')
    ylabel('CTF,MTF')
    legend('show','CTF','MTF')
    grid on
    box on
    % limit to reasonable range
    [~,idx]=min(abs(ctf-0.01));
    xlim([0 lpmm(idx)])
    ylim([0 1])
    set(gca,'YTick',0:0.25:1);
end

end


% POLARTRANS - Transforms image to polar coordinates
%
% Usage:    pim = polartrans(im, nrad, ntheta, cx, cy, linlog, shape)
%
% Arguments:
%           im     - image to be transformed.
%           nrad   - number of radius values.
%           ntheta - number of theta values.
%           cx, cy - optional specification of origin.  If this is not
%                    specified it defaults to the centre of the image.
%           linlog - optional string 'linear' or 'log' to obtain a
%                    transformation with linear or logarithmic radius
%                    values. linear is the default.
%           shape  - optional string 'full' or 'valid'
%                    'full' results in the full polar transform being
%                    returned (the circle that fully encloses the original
%                    image). This is the default.
%                    'valid' returns the polar transform of the largest
%                    circle that can fit within the image.
%
% Returns   pim    - image in polar coordinates with radius increasing
%                    down the rows and theta along the columns. The size
%                    of the image is nrad x ntheta.  Note that theta is
%                    +ve clockwise as x is considered +ve along the
%                    columns and y +ve down the rows.
%
% When specifying the origin it is assumed that the top left pixel has
% coordinates (1,1).

% Copyright (c) 2002 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% December 2002
% November 2006  Correction to calculation of maxlogr (thanks to Chang Lei)

function [pim,rmax] = polartrans(im, nrad, ntheta, cx, cy, linlog, shape)

[rows, cols] = size(im);

if nargin==3         % Set origin to centre.
    cx = cols/2+.5;  % Add 0.5 because indexing starts at 1
    cy = rows/2+.5;
end

if nargin < 7, shape = 'full'; end
if nargin < 6, linlog = 'linear'; end

if strcmp(shape,'full')         % Find maximum radius value
    dx = max([cx-1, cols-cx]);
    dy = max([cy-1, rows-cy]);
    rmax = sqrt(dx^2+dy^2);
elseif strcmp(shape,'valid')    % Find minimum radius value
    rmax = min([cx-1, cols-cx, cy-1, rows-cy]);
else
    error('Invalid shape specification');
end

% Increments in radius and theta

deltatheta = 2*pi/ntheta;

if strcmp(linlog,'linear')
    deltarad = rmax/(nrad-1);
    [theta, radius] = meshgrid([0:ntheta-1]*deltatheta, [0:nrad-1]*deltarad);
elseif strcmp(linlog,'log')
    maxlogr = log(rmax);
    deltalogr = maxlogr/(nrad-1);
    [theta, radius] = meshgrid([0:ntheta-1]*deltatheta, exp([0:nrad-1]*deltalogr));
else
    error('Invalid radial transformtion (must be linear or log)');
end

xi = radius.*cos(theta) + cx;  % Locations in image to interpolate data
yi = radius.*sin(theta) + cy;  % from.

[x,y] = meshgrid([1:cols],[1:rows]);
pim = interp2(x, y, double(im), xi, yi);

end

function mtf = createMTFfromFHWM(lpmm,fwhm)
mtf = gaussian(lpmm,0,fwhm);
mtf=mtf/max(mtf);
end

function psf = createPSFfromFHWM(x,fwhm)
psf = gaussian(x,0,fwhm);
psf=psf/max(psf);
end

function [x,psf] = createPSFfromMTF(lpmm,mtf)
% PSF is inverse FFT of PSF
% but attention, if lpmm has only values >= 0, make it symmetric for the
% IFFT
if all(lpmm>=0)
    lpmm = [fliplr(-lpmm(2:end)),lpmm];
    mtf = [fliplr(mtf(2:end)),mtf];
end

n = 2^19;
Y = fftshift(ifft(mtf,n));  % shift to center
Fs = 1/mean(diff(lpmm));    % get the sampling interval
P = abs(Y/n);            % magnitude
psf = P/max(P);
x = (-n/2:n/2-1)*(Fs/n);
end


function [lpmm,mtf] = createMTFfromPSF(x,psf)
% MTF is FFT of PSF
% actually, this is the same as createPSFfromMTF with the difference that
% FFT is used instead of iFFT
if all(x>=0)
    x = [fliplr(-x(2:end)),x];
    psf = [fliplr(psf(2:end)),psf];
end

n = 2^19;
Y = fftshift(fft(psf,n));
Fs = 1/mean(diff(x));
P = abs(Y/n);
mtf = P;
mtf = mtf/max(mtf);
lpmm = (-n/2:n/2-1)*(Fs/n);
end

function FWHM = getFWHM(xlpmm,mtfpsf)

% if xlpmm is symmetrical, take only the right half of the bell curve
% if not, use the entire set
if all(xlpmm>=0)
    idx = 1;
else
    % take right half of bell curve
    idx = (length(mtfpsf)-1)/2;
end
righthalf = mtfpsf(idx:end);
posxlpmm = xlpmm(idx:end);

[~,idx] = unique(righthalf);

FWHM = 2*interp1(righthalf(idx),posxlpmm(idx),0.5);
end

function [mtf_filled,mtf_scaled,fobj]=extrapolateLowLPMM(lpmm,mtf,method)

% if method == gaussian, fit a gaussian to the data
if strcmp(method,'gaussian')
    x=lpmm(~isnan(mtf));
    y=mtf(~isnan(mtf));
    % gaussian that is centered at 0
    fobj = fit(x(:),y(:),'gauss1','Lower',[0,0,0],'Upper',[Inf,0,Inf]);
    mtf_filled = fobj(lpmm)'; % transpose to be again of same dimensions
elseif strcmp(method,'gaussian2')
    x=lpmm(~isnan(mtf));
    y=mtf(~isnan(mtf));
    % gaussian that is centered at 0
    fobj = fit(x(:),y(:),'gauss2','Lower',[0,0,0,0,0,0],'Upper',[Inf,0,Inf,Inf,0,Inf]);
    mtf_filled = fobj(lpmm)'; % transpose to be again of same dimensions
elseif strcmp(method,'lorentzian')
    x=lpmm(~isnan(mtf));
    y=mtf(~isnan(mtf));
    % gaussian that is centered at 0
    eqn = 'a/(x^2+a^2)*b';
    fobj = fit(x(:),y(:),eqn,'Lower',[0,0],'Upper',[Inf,Inf]);
    mtf_filled = fobj(lpmm)'; % transpose to be again of same dimensions
else
    % use fillmissing
    mtf_filled = fillmissing(mtf,method)';
    fobj=[];
end
% renormalize! this is important because if low lpmm values are not
% available, DC contrast is probably overestimated.
mtf_filled = mtf_filled/max(mtf_filled);
% scale the original mtf to fit the curve
notnanidx = find(~isnan(mtf));
mtf_scaled = mtf/mtf(notnanidx(1))*mtf_filled(notnanidx(1));
end

function kernel = create2DPSFfromFWHM(imsize,fwhm,pxsize,magnification)
% create a 2D gaussian that is centered with subpixel accuracy on the grid
% do so by ensuring the kernel has odd size in both dimensions such that
% the center is on a pixel and not on a corner
nx = imsize(1);
ny = imsize(2);

if (-1)^nx == 1
    nx = nx+1;
end
if (-1)^ny == 1
    ny = ny+1;
end

% get the center pixel
cx = (nx+1)/2;
cy = (ny+1)/2;

% convert FWHM to sigma in px
% attention: the FWHM is in mm, the pxsize in µm.
% also: FWHM of PSF is object side, we need image side, so multiply by
% magnification
fwhm_px = fwhm*1000/pxsize*magnification;
sigma = fwhm_px/(2*sqrt(2*log(2))); % convert FWHM to RMS of gauss
x = (cx-nx:1:nx-cx)';
y = (cy-ny:1:ny-cy);

kernel = exp( - (x.^2)./(2 .* sigma.^2) - (y.^2)./(2 .* sigma.^2));

% now interpolate back on the original grid
kernel = imresize(kernel,[imsize(1),imsize(2)],'bicubic');

% normalize to the peak
kernel = kernel./max(kernel(:));

end

function g = gaussian(x,pos,wid)
g = exp(-((x-pos)./(0.6006.*wid)) .^2);
end
