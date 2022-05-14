clc; close all; clear all; 

%% Loading image and display 
[I_rgb, map]= imread ('blister.jpg'); 
figure (1), imshow (I_rgb, []), title ('Blister')

% Transformation from RGB to intensity image
I= im2double(rgb2gray (I_rgb));

figure (2),subplot (1, 2, 1), imshow (I, []), title ('Blister Intensity Image'), ...
    subplot (1, 2, 2), imhist (I), ylim ([0 2.6e05]),ylabel ('# pixels'),title ('Blister Histogram');  

%% Different illumination effects removal with homomorphic filtering
I= im2double (I); 

% Log domain conversion
I_log = log(1 + I); 

% Gaussian high-pass filter construction

% Filter size: because we have to pad the image
M = 2*size(I,1) + 1; 
N = 2*size(I,2) + 1;

sigma = 8;

[X, Y] = meshgrid(1:N,1:M);
centerX = ceil(N/2);
centerY = ceil(M/2);
gaussianNumerator = (X - centerX).^2 + (Y - centerY).^2;
H = exp(-gaussianNumerator./(2*sigma.^2));
H = 1 - H;

figure (8), imshow(H,'InitialMagnification',25), title ('Gaussian filter')

H = fftshift(H); 

%  Log-transformed image in the frequency domain with zero-padding
If = fft2(I_log, M, N);

% High-pass filtering and compututation of the inverse-FFT
Iout = real(ifft2(H.*If));

% Crop the image back to the original unpadded size
Iout = Iout(1:size(I_rgb,1),1:size(I_rgb,2)); 

% Inverse log function to get get the homomorphic filtered image
Ihmf = exp(Iout) - 1;

figure (9), imshow (Ihmf, []), title ('Homomorphic filtered image')
% figure (10), imhist (Ihmf)

%% Blister external contour detection 

% Mask construction
im_height= size (I_rgb, 1);
im_width=  size (I_rgb, 2);
mask= zeros (im_height, im_width); 
BW1= imcomplement(imextendedmin (Ihmf, 0.43)); 
BW2= imcomplement(imextendedmax (Ihmf, 0.56));
mask= imadd (BW1, BW2); 
figure (11), imshow (mask), title ('Blister mask')

imcomp= imcomplement (mask); 
se= strel ('disk', 9); 
imcomp= imclose (imcomp,se); 

% Boundary extraction
[B,L,N] = bwboundaries(imcomp, 'noholes');
figure (12), imshow(I_rgb, []), hold on, visboundaries(B,'Color','r'), title ('Blister external contour'), hold off

% Aspect ratio external contour computation 
STATS = regionprops(L,'all'); 
asp_ratio= STATS.Eccentricity; 
 
sprintf('External contour aspect ratio: %0.3f', asp_ratio)

%% Search in each side for the expiration date and locate the side that contains it

% 4 ROI sides extraction 
b_sideP=[996.5000 2.9565e+03 1184 200]; ARb= 0.1689; 
l_sideP=[588.5000 916.5000 208 1640]; ARl=7.8846; 
t_sideP=[844.5000 492.5000 1520 176]; ARt= 0.1158; 
r_sideP=[2.3645e+03 1.0605e+03 264 1488]; ARr= 5.6364; 

side(1).rect= b_sideP; side(1).aspectRatio= ARb; 
side(2).rect= l_sideP; side(2).aspectRatio= ARl; 
side(3).rect= t_sideP; side(3).aspectRatio= ARt; 
side(4).rect= r_sideP; side(4).aspectRatio= ARr; 

I_frame=zeros (im_height, im_width); 
M_frame=zeros (im_height, im_width);
for i=1:length (side)
    Itemp_i= Ihmf; 
    figure (13), imshow (Itemp_i, []);
    dr= drawrectangle ('AspectRatio',side(i).aspectRatio, 'Position',side(i).rect); 
    BW=createMask(dr); 
    I_frame= imadd(I_frame,BW.*Itemp_i); 
    
    Itemp_m= mask; 
    figure (14), imshow (Itemp_m, []);
    drm= drawrectangle ('AspectRatio',side(i).aspectRatio, 'Position',side(i).rect); 
    BW=createMask(drm); 
    M_frame= imadd(M_frame,BW.*Itemp_m);
    
end 
figure (15), subplot (1,2,1),imshow (I_frame, []),title ('Sides of the HF blister'), subplot (1,2,2),imshow (M_frame, []), title ('Sides of the blister mask'); 

% Properties extraction for date localization
label_date = bwlabel(M_frame);
meas_date= regionprops(label_date, I_frame, 'MeanIntensity', 'Centroid');
int_date = [meas_date.MeanIntensity];
centr_date = [meas_date.Centroid];
x_centr = centr_date(1:2:end); y_centr = centr_date(2:2:end);
date = find(int_date > 0.17);
date_mask= ismember (label_date, date); 
for k = 1 : length(date)
    x(k) = x_centr(date(k));
    y(k) = y_centr(date(k));
end

figure (16), imshow (date_mask)
figure (17),imshow (I_rgb, []), hold on, plot (x,y, '.r'), legend ('Date digits'), title ('Date localization') 
hold off

%% Localization of missing and present pills 

mask=bwmorph (mask,'majority'); 
figure (18), imshow (mask)
label_pill = bwlabel(mask);
meas_pill= regionprops(label_pill, I_frame, 'Area', 'Centroid');

% Properties extraction for missing pills localization
area_miss=[meas_pill.Area]; 
centr_miss= [meas_pill.Centroid];
xM_centr = centr_miss(1:2:end); yM_centr = centr_miss(2:2:end);
miss_pill= find (area_miss<200000 & area_miss>10000);
MP_mask= ismember (label_pill, miss_pill); 
for k = 1 : length(miss_pill)
    x_M(k) = xM_centr(miss_pill(k));
    y_M(k) = yM_centr(miss_pill(k));
end

% Properties extraction for present pills localization
area_full=[meas_pill.Area]; 
centr_full= [meas_pill.Centroid];
xF_centr = centr_full(1:2:end); yF_centr = centr_full(2:2:end);
full_blis= find (area_full>214000 & area_full<285000); 
FP_mask= ismember (label_pill, full_blis);
for k = 1 : length(full_blis)
    x_F(k) = xF_centr(full_blis(k));
    y_F(k) = yF_centr(full_blis(k));
end

figure (19),subplot (1,2,1), imshow (MP_mask), title ('Empty blister'), ...
    subplot (1,2,2), imshow (FP_mask), title ('Full blister'); 

figure (20),imshow (I_rgb, []), hold on, ...
    plot(x_M, y_M, 'xr'), hold on, plot(x_F, y_F, 'og'),...
        legend ('Pill missing','Pill present'), title ('Localization of the precence of pills')
    



