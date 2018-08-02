% This script reads in an image and turns it into an array of 1's and 0s,
% where 1 is where the diffusion occurs and 0 is where the diffusion will
% not occur (This is called a domain parameter). Then writes it to a data
% file to be read into FORTRAN.

clear all;
clc;

% Reads in image, converts to black and white and makes it so 1 is
% non restricted area for diffusion and 0 is restricted.
ImageA = imread('MSU.jpg');
imshow(ImageA)
Inte = 0.2989 * ImageA(:,:,1) + 0.5870 * ImageA(:,:,2) + 0.1140 * ImageA(:,:,3);
A = (1-Inte/255);
pcolor(A);
WriteArray_FortranBinary('MSU_DP.dat',A)