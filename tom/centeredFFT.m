function [X,freq]=centeredFFT(x,dt)
%function [X,freq]=centeredFFT(x,dt)
%
% Computes FFT, with zero frequency in the center, and returns 
%  dimensional frequency vector.
%
% Inputs:
%    x, vector to be transformed
%    dt, time increment
%
% Adapted from a function written by Quan Quach of blinkdagger.com 
% Tom Farrar, 2016, jfarrar@whoi.edu



N=length(x);
%Generate frequency index
if mod(N,2)==0
    m=-N/2:N/2-1; % N even
else
    m=-(N-1)/2:(N-1)/2; % N odd
end
freq=m/(N*dt);  %the dimensional frequency scale
X=fft(x);
X=fftshift(X); %swaps the halves of the FFT vector so that the zero frequency is in the center\\
%If you are going to compute an IFFT, first use X=ifftshift(X) to undo the shift}
