function [ f , x , y ] = rsgene2D_anisotrop( N , M, rL , rW, h , Lx , Ly , LN)
%
% [f,x,y] = rsgene2D(N,rL,h,clx,cly) 
%
% generates a square 2-dimensional random rough surface f(x,y) with NxN 
% surface points. The surface has a Gaussian height distribution and 
% exponential autocovariance functions (in both x and y), where rL is the 
% length of the surface side, h is the RMS height and clx and cly are the 
% correlation lengths in x and y. 
%
% Input:    N   - number of surface points (along square side)
%           rL  - length of surface (along square side)
%           h   - rms height
%           L - correlation length (in x and y)
%
% Output:   f  - surface heights
%           x  - surface points
%           y  - surface points
%

% (C) Jack Landy & Alex Komarov, 2014 (adapted from code originally
% developed by David Bergström)

format long;

x = linspace(-rL/2,rL/2,N); y = linspace(-rW/2,rW/2,M);
[X,Y] = meshgrid(x,y);

if LN==1
    
    Z0 = sqrt(log(1 + h^2/1^2))*randn(M,N);
    fft2_Z = fft2(Z0);
    clear Z0

    C = exp(-sqrt((X/Lx).^2 + (Y/Ly).^2)); % Exponential correlation function
    fft2_F = sqrt(fft2(C));
    clear C

    f = exp(real(ifft2(fft2_Z.*fft2_F)) + log(1^2/sqrt(1^2+h^2))) - 1;
    f = f - mean(f(:));
    clear fft2_F fft2_Z
    
else
    
    Z0 = h*randn(M,N);
    fft2_Z = fft2(Z0);
    clear Z0

    C = exp(-sqrt((X/Lx).^2 + (Y/Ly).^2)); % Exponential correlation function
    fft2_F = sqrt(fft2(C));
    clear C

    f = real(ifft2(fft2_Z.*fft2_F));
    f = f - mean(f(:));
    clear fft2_F fft2_Z

end
