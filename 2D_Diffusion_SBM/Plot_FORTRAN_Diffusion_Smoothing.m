clc
clear

% Reads in concentration dat file created in FORTRAN
A = ReadArray_FortranBinary('DP_Smoothed_Animate.dat',3);

[time, M, N] = size(A);

% Makes gif
gif('DP_Smoothing.gif')

for i = 1:time
    rshape = reshape(A(i,:,:),[M N]);
    pcolor(rshape)
    hold on;
    colormap jet;
    shading flat;
    axis equal;
    hold off;
    pause(0.001);
    gif
end

web('DP_Smoothing.gif');

A_end = reshape(A(19,:,:),[M N]);
A_end = A_end';
WriteArray_FortranBinary('DP_Smoothed_MATLAB.dat',A_end);