clc
clear

% Reads in concentration dat file created in FORTRAN
A = ReadArray_FortranBinary('MSU_Conc_Data.dat',3);
[time, M, N] = size(A);

% Makes gif
gif('MSU_Diffusion.gif')

for i = 1:100:time
    rshape = reshape(A(i,:,:),[M N]);
    rshape = rshape';
    pcolor(rshape)
    hold on;
    colormap jet;
    shading flat;
    axis equal;
    hold off;
    pause(0.001);
    gif
end

web('MSU_Diffusion.gif');

