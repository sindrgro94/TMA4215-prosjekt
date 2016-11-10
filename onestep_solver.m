%% Introduction
% This script answers task 8
%' @jac jacobi matrix

jac = jacobi('test')

[tnext, ynext, le, iflag] = onestep(f,jac,tn,yn,h,Tolit)