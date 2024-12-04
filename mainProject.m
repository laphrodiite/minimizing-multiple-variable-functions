clc, clearvars, close all

% Initialising function f.
sympref('FloatingPointOutput',true);
syms x y
f = x^5 * exp(-x^2-y^2);

% Calculation needed functions
v= [x y];
fGrad = gradient(f, v);
fHess = hessian(f, v);

% Plotting to get a clear image.
figure(1);
ezsurf(f);
figure(2);
ezcontour(f);
%%

% Testing 

gDes(f, fGrad, -1, 1, 0.5, false,false);

newton(f, fGrad, fHess, 1, -1, 0.5, false, false);

leven_marq(f, fGrad, fHess, 1, -1, 0.5, false, false)
