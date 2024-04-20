clc, clearvars, close all

% Initialising function f.
sympref('FloatingPointOutput',true);
syms x y
f = x^3 * exp(-x^2-y^4);

% Calculation needed functions
v= [x y];
fGrad = gradient(f, v);
fHess = hessian(f, v);

% Plotting to get a clear image.
figure(1);
ezsurf(f);
figure(2);
ezcontour(f);

% Minimum found should be about -0.409916 @ [-1.22472, 0]

% Testing 

gDes(f, fGrad, -1, -1, 0.5, true);

newton(f, fGrad, fHess, -1, -1, 0.5, false);

leven_marq(f, fGrad, fHess, -1, -1, 0.5, false)
