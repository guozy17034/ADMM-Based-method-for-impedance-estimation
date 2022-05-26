function [H,F] = getHession_F(fn)
% fn : function name
% H :hessian matrix
% F :一次项系数
syms x y lambda rho;
if strcmp(fn,'f1')
    f = 100*(2*x^2+2)^2 + lambda*(2*x + 3*y -5) + rho/2*(2*x + 3*y -5)^2;
    H = hessian(f,x);
    F = (2*lambda + (rho*(12*y - 20))/2 - 2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fcol = collect(f,{'x'});
    disp(fcol);
elseif strcmp(fn,'f2')
    f = (y-2)^2 + lambda*(2*x + 3*y -5) + rho/2*(2*x + 3*y -5)^2;
    H = hessian(f,y);
    F = (3*lambda + (rho*(12*x - 30))/2 - 4);
    fcol = collect(f,{'y'});
    disp(fcol);
end
% fcol = collect(f,{'x'});
% fcol = collect(f,{'y'});
% disp(fcol);
end

