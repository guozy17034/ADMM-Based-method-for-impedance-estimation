function [x,y] = solve_admm(param)

x = param.x0;
y = param.y0;
lambda = param.lambda;
beta = param.beta;
rho = param.rho;
rhomax = param.rhomax;
Hx = param.Hx;
Fx = param.Fx;
Hy = param.Hy;
Fy = param.Fy;
%%
options = optimoptions('quadprog','Algorithm','interior-point-convex');
xlb = 0;
xub = 3;
ylb = 1;
yub = 4;
maxIter = param.maxIter;
i = 1;
% for plot
funval = zeros(maxIter-1,1);
iterNum = zeros(maxIter-1,1);
while 1
    if i == maxIter
        break;
    end
    % solve x
    Hxx = eval(Hx);
    Fxx = eval(Fx);
    x = quadprog(Hxx,Fxx,[],[],[],[],xlb,xub,[],options);% descend
    % solve y
    Hyy = eval(Hy);
    Fyy = eval(Fy);
    y = quadprog(Hyy,Fyy,[],[],[],[],ylb,yub,[],options);%descend
    % update lambda
    lambda = lambda + rho*(2*x + 3*y -5); % ascend
    % rho = min(rhomax,beta*rho);
    funval(i) = compute_fval(x,y);
    iterNum(i) = i;
    i = i + 1;
end
figure;
err=abs(400.1123-funval);
plot(iterNum,err,'-r');

end

function  fval  = compute_fval(x,y)
fval = 100*(2*x^2+2)^2 + (y-2)^2;
end