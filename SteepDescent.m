function[inform,x] = SteepDescent(fun, x, sdparams)
%sdparams = struct('maxit', 1000, 'toler', 1.0e-4);
%inform.status inform.iter
global numf numg
numf = 0;
numg = 0;
step_count = 0;
j = 1;
%X = zeros(sdparams.maxit);
%alfa = ones(sdparams.maxit);
alfa(1) = 1;

X(:,1) = x.p;
params = struct('ftol', 1e-4, 'gtol', 0.9, 'xtol', 1e-6,'stpmin',0,'stpmax',1e20,'maxfev',10000);
while step_count <= sdparams.maxit
    if abs(fun(X(:,j),2)) < sdparams.toler
        %gradient tolerance is achieved
        inform.status = 1;
        inform.iter = step_count;
        x.p = X(:,j);
        x.f = fun(X(:,j),1);
        x.g = fun(X(:,j),2);
       return;
    else
        if j > 1 
            alfa(j) = max(10 * params.xtol, fun(X(:,j-1),2)'*fun(X(:,j-1),2) * alfa(j-1)/(fun(X(:,j),2)'*fun(X(:,j),2)));
        end
        x_stepsize = struct('f',fun(X(:,j),1),'p',X(:,j),'g',fun(X(:,j),2));
        X(:,j+1) = X(:,j) - fun(X(:,j),2) * cvsrch(fun, x_stepsize, -fun(X(:,j),2), alfa(j),params);
        j = j + 1;
        step_count = step_count + 1;
    end
end
x.p = X(:,j);
x.f = fun(X(:,j),1);
x.g = fun(X(:,j),2);
inform.status = 0;
inform.iter = step_count;
return;
      