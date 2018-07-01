function grdd = fn_hesscd(Theta,grdh, MOM,W,T,mult)
% Computing numerical hessian using a central difference with
%                    function grdd = hesscd(fcn,x0,grdh,Passed variables1)
%
%   fcn: a string naming the objective function.
%   x0: a column vector n*1, at which point the hessian is evaluated.
%   grdh: step size.
%   grdd: hessian matrix (second derivative), n*n.

stps = eps^(1/3);
% eps: floating point relative accuracy or machine precision: 2.22e-16
% stps: step size recommended by Dennis and Schnabel: 6.006e-6

x0 = Theta(:);
f0 = object(Theta,MOM,W,T,mult);

% ** initializations
k    = length(x0);
grdd = zeros(k);

% ** Computation of stepsize (dh)
if all(grdh)
    dh = grdh;
else
    ax0 = abs(x0);
    if all(x0)
        dax0 = x0./ax0;
    else
        dax0 = 1;
    end
    dh = stps * (max([ax0 (1e-2)*ones(k,1)]'))' .* dax0;
end

xdh = x0 + dh;
dh  = xdh - x0;    % This increases precision slightly
dhm = dh(:,ones(k,1));
ee  = eye(k) .* dhm;

i = 1;
while i <= k
    j = i;
    while j <= k
        %fune1 = Lnx( x0 + ee(:,i) + ee(:,j), data, invSIG, f,fx,fy,fxp,fyp, sth, seps, sgam );
        fune1 = object(x0 + ee(:,i) + ee(:,j),MOM,W,T,mult);

        %        logpri=priordeephabits(x0 + ee(:,i) + ee(:,j));
%        fune1 = fune1 + logpri;
        
        
        fune2 = object(x0 - ee(:,i) + ee(:,j),MOM,W,T,mult);
%        logpri=priordeephabits(x0 - ee(:,i) + ee(:,j));
%        fune2 = fune2 + logpri;
        
        fune3 = object(x0 + ee(:,i) - ee(:,j),MOM,W,T,mult);
%        logpri=priordeephabits(x0 + ee(:,i) - ee(:,j));
%       fune3 = fune3 + logpri;
        
        fune4 = object(x0 - ee(:,i) - ee(:,j),MOM,W,T,mult);
        %        logpri=priordeephabits(x0 - ee(:,i) - ee(:,j));
%        fune4 = fune4 + logpri;
        
        grdd(i,j) = (fune1 - fune2 - fune3 + fune4)  / (4 * dh(i) * dh(j));
        
        if i ~= j
            grdd(j,i) = grdd(i,j);
        end
    j = j+1;
    end
    i = i+1;
end