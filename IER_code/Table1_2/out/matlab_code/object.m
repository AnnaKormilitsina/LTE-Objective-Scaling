function Ln = object(Theta,MOM, W,T,mult)
% f(1,1) = mean(data) - Theta(1);
% f(2,1) = sum(data - mean(data)).^2 /T - Theta(2);

f(1,1) = MOM(1) - Theta(1);
f(2,1) = MOM(2) - Theta(2);

Ln = - 0.5*T* f' * W * f;
%Ln = - f' * W * f; 