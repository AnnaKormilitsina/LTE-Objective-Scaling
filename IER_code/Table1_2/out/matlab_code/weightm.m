function W = weightm(data, n, mu, sig2);

% m1 = mu;%mean(data);
% m2 = sig2;%var(data);
% m4 = 3*m2^2;
% W(1,1) = m2/n;
% W(2,2) = (n-1)^2/n^3*m4 - (n-1)*(n-3)/n^3*m2^2;
% W(1,2) = -m1^3 + (-6/n^3-3/n+6/n^2)*m1*m2;%3/n *m1*m2;
% W(2,1) = W(1,2);
%W(2,2) = ( (data-mu*ones(size(data))).^4'*(data-mu*ones(size(data))).^4 )/n-sig2^2;
%W(1,2) = ( (data-mu*ones(size(data))).^3'*(data-mu*ones(size(data))).^3 )/n;

W(1,1) = sig2;
W(1,2) = sum((data - mu*ones(size(data))).^3)/n;
W(2,1) = W(1,2);
W(2,2) = sum((data - mu*ones(size(data))).^4)/n -sig2^2;

W = inv(W);