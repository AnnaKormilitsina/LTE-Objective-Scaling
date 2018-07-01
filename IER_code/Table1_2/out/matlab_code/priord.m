function p = priord(Theta)

Theta = Theta(:);
p = 0;

% for a beta distribution, betapdf(x,a,b) gives draw from distribution with
% mu= a/(a+b), and sigmasq= ab/(a+b)^2(a+b+1), so we can solve 
% alfa= mu*(mu*(1-mu)/var -1)
% betta =(1-mu)*(mu*(1-mu)/var -1)
% but notice its variance, not Pstdev in the formula
%
% for a gamma distribution gammapdf(x, k, theta) give sdraw from
% distribution with mu = k*theta and sigmasq = k*theta^2, so we can solve,
% theta = sigmasq/mu and k = mu^2/sigmasq

%         mean     std   
%           1       2  %
DTheta = [ 'U',   'U']';
Lo     = [ -1000.0      0.0]';
Up     = [  1000.0,  1000.0]';
Pmean  = [ 0.5,   0.25]';
Pstd   = [0.2,    0.1]';  

if sum( Theta < Lo) > 0 | sum(Theta > Up) > 0
    p = -inf;
    return
else
    index_IG = find(DTheta == 'I');
    index_B  = find(DTheta == 'B');
    index_G  = find(DTheta == 'G');
    index_U  = find(DTheta == 'U');
    index_N  = find(DTheta == 'N');

    a(index_B,1) = Pmean(index_B).*( Pmean(index_B).*(1-Pmean(index_B))./(Pstd(index_B).^2)-1);
    b(index_B,1) = (1-Pmean(index_B)).*( Pmean(index_B).*(1-Pmean(index_B))./(Pstd(index_B).^2)-1);
    
    a(index_N,1) = Pmean(index_N);
    b(index_N,1) = Pstd(index_N);

    a(index_G,1) = (Pmean(index_G)./Pstd(index_G)).^2; % kapa
    b(index_G,1) = Pstd(index_G).^2./Pmean(index_G);   % theta

    a(index_IG,1) = (Pmean(index_IG)./Pstd(index_IG)).^2 + 2; % alpha
    b(index_IG,1) = Pmean(index_IG).*(a(index_IG)-1);         % beta

    a(index_U,1) = Lo(index_U);
    b(index_U,1) = Up(index_U);
        
    p  = 0;
    % UNIFORM DISTRIBUTIONS
    p  = sum( log( 1./(b(index_U) - a(index_U)) ) ) + p;
    % BETA DISTRIBUTIONS
    p  = sum( log( Theta(index_B).^ (a(index_B) - 1).*(1 - Theta(index_B)).^(b(index_B) - 1)./ beta(a(index_B),b(index_B)) ) ) + p;
    % NORMAL DISTRIBUTIONS
    % p = log(normpdf(Theta(i), 1, 1))+p; %investment adjustment cost
    xn = (Theta(index_B)-a(index_B))./b(index_B);
    p  = sum( -0.5*xn.^2 - 0.5*log(2*pi) - log(b(index_B) )) + p;
    % GAMA DISTRIBUTIONS
    p = sum((a(index_G) - 1).*log(Theta(index_G)) - (Theta(index_G) ./ b(index_G)) - gammaln(a(index_G)) - a(index_G).*log(b(index_G)) )+ p;
    % INV GAMA DISTRIBUTIONS
    p = sum( log(  b(index_IG).^a(index_IG)./gamma(a(index_IG)).*Theta(index_IG).^(-a(index_IG)-1).*exp(-b(index_IG)./Theta(index_IG)) )) + p;
end