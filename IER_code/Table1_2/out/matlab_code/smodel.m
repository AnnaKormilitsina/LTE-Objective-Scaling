clear all
close all
clc

format short g
mult   =   500;
sig_th =  0.05;
% 1    5    10   50    100  500  1000
% 1.5  1    0.5  0.2  0.15  0.05 0.01
T = 200;
mu = 0.5;
if (1 == 2)
    % Normal with mean mu and std = mu
    data = mu + mu*randn(T,1);

    MOMdata(1) = mean(data);
    MOMdata(2) = var(data);

    fid = fopen(['MOMdata_',num2str(T),'.txt'],'w');
    for i = 1:length(MOMdata(:))
        fprintf(fid,[ num2str(MOMdata(i)) , 'D+0 \r\n']);
    end
    fclose(fid);
   
    % Efficient weighting
    W = weightm(data,T,0.5,0.5^2);
    fid = fopen(['invSIGdata_',num2str(T),'.txt'],'w');
    for i = 1:length(W(:))
        fprintf(fid,[ num2str(W(i)) , 'D+0 \r\n']);
    end
    fclose(fid);
    eval(['save MOMdata_',num2str(T),'.mat MOMdata W data'])
end

load MOMData_200

count = 0;
acc   = 0;

%Gmm
Theta  = [0.5 0.25];
ntheta = size(Theta,2);

loglike = object(Theta,MOMdata,W,T,mult)
logpri  = priord(Theta)

hessian = 1;
if hessian == 1 
         hessian_x     = fn_hesscd(Theta,1.0e-5, MOMdata,W,T,mult);
         H0            = -inv(hessian_x);
         [V,D]         = eig(H0);
         D             = abs(D);
         H0            = V*D*V';
         D             = chol(H0);
else
      D = 1;
end

tic;
Niter = 200000;
for j = 1:Niter
    if rem(j,10000) == 0  
       [ j  loglike 100*acc/j ]
       Theta
    end
        newTheta  = Theta + (sig_th*randn(size(Theta)))*D;
        newlogpri = priord(newTheta);
        while newlogpri == -inf 
            newTheta  = Theta + (sig_th*randn(size(Theta)))*D;%mvnrnd(zeros(1,length(Theta)), variance_normal); 
            newlogpri = priord(newTheta);
        end
        newloglike = object(newTheta,MOMdata,W, T, mult);
        if newloglike == -Inf;
                ratio = 0;
        else
                ratio = exp(mult*(newloglike + newlogpri - (loglike + logpri)));
                count = count + 1;
                if rand <= ratio % we accept with prob ratio
                        loglike = newloglike;
                        logpri  = newlogpri;
                        Theta   = newTheta;
                        acc     = acc + 1; % keeping track of the acceptances
                end
        end
        Theta_s(j,:)   = Theta;
        loglike_s(j,1) = loglike;
        logpri_s(j,1)  = logpri;
end
Theta           = mean(Theta_s);
acceptance_rate = acc/Niter*100;

toc;
hours = toc/3600;

Vg = zeros(ntheta,ntheta);
for i = 1:Niter
   Vg = Vg +  (Theta_s(i,:)-Theta)'*(Theta_s(i,:)-Theta);
end
Vg = T*Vg/Niter
Vtheta = mult^2 * Vg*W*Vg

fid = fopen(['xVg_',num2str(mult),'_200.txt'],'w');
for i = 1:length(Theta)^2
       fprintf(fid,[ num2str(Vg(i)) , ' \r\n']);
end
fclose(fid);

fid = fopen(['xVtheta_',num2str(mult),'_200.txt'],'w');
for i = 1:length(Theta)^2
       fprintf(fid,[ num2str(Vtheta(i)) , ' \r\n']);
end
fclose(fid);

eval(['save results_mu',num2str(mult),'.mat'])