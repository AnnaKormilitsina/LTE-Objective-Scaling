clear all
close all
clc

num_theta = 7;
names = char('\kappa ','\rho_{\theta} ',' \rho_{\epsilon} ',' \rho_{\gamma} ',' \sigma^2_{\theta} ','\sigma^2_{\epsilon}', '\sigma^2_{\gamma}' ); 
mult = 1000;
y = [];
for i = 1:num_theta
    eval(['load x',num2str(i),'_',num2str(mult),'_mom.txt'])
    eval(['x =  x',num2str(i),'_',num2str(mult),'_mom;'])
    y = [y x];
end
start = ceil(0.1*length(y));

h = figure;
for i = 1:num_theta
    subplot(2,4,i)
    hist(y(start:end,i,1),20)
    axis tight
    title(names(i,:))
end
print(h,'-dps','hist1000.eps')
print(h,'-dpdf','hist1000.pdf')
h = figure;
for i = 1:num_theta
    subplot(2,4,i)
    plot(y(1:end,i,1))
    axis tight
    title(names(i,:))
end
