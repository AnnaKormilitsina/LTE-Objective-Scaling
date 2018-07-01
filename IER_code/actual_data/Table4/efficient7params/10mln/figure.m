clear all
%close all
%clc

num_theta = 7;
names = char('\kappa ','\rho_{\theta} ',' \rho_{\epsilon} ',' \rho_{\gamma} ',' \sigma^2_{\theta} ','\sigma^2_{\epsilon}', '\sigma^2_{\gamma}' ); 
mult = 1;
y = [];
for i = 1:num_theta
    eval(['load x',num2str(i),'_',num2str(mult),'_mom.txt'])
    eval(['x =  x',num2str(i),'_',num2str(mult),'_mom;'])
    y = [y x];
end

if(mult == 1000)

end

h = figure;
for i = 1:num_theta
    subplot(2,4,i)
    hist(y(:,i,1),20)
    axis tight
    title(names(i,:))
end
print(h,'-dps','hist1.eps')
print(h,'-dpdf','hist1.pdf')
h = figure;
for i = 1:num_theta
    subplot(2,4,i)
    plot(y(:,i,1))
    axis tight
    title(names(i,:))
end
