clear all
close all
clc

num_theta = 7;
vector = [1 10 100 1000];
a = []; b = []; s = [];
for i = 1:length(vector)
    eval([ 'load xVtheta_',num2str(vector(i)),'.txt'])
    vtheta(i,:) = diag(reshape(eval(['xVtheta_',num2str(vector(i))]),num_theta,num_theta));
    eval([ 'load xVg_',num2str(vector(i)),'.txt'])
    vg(i,:) = diag(reshape(eval(['xVg_',num2str(vector(i))]),num_theta,num_theta));
    a = [ char('&');  a]; %char('&','&','&','&');
    b = [ char('\\'); b]; %,'\\','\\','\\');
    s = [ char('$');  s]; %,'$','$','$');
end
format short g
table_vg = [s num2str(vector')  s];
table_vtheta = table_vg;
for i = 1:num_theta 
    table_vg = [ table_vg  a num2str(vg(:,i),'% 10.3g')];%  a num2str(vg(:,2),'% 10.2g') a num2str(vg(:,3),'% 10.2g') a num2str(vg(:,4),'% 10.2g') a num2str(0.0001*vg(:,5),'% 10.2g') a num2str(0.0001*vg(:,6),'% 10.2g')  b ]
    table_vtheta = [table_vtheta a num2str(vtheta(:,i),'% 10.3g') ];% a num2str(vtheta(:,2),'% 10.2g') a num2str(vtheta(:,3),'% 10.2g') a num2str(vtheta(:,4),'% 10.2g') a num2str(0.0001*vtheta(:,5),'% 10.2g') a num2str(0.0001*vtheta(:,6),'% 10.2g')  b ]
end
disp('v lte')
table_vg = [table_vg b]
disp('vtheta')
table_vtheta = [table_vtheta b]