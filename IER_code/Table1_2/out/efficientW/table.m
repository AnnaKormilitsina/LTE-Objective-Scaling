clear all
close all
clc

vector = [1 5 10  50 100 500 1000]';
a = []; b = [];
for i = 1:length(vector)
    eval(['load xVtheta_',num2str(vector(i)),'_200.txt'])
    eval(['load xVg_',num2str(vector(i)),'_200.txt'])
    eval(['vtheta(i,:) = xVtheta_',num2str(vector(i)),'_200;']);
    eval(['vg(i,:) = xVg_',num2str(vector(i)),'_200;']);
    a = [a; '&'];
    b = [b; '\\'];
end
format short g
Table1 = [ num2str(vector) a num2str(vg(:,1),'% 10.2g') a num2str(vtheta(:,1),'% 10.2g') a num2str(vg(:,2),'% 10.2g') a num2str(vtheta(:,2),'% 10.2g') a num2str(vg(:,4),'% 10.2g') a num2str(vtheta(:,4),'% 10.2g') b ]