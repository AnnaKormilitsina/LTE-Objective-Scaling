close all
clear all
clc

num_theta = 2;
vector = [1 5 10 50 100 500 1000];
a = []; b = []; s = [];
for i = 1:length(vector)
    eval([ 'load xcovprob_',num2str(vector(i)),'.txt'])
    eval([ 'load xcovprobg_',num2str(vector(i)),'.txt'])
    y(i,1:num_theta)             = eval(['xcovprobg_',num2str(vector(i)),';']);
    y(i,num_theta+1:2*num_theta) = eval(['xcovprob_', num2str(vector(i)),';']);
    a = [ char('&');  a];
    b = [ char('\\'); b];
    s = [ char('$');  s];
end
%y = [vector' y];
table2 = [s num2str(vector')  s];
for i = 1:num_theta*2
    table2 = [ table2  a num2str(y(:,i),'% 10.3g')];
end
save table2_eff.mat table2
table2 = [table2 b]
