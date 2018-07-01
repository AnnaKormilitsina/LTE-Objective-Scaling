close all
clear all
clc

num_theta = 7;
vector = [1  10  100  1000];
a = []; b = []; s = [];
for i = 1:length(vector)
    eval([ 'load xcount_',num2str(vector(i)),'.txt'])
    y(i,1:num_theta) = eval(['xcount_', num2str(vector(i)),';']);
    a = [ char('&');  a];
    b = [ char('\\'); b];
    s = [ char('$');  s];
end
%y = [vector' y];
table2 = [s num2str(vector')  s];
for i = 1:num_theta
    table2 = [ table2  a num2str(y(:,i),'% 10.3g')];
end
table2 = [table2 b]
