clear
close all
clc

L=1;
N=100;
x=linspace(0,1,N);
Dx=x(2)-x(1);
t_end=1;
h=0.001;
Ub=1;

Uin=exp(-10*x)';

U=BDF(Uin,Ub,N,Dx,h,t_end);
plot(x,U(:,100));