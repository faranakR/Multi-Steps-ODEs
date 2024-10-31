function [Y,h,Conv]=AB_AM_1(Yn,fun,h,tn,ETol,Frac)
%AB-1 predictor step
Y0=Yn(:,end)+h*fun(Yn(:,end),tn(end-1));
%AM-1 corrector step
Y=Yn(:,end)+h*fun(Y0,tn(end));
%Error estimation
EST=abs((-0.5/(-0.5-0.5))*norm(Y-Y0));
%Stepsize prediction
r=(Frac*ETol/EST)^(1/2);
h=r*h;
%Check whether current step succeeded
if EST>ETol
    Conv=0; %Unsuccessful
else
    Conv=1; %Successful
end
end