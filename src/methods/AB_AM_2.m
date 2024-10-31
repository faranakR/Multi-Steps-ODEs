function [Y,h,Conv]=AB_AM_2(Yn,fun,h,t,ETol,Frac)
%AB-2 predictor step
Y0=Yn(:,end)+(h(2)+0.5*h(2)^2/h(1))*fun(Yn(:,end),t(end-1))-...
    (0.5*h(2)^2/h(1))*fun(Yn(:,end-1),t(end-2));
%AM-2 corrector step
Y=Yn(:,end)+0.5*h(2)*(fun(Y0,t(end))+fun(Yn(:,end),t(end-1)));
%Error estimation
EST=abs((-1/12/(-1/12-5/12))*norm(Y-Y0,2));
%Stepsize prediction
r=(Frac*ETol/EST)^(1/3);
h(1)=h(2);
h(2)=r*h(1);
%Check whether current step succeeded
if EST>ETol
    Conv=0; %Unsuccessful
else
    Conv=1; %Successful
end
end