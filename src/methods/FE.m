function Y=FE(Yn,fun,h,tn)
Y=Yn(:,end)+h*fun(Yn(:,end),tn(end-1));
end