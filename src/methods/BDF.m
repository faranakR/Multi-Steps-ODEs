function [U,t]=BDF(Uin,Ub,N,Dx,h,t_end)
t=[0,h];
U=Uin;
i=1;
while t(i+1)<=t_end
    Sol(:,1)=U(2:end,i);
    while 1
        A=diag(-1/Dx*ones(1,N-1))+diag(1/Dx*ones(1,N-2),-1);
        B=[Ub/Dx;zeros(N-2,1)];
        Sol(:,2)=U(2:end,i)+h*(A*Sol(:,1)+B);
        Err=sum(abs(Sol(:,2)-Sol(:,1)))/(N-1);
        if Err<1e-6
            U(:,i+1)=[Ub;Sol(:,2)];
            i=i+1;
            t(i+1)=t(i)+h;
            break;
        end
        Sol(:,1)=Sol(:,2);
    end
end

end