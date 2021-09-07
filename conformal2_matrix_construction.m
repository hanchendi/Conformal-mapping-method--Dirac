clear
clc

load('number.mat');
load('zero_point.mat');
load('parameter.mat');
load('normalized.mat')

l=parameter(1);

M1=zeros(l,l);
M2=zeros(l,l);
M3=zeros(l,l);
M4=zeros(l,l);
M5=zeros(l,l);

for i=1:l
    for j=1:i
        
        if number(i)==number(j)
            M1(i,j)=normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i), r*zero_point(i)).*besselj(number(j), r*zero_point(j))+besselj(number(i)+1, r*zero_point(i)).*besselj(number(j)+1, r*zero_point(j)))...
            .*r.^3,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
        end
        M1(i,j)=M1(i,j)*2*pi;
        
        if abs(number(i)-number(j))==1
            M2(i,j)=normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i), r*zero_point(i)).*besselj(number(j), r*zero_point(j))+besselj(number(i)+1, r*zero_point(i)).*besselj(number(j)+1, r*zero_point(j)))...
            .*r.^2,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
        end
        M2(i,j)=M2(i,j)*pi;
        
        if number(i)==number(j)
            M3(i,j)=normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i), r*zero_point(i)).*besselj(number(j), r*zero_point(j))+besselj(number(i)+1, r*zero_point(i)).*besselj(number(j)+1, r*zero_point(j)))...
            .*r.^5,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
        end
        M3(i,j)=M3(i,j)*2*pi;
        
        if abs(number(i)-number(j))==1
            M4(i,j)=normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i), r*zero_point(i)).*besselj(number(j), r*zero_point(j))+besselj(number(i)+1, r*zero_point(i)).*besselj(number(j)+1, r*zero_point(j)))...
            .*r.^4,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
        end
        M4(i,j)=M4(i,j)*pi;
        
        if abs(number(i)-number(j))==2
            M5(i,j)=normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i), r*zero_point(i)).*besselj(number(j), r*zero_point(j))+besselj(number(i)+1, r*zero_point(i)).*besselj(number(j)+1, r*zero_point(j)))...
            .*r.^3,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
        end
        M5(i,j)=M5(i,j)*pi;
    end
    disp(i/l)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%conformal transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=0.2;
c=0.2;
delta=pi/3;
B4=zeros(l,l);
B5=zeros(l,l);

for i=1:l
    for j=1:i
        if number(i)-number(j)==1
            B4(i,j)=1/sqrt(-1)*M4(i,j);
        elseif number(i)-number(j)==-1
            B4(i,j)=-1/sqrt(-1)*M4(i,j);
        elseif number(i)-number(j)==2
            B5(i,j)=1/sqrt(-1)*M5(i,j);
        elseif number(i)-number(j)==-2
            B5(i,j)=-1/sqrt(-1)*M5(i,j);
        end
    end
end

for i=1:l
    for j=i:l
        M1(i,j)=M1(j,i);
        M2(i,j)=M2(j,i);
        M3(i,j)=M3(j,i);
        M4(i,j)=M4(j,i);
        M5(i,j)=M5(j,i);
        B4(i,j)=conj(B4(j,i));
        B5(i,j)=conj(B5(j,i));
    end
end

J=1/(1+2*b^2+3*c^2).*(eye(l)+4*b^2.*M1+4*b.*M2+9*c^2.*M3+12*b*c*(cos(delta).*M4-sin(delta).*B4)+6*c*(cos(delta).*M5-sin(delta).*B5));

clear M1 M2 M3 M4 M5 B4 B5
for i=1:l
    U(i,i)=1/zero_point(i);
end
J1=U*J*U;
disp(1)
[D T]=eig(J1);
disp(2)
for i=1:l
    eigv_statical(i)=sqrt(1/real(T(i,i)));
    eigs_statical(:,i)=U*D(:,i);
end
 
eigs_statical=eigs_statical(1:l,1:5000);           
save([pwd,'/eigv_statical.mat'],'eigv_statical','-v7.3');
save([pwd,'/eigs_statical.mat'],'eigs_statical','-v7.3');
