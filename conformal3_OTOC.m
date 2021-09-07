clear
clc

load('eigs_statical.mat');
load('zero_point.mat');
load('normalized.mat');
load('number.mat');
load('parameter.mat');

M=8000;
N=2000;
b=0.2;
c=0.2;
delta=pi/3;
l=parameter(1);
dr=1/N;
rr=dr:dr:1;
rr=rr';
dtheta=2*pi/N;
theta=dtheta:dtheta:2*pi;


xx1=rr*cos(theta);
yy1=rr*sin(theta);
z=xx1+sqrt(-1)*yy1;
w=(z+b*z.^2+c*exp(sqrt(-1)*delta)*z.^3)./sqrt(1+2*b^2+3*c^2);
g=(1+2*b*z+3*c*exp(sqrt(-1)*delta)*z.^2)./sqrt(1+2*b^2+3*c^2);
xx2=real(w);
yy2=imag(w);

psi1=zeros(N,N);
psi2=zeros(N,N);
psi=zeros(N,N);

for i=1:N
    NN(i,:)=rr(i)*dr*dtheta*abs(g(i,:)).^2;
end

for mode=1601:2000
    psi1=zeros(N,N);
    psi2=zeros(N,N);
A(1:l)=abs(eigs_statical(1:l,mode));
for i=1:M
    a(i)=find(A==max(A));
    A(a(i))=0;
end

for i=1:M
    psi1=psi1+eigs_statical(a(i),mode)*normalized(a(i))*(besselj(number(a(i)),rr.*zero_point(a(i)))*exp(sqrt(-1)*(number(a(i)))*theta));
    psi2=psi2+eigs_statical(a(i),mode)*normalized(a(i))*(besselj(number(a(i))+1,rr.*zero_point(a(i)))*sqrt(-1)*exp(sqrt(-1)*(number(a(i))+1)*theta));
end
psi=abs(psi1).^2+abs(psi2).^2;
wave_normalized=1/sqrt(sum(sum(psi.*NN)));
psi1=psi1*wave_normalized;
psi2=psi2*wave_normalized;
%figure()
%mesh(xx2,yy2,abs(psi1).^2+abs(psi2).^2)

save([pwd,'/OTOC/Comformal_1_',num2str(mode),'.mat'],'psi1');
save([pwd,'/OTOC/Comformal_2_',num2str(mode),'.mat'],'psi2');
%close()
disp(mode)

clear A a
end
save([pwd,'/OTOC/NN.mat'],'NN');
save([pwd,'/OTOC/xx2.mat'],'xx2');
save([pwd,'/OTOC/yy2.mat'],'yy2');
