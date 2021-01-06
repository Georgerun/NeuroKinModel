% This is a function to perform 2CXM reversible model using PET data.
% Author Dr Giorgos Papanastasiou
% 23/02/2018

function F = Model2CXM_PET(params,t,AIF,m)

Timeoff=zeros(1,m);
AIF_New=[Timeoff AIF];
AIF=AIF_New(1:length(AIF));

% Model equation
K1=params(1);k2=params(2);k3=params(3);k4=params(4);
R=(K1.*((((k2+k3+k4+sqrt((k2+k3+k4).^2-4.*k2.*k4)))./2)-k3-k4)./(sqrt((k2+k3+k4).^2-4.*k2.*k4))).*exp((-(k2+k3+k4+sqrt((k2+k3+k4).^2-4.*k2.*k4))./2).*t)+(K1.*((((k2+k3+k4-sqrt((k2+k3+k4).^2-4.*k2.*k4)))./2)-k3-k4)./-(sqrt((k2+k3+k4).^2-4.*k2.*k4))).*exp((-(k2+k3+k4-sqrt((k2+k3+k4).^2-4.*k2.*k4))./2).*t);

F=conv(AIF,R)*t(2);
F=F(1:7000);

% Standard (initial) length
% F=F(1:7550);
end

% [x2c,fval2c]=fmincon(@residual2CXMPET,[0.001 0.0001 0.0001 0.0001],[],[],[],[],[0.0005 0.00005 0.000005 0.000005],[0.01 0.01 0.001 0.001],[],options,Int_time(1:90),Int_Caif(1:90),0,Int_Csupmot(1:90));
% figure,plot(Int_time(1:90),Int_Csupmot(1,1:90),'-ob');hold on;plot(Int_time(1:90),Model2CXM_PET(x2c,Int_time(1:90),Int_Caif(1:90),0),'r');