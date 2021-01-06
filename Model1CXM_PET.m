% This is a function to perform 1CXM model (including Extraction fraction) using PET data.
% Author Dr Giorgos Papanastasiou
% 23/02/2018

function F = Model1CXM_PET(params,t,AIF,m)

Timeoff=zeros(1,m);
AIF_New=[Timeoff AIF];
AIF=AIF_New(1:length(AIF));

% Model equation (1:1750 might be the optimum time window for this model) 
F=params(1);E=params(2);k2=params(3);
R=F.*E.*exp(-k2.*t);
F=conv(AIF,R)*t(2);
F=F(1:7100);
end

% [x1c,fval1c]=fmincon(@residual1CXMPET,[0.001 0.0001 0.0001],[],[],[],[],[0.0005 0.00005 0.000005],[0.01 0.01 0.001],[],options,Int_time(1:90),Int_Caif(1:90),0,Int_Csupmot(1:90));
% figure,plot(Int_time(1:90),Int_Csupmot(1,1:90),'-ob');hold on;plot(Int_time(1:90),Model1CXM_PET(x1c,Int_time(1:90),Int_Caif(1:90),0),'r');