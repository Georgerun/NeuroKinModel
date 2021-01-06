% This is a function to perform Full Reference Tissue modelling using PET data.
% Author Dr Giorgos Papanastasiou
% 23/02/2018

function F = ModelFullRTM_PET(params,t,Ref,m)

Timeoff=zeros(1,m);
Ref_New=[Timeoff Ref];
Ref=Ref_New(1:length(Ref));

% Model equation
R1=params(1);k3=params(2);k4=params(3);k2=params(4);

s=k2+k3+k4;
r=k2/R1;
q=4.*k2.*k4;
p=sqrt(s^2-q);

c=(s+p)./2;
d=(s-p)./2;

a=((k3+k4-c).*(c-r))./p;
b=((d-k3-k4).*(d-r))./p;

Rr1=a.*exp(-c.*t);

Rr2=b.*exp(-d.*t);

Component1=R1.*Ref;

Component2=conv(Ref,Rr1)*t(2);
Component2=Component2(1:7000);

Component3=conv(Ref,Rr2)*t(2);
Component3=Component3(1:7000);

F=Component1+Component2+Component3;

end

