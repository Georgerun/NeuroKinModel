% This is a function to perform Simplified Reference Tissue modelling using PET data.
% Author Dr Giorgos Papanastasiou
% 23/02/2018

function F = ModelSRTM_PET(params,t,Ref,m)

Timeoff=zeros(1,m);
Ref_New=[Timeoff Ref];
Ref=Ref_New(1:length(Ref));

% Model equation
R1=params(1);k2=params(2);BP=params(3);

% F1=(R1.*Ref+(k2-(R1.*k2)./(1+BP)).*Ref);
% F2=exp(-(k2.*t)./(1+BP));
% F=conv(F1,F2)*t(2);

F1=((k2-(R1.*k2)./(1+BP)));
F2=k2./(1+BP);
R=F1.*exp(-F2.*t);

Component1=R1.*Ref;
Component2=conv(Ref,R)*t(2);
Component2=Component2(1:7100);

F=Component1+Component2;

end

