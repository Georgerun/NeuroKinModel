function [G]=residualFullRTMPET(params,t,Ref,m,Int_Gd_Conc)

%residual at Fourier Fermi function
G=sum((ModelFullRTM_PET(params,t,Ref,m)-Int_Gd_Conc).^2);

end
