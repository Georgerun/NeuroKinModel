function [G]=residualSRTMPET(params,t,Ref,m,Int_Gd_Conc)

%residual at Fourier Fermi function
G=sum((ModelSRTM_PET(params,t,Ref,m)-Int_Gd_Conc).^2);

end
