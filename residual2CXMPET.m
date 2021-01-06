function [G]=residual2CXMPET(params,t,AIF,m,Int_Gd_Conc)

%residual at Fourier Fermi function
G=sum((Model2CXM_PET(params,t,AIF,m)-Int_Gd_Conc).^2);

end
