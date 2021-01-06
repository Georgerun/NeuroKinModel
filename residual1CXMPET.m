function G=residual1CXMPET(params,t,AIF,m,Int_Gd_Conc)

%residual at Fourier Fermi function
G=sum((Model1CXM_PET(params,t,AIF,m)-Int_Gd_Conc).^2);

end
 