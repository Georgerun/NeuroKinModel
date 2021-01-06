function [ComAif] = PrepareComAIF(Int_Ref, K1r, k2r)

% This function calculates an AIF based on the Ref data, by solving a 1-TC
% model, and assuming that free tracer is the same in the plasma and in the reference tissue. 

% It uses average values for K1r and k2r in this cohort, from previous
% measurements using the population-based AIF. 

% It also uses the Int_Ref derived from pmod

% Convert K1r, k2r into sec

K1r=K1r/60;
k2r=k2r/60;

%First, differentiate Int_Ref
DerRef=diff(Int_Ref);

% Then, compute ComnAif 
ComAif=(1./K1r).*(DerRef+k2r.*Int_Ref(1:end-1));

% Replace negative noisy points in the beginning, coming up due to differentiation
ComAif(ComAif<0)=0;

% Replace very high or low discontinuities after the first-pass, coming up due to differentiation

NMatrix=(ComAif(1350:end));

idx = find(NMatrix > 2.5);
idx2 = find(NMatrix <= 0);

NMatrix(idx)=NMatrix(idx-1);
NMatrix(idx2)=NMatrix(idx2-1);

ComAif=[ComAif(1:1349) NMatrix];

figure,plot(ComAif,'b-')
%ComTaif=Int_Taif(1:length(ComAif));

end