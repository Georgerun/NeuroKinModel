function [TimeTotalInSec, Int_time, Put_all_data1, Put_all_data2] = CreateExponDecay(Tissue_1, Tissue_2, Time_1, Time_2, Phase1_time, Phase2_time)

% This function creates exponential decay phase in between the two disconnected phases (phase 1 and 2). 
% The duration of this exponential decay is about 60 min, whilst phases 1 and 2 is 30 min each. but accurate information 
% is here incorporated in the processing pipeline.

% Author Dr. Giorgos Papanastasiou
% 07/08/2018

% close all 

% Perform interpolation to appropriately extract data with 1 sec step.

% Prepare Int_time vector, matching the full temporal profile of the data.
% To calculate this, the 2 phase times are each separated by colons (hh:mm:ss), and then difference between the two phases is
% calculated (in sec). Phase1_time and 2_time are the beginning of the acquisition
% times for each phase. 

Tr1 = [fix(Phase1_time(1)/1E+4) fix(Phase1_time(1)/1E+2)-1E+2*fix(Phase1_time(1)/1E+4) rem(Phase1_time(1), 1E+2)];
Phase1_time = sprintf('%02d:%02d:%02d', Tr1);

Tr2 = [fix(Phase2_time(1)/1E+4) fix(Phase2_time(1)/1E+2)-1E+2*fix(Phase2_time(1)/1E+4) rem(Phase2_time(1), 1E+2)];
Phase2_time = sprintf('%02d:%02d:%02d', Tr2);

%TimeDif=Phase2_time(1)-Phase1_time(1);
% TimeDif
%TimeDifInSec=3600.*(TimeDif./10000);

TimeDifInSec=diff(datenum([Phase1_time;Phase2_time]))*24*3600;

% The above includes the Time_1 info. So, now calculate the time difference
% between phases only (subtracting Time_1 info)
TimeDifInSecBP=TimeDifInSec-Time_1(end);
TimeDifInSecBP=round(TimeDifInSecBP);
TimeTotalInSec=Time_1(end)+TimeDifInSecBP+Time_2(end);

Int_time=0:TimeTotalInSec;

% Round up or down accordingly the Time matrices before Interpolation
Time_1=round(Time_1);
Time_1=[0 Time_1];
Time_2=round(Time_2);

Tissue_1=Tissue_1';
Tissue_2=Tissue_2';

Tissue_1=[0 Tissue_1];

Int_Tissue_1=interp1(Time_1, Tissue_1, Int_time(1:Time_1(end)),'spline');
Int_Tissue_2=interp1(Time_2, Tissue_2, Int_time(Time_2(1):Time_2(end)),'spline');

% Perform linear fitting between the ending point of phase 1 and starting
% point of phase 2:
coefficients = polyfit([Int_time(Time_1(end)),Int_time(Time_2(1)+Time_1(end)+TimeDifInSecBP)],[Int_Tissue_1(Time_1(end)), Int_Tissue_2(1)], 1);

% As part of it, create time vector x which is the actual time difference
% between the two phases. Then, you calculate the linear part y. 
x=Int_time(Time_1(end))+1:Time_1(end)+TimeDifInSecBP+Time_2(1)-1;
y=coefficients(1)*x+coefficients(2);

% And create a new matrix incorporating this 60 min linear fit between
% phases 1 and 2. 
New_Matrix=[Int_Tissue_1 y Int_Tissue_2];

% Sanity graphs

% 
% figure,plot(Int_time(1:length(Int_Tissue_1)),Int_Tissue_1,'b-o')
Uptophase2=Time_1(end)+TimeDifInSecBP+Time_2(1)-1;
% hold on;plot(Int_time(Uptophase2:length(Int_time)-1),Int_Tissue_2,'b-o')
% 

%figure,plot(Int_time, New_Matrix)

% Create exponential decay(s) at the time period corresponding to the
% linear phase. 

Matrix_to_fit=New_Matrix(Int_time(Time_1(end)+1:Time_1(end)+TimeDifInSecBP+Time_2(1)-1));

% Define exponential decays (2CM-based) 
% g = fittype('A*exp(-b*x)');
g2 = fittype('(K1.*((((k2+k3+k4+sqrt((k2+k3+k4).^2-4.*k2.*k4)))./2)-k3-k4)./(sqrt((k2+k3+k4).^2-4.*k2.*k4))).*exp((-(k2+k3+k4+sqrt((k2+k3+k4).^2-4.*k2.*k4))./2).*x)+(K1.*((((k2+k3+k4-sqrt((k2+k3+k4).^2-4.*k2.*k4)))./2)-k3-k4)./-(sqrt((k2+k3+k4).^2-4.*k2.*k4))).*exp((-(k2+k3+k4-sqrt((k2+k3+k4).^2-4.*k2.*k4))./2).*x);');

%2 different exponential decays in linear part (dividing the linear part by two)
Gethalfmatrix=round(length(Matrix_to_fit)/2);

% Investigate whether matrix lengths are correct and account for any 1-2-3
% point inconsistencies due to differences between input (odd or even) data lengths. 
if Gethalfmatrix*2==length(Matrix_to_fit)
Matrix_to_fit=Matrix_to_fit;
end

if Gethalfmatrix*2>length(Matrix_to_fit)
Matrix_to_fit=[Matrix_to_fit Matrix_to_fit(end)];
end

if Gethalfmatrix*2<length(Matrix_to_fit)
Matrix_to_fit=Matrix_to_fit(1:end-1);
end

Time_to_fit=1:length(Matrix_to_fit);

% Fit a 2CM-based exponential decay in normalised linear part y2
% First, fit exp decay in the standard linear part.  

f01 = fit(Time_to_fit',Matrix_to_fit',g2,'StartPoint',[8 0.005 0.001 0.001]);
%figure,plot(Time_to_fit,Matrix_to_fit,'o',Time_to_fit,f01(Time_to_fit),'r-')

Exp_Matrix1=[f01(Time_to_fit(1,1:end))];
Exp_Matrix1=reshape(Exp_Matrix1,[],1);

% Normalise the fitted area, by subtracting the overestimated values at the
% ending points from the ending points themselves. 
S1=f01(Time_to_fit(1))-Matrix_to_fit(1); S2=f01(Time_to_fit(end))-Matrix_to_fit(end);
P1=Matrix_to_fit(1)-S1; P2=Matrix_to_fit(end)-S2;

% Re-calculate a linear fit in the normalised points
coefficients2 = polyfit([Int_time(Time_1(end)),Int_time(Time_2(1)+Time_1(end)+TimeDifInSecBP)],[P1, P2], 1);
y2=coefficients2(1)*x+coefficients2(2);

% Add this normalised linear part in the Matrix
NNew_Matrix=[Int_Tissue_1 y2 Int_Tissue_2];
New_Matrix_to_fit=NNew_Matrix(Int_time(Time_1(end)+1:Time_1(end)+TimeDifInSecBP+Time_2(1)-1));

% Investigate whether matrix lengths are correct and account for any 1-2-3
% point inconsistencies due to differences between input (odd or even) data lengths. 
if length(New_Matrix_to_fit)==length(Matrix_to_fit)
New_Matrix_to_fit=New_Matrix_to_fit;
end

if length(New_Matrix_to_fit)<length(Matrix_to_fit)
New_Matrix_to_fit=[New_Matrix_to_fit New_Matrix_to_fit(end)];
end

if length(New_Matrix_to_fit)>length(Matrix_to_fit)
New_Matrix_to_fit=New_Matrix_to_fit(1:end-1);
end

% Perform 2CM model fit in the normalised Matrix
fn = fit(Time_to_fit',New_Matrix_to_fit',g2,'StartPoint',[8 0.0005 0.0001 0.0001]);
%figure,plot(Time_to_fit,New_Matrix_to_fit,'o',Time_to_fit,fn(Time_to_fit),'r-')

NExp_Matrix=[fn(Time_to_fit(1,1:end))];
NExp_Matrix=reshape(NExp_Matrix,[],1);

Put_all_data1=[New_Matrix(1:Int_time(Time_1(end))) NExp_Matrix' New_Matrix(Int_time(Uptophase2+1:end))];
% figure,plot(Put_all_data1,'b-o')

% Fit 2CM-based exp decay across the NNew_Matrix (Matrix_to_fit_NNM, starting
% from x amount of data after the peak down to x amount of data after the beginning of phase 2)
% [MaxV, Pos]=max(New_Matrix,[],2);
[NMaxV, NPos]=max(NNew_Matrix,[],2);

NPostostart=NPos+500;

NTimepointsonphases=Time_1(end)-NPostostart;

Matrix_to_fit_NNM=NNew_Matrix(Int_time(NPostostart+1):Time_1(end)+TimeDifInSecBP+Time_2(1)-1+NTimepointsonphases);
Time_to_fit_NNM=1:length(Matrix_to_fit_NNM);

fnnm = fit(Time_to_fit_NNM',Matrix_to_fit_NNM',g2,'StartPoint',[NMaxV 0.0005 0.0001 0.0001]);
%
%
% figure,plot(Time_to_fit_NNM,Matrix_to_fit_NNM,'o',Time_to_fit_NNM, fnnm(Time_to_fit_NNM),'r-')
%
%

Exp_Matrix_NNM=[fnnm(Time_to_fit_NNM(NTimepointsonphases+1):Time_to_fit_NNM(end-NTimepointsonphases))];
Exp_Matrix_NNM=reshape(Exp_Matrix_NNM,[],1);

Put_all_data2=[NNew_Matrix(1:Int_time(Time_1(end))) Exp_Matrix_NNM' NNew_Matrix(Int_time(Uptophase2+3:end))];
% figure,plot(Put_all_data2,'b-o')

Put_all_data2=[0 Put_all_data2];

Int_time=Int_time(1:length(Put_all_data2));

end

  