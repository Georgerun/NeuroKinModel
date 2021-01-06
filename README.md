# NeuroKinModel. Please cite the relevant manuscript: https://doi.org/10.1016/j.neuroimage.2020.117482: 

Giorgos Papanastasiou, Mark A. Rodrigues, Chengjia Wang, Kerstin Heurling, Christophe Lucatelli, Rustam Al-Shahi Salman, Joanna M. Wardlaw, Edwin J.R. van Beek, Gerard Thompson,
Pharmacokinetic modelling for the simultaneous assessment of perfusion and 18F-flutemetamol uptake in cerebral amyloid angiopathy using a reduced PET-MR acquisition time: Proof of concept,
NeuroImage (225), 2021, 117482.
ISSN 1053-8119,

(http://www.sciencedirect.com/science/article/pii/S1053811920309678)

The basic code to perform kinetic modelling and to generate model-based input functions is provided. Its basic components are:

The Model1CXM_PET, Model2CXM_PETM, ModelSRTM_PET and ModelFullRTM_PET codes are the functions responsible to perform 
1-compartmental, 2-compartmental, simplified reference tissue and full reference tissue kinetic modelling, respectively.

The "residual" components implement the nonlinear model fitting for each model. 

The PrepareComAIF calculated the model-based input function, by using standard (estimated) values for K1ʹ (=0.29 ml/min/ml) and 
k2ʹ (=0.11 1/min). 
These values can be modified accordingly for different studies, by adding your estimated or literature values for K1ʹ and k2ʹ.  

Finally, the CreateExponDecay code performs exponential decay between phases 1 and 2, as mathematically described in our manuscript. 
