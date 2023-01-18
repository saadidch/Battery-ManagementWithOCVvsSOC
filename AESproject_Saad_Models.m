%Battery C1202

SOC_b1= [0 0.0158 0.0315 0.0473 0.063 0.0788 0.0945 0.1062 0.1179 0.1296 0.1413 0.153 0.1885 0.2239 0.2594 0.2949 0.3303 0.384 0.4376 0.4912 0.5448 0.5985 0.6547 0.711 0.7673 0.8235 0.8798 0.9039 0.9279 0.9519 0.976 1 ];
OCV_b1= [2.6929 3.0649 3.2384 3.3177 3.3552 3.3762 3.3923 3.4039 3.4162 3.429 3.4425 3.4561 3.4967 3.5324 3.5619 3.5868 3.6094 3.644 3.6836 3.7298 3.7817 3.8368 3.8946 3.949 3.9974 4.0392 4.0759 4.0912 4.1073 4.125 4.1457 4.171 ];

%Battery C1203

SOC_b2=[0 0.0159 0.0317	0.0476 0.0634 0.0793 0.0951	0.1069 0.1186 0.1304 0.1421	0.1539 0.1895 0.225	0.2606 0.2962 0.3318 0.3855	0.4392 0.4929 0.5466 0.6003	0.6569 0.7135 0.7701 0.8267	0.8833 0.9067 0.93 0.9533 0.9767 1 ];
OCV_b2= [2.6902 3.064 3.2381 3.3175 3.3551 3.3761 3.3923 3.404 3.4162 3.4291 3.4425 3.4562 3.4968 3.5325 3.5622 3.5874 3.6103 3.6453 3.6852 3.7315 3.7834 3.8384 3.8963 3.9508 3.9995 4.0415 4.0784 4.0933 4.1089 4.1259 4.1456 4.1693 ];

%Battery C1204


SOC_b3= [0 0.0157 0.0315 0.0472 0.0629 0.0786 0.0944 0.1063 0.1182 0.1301 0.142 0.1539 0.1896 0.2254 0.2611 0.2968 0.3325 0.3864 0.4403 0.4943 0.5482 0.6021 0.6582 0.7144 0.7706 0.8267 0.8829 0.9063 0.9297 0.9532 0.9766 1 ];
OCV_b3= [2.7003 3.0671 3.2385 3.3171 3.3543 3.3751 3.391 3.4028 3.4151 3.4282 3.4418 3.4557 3.4965 3.5323 3.5621 3.5873 3.6101 3.6451 3.6849 3.7313 3.7833 3.8384 3.8957 3.9498 3.9982 4.0401 4.0773 4.0925 4.1083 4.1256 4.1456 4.1696 ];

%Battery C1205

SOC_b4= [0 0.0164 0.0328 0.0492 0.0656 0.082 0.0984 0.1092 0.12 0.1308 0.1416 0.1523 0.188 0.2237 0.2594 0.2951 0.3308 0.3867 0.4426 0.4985 0.5544 0.6103 0.6709 0.7315 0.792 0.8526 0.9132 0.9306 0.9479 0.9653 0.9826 1 ];
OCV_b4= [2.7296 3.0857 3.2497 3.3247 3.3609 3.3821 3.3991 3.41 3.4212 3.4329 3.4448 3.457 3.4963 3.5314 3.5611 3.5865 3.6099 3.6471 3.6893 3.7382 3.7931 3.8511 3.9139 3.9728 4.025 4.0694 4.108 4.1187 4.1297 4.1412 4.1537 4.1676 ];
%---------------------------------------------------------------------------------------------

E = 0.175; %scaling
One= [ones(32,1)];
SOC = 0:.0001:1;
z=SOC*(1-2*E) + E;


%================================================================ For Battery C1202 =============================================%
figure(2)
%Unnewehr Universal Model--------------------------------------------------
V=OCV_b1'; s=SOC_b1'; P=[One s]; k= inv(P'*P)*P'*V; % 'k' parameter calculation 
kp_unb1=k; OCV= kp_unb1(1) + kp_unb1(2)*(SOC); 
plot(SOC,OCV) %Plotting
hold on
%Shepherd model------------------------------------------------------------
V=OCV_b1'; ScaledS= SOC_b1*(1-2*E) + E; s=ScaledS'; P=[One 1./s]; 
k= inv(P'*P)*P'*V; % 'k' parameter calculation 
kp_shb1=k; OCV= kp_shb1(1) + kp_shb1(2)*(1./z); plot(SOC,OCV) %Plotting
% Nernst Model-------------------------------------------------------------
V=OCV_b1'; ScaledS= SOC_b1*(1-2*E) + E; s=ScaledS';
P=[One log(s) log(1-s)]; k= inv(P'*P)*P'*V; kp_NNb1=k;  % 'k' parameter calculation
OCV= kp_NNb1(1)  + kp_NNb1(2)*log(z) + kp_NNb1(3)*log(1-z); 
plot(SOC,OCV)%Plotting
% Combined Model-----------------------------------------------------------
V=OCV_b1'; ScaledS= SOC_b1*(1-2*E) + E; s=ScaledS';
P=[One 1./s s log(s) log(1-s)]; 
k= inv(P'*P)*P'*V; kp_cmb1=k; % 'k' parameter calculation 
OCV= kp_cmb1(1) + kp_cmb1(2)*(1./z) +kp_cmb1(3)*z + kp_cmb1(4)*log(z)...
    + kp_cmb1(5)*log(1-z);
plot(SOC,OCV) %Plotting
% Combined +3 Model--------------------------------------------------------
V=OCV_b1'; ScaledS= SOC_b1*(1-2*E) + E; s=ScaledS';
P=[One 1./s 1./(s.^2) 1./(s.^3) 1./(s.^4) s log(s) log(1-s)]; 
k= inv(P'*P)*P'*V; kp_cm3b1=k; % 'k' parameter calculation 
OCV = kp_cm3b1(1)+ kp_cm3b1(2)*(1./z) + kp_cm3b1(3)*(1./(z.^2))+ ...
    kp_cm3b1(4)*(1./(z.^3)) + kp_cm3b1(5)*(1./(z.^4))+ kp_cm3b1(6)*(z) ...
    + kp_cm3b1(7)*(log(z))+ kp_cm3b1(8)*(log(1-z));
plot(SOC,OCV)%Plotting
% Polynomial Model taking m=3 and n=2--------------------------------------
V=OCV_b1'; ScaledS= SOC_b1*(1-2*E) + E; s=ScaledS';
P=[One s s.^2 s.^3 1./(s) 1./(s.^2)]; k= inv(P'*P)*P'*V; kp_pob1=k; % 'k' parameter calculation
OCV= kp_pob1(1)+ kp_pob1(2)*(z) +kp_pob1(3)*(z.^2) + kp_pob1(4)*(z.^3)+ ...
    kp_pob1(5)*(1./z) + kp_pob1(6)*(1./(z.^2)); 
plot(SOC,OCV)%Plotting
% Exponential Model taking m=3 and n=2-------------------------------------
V=OCV_b1'; s=SOC_b1'; P=[One exp(s) exp(s.^2) exp(s.^3) exp(-s) exp(-(s.^2))]; 
k= inv(P'*P)*P'*V; kp_exb1=k; % 'k' parameter calculation
OCV= kp_exb1(1)+ kp_exb1(2)*exp(SOC) +kp_exb1(3)*exp(SOC.^2)...
    + kp_exb1(4)*exp(SOC.^3)+ kp_exb1(5)*exp(-SOC) + kp_exb1(6)*exp(-(SOC.^2));
plot(SOC,OCV) %Plotting
%Combined Plots------------------------------------------------------------
plot(SOC_b1,OCV_b1)% Plotting
legend("Unnewehr Universal Model", "Shepherd Model","Nernst Model",...
    "Combine Model","Combine +3 Model","Polynomial Model","Exponential Model","Given Data Model")
title("Results for Battery C1202");xlabel('(SOC)'); ylabel('(OCV)');

%================================================================ For Battery C1203 =============================================%
SOC = 0:.0001:1;
z=SOC*(1-2*E) + E;
figure(3)

%Unnewehr Universal Model
V=OCV_b2';
s=SOC_b2';
P=[One s];
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_unb2=k; %
OCV= kp_unb2(1) + kp_unb2(2)*(SOC);
%Plotting
plot(SOC,OCV)
hold on

%Sheperd Model
V=OCV_b2';
ScaledS= SOC_b2*(1-2*E) + E;
s=ScaledS';
P=[One 1./s]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_shb2=k;
OCV= kp_shb2(1) + kp_shb2(2)*(1./z);
%Plotting
plot(SOC,OCV)

%Nernst Model
V=OCV_b2';
ScaledS= SOC_b2*(1-2*E) + E;
s=ScaledS';
P=[One log(s) log(1-s)]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_NNb2=k;
OCV= kp_NNb2(1)  + kp_NNb2(2)*log(z) + kp_NNb2(3)*log(1-z);
%Plotting
plot(SOC,OCV)

%Combined Model
V=OCV_b2';
ScaledS= SOC_b2*(1-2*E) + E;
s=ScaledS';
P=[One 1./s s log(s) log(1-s)]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_cmb2=k;
OCV= kp_cmb2(1) + kp_cmb2(2)*(1./z) +kp_cmb2(3)*z + kp_cmb2(4)*log(z) + kp_cmb2(5)*log(1-z);
%Plotting
plot(SOC,OCV)

%Combined +3 Model 
V=OCV_b2';
ScaledS= SOC_b2*(1-2*E) + E;
s=ScaledS';
P=[One 1./s 1./(s.^2) 1./(s.^3) 1./(s.^4) s log(s) log(1-s)]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_cm3b2=k;
OCV = kp_cm3b2(1)+ kp_cm3b2(2)*(1./z) + kp_cm3b2(3)*(1./(z.^2))+ kp_cm3b2(4)*(1./(z.^3))...
    + kp_cm3b2(5)*(1./(z.^4))+ kp_cm3b2(6)*(z) + kp_cm3b2(7)*(log(z))+ kp_cm3b2(8)*(log(1-z));
%Plotting
plot(SOC,OCV)

%Polynomial Model
V=OCV_b2';
ScaledS= SOC_b2*(1-2*E) + E;
s=ScaledS';
P=[One s s.^2 s.^3 1./(s) 1./(s.^2)]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_pob2=k;
OCV= kp_pob2(1)+ kp_pob2(2)*(z) +kp_pob2(3)*(z.^2) + kp_pob2(4)*(z.^3)+ ...
    kp_pob2(5)*(1./z) + kp_pob2(6)*(1./(z.^2));
%Plotting
plot(SOC,OCV)

%Exponential Model
V=OCV_b2';
s=SOC_b2';
P=[One exp(s) exp(s.^2) exp(s.^3) exp(-s) exp(-(s.^2))]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_exb2=k;
OCV= kp_exb2(1)+ kp_exb2(2)*exp(SOC) +kp_exb2(3)*exp(SOC.^2) ...
    + kp_exb2(4)*exp(SOC.^3)+ kp_exb2(5)*exp(-SOC) + kp_exb2(6)*exp(-(SOC.^2));
%Plotting
plot(SOC,OCV)
%Combined Plot
plot(SOC_b2,OCV_b2) % Plotting
legend("Unnewehr Universal Model", "Shepherd Model","Nernst Model",...
    "Combine Model","Combine +3 Model","Polynomial Model","Exponential Model","Given Data Model")
title("Results for Battery C1203");
xlabel('(SOC)');
ylabel('(OCV)');

%================================================================ For Battery C1204 =============================================%
SOC = 0:.0001:1;
z=SOC*(1-2*E) + E;
figure(4)

%Unnewehr Universal Model
V=OCV_b3';
s=SOC_b3';
P=[One s];
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_unb3=k;
OCV= kp_unb3(1) + kp_unb3(2)*(SOC);
%Plotting
plot(SOC,OCV)
hold on

%Sheperd Model
V=OCV_b3';
ScaledS= SOC_b3*(1-2*E) + E;
s=ScaledS';
P=[One 1./s]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_shb3=k;
OCV= kp_shb3(1) + kp_shb3(2)*(1./z);
%Plotting
plot(SOC,OCV)

%Nernst Model
V=OCV_b3';
ScaledS= SOC_b3*(1-2*E) + E;
s=ScaledS';
P=[One log(s) log(1-s)]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_NNb3=k;
OCV= kp_NNb3(1)  + kp_NNb3(2)*log(z) + kp_NNb3(3)*log(1-z);
%Plotting
plot(SOC,OCV)

%Combined Model
V=OCV_b3';
ScaledS= SOC_b3*(1-2*E) + E;
s=ScaledS';
P=[One 1./s s log(s) log(1-s)]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_cmb3=k;
OCV= kp_cmb3(1) + kp_cmb3(2)*(1./z) +kp_cmb3(3)*z + kp_cmb3(4)*log(z)...
    + kp_cmb3(5)*log(1-z);
%Plotting
plot(SOC,OCV)

%Combined +3 Model
V=OCV_b3';
ScaledS= SOC_b3*(1-2*E) + E;
s=ScaledS';
P=[One 1./s 1./(s.^2) 1./(s.^3) 1./(s.^4) s log(s) log(1-s)]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_cm3b3=k;
OCV = kp_cm3b3(1)+ kp_cm3b3(2)*(1./z) + kp_cm3b3(3)*(1./(z.^2))+...
    kp_cm3b3(4)*(1./(z.^3)) + kp_cm3b3(5)*(1./(z.^4))+ kp_cm3b3(6)*(z)...
    + kp_cm3b3(7)*(log(z))+ kp_cm3b3(8)*(log(1-z));
%Plotting
plot(SOC,OCV)

%Polynomial Model
V=OCV_b3';
ScaledS= SOC_b3*(1-2*E) + E;
s=ScaledS';
P=[One s s.^2 s.^3 1./(s) 1./(s.^2)]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_pob3=k;
OCV= kp_pob3(1)+ kp_pob3(2)*(z) +kp_pob3(3)*(z.^2) + kp_pob3(4)*(z.^3)...
    + kp_pob3(5)*(1./z) + kp_pob3(6)*(1./(z.^2));
%Plotting
plot(SOC,OCV)

%Exponential Model
V=OCV_b3';
s=SOC_b3';
P=[One exp(s) exp(s.^2) exp(s.^3) exp(-s) exp(-(s.^2))]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_exb3=k;
OCV= kp_exb3(1)+ kp_exb3(2)*exp(SOC) +kp_exb3(3)*exp(SOC.^2)...
    + kp_exb3(4)*exp(SOC.^3)+ kp_exb3(5)*exp(-SOC) +...
    kp_exb3(6)*exp(-(SOC.^2));
%Plotting
plot(SOC,OCV)
%Combined Plot
plot(SOC_b3,OCV_b3) % Plotting
legend("Unnewehr Universal Model", " Shepherd Model "," Nernst Model ",...
    " Combine Model "," Combine +3 Model "," Polynomial Model ",...
    " Exponential Model "," Given Data Model ")
title("Results for Battery C1204");
xlabel('(SOC)');
ylabel('(OCV)');

%================================================================ For Battery C1205 =============================================%
SOC = 0:.0001:1;
z=SOC*(1-2*E) + E;
figure(5)

%Unnewehr Universal Model
V=OCV_b4';
s=SOC_b4';
P=[One s];
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_unb4=k;
OCV= kp_unb4(1) + kp_unb4(2)*(SOC); 
%Plotting
plot(SOC,OCV)
hold on

%Sheperd Model
V=OCV_b4';
ScaledS= SOC_b4*(1-2*E) + E;
s=ScaledS';
P=[One 1./s]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_shb4=k;
OCV= kp_shb4(1) + kp_shb4(2)*(1./z);
%Plotting
plot(SOC,OCV)

%Nernst Model
V=OCV_b4';
ScaledS= SOC_b4*(1-2*E) + E;
s=ScaledS';
P=[One log(s) log(1-s)]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_NNb4=k;
OCV= kp_NNb4(1)  + kp_NNb4(2)*log(z) + kp_NNb4(3)*log(1-z);
%Plotting
plot(SOC,OCV)

%Combined Model
V=OCV_b4';
ScaledS= SOC_b4*(1-2*E) + E;
s=ScaledS';
P=[One 1./s s log(s) log(1-s)]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_cmb4=k;
OCV= kp_cmb4(1) + kp_cmb4(2)*(1./z) +kp_cmb4(3)*z ...
    + kp_cmb4(4)*log(z) + kp_cmb4(5)*log(1-z);
%Plotting
plot(SOC,OCV)

%Combined 3+ Model
V=OCV_b4';
ScaledS= SOC_b4*(1-2*E) + E;
s=ScaledS';
P=[One 1./s 1./(s.^2) 1./(s.^3) 1./(s.^4) s log(s) log(1-s)]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_cm3b4=k;
OCV = kp_cm3b4(1)+ kp_cm3b4(2)*(1./z) + kp_cm3b4(3)*(1./(z.^2))...
    + kp_cm3b4(4)*(1./(z.^3)) + kp_cm3b4(5)*(1./(z.^4))+ kp_cm3b4(6)*(z)...
    + kp_cm3b4(7)*(log(z))+ kp_cm3b4(8)*(log(1-z));
%Plotting
plot(SOC,OCV)

%Polynomial Model
V=OCV_b4';
ScaledS= SOC_b4*(1-2*E) + E;
s=ScaledS';
P=[One s s.^2 s.^3 1./(s) 1./(s.^2)]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_pob4=k;
OCV= kp_pob4(1)+ kp_pob4(2)*(z) +kp_pob4(3)*(z.^2) + ...
    kp_pob4(4)*(z.^3)+ kp_pob4(5)*(1./z) + kp_pob4(6)*(1./(z.^2));
%Plotting
plot(SOC,OCV)

%Exponential Model
V=OCV_b4';
s=SOC_b4';
P=[One exp(s) exp(s.^2) exp(s.^3) exp(-s) exp(-(s.^2))]; 
k= inv(P'*P)*P'*V;
% 'k' parameter calculation
kp_exb4=k;
OCV= kp_exb4(1)+ kp_exb4(2)*exp(SOC) +kp_exb4(3)*exp(SOC.^2)...
    + kp_exb4(4)*exp(SOC.^3)+ kp_exb4(5)*exp(-SOC) + kp_exb4(6)*exp(-(SOC.^2));
%Plotting
plot(SOC,OCV)
%Combined Plot
plot(SOC_b4,OCV_b4)% Plotting
legend("Unnewehr Universal Model", "Shepherd Model","Nernst Model",...
    "Combine Model","Combine +3 Model","Polynomial Model","Exponential Model"...
    ,"Given Data Model")
title("Results for Battery C1205");
xlabel('(SOC)');
ylabel('(OCV)');

%================================================================================================================================%

%----------------------------------------------------------OCV Modelling Errors--------------------------------------------------%

%========================================================== For Battery C1202 ===================================================%
z=SOC_b1*(1-2*E) + E;
figure(6)
%Unnewehr Universal Model
subplot(3,3,1);
OCV= kp_unb1(1) + kp_unb1(2)*(SOC_b1);
T1= OCV_b1-OCV;plot(SOC_b1,T1)
title("Unnewhr");xlabel('(SOC)');ylabel('(OCV)');
%Shepherd Model
subplot(3,3,2);
OCV= kp_shb1(1) + kp_shb1(2)*(1./z);
T1= OCV_b1-OCV;plot(SOC_b1,T1)
title("Shepherd");xlabel('(SOC)');ylabel('(OCV)');
%Nernst Model
subplot(3,3,3);
OCV= kp_NNb1(1)  + kp_NNb1(2)*log(z) + kp_NNb1(3)*log(1-z);
T1= OCV_b1-OCV;plot(SOC_b1,T1)
title("Nernst");xlabel('(SOC)');ylabel('(OCV)');
%Combine Model
subplot(3,3,4);
OCV= kp_cmb1(1) + kp_cmb1(2)*(1./z) +kp_cmb1(3)*z + kp_cmb1(4)*log(z) + kp_cmb1(5)*log(1-z);
T1= OCV_b1-OCV;plot(SOC_b1,T1)
title("Combine");xlabel('(SOC)');ylabel('(OCV)');
%Combine +3 Model
subplot(3,3,5);
OCV = kp_cm3b1(1)+ kp_cm3b1(2)*(1./z) + kp_cm3b1(3)*(1./(z.^2))+...
    kp_cm3b1(4)*(1./(z.^3)) + kp_cm3b1(5)*(1./(z.^4))+ kp_cm3b1(6)*(z)...
    + kp_cm3b1(7)*(log(z))+ kp_cm3b1(8)*(log(1-z)); 
T1= OCV_b1-OCV;plot(SOC_b1,T1)
title("Combine +3");xlabel('(SOC)');
ylabel('(OCV)');
%Polynomial Model
subplot(3,3,6);
OCV= kp_pob1(1)+ kp_pob1(2)*(z) +kp_pob1(3)*(z.^2) + kp_pob1(4)*(z.^3)...
    + kp_pob1(5)*(1./z) + kp_pob1(6)*(1./(z.^2));
T1= OCV_b1-OCV;plot(SOC_b1,T1)
title("Polynomial");xlabel('(SOC)');ylabel('(OCV)');
%Exponential Model
subplot(3,3,8);
OCV= kp_exb1(1)+ kp_exb1(2)*exp(SOC_b1) +kp_exb1(3)*exp(SOC_b1.^2)...
    + kp_exb1(4)*exp(SOC_b1.^3)+ kp_exb1(5)*exp(-SOC_b1) + kp_exb1(6)*exp(-(SOC_b1.^2));
T1= OCV_b1-OCV;plot(SOC_b1,T1)
title("Exponential");xlabel('(SOC)');ylabel('(OCV)');
sgtitle("OCV Modelling Errors for Battery C1202")
%========================================================== For Battery C1203 ===================================================%
z=SOC_b2*(1-2*E) + E;
figure(7)
%Unnewehr Universal Model
subplot(3,3,1);
OCV= kp_unb2(1) + kp_unb2(2)*(SOC_b2);
T1= OCV_b2-OCV;plot(SOC_b2,T1)
title("Unnewhr");xlabel('(SOC)');ylabel('(OCV)');
%Shepherd Model
subplot(3,3,2);
OCV= kp_shb2(1) + kp_shb2(2)*(1./z);
T1= OCV_b2-OCV;plot(SOC_b2,T1)
title("Shepherd");xlabel('(SOC)');ylabel('(OCV)');
%Nernst Model
subplot(3,3,3);
OCV= kp_NNb2(1)  + kp_NNb2(2)*log(z) + kp_NNb2(3)*log(1-z);
T1= OCV_b2-OCV;plot(SOC_b2,T1)
title("Nernst");xlabel('(SOC)');ylabel('(OCV)');
%Combine Model
subplot(3,3,4);
OCV= kp_cmb2(1) + kp_cmb2(2)*(1./z) +kp_cmb2(3)*z + kp_cmb2(4)*log(z)...
    + kp_cmb2(5)*log(1-z);
T1= OCV_b2-OCV;plot(SOC_b2,T1)
title("Combine");xlabel('(SOC)');ylabel('(OCV)');
%Combine +3 Model
subplot(3,3,5);
OCV = kp_cm3b2(1)+ kp_cm3b2(2)*(1./z) + kp_cm3b2(3)*(1./(z.^2))...
    + kp_cm3b2(4)*(1./(z.^3)) + kp_cm3b2(5)*(1./(z.^4))+ kp_cm3b2(6)*(z)...
    + kp_cm3b2(7)*(log(z))+ kp_cm3b2(8)*(log(1-z)); 
T1= OCV_b2-OCV;
plot(SOC_b2,T1)
title("Combine +3");xlabel('(SOC)');ylabel('(OCV)');
%Polynomial Model
subplot(3,3,6);
OCV= kp_pob2(1)+ kp_pob2(2)*(z) +kp_pob2(3)*(z.^2) +...,
    kp_pob2(4)*(z.^3)+ kp_pob2(5)*(1./z) + kp_pob2(6)*(1./(z.^2));
T1= OCV_b2-OCV;plot(SOC_b2,T1)
title("Polynomial");xlabel('(SOC)');ylabel('(OCV)');
%Exponential Model
subplot(3,3,8);
OCV= kp_exb2(1)+ kp_exb2(2)*exp(SOC_b2) +kp_exb2(3)*exp(SOC_b2.^2)...
    + kp_exb2(4)*exp(SOC_b2.^3)+ kp_exb2(5)*exp(-SOC_b2) + kp_exb2(6)*exp(-(SOC_b2.^2));
T1= OCV_b2-OCV;plot(SOC_b2,T1)
title("Exponential");xlabel('(SOC)');ylabel('(OCV)');sgtitle("OCV Modelling Errors for Battery C1203")
%========================================================== For Battery C1204 ===================================================%
z=SOC_b3*(1-2*E) + E;
figure(8)
%Unnewehr Universal Model
subplot(3,3,1);
OCV= kp_unb3(1) + kp_unb3(2)*(SOC_b3);
T1= OCV_b3-OCV;plot(SOC_b3,T1)
title("Unnewhr");xlabel('(SOC)');ylabel('(OCV)');
%Shepherd Model
subplot(3,3,2);
OCV= kp_shb3(1) + kp_shb3(2)*(1./z);
T1= OCV_b3-OCV;plot(SOC_b3,T1)
title("Shepherd");xlabel('(SOC)');ylabel('(OCV)');
%Nernst Model
subplot(3,3,3);
OCV= kp_NNb3(1)  + kp_NNb3(2)*log(z) + kp_NNb3(3)*log(1-z);
T1= OCV_b3-OCV;plot(SOC_b3,T1)
title("Nernst");xlabel('(SOC)');ylabel('(OCV)');
%Combine Model
subplot(3,3,4);
OCV= kp_cmb3(1) + kp_cmb3(2)*(1./z) +kp_cmb3(3)*z + kp_cmb3(4)*log(z)...
    + kp_cmb3(5)*log(1-z);
T1= OCV_b3-OCV;plot(SOC_b3,T1)
title("Combine");xlabel('(SOC)');ylabel('(OCV)');
%Combine +3 Model
subplot(3,3,5);
OCV = kp_cm3b3(1)+ kp_cm3b3(2)*(1./z) + kp_cm3b3(3)*(1./(z.^2))...
    + kp_cm3b3(4)*(1./(z.^3)) + kp_cm3b3(5)*(1./(z.^4))+ kp_cm3b3(6)*(z)...
    + kp_cm3b3(7)*(log(z))+ kp_cm3b3(8)*(log(1-z)); 
T1= OCV_b3-OCV;plot(SOC_b3,T1)
title("Combine +3");xlabel('(SOC)');ylabel('(OCV)');
%Polynomial Model
subplot(3,3,6);
OCV= kp_pob3(1)+ kp_pob3(2)*(z) +kp_pob3(3)*(z.^2) + kp_pob3(4)*(z.^3)...
    + kp_pob3(5)*(1./z) + kp_pob3(6)*(1./(z.^2));
T1= OCV_b3-OCV;plot(SOC_b3,T1)
title("Polynomial");xlabel('(SOC)');ylabel('(OCV)');
%Exponential Model
subplot(3,3,8);
OCV= kp_exb3(1)+ kp_exb3(2)*exp(SOC_b3) +kp_exb3(3)*exp(SOC_b3.^2)...
    + kp_exb3(4)*exp(SOC_b3.^3)+ kp_exb3(5)*exp(-SOC_b3) + kp_exb3(6)*exp(-(SOC_b3.^2));
T1= OCV_b3-OCV;plot(SOC_b3,T1)
title("Exponential");xlabel('(SOC)');ylabel('(OCV)');
sgtitle("OCV Modelling Errors for Battery C1204")
%========================================================== For Battery C1205 ===================================================%

z=SOC_b4*(1-2*E) + E;
figure(9)
%Unnewehr Universal Model
subplot(3,3,1);
OCV= kp_unb4(1) + kp_unb4(2)*(SOC_b4);
T1= OCV_b4-OCV;plot(SOC_b4,T1)
title("Unnewhr");xlabel('(SOC)');ylabel('(OCV)');
%Shepherd Model
subplot(3,3,2);
OCV= kp_shb4(1) + kp_shb4(2)*(1./z);
T1= OCV_b4-OCV;plot(SOC_b4,T1)
title("Shepherd");xlabel('(SOC)');ylabel('(OCV)');
%Nernst Model
subplot(3,3,3);
OCV= kp_NNb4(1)  + kp_NNb4(2)*log(z) + kp_NNb4(3)*log(1-z);
T1= OCV_b4-OCV;plot(SOC_b4,T1)
title("Nernst");xlabel('(SOC)');ylabel('(OCV)');
%Combine Model
subplot(3,3,4);
OCV= kp_cmb4(1) + kp_cmb4(2)*(1./z) +kp_cmb4(3)*z + kp_cmb4(4)*log(z)...
    + kp_cmb4(5)*log(1-z);
T1= OCV_b4-OCV;plot(SOC_b4,T1)
title("Combine");xlabel('(SOC)');ylabel('(OCV)');
%Combine +3 Model
subplot(3,3,5);
OCV = kp_cm3b4(1)+ kp_cm3b4(2)*(1./z) + kp_cm3b4(3)*(1./(z.^2))+...
    kp_cm3b4(4)*(1./(z.^3)) + kp_cm3b4(5)*(1./(z.^4))+ kp_cm3b4(6)*(z)...
    + kp_cm3b4(7)*(log(z))+ kp_cm3b4(8)*(log(1-z)); 
T1= OCV_b4-OCV;plot(SOC_b4,T1)
title("Combine +3");xlabel('(SOC)');ylabel('(OCV)');
%Polynomial Model
subplot(3,3,6);
OCV= kp_pob4(1)+ kp_pob4(2)*(z) +kp_pob4(3)*(z.^2) + kp_pob4(4)*(z.^3)...
    + kp_pob4(5)*(1./z) + kp_pob4(6)*(1./(z.^2));
T1= OCV_b4-OCV;plot(SOC_b4,T1)
title("Polynomial");xlabel('(SOC)');ylabel('(OCV)');
%Exponential Model
subplot(3,3,8);
OCV= kp_exb4(1)+ kp_exb4(2)*exp(SOC_b4) +kp_exb4(3)*exp(SOC_b4.^2)...
    + kp_exb4(4)*exp(SOC_b4.^3)+ kp_exb4(5)*exp(-SOC_b4) + kp_exb4(6)*exp(-(SOC_b4.^2));
T1= OCV_b4-OCV;plot(SOC_b4,T1) % Plotting
title("Exponential");xlabel('(SOC)');ylabel('(OCV)');
sgtitle("OCV Modelling Errors for Battery C1205")

%================================================================================================================================%

%---------------------------------------------------------  Error Calculation  --------------------------------------------------%

%========================================================== For Battery C1202 ===================================================%

%Best Fit
%Unnewehr Universal Model
Vmean = 0;sum_a=0;sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b1(i);
end
Vmean=Vmean/32;
Vhat= kp_unb1(1) + kp_unb1(2)*(SOC_b1);
a=abs(Vhat-OCV_b1);b=abs(OCV_b1-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_U_1202=1-(sum_a/sum_b);
%Shepherd Model
Vmean = 0;sum_a=0;sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b1(i);
end
Vmean=Vmean/32; z=SOC_b1*(1-2*E) + E;
Vhat= kp_shb1(1) + kp_shb1(2)*(1./z);
a=abs(Vhat-OCV_b1); b=abs(OCV_b1-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_S_1202=1-(sum_a/sum_b);
%Nernst Model
Vmean = 0;sum_a=0;sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b1(i);
end
Vmean=Vmean/32;z=SOC_b1*(1-2*E) + E;
Vhat= kp_NNb1(1) + kp_NNb1(2)*log(z) + kp_NNb1(3)*log(1-z);
a=abs(Vhat-OCV_b1);b=abs(OCV_b1-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_nn_1202=1-(sum_a/sum_b);
%Combine Model
Vmean = 0;sum_a=0;sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b1(i);
end
Vmean=Vmean/32; z=SOC_b1*(1-2*E) + E;
Vhat= kp_cmb1(1) + kp_cmb1(2)*(1./z) +kp_cmb1(3)*z...
    + kp_cmb1(4)*log(z) + kp_cmb1(5)*log(1-z);
a=abs(Vhat-OCV_b1); b=abs(OCV_b1-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_cm_1202=1-(sum_a/sum_b);
%Combine +3 Model
Vmean = 0;sum_a=0;sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b1(i);
end
Vmean=Vmean/32;z=SOC_b1*(1-2*E) + E;
Vhat= kp_cm3b1(1)+ kp_cm3b1(2)*(1./z) + kp_cm3b1(3)*(1./(z.^2))...
    + kp_cm3b1(4)*(1./(z.^3)) + kp_cm3b1(5)*(1./(z.^4))+...
    kp_cm3b1(6)*(z) + kp_cm3b1(7)*(log(z))+ kp_cm3b1(8)*(log(1-z));
a=abs(Vhat-OCV_b1); b=abs(OCV_b1-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_cm3_1202=1-(sum_a/sum_b);
%Polynomial Model
Vmean = 0;sum_a=0;sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b1(i);
end
Vmean=Vmean/32; z=SOC_b1*(1-2*E) + E;
Vhat= kp_pob1(1)+ kp_pob1(2)*(z) +kp_pob1(3)*(z.^2)...
    + kp_pob1(4)*(z.^3)+ kp_pob1(5)*(1./z) + kp_pob1(6)*(1./(z.^2));
a=abs(Vhat-OCV_b1); b=abs(OCV_b1-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_pol_1202=1-(sum_a/sum_b);
%Exponential Model
Vmean = 0; sum_a=0; sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b1(i);
end
Vmean=Vmean/32; z=SOC_b1*(1-2*E) + E;
Vhat= kp_exb1(1)+ kp_exb1(2)*exp(SOC_b1)...
    +kp_exb1(3)*exp(SOC_b1.^2) + kp_exb1(4)*exp(SOC_b1.^3)...
    + kp_exb1(5)*exp(-SOC_b1) + kp_exb1(6)*exp(-(SOC_b1.^2));
a=abs(Vhat-OCV_b1); b=abs(OCV_b1-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_ex_1202=1-(sum_a/sum_b);
%---------------------------------------------------------------------------------------------
%R-sqaure
%Unnewehr Universal Model
Vmean = 0;sum_a=0;sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b1(i);
end
Vmean=Vmean/32;
Vhat= kp_unb1(1) + kp_unb1(2)*(SOC_b1);
a=(abs(Vhat-OCV_b1)).^2; b=(abs(OCV_b1-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_U_1202=1-(sum_a/sum_b);
%Shepherd Model
Vmean = 0;sum_a=0;sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b1(i);
end
Vmean=Vmean/32; z=SOC_b1*(1-2*E) + E;
Vhat= kp_shb1(1) + kp_shb1(2)*(1./z);
a=(abs(Vhat-OCV_b1)).^2; b=(abs(OCV_b1-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_S_1202=1-(sum_a/sum_b);
%Nernst Model
Vmean = 0; sum_a=0; sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b1(i);
end
Vmean=Vmean/32; z=SOC_b1*(1-2*E) + E;
Vhat= kp_NNb1(1)  + kp_NNb1(2)*log(z) + kp_NNb1(3)*log(1-z);
a=(abs(Vhat-OCV_b1)).^2; b=(abs(OCV_b1-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_nn_1202=1-(sum_a/sum_b);
%Combine Model
Vmean = 0; sum_a=0; sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b1(i);
end
Vmean=Vmean/32; z=SOC_b1*(1-2*E) + E;
Vhat= kp_cmb1(1) + kp_cmb1(2)*(1./z) +kp_cmb1(3)*z + kp_cmb1(4)*log(z) + kp_cmb1(5)*log(1-z);
a=(abs(Vhat-OCV_b1)).^2; b=(abs(OCV_b1-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_cm_1202=1-(sum_a/sum_b);
%Combine +3 Model

Vmean = 0; sum_a=0; sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b1(i);
end
Vmean=Vmean/32; z=SOC_b1*(1-2*E) + E;
Vhat= kp_cm3b1(1)+ kp_cm3b1(2)*(1./z) + kp_cm3b1(3)*(1./(z.^2))...
    + kp_cm3b1(4)*(1./(z.^3)) + kp_cm3b1(5)*(1./(z.^4))+ kp_cm3b1(6)*(z) ...
    + kp_cm3b1(7)*(log(z))+ kp_cm3b1(8)*(log(1-z));
a=(abs(Vhat-OCV_b1)).^2; b=(abs(OCV_b1-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_cm3_1202=1-(sum_a/sum_b);
%Polynomial Model
Vmean = 0; sum_a=0; sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b1(i);
end
Vmean=Vmean/32; z=SOC_b1*(1-2*E) + E;
Vhat= kp_pob1(1)+ kp_pob1(2)*(z) +kp_pob1(3)*(z.^2) +...
    kp_pob1(4)*(z.^3)+ kp_pob1(5)*(1./z) + kp_pob1(6)*(1./(z.^2));
a=(abs(Vhat-OCV_b1)).^2; b=(abs(OCV_b1-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_pol_1202=1-(sum_a/sum_b);
%Exponential Model
Vmean = 0; sum_a=0; sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b1(i);
end
Vmean=Vmean/32; z=SOC_b1*(1-2*E) + E;
Vhat= kp_exb1(1)+ kp_exb1(2)*exp(SOC_b1) +kp_exb1(3)*exp(SOC_b1.^2) ...
    + kp_exb1(4)*exp(SOC_b1.^3)+ kp_exb1(5)*exp(-SOC_b1) + kp_exb1(6)*exp(-(SOC_b1.^2));
a=(abs(Vhat-OCV_b1)).^2; b=(abs(OCV_b1-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_ex_1202=1-(sum_a/sum_b);
%---------------------------------------------------------------------------------------------
%Max Error
z=SOC_b1*(1-2*E) + E;
%Unnewehr Universal Model
Vhat= kp_unb1(1) + kp_unb1(2)*(SOC_b1);
ME=abs(OCV_b1-Vhat);
max_error_u_1202=max(ME);
%Shepherd Model
Vhat= kp_shb1(1) + kp_shb1(2)*(1./z);
ME=abs(OCV_b1-Vhat);
max_error_s_1202=max(ME);
%Nernst Model
Vhat= kp_NNb1(1)  + kp_NNb1(2)*log(z) + kp_NNb1(3)*log(1-z);
ME=abs(OCV_b1-Vhat);
max_error_nn_1202=max(ME);
%Combine Model
Vhat= kp_cmb1(1) + kp_cmb1(2)*(1./z) +kp_cmb1(3)*z + kp_cmb1(4)*log(z)...
    + kp_cmb1(5)*log(1-z);
ME=abs(OCV_b1-Vhat);
max_error_cm_1202=max(ME);
%Combine +3 Model
Vhat= kp_cm3b1(1)+ kp_cm3b1(2)*(1./z) + kp_cm3b1(3)*(1./(z.^2))...
    + kp_cm3b1(4)*(1./(z.^3)) + kp_cm3b1(5)*(1./(z.^4))+ kp_cm3b1(6)*(z)...
    + kp_cm3b1(7)*(log(z))+ kp_cm3b1(8)*(log(1-z));
ME=abs(OCV_b1-Vhat);
max_error_cm3_1202=max(ME);
%Polynomial Model
Vhat= kp_pob1(1)+ kp_pob1(2)*(z) +kp_pob1(3)*(z.^2) +...
    kp_pob1(4)*(z.^3)+ kp_pob1(5)*(1./z) + kp_pob1(6)*(1./(z.^2));
ME=abs(OCV_b1-Vhat);
max_error_pol_1202=max(ME);
%Exponential Model
Vhat= kp_exb1(1)+ kp_exb1(2)*exp(SOC_b1) +kp_exb1(3)*exp(SOC_b1.^2) ...
    + kp_exb1(4)*exp(SOC_b1.^3)+ kp_exb1(5)*exp(-SOC_b1) + kp_exb1(6)*exp(-(SOC_b1.^2));
ME=abs(OCV_b1-Vhat);
max_error_ex_1202=max(ME);

%============================================================================================
%Root Mean Square Error
N=32;
z=SOC_b1*(1-2*E) + E;
%Unnewehr Universal Model
sum_a=0;
Vhat= kp_unb1(1) + kp_unb1(2)*(SOC_b1);
a=abs(OCV_b1-Vhat);
M=2;
for i=1:32
    sum_a= sum_a+ a(i);  
end
rmse_u_1202=sum_a/sqrt(N-M);

%Shepherd Model
sum_a=0;
Vhat= kp_shb1(1) + kp_shb1(2)*(1./z);
a=abs(OCV_b1-Vhat);
M=2;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_sh_1202=sum_a/sqrt(N-M);

%Nernst Model
sum_a=0;
Vhat= kp_NNb1(1)  + kp_NNb1(2)*log(z) + kp_NNb1(3)*log(1-z);
a=abs(OCV_b1-Vhat);
M=3;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_nn_1202=sum_a/sqrt(N-M);
%Combine Model
sum_a=0;
Vhat= kp_cmb1(1) + kp_cmb1(2)*(1./z) +kp_cmb1(3)*z +...
    kp_cmb1(4)*log(z) + kp_cmb1(5)*log(1-z);
a=abs(OCV_b1-Vhat);
M=5;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_cm_1202=sum_a/sqrt(N-M);

%Combine +3 Model
sum_a=0;
Vhat= kp_cm3b1(1)+ kp_cm3b1(2)*(1./z) + kp_cm3b1(3)*(1./(z.^2))+...
    kp_cm3b1(4)*(1./(z.^3)) + kp_cm3b1(5)*(1./(z.^4))+ kp_cm3b1(6)*(z)...
    + kp_cm3b1(7)*(log(z))+ kp_cm3b1(8)*(log(1-z));
a=abs(OCV_b1-Vhat);
M=8;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_cm3_1202=sum_a/sqrt(N-M);

%Polynomial Model
sum_a=0;
Vhat= kp_pob1(1)+ kp_pob1(2)*(z) +kp_pob1(3)*(z.^2) + kp_pob1(4)*(z.^3)...
    + kp_pob1(5)*(1./z) + kp_pob1(6)*(1./(z.^2));
a=abs(OCV_b1-Vhat);
M=6;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_pol_1202=sum_a/sqrt(N-M);

%Exponential Model
sum_a=0;
Vhat= kp_exb1(1)+ kp_exb1(2)*exp(SOC_b1) +kp_exb1(3)*exp(SOC_b1.^2) ...
    + kp_exb1(4)*exp(SOC_b1.^3)+ kp_exb1(5)*exp(-SOC_b1) + kp_exb1(6)*exp(-(SOC_b1.^2));
a=abs(OCV_b1-Vhat);
M=6;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_ex_1202=sum_a/sqrt(N-M);
%=========================================================================================================
%AIC
z=SOC_b1*(1-2*E) + E;
%Unnewehr Universal Model
Vhat= kp_unb1(1) + kp_unb1(2)*(SOC_b1);
S2=0;
e=OCV_b1-Vhat;
M=2;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_u_1202=32*log(S2./32) + 2*(M+1);
%Shepherd Model
Vhat= kp_shb1(1) + kp_shb1(2)*(1./z);
S2=0;
e=OCV_b1-Vhat;
M=2;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_s_1202=32*log(S2./32) + 2*(M+1);
%Nernst Model
Vhat= kp_NNb1(1)  + kp_NNb1(2)*log(z) + kp_NNb1(3)*log(1-z);
S2=0;
e=OCV_b1-Vhat;
M=3;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_nn_1202=32*log(S2./32) + 2*(M+1);

%Combine Model
Vhat= kp_cmb1(1) + kp_cmb1(2)*(1./z) +kp_cmb1(3)*z +...
    kp_cmb1(4)*log(z) + kp_cmb1(5)*log(1-z);
S2=0;
e=OCV_b1-Vhat;
M=5;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_cm_1202=32*log(S2./32) + 2*(M+1);

%Combine +3 Model
Vhat= kp_cm3b1(1)+ kp_cm3b1(2)*(1./z) + kp_cm3b1(3)*(1./(z.^2))+ kp_cm3b1(4)*(1./(z.^3))...
    + kp_cm3b1(5)*(1./(z.^4))+ kp_cm3b1(6)*(z) + kp_cm3b1(7)*(log(z))+ kp_cm3b1(8)*(log(1-z));
S2=0;
e=OCV_b1-Vhat;
M=8;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_cm3_1202=32*log(S2./32) + 2*(M+1);

%Polynomial Model
Vhat= kp_pob1(1)+ kp_pob1(2)*(z) +kp_pob1(3)*(z.^2) +...
    kp_pob1(4)*(z.^3)+ kp_pob1(5)*(1./z) + kp_pob1(6)*(1./(z.^2));
S2=0;
e=OCV_b1-Vhat;
M=6;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_pol_1202=32*log(S2./32) + 2*(M+1);

%Exponential Model
Vhat= kp_exb1(1)+ kp_exb1(2)*exp(SOC_b1) +...
    kp_exb1(3)*exp(SOC_b1.^2) + kp_exb1(4)*exp(SOC_b1.^3)+ kp_exb1(5)*exp(-SOC_b1)...
    + kp_exb1(6)*exp(-(SOC_b1.^2));
S2=0;
e=OCV_b1-Vhat;
M=6;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_ex_1202=32*log(S2./32) + 2*(M+1);

%=========================================================================================================

%Creating Metrics Table & Ranking Table - Battery C1202

AIC=[aic_u_1202;aic_s_1202;aic_nn_1202;aic_cm_1202;aic_cm3_1202;aic_pol_1202;aic_ex_1202];
RMSE=[rmse_u_1202;rmse_sh_1202;rmse_nn_1202;rmse_cm_1202;rmse_cm3_1202;rmse_pol_1202;rmse_ex_1202];
RSquare=[R2_U_1202;R2_S_1202;R2_nn_1202;R2_cm_1202;R2_cm3_1202;R2_pol_1202;R2_ex_1202];
BestFit=[BF_U_1202;BF_S_1202;BF_nn_1202;BF_cm_1202;BF_cm3_1202;BF_pol_1202;BF_ex_1202];
MaxError=[max_error_u_1202;max_error_s_1202;max_error_nn_1202;max_error_cm_1202;max_error_cm3_1202;max_error_pol_1202;max_error_ex_1202];
OCVModel=["Unnewhr";"Shepherd";"Nernst";"Combine";"Combine +3";"Polynomial";"Exponential"];
Metrics_table_1202 = table(OCVModel,AIC,RMSE,RSquare,BestFit,MaxError)
AIC=[1;2;3;4;5;6;7];
RMSE=[1;2;3;4;5;7;6];
RSquare=[1;2;3;4;5;6;7];
BestFit=[1;2;3;4;5;7;6];
MaxError=[1;2;3;4;6;5;7];
BordaRanking=[1;2;3;4;5;6;7];
OCVModel=["Combine +3";"Polynomial";"Combine";"Exponential";"Nernst";"Shepherd";"Unnewhr"];
ranking_table_1202 = table(OCVModel,AIC,RMSE,RSquare,BestFit,MaxError,BordaRanking)

%========================================================== For Battery C1203 ====================================================%

%Best Fit

%Unnewehr Universal Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean + OCV_b2(i);
end
Vmean=Vmean/32;
Vhat= kp_unb2(1) + kp_unb2(2)*(SOC_b2);
a=abs(Vhat-OCV_b2);
b=abs(OCV_b2-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_U_1203=1-(sum_a/sum_b);

%Shepherd Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b2(i);
end
Vmean=Vmean/32;
z=SOC_b2*(1-2*E) + E;
Vhat= kp_shb2(1) + kp_shb2(2)*(1./z);
a=abs(Vhat-OCV_b2);
b=abs(OCV_b2-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_S_1203=1-(sum_a/sum_b);

%Nernst Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b2(i);
end
Vmean=Vmean/32;
z=SOC_b2*(1-2*E) + E;
Vhat= kp_NNb2(1)  + kp_NNb2(2)*log(z) + kp_NNb2(3)*log(1-z);
a=abs(Vhat-OCV_b2);
b=abs(OCV_b2-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_nn_1203=1-(sum_a/sum_b);

%Combine Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b2(i);
end
Vmean=Vmean/32;
z=SOC_b2*(1-2*E) + E;
Vhat= kp_cmb2(1) + kp_cmb2(2)*(1./z) +kp_cmb2(3)*z + kp_cmb2(4)*log(z) + kp_cmb2(5)*log(1-z);
a=abs(Vhat-OCV_b2);
b=abs(OCV_b2-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_cm_1203=1-(sum_a/sum_b);

%Combine +3 Model

Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b2(i);
end
Vmean=Vmean/32;
z=SOC_b2*(1-2*E) + E;
Vhat= kp_cm3b2(1)+ kp_cm3b2(2)*(1./z) + kp_cm3b2(3)*(1./(z.^2))+...
    kp_cm3b2(4)*(1./(z.^3)) + kp_cm3b2(5)*(1./(z.^4))+ kp_cm3b2(6)*(z)...
    + kp_cm3b2(7)*(log(z))+ kp_cm3b2(8)*(log(1-z));
a=abs(Vhat-OCV_b2);
b=abs(OCV_b2-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_cm3_1203=1-(sum_a/sum_b);

%Polynomial Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b2(i);
end
Vmean=Vmean/32;
z=SOC_b2*(1-2*E) + E;
Vhat= kp_pob2(1)+ kp_pob2(2)*(z) +kp_pob2(3)*(z.^2) +...
    kp_pob2(4)*(z.^3)+ kp_pob2(5)*(1./z) + kp_pob2(6)*(1./(z.^2));
a=abs(Vhat-OCV_b2);
b=abs(OCV_b2-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_pol_1203=1-(sum_a/sum_b);

%Exponential Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b2(i);
end
Vmean=Vmean/32;
z=SOC_b2*(1-2*E) + E;
Vhat= kp_exb2(1)+ kp_exb2(2)*exp(SOC_b2) +kp_exb2(3)*exp(SOC_b2.^2)...
    + kp_exb2(4)*exp(SOC_b2.^3)+ kp_exb2(5)*exp(-SOC_b2) + kp_exb2(6)*exp(-(SOC_b2.^2));
a=abs(Vhat-OCV_b2);
b=abs(OCV_b2-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_ex_1203=1-(sum_a/sum_b);
%=========================================================================================================
%R-sqaure
%Unnewehr Universal Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b2(i);
end
Vmean=Vmean/32;
Vhat= kp_unb2(1) + kp_unb2(2)*(SOC_b2);
a=(abs(Vhat-OCV_b2)).^2;
b=(abs(OCV_b2-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_U_1203=1-(sum_a/sum_b);

%Shepherd Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b2(i);
end
Vmean=Vmean/32;
z=SOC_b2*(1-2*E) + E;
Vhat= kp_shb2(1) + kp_shb2(2)*(1./z);
a=(abs(Vhat-OCV_b2)).^2;
b=(abs(OCV_b2-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_S_1203=1-(sum_a/sum_b);

%Nernst Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b2(i);
end
Vmean=Vmean/32;
z=SOC_b2*(1-2*E) + E;
Vhat= kp_NNb2(1)  + kp_NNb2(2)*log(z) + kp_NNb2(3)*log(1-z);
a=(abs(Vhat-OCV_b2)).^2;
b=(abs(OCV_b2-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_nn_1203=1-(sum_a/sum_b);

%Combine Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b2(i);
end
Vmean=Vmean/32;
z=SOC_b2*(1-2*E) + E;
Vhat= kp_cmb2(1) + kp_cmb2(2)*(1./z) +kp_cmb2(3)*z +...
    kp_cmb2(4)*log(z) + kp_cmb2(5)*log(1-z);
a=(abs(Vhat-OCV_b2)).^2;
b=(abs(OCV_b2-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_cm_1203=1-(sum_a/sum_b);

%Combine +3 Model

Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b2(i);
end
Vmean=Vmean/32;
z=SOC_b2*(1-2*E) + E;
Vhat= kp_cm3b2(1)+ kp_cm3b2(2)*(1./z) + kp_cm3b2(3)*(1./(z.^2))...
    + kp_cm3b2(4)*(1./(z.^3)) + kp_cm3b2(5)*(1./(z.^4))+ kp_cm3b2(6)*(z) ...
    + kp_cm3b2(7)*(log(z))+ kp_cm3b2(8)*(log(1-z));
a=(abs(Vhat-OCV_b2)).^2;
b=(abs(OCV_b2-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_cm3_1203=1-(sum_a/sum_b);

%Polynomial Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b2(i);
end
Vmean=Vmean/32;
z=SOC_b2*(1-2*E) + E;
Vhat= kp_pob2(1)+ kp_pob2(2)*(z) +kp_pob2(3)*(z.^2) +...
    kp_pob2(4)*(z.^3)+ kp_pob2(5)*(1./z) + kp_pob2(6)*(1./(z.^2));
a=(abs(Vhat-OCV_b2)).^2;
b=(abs(OCV_b2-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_pol_1203=1-(sum_a/sum_b);

%Exponential Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b2(i);
end
Vmean=Vmean/32;
z=SOC_b2*(1-2*E) + E;
Vhat= kp_exb2(1)+ kp_exb2(2)*exp(SOC_b2) +...
    kp_exb2(3)*exp(SOC_b2.^2) + kp_exb2(4)*exp(SOC_b2.^3)+...
    kp_exb2(5)*exp(-SOC_b2) + kp_exb2(6)*exp(-(SOC_b2.^2));
a=(abs(Vhat-OCV_b2)).^2;
b=(abs(OCV_b2-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_ex_1203=1-(sum_a/sum_b);
%=========================================================================================================
%Max Error

z=SOC_b2*(1-2*E) + E;
%Unnewehr Universal Model
Vhat= kp_unb2(1) + kp_unb2(2)*(SOC_b2);
ME=abs(OCV_b2-Vhat);
max_error_u_1203=max(ME);
%Shepherd Model
Vhat= kp_shb2(1) + kp_shb2(2)*(1./z);
ME=abs(OCV_b2-Vhat);
max_error_s_1203=max(ME);
%Nernst Model
Vhat= kp_NNb2(1)  + kp_NNb2(2)*log(z) + kp_NNb2(3)*log(1-z);
ME=abs(OCV_b2-Vhat);
max_error_nn_1203=max(ME);
%Combine Model
Vhat= kp_cmb2(1) + kp_cmb2(2)*(1./z) +kp_cmb2(3)*z +...
    kp_cmb2(4)*log(z) + kp_cmb2(5)*log(1-z);
ME=abs(OCV_b2-Vhat);
max_error_cm_1203=max(ME);
%Combine +3 Model
Vhat= kp_cm3b2(1)+ kp_cm3b2(2)*(1./z) + kp_cm3b2(3)*(1./(z.^2))+...
    kp_cm3b2(4)*(1./(z.^3)) + kp_cm3b2(5)*(1./(z.^4))+ kp_cm3b2(6)*(z)...
    + kp_cm3b2(7)*(log(z))+ kp_cm3b2(8)*(log(1-z));
ME=abs(OCV_b2-Vhat);
max_error_cm3_1203=max(ME);
%Polynomial Model
Vhat= kp_pob2(1)+ kp_pob2(2)*(z) +kp_pob2(3)*(z.^2) + kp_pob2(4)*(z.^3)...
    + kp_pob2(5)*(1./z) + kp_pob2(6)*(1./(z.^2));
ME=abs(OCV_b2-Vhat);
max_error_pol_1203=max(ME);
%Exponential Model
Vhat= kp_exb2(1)+ kp_exb2(2)*exp(SOC_b2) +kp_exb2(3)*exp(SOC_b2.^2)...
    + kp_exb2(4)*exp(SOC_b2.^3)+ kp_exb2(5)*exp(-SOC_b2) + kp_exb2(6)*exp(-(SOC_b2.^2));
ME=abs(OCV_b2-Vhat);
max_error_ex_1203=max(ME);

%=========================================================================================================
%Root Mean Square Error

N=32;
z=SOC_b2*(1-2*E) + E;
%Unnewehr Universal Model
sum_a=0;
Vhat= kp_unb2(1) + kp_unb2(2)*(SOC_b2);
a=abs(OCV_b2-Vhat);
M=2;
for i=1:32
    sum_a= sum_a+ a(i);  
end
rmse_u_1203=sum_a/sqrt(N-M);

%Shepherd Model
sum_a=0;
Vhat= kp_shb2(1) + kp_shb2(2)*(1./z);
a=abs(OCV_b2-Vhat);
M=2;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_s_1203=sum_a/sqrt(N-M);

%Nernst Model
sum_a=0;
Vhat= kp_NNb2(1)  + kp_NNb2(2)*log(z) + kp_NNb2(3)*log(1-z);
a=abs(OCV_b2-Vhat);
M=3;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_nn_1203=sum_a/sqrt(N-M);
%Combine Model
sum_a=0;
Vhat= kp_cmb2(1) + kp_cmb2(2)*(1./z) +kp_cmb2(3)*z + kp_cmb2(4)*log(z) + kp_cmb2(5)*log(1-z);
a=abs(OCV_b2-Vhat);
M=5;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_cm_1203=sum_a/sqrt(N-M);

%Combine +3 Model
sum_a=0;
Vhat= kp_cm3b2(1)+ kp_cm3b2(2)*(1./z) + kp_cm3b2(3)*(1./(z.^2))...
    + kp_cm3b2(4)*(1./(z.^3)) + kp_cm3b2(5)*(1./(z.^4))+ kp_cm3b2(6)*(z)...
    + kp_cm3b2(7)*(log(z))+ kp_cm3b2(8)*(log(1-z));
a=abs(OCV_b2-Vhat);
M=8;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_cm3_1203=sum_a/sqrt(N-M);

%Polynomial Model
sum_a=0;
Vhat= kp_pob2(1)+ kp_pob2(2)*(z) +kp_pob2(3)*(z.^2)...
    + kp_pob2(4)*(z.^3)+ kp_pob2(5)*(1./z) + kp_pob2(6)*(1./(z.^2));
a=abs(OCV_b2-Vhat);
M=6;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_pol_1203=sum_a/sqrt(N-M);

%Exponential Model
sum_a=0;
Vhat= kp_exb2(1)+ kp_exb2(2)*exp(SOC_b2) +kp_exb2(3)*exp(SOC_b2.^2)...
    + kp_exb2(4)*exp(SOC_b2.^3)+ kp_exb2(5)*exp(-SOC_b2) + kp_exb2(6)*exp(-(SOC_b2.^2));
a=abs(OCV_b2-Vhat);
M=6;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_ex_1203=sum_a/sqrt(N-M);
%=========================================================================================================
%AIC
z=SOC_b2*(1-2*E) + E;
%Unnewehr Universal Model
Vhat= kp_unb2(1) + kp_unb2(2)*(SOC_b2);
S2=0;
e=OCV_b2-Vhat;
M=2;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_u_1203=32*log(S2./32) + 2*(M+1);

%Shepherd Model
Vhat= kp_shb2(1) + kp_shb2(2)*(1./z);
S2=0;
e=OCV_b2-Vhat;
M=2;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_sh_1203=32*log(S2./32) + 2*(M+1);

%Nernst Model
Vhat= kp_NNb2(1)  + kp_NNb2(2)*log(z) + kp_NNb2(3)*log(1-z);
S2=0;
e=OCV_b2-Vhat;
M=3;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_nn_1203=32*log(S2./32) + 2*(M+1);

%Combine Model
Vhat= kp_cmb2(1) + kp_cmb2(2)*(1./z) +kp_cmb2(3)*z +...
    kp_cmb2(4)*log(z) + kp_cmb2(5)*log(1-z);
S2=0;
e=OCV_b2-Vhat;
M=5;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_cm_1203=32*log(S2./32) + 2*(M+1);

%Combine +3 Model
Vhat= kp_cm3b2(1)+ kp_cm3b2(2)*(1./z) + kp_cm3b2(3)*(1./(z.^2))+...
    kp_cm3b2(4)*(1./(z.^3)) + kp_cm3b2(5)*(1./(z.^4))+ kp_cm3b2(6)*(z)...
    + kp_cm3b2(7)*(log(z))+ kp_cm3b2(8)*(log(1-z));
S2=0;
e=OCV_b2-Vhat;
M=8;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_cm3_1203=32*log(S2./32) + 2*(M+1);

%Polynomial Model
Vhat= kp_pob2(1)+ kp_pob2(2)*(z) +kp_pob2(3)*(z.^2)...
    + kp_pob2(4)*(z.^3)+ kp_pob2(5)*(1./z) + kp_pob2(6)*(1./(z.^2));
S2=0;
e=OCV_b2-Vhat;
M=6;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_pol_1203=32*log(S2./32) + 2*(M+1);

%Exponential Model
Vhat= kp_exb2(1)+ kp_exb2(2)*exp(SOC_b2) +kp_exb2(3)*exp(SOC_b2.^2)...
    + kp_exb2(4)*exp(SOC_b2.^3)+ kp_exb2(5)*exp(-SOC_b2) + kp_exb2(6)*exp(-(SOC_b2.^2));
S2=0;
e=OCV_b2-Vhat;
M=6;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_ex_1203=32*log(S2./32) + 2*(M+1);
%=========================================================================================================

%Creating Metrics Table & Ranking Table - Battery C1203

AIC=[aic_u_1203;aic_sh_1203;aic_nn_1203;aic_cm_1203;aic_cm3_1203;aic_pol_1203;aic_ex_1203];
RMSE=[rmse_u_1203;rmse_s_1203;rmse_nn_1203;rmse_cm_1203;rmse_cm3_1203;rmse_pol_1203;rmse_ex_1203];
RSquare=[R2_U_1203;R2_S_1203;R2_nn_1203;R2_cm_1203;R2_cm3_1203;R2_pol_1203;R2_ex_1203];
BestFit=[BF_U_1203;BF_S_1203;BF_nn_1203;BF_cm_1203;BF_cm3_1203;BF_pol_1203;BF_ex_1203];
MaxError=[max_error_u_1203;max_error_s_1203;max_error_nn_1203;max_error_cm_1203;max_error_cm3_1203;max_error_pol_1203;max_error_ex_1203];
OCVModel=["Unnewhr";"Shepherd";"Nernst";"Combine";"Combine +3";"Polynomial";"Exponential"];
Metrics_table_1203 = table(OCVModel,AIC,RMSE,RSquare,BestFit,MaxError)
AIC=[1;2;3;4;5;6;7];
RMSE=[1;2;3;4;5;7;6];
RSquare=[1;2;3;4;5;6;7];
BestFit=[1;2;3;4;5;7;6];
MaxError=[1;2;3;4;6;5;7];
BordaRanking=[1;2;3;4;5;6;7];
OCVModel=["Combine +3";"Polynomial";"Combine";"Exponential";"Nernst";"Shepherd";"Unnewhr"];
ranking_table_1203 = table(OCVModel,AIC,RMSE,RSquare,BestFit,MaxError,BordaRanking)

%========================================================== For Battery C1204 ===================================================%

%Best Fit

%Unnewehr Universal Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b3(i);
end
Vmean=Vmean/32;
Vhat= kp_unb3(1) + kp_unb3(2)*(SOC_b3);
a=abs(Vhat-OCV_b3);
b=abs(OCV_b3-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_U_1204=1-(sum_a/sum_b);

%Shepherd Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b3(i);
end
Vmean=Vmean/32;
z=SOC_b3*(1-2*E) + E;
Vhat= kp_shb3(1) + kp_shb3(2)*(1./z);
a=abs(Vhat-OCV_b3);
b=abs(OCV_b3-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_S_1204=1-(sum_a/sum_b);

%Nernst Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b3(i);
end
Vmean=Vmean/32;
z=SOC_b3*(1-2*E) + E;
Vhat= kp_NNb3(1)  + kp_NNb3(2)*log(z) + kp_NNb3(3)*log(1-z);
a=abs(Vhat-OCV_b3);
b=abs(OCV_b3-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_nn_1204=1-(sum_a/sum_b);

%Combine Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b3(i);
end
Vmean=Vmean/32;
z=SOC_b3*(1-2*E) + E;
Vhat= kp_cmb3(1) + kp_cmb3(2)*(1./z) +kp_cmb3(3)*z +...
    kp_cmb3(4)*log(z) + kp_cmb3(5)*log(1-z);
a=abs(Vhat-OCV_b3);
b=abs(OCV_b3-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_cm_1204=1-(sum_a/sum_b);

%Combine +3 Model

Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b3(i);
end
Vmean=Vmean/32;
z=SOC_b3*(1-2*E) + E;
Vhat= kp_cm3b3(1)+ kp_cm3b3(2)*(1./z) + kp_cm3b3(3)*(1./(z.^2))+...
    kp_cm3b3(4)*(1./(z.^3)) + kp_cm3b3(5)*(1./(z.^4))+ kp_cm3b3(6)*(z)...
    + kp_cm3b3(7)*(log(z))+ kp_cm3b3(8)*(log(1-z));
a=abs(Vhat-OCV_b3);
b=abs(OCV_b3-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_cm3_1204=1-(sum_a/sum_b);

%Polynomial Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b3(i);
end
Vmean=Vmean/32;
z=SOC_b3*(1-2*E) + E;
Vhat= kp_pob3(1)+ kp_pob3(2)*(z) +kp_pob3(3)*(z.^2) +...
    kp_pob3(4)*(z.^3)+ kp_pob3(5)*(1./z) + kp_pob3(6)*(1./(z.^2));
a=abs(Vhat-OCV_b3);
b=abs(OCV_b3-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_pol_1204=1-(sum_a/sum_b);

%Exponential Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b3(i);
end
Vmean=Vmean/32;
z=SOC_b3*(1-2*E) + E;
Vhat= kp_exb3(1)+ kp_exb3(2)*exp(SOC_b3) +kp_exb3(3)*exp(SOC_b3.^2) ...
    + kp_exb3(4)*exp(SOC_b3.^3)+ kp_exb3(5)*exp(-SOC_b3) + kp_exb3(6)*exp(-(SOC_b3.^2));
a=abs(Vhat-OCV_b3);
b=abs(OCV_b3-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_ex_1204=1-(sum_a/sum_b);
%=========================================================================================================
%R-sqaure

%Unnewehr Universal Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b3(i);
end
Vmean=Vmean/32;
Vhat= kp_unb3(1) + kp_unb3(2)*(SOC_b3);
a=(abs(Vhat-OCV_b3)).^2;
b=(abs(OCV_b3-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_U_1204=1-(sum_a/sum_b);

%Shepherd Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b3(i);
end
Vmean=Vmean/32;
z=SOC_b3*(1-2*E) + E;
Vhat= kp_shb3(1) + kp_shb3(2)*(1./z);
a=(abs(Vhat-OCV_b3)).^2;
b=(abs(OCV_b3-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_S_1204=1-(sum_a/sum_b);

%Nernst Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b3(i);
end
Vmean=Vmean/32;
z=SOC_b3*(1-2*E) + E;
Vhat= kp_NNb3(1)  + kp_NNb3(2)*log(z) + kp_NNb3(3)*log(1-z);
a=(abs(Vhat-OCV_b3)).^2;
b=(abs(OCV_b3-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_nn_1204=1-(sum_a/sum_b);

%Combine Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b3(i);
end
Vmean=Vmean/32;
z=SOC_b3*(1-2*E) + E;
Vhat= kp_cmb3(1) + kp_cmb3(2)*(1./z) +kp_cmb3(3)*z +...
    kp_cmb3(4)*log(z) + kp_cmb3(5)*log(1-z);
a=(abs(Vhat-OCV_b3)).^2;
b=(abs(OCV_b3-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_cm_1204=1-(sum_a/sum_b);

%Combine +3 Model

Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b3(i);
end
Vmean=Vmean/32;
z=SOC_b3*(1-2*E) + E;
Vhat= kp_cm3b3(1)+ kp_cm3b3(2)*(1./z) + kp_cm3b3(3)*(1./(z.^2))...
    + kp_cm3b3(4)*(1./(z.^3)) + kp_cm3b3(5)*(1./(z.^4))+ kp_cm3b3(6)*(z)...
    + kp_cm3b3(7)*(log(z))+ kp_cm3b3(8)*(log(1-z));
a=(abs(Vhat-OCV_b3)).^2;
b=(abs(OCV_b3-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_cm3_1204=1-(sum_a/sum_b);

%Polynomial Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b3(i);
end
Vmean=Vmean/32;
z=SOC_b3*(1-2*E) + E;
Vhat= kp_pob3(1)+ kp_pob3(2)*(z) +kp_pob3(3)*(z.^2) ...
    + kp_pob3(4)*(z.^3)+ kp_pob3(5)*(1./z) + kp_pob3(6)*(1./(z.^2));
a=(abs(Vhat-OCV_b3)).^2;
b=(abs(OCV_b3-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_pol_1204=1-(sum_a/sum_b);

%Exponential Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b3(i);
end
Vmean=Vmean/32;
z=SOC_b3*(1-2*E) + E;
Vhat= kp_exb3(1)+ kp_exb3(2)*exp(SOC_b3) +kp_exb3(3)*exp(SOC_b3.^2)...
    + kp_exb3(4)*exp(SOC_b3.^3)+ kp_exb3(5)*exp(-SOC_b3) + kp_exb3(6)*exp(-(SOC_b3.^2));
a=(abs(Vhat-OCV_b3)).^2;
b=(abs(OCV_b3-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_ex_1204=1-(sum_a/sum_b);
%=========================================================================================================
%Max Error

z=SOC_b3*(1-2*E) + E;
%Unnewehr Universal Model
Vhat= kp_unb3(1) + kp_unb3(2)*(SOC_b3);
ME=abs(OCV_b3-Vhat);
max_error_u_1204=max(ME);
%Shepherd Model
Vhat= kp_shb3(1) + kp_shb3(2)*(1./z);
ME=abs(OCV_b3-Vhat);
max_error_s_1204=max(ME);
%Nernst Model
Vhat= kp_NNb3(1)  + kp_NNb3(2)*log(z) + kp_NNb3(3)*log(1-z);
ME=abs(OCV_b3-Vhat);
max_error_nn_1204=max(ME);
%Combine Model
Vhat= kp_cmb3(1) + kp_cmb3(2)*(1./z) +kp_cmb3(3)*z ...
    + kp_cmb3(4)*log(z) + kp_cmb3(5)*log(1-z);
ME=abs(OCV_b3-Vhat);
max_error_cm_1204=max(ME);
%Combine +3 Model
Vhat= kp_cm3b3(1)+ kp_cm3b3(2)*(1./z) + kp_cm3b3(3)*(1./(z.^2))...
    + kp_cm3b3(4)*(1./(z.^3)) + kp_cm3b3(5)*(1./(z.^4))+ kp_cm3b3(6)*(z)...
    + kp_cm3b3(7)*(log(z))+ kp_cm3b3(8)*(log(1-z));
ME=abs(OCV_b3-Vhat);
max_error_cm3_1204=max(ME);
%Polynomial Model
Vhat= kp_pob3(1)+ kp_pob3(2)*(z) +kp_pob3(3)*(z.^2) + kp_pob3(4)*(z.^3)...
    + kp_pob3(5)*(1./z) + kp_pob3(6)*(1./(z.^2));
ME=abs(OCV_b3-Vhat);
max_error_pol_1204=max(ME);
%Exponential Model
Vhat= kp_exb3(1)+ kp_exb3(2)*exp(SOC_b3) +kp_exb3(3)*exp(SOC_b3.^2)...
    + kp_exb3(4)*exp(SOC_b3.^3)+ kp_exb3(5)*exp(-SOC_b3) + kp_exb3(6)*exp(-(SOC_b3.^2));
ME=abs(OCV_b3-Vhat);
max_error_ex_1204=max(ME);

%=========================================================================================================
%Root Mean Square Error

N=32;
z=SOC_b3*(1-2*E) + E;
%Unnewehr Universal Model
sum_a=0;
Vhat= kp_unb3(1) + kp_unb3(2)*(SOC_b3);
a=abs(OCV_b3-Vhat);
M=2;
for i=1:32
    sum_a= sum_a+ a(i);  
end
rmse_u_1204=sum_a/sqrt(N-M);

%Shepherd Model
sum_a=0;
Vhat= kp_shb3(1) + kp_shb3(2)*(1./z);
a=abs(OCV_b3-Vhat);
M=2;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_s_1204=sum_a/sqrt(N-M);

%Nernst Model
sum_a=0;
Vhat= kp_NNb3(1)  + kp_NNb3(2)*log(z) + kp_NNb3(3)*log(1-z);
a=abs(OCV_b3-Vhat);
M=3;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_nn_1204=sum_a/sqrt(N-M);
%Combine Model
sum_a=0;
Vhat= kp_cmb3(1) + kp_cmb3(2)*(1./z) +kp_cmb3(3)*z ...
    + kp_cmb3(4)*log(z) + kp_cmb3(5)*log(1-z);
a=abs(OCV_b3-Vhat);
M=5;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_cm_1204=sum_a/sqrt(N-M);

%Combine +3 Model
sum_a=0;
Vhat= kp_cm3b3(1)+ kp_cm3b3(2)*(1./z) + kp_cm3b3(3)*(1./(z.^2))...
    + kp_cm3b3(4)*(1./(z.^3)) + kp_cm3b3(5)*(1./(z.^4))+ kp_cm3b3(6)*(z) ...
    + kp_cm3b3(7)*(log(z))+ kp_cm3b3(8)*(log(1-z));
a=abs(OCV_b3-Vhat);
M=8;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_cm3_1204=sum_a/sqrt(N-M);

%Polynomial Model
sum_a=0;
Vhat= kp_pob3(1)+ kp_pob3(2)*(z) +kp_pob3(3)*(z.^2)...
    + kp_pob3(4)*(z.^3)+ kp_pob3(5)*(1./z) + kp_pob3(6)*(1./(z.^2));
a=abs(OCV_b3-Vhat);
M=6;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_pol_1204=sum_a/sqrt(N-M);

%Exponential Model
sum_a=0;
Vhat= kp_exb3(1)+ kp_exb3(2)*exp(SOC_b3) +...
    kp_exb3(3)*exp(SOC_b3.^2) + kp_exb3(4)*exp(SOC_b3.^3)+ ...
    kp_exb3(5)*exp(-SOC_b3) + kp_exb3(6)*exp(-(SOC_b3.^2));
a=abs(OCV_b3-Vhat);
M=6;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_ex_1204=sum_a/sqrt(N-M);
%=========================================================================================================
%AIC

z=SOC_b3*(1-2*E) + E;
%Unnewehr Universal Model
Vhat= kp_unb3(1) + kp_unb3(2)*(SOC_b3);
S2=0;
e=OCV_b3-Vhat;
M=2;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_u_1204=32*log(S2./32) + 2*(M+1);

%Shepherd Model
Vhat= kp_shb3(1) + kp_shb3(2)*(1./z);
S2=0;
e=OCV_b3-Vhat;
M=2;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_sh_1204=32*log(S2./32) + 2*(M+1);

%Nernst Model
Vhat= kp_NNb3(1)  + kp_NNb3(2)*log(z) + kp_NNb3(3)*log(1-z);
S2=0;
e=OCV_b3-Vhat;
M=3;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_nn_1204=32*log(S2./32) + 2*(M+1);

%Combine Model
Vhat= kp_cmb3(1) + kp_cmb3(2)*(1./z) +kp_cmb3(3)*z ...
    + kp_cmb3(4)*log(z) + kp_cmb3(5)*log(1-z);
S2=0;
e=OCV_b3-Vhat;
M=5;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_cm_1204=32*log(S2./32) + 2*(M+1);

%Combine +3 Model
Vhat= kp_cm3b3(1)+ kp_cm3b3(2)*(1./z) + kp_cm3b3(3)*(1./(z.^2))...
    + kp_cm3b3(4)*(1./(z.^3)) + kp_cm3b3(5)*(1./(z.^4))+ kp_cm3b3(6)*(z)...
    + kp_cm3b3(7)*(log(z))+ kp_cm3b3(8)*(log(1-z));
S2=0;
e=OCV_b3-Vhat;
M=8;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_cm3_1204=32*log(S2./32) + 2*(M+1);

%Polynomial Model
Vhat= kp_pob3(1)+ kp_pob3(2)*(z) +kp_pob3(3)*(z.^2) + ...
    kp_pob3(4)*(z.^3)+ kp_pob3(5)*(1./z) + kp_pob3(6)*(1./(z.^2));
S2=0;
e=OCV_b3-Vhat;
M=6;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_pol_1204=32*log(S2./32) + 2*(M+1);

%Exponential Model
Vhat= kp_exb3(1)+ kp_exb3(2)*exp(SOC_b3) +kp_exb3(3)*exp(SOC_b3.^2)...
    + kp_exb3(4)*exp(SOC_b3.^3)+ kp_exb3(5)*exp(-SOC_b3) + kp_exb3(6)*exp(-(SOC_b3.^2));
S2=0;
e=OCV_b3-Vhat;
M=6;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_ex_1204=32*log(S2./32) + 2*(M+1);
%=========================================================================================================
%Creating Metrics Table & Ranking Table - Battery C1204

AIC=[aic_u_1204;aic_sh_1204;aic_nn_1204;aic_cm_1204;aic_cm3_1204;aic_pol_1204;aic_ex_1204];
RMSE=[rmse_u_1204;rmse_s_1204;rmse_nn_1204;rmse_cm_1204;rmse_cm3_1204;rmse_pol_1204;rmse_ex_1204];
RSquare=[R2_U_1204;R2_S_1204;R2_nn_1204;R2_cm_1204;R2_cm3_1204;R2_pol_1204;R2_ex_1204];
BestFit=[BF_U_1204;BF_S_1204;BF_nn_1204;BF_cm_1204;BF_cm3_1204;BF_pol_1204;BF_ex_1204];
MaxError=[max_error_u_1204;max_error_s_1204;max_error_nn_1204;max_error_cm_1204;max_error_cm3_1204;max_error_pol_1204;max_error_ex_1204];
OCVModel=["Unnewhr";"Shepherd";"Nernst";"Combine";"Combine +3";"Polynomial";"Exponential"];
Metrics_table_1204 = table(OCVModel,AIC,RMSE,RSquare,BestFit,MaxError)
AIC=[1;2;3;4;5;6;7];
RMSE=[1;2;3;4;5;7;6];
RSquare=[1;2;3;4;5;6;7];
BestFit=[1;2;3;4;5;7;6];
MaxError=[1;2;3;4;6;5;7];
BordaRanking=[1;2;3;4;5;6;7];
OCVModel=["Combine +3";"Polynomial";"Combine";"Exponential";"Nernst";"Shepherd";"Unnewhr"];
ranking_table_1204 = table(OCVModel,AIC,RMSE,RSquare,BestFit,MaxError,BordaRanking)

%========================================================== For Battery C1205 ===================================================%%

%Best Fit

%Unnewehr Universal Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b4(i);
end
Vmean=Vmean/32;
Vhat= kp_unb4(1) + kp_unb4(2)*(SOC_b4);
a=abs(Vhat-OCV_b4);
b=abs(OCV_b4-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BestFitLinearBatt4=1-(sum_a/sum_b);

%Shepherd Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b4(i);
end
Vmean=Vmean/32;
z=SOC_b4*(1-2*E) + E;
Vhat= kp_shb4(1) + kp_shb4(2)*(1./z);
a=abs(Vhat-OCV_b4);
b=abs(OCV_b4-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BestFitShepherdBatt4=1-(sum_a/sum_b);

%Nernst Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b4(i);
end
Vmean=Vmean/32;
z=SOC_b4*(1-2*E) + E;
Vhat= kp_NNb4(1)  + kp_NNb4(2)*log(z) + kp_NNb4(3)*log(1-z);
a=abs(Vhat-OCV_b4);
b=abs(OCV_b4-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_nn_1205=1-(sum_a/sum_b);

%Combine Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b4(i);
end
Vmean=Vmean/32;
z=SOC_b4*(1-2*E) + E;
Vhat= kp_cmb4(1) + kp_cmb4(2)*(1./z) +kp_cmb4(3)*z +...
    kp_cmb4(4)*log(z) + kp_cmb4(5)*log(1-z);
a=abs(Vhat-OCV_b4);
b=abs(OCV_b4-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_cm_1205=1-(sum_a/sum_b);

%Combine +3 Model

Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b4(i);
end
Vmean=Vmean/32;
z=SOC_b4*(1-2*E) + E;
Vhat= kp_cm3b4(1)+ kp_cm3b4(2)*(1./z) + kp_cm3b4(3)*(1./(z.^2))...
    + kp_cm3b4(4)*(1./(z.^3)) + kp_cm3b4(5)*(1./(z.^4))+ kp_cm3b4(6)*(z)...
    + kp_cm3b4(7)*(log(z))+ kp_cm3b4(8)*(log(1-z));
a=abs(Vhat-OCV_b4);
b=abs(OCV_b4-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_cm3_1205=1-(sum_a/sum_b);

%Polynomial Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b4(i);
end
Vmean=Vmean/32;
z=SOC_b4*(1-2*E) + E;
Vhat= kp_pob4(1)+ kp_pob4(2)*(z) +kp_pob4(3)*(z.^2) + ...
    kp_pob4(4)*(z.^3)+ kp_pob4(5)*(1./z) + kp_pob4(6)*(1./(z.^2));
a=abs(Vhat-OCV_b4);
b=abs(OCV_b4-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_pol_1205=1-(sum_a/sum_b);

%Exponential Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b4(i);
end
Vmean=Vmean/32;
z=SOC_b4*(1-2*E) + E;
Vhat= kp_exb4(1)+ kp_exb4(2)*exp(SOC_b4) +kp_exb4(3)*exp(SOC_b4.^2)...
    + kp_exb4(4)*exp(SOC_b4.^3)+ kp_exb4(5)*exp(-SOC_b4) + kp_exb4(6)*exp(-(SOC_b4.^2));
a=abs(Vhat-OCV_b4);
b=abs(OCV_b4-Vmean);
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
BF_ex_1205=1-(sum_a/sum_b);
%=========================================================================================================
%R-sqaure

%Unnewehr universal model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b4(i);
end
Vmean=Vmean/32;
Vhat= kp_unb4(1) + kp_unb4(2)*(SOC_b4);
a=(abs(Vhat-OCV_b4)).^2;
b=(abs(OCV_b4-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_U_1205=1-(sum_a/sum_b);

%Shepherd Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b4(i);
end
Vmean=Vmean/32;
z=SOC_b4*(1-2*E) + E;
Vhat= kp_shb4(1) + kp_shb4(2)*(1./z);
a=(abs(Vhat-OCV_b4)).^2;
b=(abs(OCV_b4-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_S_1205=1-(sum_a/sum_b);

%Nernst Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b4(i);
end
Vmean=Vmean/32;
z=SOC_b4*(1-2*E) + E;
Vhat= kp_NNb4(1)  + kp_NNb4(2)*log(z) + kp_NNb4(3)*log(1-z);
a=(abs(Vhat-OCV_b4)).^2;
b=(abs(OCV_b4-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_nn_1205=1-(sum_a/sum_b);

%Combine Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b4(i);
end
Vmean=Vmean/32;
z=SOC_b4*(1-2*E) + E;
Vhat= kp_cmb4(1) + kp_cmb4(2)*(1./z) +kp_cmb4(3)*z +...
    kp_cmb4(4)*log(z) + kp_cmb4(5)*log(1-z);
a=(abs(Vhat-OCV_b4)).^2;
b=(abs(OCV_b4-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_cm_1205=1-(sum_a/sum_b);

%Combine +3 Model

Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b4(i);
end
Vmean=Vmean/32;
z=SOC_b4*(1-2*E) + E;
Vhat= kp_cm3b4(1)+ kp_cm3b4(2)*(1./z) + kp_cm3b4(3)*(1./(z.^2))+ ...
    kp_cm3b4(4)*(1./(z.^3)) + kp_cm3b4(5)*(1./(z.^4))+ kp_cm3b4(6)*(z)...
    + kp_cm3b4(7)*(log(z))+ kp_cm3b4(8)*(log(1-z));
a=(abs(Vhat-OCV_b4)).^2;
b=(abs(OCV_b4-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_cm3_1205=1-(sum_a/sum_b);

%Polynomial Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b4(i);
end
Vmean=Vmean/32;
z=SOC_b4*(1-2*E) + E;
Vhat= kp_pob4(1)+ kp_pob4(2)*(z) +kp_pob4(3)*(z.^2) ...
    + kp_pob4(4)*(z.^3)+ kp_pob4(5)*(1./z) + kp_pob4(6)*(1./(z.^2));
a=(abs(Vhat-OCV_b4)).^2;
b=(abs(OCV_b4-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_pol_1205=1-(sum_a/sum_b);

%Exponential Model
Vmean = 0;
sum_a=0;
sum_b=0;
for i= 1:32
    Vmean = Vmean+ OCV_b4(i);
end
Vmean=Vmean/32;
z=SOC_b4*(1-2*E) + E;
Vhat= kp_exb4(1)+ kp_exb4(2)*exp(SOC_b4) +...
    kp_exb4(3)*exp(SOC_b4.^2) + kp_exb4(4)*exp(SOC_b4.^3)+...
    kp_exb4(5)*exp(-SOC_b4) + kp_exb4(6)*exp(-(SOC_b4.^2));
a=(abs(Vhat-OCV_b4)).^2;
b=(abs(OCV_b4-Vmean)).^2;
for i=1:32
    sum_a= sum_a+ a(i);
    sum_b= sum_b+ b(i);
end
R2_ex_1205=1-(sum_a/sum_b);
%=========================================================================================================
%Max Error

z=SOC_b4*(1-2*E) + E;
%Unnewehr Universal Model
Vhat= kp_unb4(1) + kp_unb4(2)*(SOC_b4);
ME=abs(OCV_b4-Vhat);
max_error_u_1205=max(ME);
%Shepherd Model
Vhat= kp_shb4(1) + kp_shb4(2)*(1./z);
ME=abs(OCV_b4-Vhat);
max_error_s_1205=max(ME);
%Nernst Model
Vhat= kp_NNb4(1)  + kp_NNb4(2)*log(z) + kp_NNb4(3)*log(1-z);
ME=abs(OCV_b4-Vhat);
max_error_nn_1205=max(ME);
%Combine Model
Vhat= kp_cmb4(1) + kp_cmb4(2)*(1./z) +kp_cmb4(3)*z ...
    + kp_cmb4(4)*log(z) + kp_cmb4(5)*log(1-z);
ME=abs(OCV_b4-Vhat);
max_error_cm_1205=max(ME);
%Combine +3 Model
Vhat= kp_cm3b4(1)+ kp_cm3b4(2)*(1./z) + ...
    kp_cm3b4(3)*(1./(z.^2))+ kp_cm3b4(4)*(1./(z.^3)) ...
    + kp_cm3b4(5)*(1./(z.^4))+ kp_cm3b4(6)*(z) + kp_cm3b4(7)*(log(z))+...
    kp_cm3b4(8)*(log(1-z));
ME=abs(OCV_b4-Vhat);
max_error_cm3_1205=max(ME);
%Polynomial Model
Vhat= kp_pob4(1)+ kp_pob4(2)*(z) +kp_pob4(3)*(z.^2)...
    + kp_pob4(4)*(z.^3)+ kp_pob4(5)*(1./z) + kp_pob4(6)*(1./(z.^2));
ME=abs(OCV_b4-Vhat);
max_error_pol_1205=max(ME);
%Exponential Model
Vhat= kp_exb4(1)+ kp_exb4(2)*exp(SOC_b4) +...
    kp_exb4(3)*exp(SOC_b4.^2) + kp_exb4(4)*exp(SOC_b4.^3)+...
    kp_exb4(5)*exp(-SOC_b4) + kp_exb4(6)*exp(-(SOC_b4.^2));
ME=abs(OCV_b4-Vhat);
max_error_ex_1205=max(ME);

%=========================================================================================================
%Root Mean Square Error

N=32;
z=SOC_b4*(1-2*E) + E;
%Unnewehr Universal Model
sum_a=0;
Vhat= kp_unb4(1) + kp_unb4(2)*(SOC_b4);
a=abs(OCV_b4-Vhat);
M=2;
for i=1:32
    sum_a= sum_a+ a(i);  
end
rmse_u_1205=sum_a/sqrt(N-M);

%Shepherd Model
sum_a=0;
Vhat= kp_shb4(1) + kp_shb4(2)*(1./z);
a=abs(OCV_b4-Vhat);
M=2;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_s_1205=sum_a/sqrt(N-M);

%Nernst Model
sum_a=0;
Vhat= kp_NNb4(1)  + kp_NNb4(2)*log(z) + kp_NNb4(3)*log(1-z);
a=abs(OCV_b4-Vhat);
M=3;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_nn_1205=sum_a/sqrt(N-M);
%Combine Model
sum_a=0;
Vhat= kp_cmb4(1) + kp_cmb4(2)*(1./z) +kp_cmb4(3)*z...
    + kp_cmb4(4)*log(z) + kp_cmb4(5)*log(1-z);
a=abs(OCV_b4-Vhat);
M=5;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_cm_1205=sum_a/sqrt(N-M);

%Combine +3 Model
sum_a=0;
Vhat= kp_cm3b4(1)+ kp_cm3b4(2)*(1./z) + kp_cm3b4(3)*(1./(z.^2))...
    + kp_cm3b4(4)*(1./(z.^3)) + kp_cm3b4(5)*(1./(z.^4))+ kp_cm3b4(6)*(z)...
    + kp_cm3b4(7)*(log(z))+ kp_cm3b4(8)*(log(1-z));
a=abs(OCV_b4-Vhat);
M=8;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_cm3_1205=sum_a/sqrt(N-M);

%Polynomial Model
sum_a=0;
Vhat= kp_pob4(1)+ kp_pob4(2)*(z) +kp_pob4(3)*(z.^2) +...
    kp_pob4(4)*(z.^3)+ kp_pob4(5)*(1./z) + kp_pob4(6)*(1./(z.^2));
a=abs(OCV_b4-Vhat);
M=6;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_pol_1205=sum_a/sqrt(N-M);

%Exponential Model
sum_a=0;
Vhat= kp_exb4(1)+ kp_exb4(2)*exp(SOC_b4) +...
    kp_exb4(3)*exp(SOC_b4.^2) + kp_exb4(4)*exp(SOC_b4.^3)+ ...
    kp_exb4(5)*exp(-SOC_b4) + kp_exb4(6)*exp(-(SOC_b4.^2));
a=abs(OCV_b4-Vhat);
M=6;
for i=1:32
    sum_a= sum_a+ a(i);
end
rmse_ex_1205=sum_a/sqrt(N-M);
%=========================================================================================================
%AIC

z=SOC_b4*(1-2*E) + E;
%Unnewehr Universal Model
Vhat= kp_unb4(1) + kp_unb4(2)*(SOC_b4);
S2=0;
e=OCV_b4-Vhat;
M=2;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_u_1205=32*log(S2./32) + 2*(M+1);

%Shepherd Model
Vhat= kp_shb4(1) + kp_shb4(2)*(1./z);
S2=0;
e=OCV_b4-Vhat;
M=2;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_sh_1205=32*log(S2./32) + 2*(M+1);

%Nernst Model
Vhat= kp_NNb4(1)  + kp_NNb4(2)*log(z) + kp_NNb4(3)*log(1-z);
S2=0;
e=OCV_b4-Vhat;
M=3;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_nn_1205=32*log(S2./32) + 2*(M+1);

%Combine Model
Vhat= kp_cmb4(1) + kp_cmb4(2)*(1./z) ...
    +kp_cmb4(3)*z + kp_cmb4(4)*log(z) + kp_cmb4(5)*log(1-z);
S2=0;
e=OCV_b4-Vhat;
M=5;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_cm_1205=32*log(S2./32) + 2*(M+1);

%Combine +3 Model
Vhat= kp_cm3b4(1)+ kp_cm3b4(2)*(1./z) +...
    kp_cm3b4(3)*(1./(z.^2))+ kp_cm3b4(4)*(1./(z.^3)) +...
    kp_cm3b4(5)*(1./(z.^4))+ kp_cm3b4(6)*(z) + kp_cm3b4(7)*(log(z))...
    + kp_cm3b4(8)*(log(1-z));
S2=0;
e=OCV_b4-Vhat;
M=8;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_cm3_1205=32*log(S2./32) + 2*(M+1);

%Polynomial Model
Vhat= kp_pob4(1)+ kp_pob4(2)*(z) +kp_pob4(3)*(z.^2)...
    + kp_pob4(4)*(z.^3)+ kp_pob4(5)*(1./z) + kp_pob4(6)*(1./(z.^2));
S2=0;
e=OCV_b4-Vhat;
M=6;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_pol_1205=32*log(S2./32) + 2*(M+1);

%Exponential Model
Vhat= kp_exb4(1)+ kp_exb4(2)*exp(SOC_b4) +...
    kp_exb4(3)*exp(SOC_b4.^2) + kp_exb4(4)*exp(SOC_b4.^3)+...
    kp_exb4(5)*exp(-SOC_b4) + kp_exb4(6)*exp(-(SOC_b4.^2));
S2=0;
e=OCV_b4-Vhat;
M=6;
for i=1:32
    S2 = S2 + (e(i).^2);
end
aic_ex_1205=32*log(S2./32) + 2*(M+1);
%=========================================================================================================

%Creating Metrics Table & Ranking Table - Battery C1205

AIC=[aic_u_1205;aic_sh_1205;aic_nn_1205;aic_cm_1205;aic_cm3_1205;aic_pol_1205;aic_ex_1205];
RMSE=[rmse_u_1205;rmse_s_1205;rmse_nn_1205;rmse_cm_1205;rmse_cm3_1205;rmse_pol_1205;rmse_ex_1205];
RSquare=[R2_U_1205;R2_S_1205;R2_nn_1205;R2_cm_1205;R2_cm3_1205;R2_pol_1205;R2_ex_1205];
BestFit=[BestFitLinearBatt4;BestFitShepherdBatt4;BF_nn_1205;BF_cm_1205;BF_cm3_1205;BF_pol_1205;BF_ex_1205];
MaxError=[max_error_u_1205;max_error_s_1205;max_error_nn_1205;max_error_cm_1205;max_error_cm3_1205;max_error_pol_1205;max_error_ex_1205];
OCVModel=["Unnewhr";"Shepherd";"Nernst";"Combine";"Combine +3";"Polynomial";"Exponential"];
Metrics_table_1205 = table(OCVModel,AIC,RMSE,RSquare,BestFit,MaxError)
AIC=[1;2;3;4;5;6;7];
RMSE=[1;2;3;4;5;7;6];
RSquare=[1;2;3;4;5;6;7];
BestFit=[1;2;3;4;5;7;6];
MaxError=[1;2;3;5;6;4;7];
BordaRanking=[1;2;3;4;5;6;7];
OCVModel=["Combine +3";"Polynomial";"Combine";"Exponential";"Nernst";"Shepherd";"Unnewhr"];
ranking_table_1205 = table(OCVModel,AIC,RMSE,RSquare,BestFit,MaxError,BordaRanking)
