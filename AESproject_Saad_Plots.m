%importing Battery data to matlab
t = readtable('Data.csv');

%------------------------------------------------------------------------
% ------------Q#1 - Plot the OCV-SOC  data-------------------------------
%For Battery C1202
figure(1)
subplot(2,2,1)
SOC_b1 = t{:,1}; OCV_b1 = t{:,2};
plot(SOC_b1,OCV_b1);
title("Battery C1202");
xlabel('State of Charge (SOC)');
ylabel('Voltage (OCV)');
%-------------------------------
%For Battery C1203
figure(1)
subplot(2,2,2)
SOC_b2 = t{:,3}; OCV_b2 = t{:,4};
plot(SOC_b2,OCV_b2);
title("Battery C1203");
xlabel('State of Charge (SOC)');
ylabel('Voltage (OCV)');
%-------------------------------
%For Battery C1204
figure(1)
subplot(2,2,3)
SOC_b3 = t{:,5}; OCV_b3 = t{:,6};
plot(SOC_b3,OCV_b3);
title("Battery C1204");
xlabel('State of Charge (SOC)');
ylabel('Voltage (OCV)');
%-------------------------------
%For Battery C1205
figure(1)
subplot(2,2,4)
SOC_b4 = t{:,7}; OCV_b4 = t{:,8};
plot(SOC_b4,OCV_b4);
title("Battery C1205");
xlabel('State of Charge (SOC)');
ylabel('Voltage (OCV)');
%------------------------------------------------------------------------
