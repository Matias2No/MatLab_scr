%% Direct light

files = dir('C:\Matias\PosDoc\Publicaciones\Reproducciones\Measurements\SPR-setup\SpectroThorlabs-SLS201\03112021\Absorbancia\1)Light\*.txt');



%prealoca los vectores mean y std para cada angulo
intensities1= zeros(2048,5); %para el Prom y STD de todo el espectro
peak_spectrum_Ns = zeros(1,length(files));
peak_spectrum = zeros(1,length(files));
NsLight = zeros(2048,5);




figure(1)
clf


for i = 1:length(files)
a=importdata([files(i).folder '\' files(i).name]);
wavelengths = a(:,1);
absorbancia = a(:,2);
%m = 1:length(wavelengths);


subplot(2,2,1)

plot(wavelengths, absorbancia, '-','linewidth',3);
hold on
xlim([450 1000]);
xlabel('Wavelength (nm)');
ylabel('Absorbancia (OD)');
title('Direct Light');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;




% %%Calculo del pico de cada medicion
% n=10       ;   % polinomion grade
% [p,~,mu] = polyfit(wavelengths,absorbancia,n); %fiteo con centering and scaling to improve the numerical properties.
% fitted = polyval(p, wavelengths,[],mu);   % Create polynom
% plot(wavelengths, fitted, '-.');
% [Max,LOCS] = max(fitted);              % Find index of peak in fitted function
% peak_spectrum(i) = wavelengths(LOCS);



%Para prom y STD de todo el espectro
intensities1(:,i) = a(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Normalizacion de los espectros

Ns = normalize (absorbancia,'range');

%%%Normalizacion para hacer promedio
NsLight(:,i) = Ns(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(2,2,3)
plot(wavelengths,Ns,'linewidth',3)
hold on
xlim([450 1000]);
ylim([0 1.1])
xlabel('Wavelength (nm)');
ylabel('Normalized Intensity');
title('Normalized Light');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;



%Calculo de los picos sobre la normalizacion

% n=15       ;   % polinomion grade
% [pNs,~,muNs] = polyfit(wavelengths,Ns,n); %fiteo con centering and scaling to improve the numerical properties.
% fittedNs = polyval(pNs, wavelengths,[],muNs);   % Create polynom
% plot(wavelengths, fittedNs, '-.');
% [MaxNs,LOCSNs] = max(fittedNs);              % Find index of peak in fitted function
% peak_spectrum_Ns(i) = wavelengths(LOCSNs);






end



%%Promedio y STD de todo el espectro
% 


Light_mean = mean(intensities1,2);
Light_STD = std(intensities1,0,2);



subplot(2,2,2)

plot(wavelengths,Light_mean, '-','linewidth',3);
xlim([450 1000]);
xlabel('Wavelength (nm)');
ylabel('Absorbancia (OD)');
%title('Microchannel - 0^{\circ}');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;


title('Average Normalized Light')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%Promedio y STD de todo el espectro normalizado
% 


Light_meanNS = mean(NsLight,2);
Light_STDNS = std(NsLight,0,2);



subplot(2,2,4)

plot(wavelengths,Light_meanNS, '-','linewidth',3);
xlim([450 1000]);
xlabel('Wavelength (nm)');
ylabel('Absorbancia (OD)');
%title('Microchannel - 0^{\circ}');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;


title('Average of Direct Light Normalized absorbance')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%Grafica los picos calculados


% spects = (1:length(files));
% subplot(2,2,2)
% plot(spects,peak_spectrum, 'ro')
% xlabel('Measurement N^{\circ}');
% ylabel('Peak (nm)');
% %title('Microchannel - 0^{\circ}');
% 
% set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
% grid on;
% 
% 
% 
% subplot(2,2,4)
% plot(spects,peak_spectrum_Ns, 'ro')
% xlabel('Measurement N^{\circ}');
% ylabel('Normalized Peak (nm)');
% %title('Microchannel - 0^{\circ}');
% set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
% grid on;
% 
% 
% 
% sgtitle('Direct light','FontSize',16)




%% Cubeta

files = dir('C:\Matias\PosDoc\Publicaciones\Reproducciones\Measurements\SPR-setup\SpectroThorlabs-SLS201\03112021\Absorbancia\2)Cubeta\*.txt');



%prealoca los vectores mean y std para cada angulo
intensities1= zeros(2048,5); %para el Prom y STD de todo el espectro
peak_spectrum_Ns = zeros(1,length(files));
peak_spectrum = zeros(1,length(files));
NsCub =  zeros(2048,5);




figure(2)
clf


for i = 1:length(files)
a=importdata([files(i).folder '\' files(i).name]);
wavelengths = a(:,1);
absorbancia = a(:,2);
%m = 1:length(wavelengths);


subplot(2,2,1)

plot(wavelengths, absorbancia, '-','linewidth',3);
hold on
xlim([450 1000]);
xlabel('Wavelength (nm)');
ylabel('Absorbancia (OD)');
title('Cuvette');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;




% %%Calculo del pico de cada medicion
% n=10       ;   % polinomion grade
% [p,~,mu] = polyfit(wavelengths,absorbancia,n); %fiteo con centering and scaling to improve the numerical properties.
% fitted = polyval(p, wavelengths,[],mu);   % Create polynom
% plot(wavelengths, fitted, '-.');
% [Max,LOCS] = max(fitted);              % Find index of peak in fitted function
% peak_spectrum(i) = wavelengths(LOCS);



%Normalizacion de los espectros

Ns = normalize (absorbancia,'range');


%%%Normalizacion para hacer promedio
NsCub(:,i) = Ns(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3)
plot(wavelengths,Ns,'linewidth',3)
hold on
xlim([450 1000]);
ylim([0 1.1])
xlabel('Wavelength (nm)');
ylabel('Normalized Absorbance-Cuvette');
title('Normalized Cuvette Absorbance');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;



%Calculo de los picos sobre la normalizacion

% n=15       ;   % polinomion grade
% [pNs,~,muNs] = polyfit(wavelengths,Ns,n); %fiteo con centering and scaling to improve the numerical properties.
% fittedNs = polyval(pNs, wavelengths,[],muNs);   % Create polynom
% plot(wavelengths, fittedNs, '-.');
% [MaxNs,LOCSNs] = max(fittedNs);              % Find index of peak in fitted function
% peak_spectrum_Ns(i) = wavelengths(LOCSNs);


%Para prom y STD de todo el espectro
intensities1(:,i) = a(:,2);




end



%%Promedio y STD de todo el espectro
% 


Cub_mean = mean(intensities1,2);
Cub_STD = std(intensities1,0,2);



subplot(2,2,2)

plot(wavelengths,Cub_mean, '-','linewidth',3);
xlim([450 1000]);
xlabel('Wavelength (nm)');
ylabel('Absorbancia (OD)');
%title('Microchannel - 0^{\circ}');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;


title('Average of Cuvette absorbance')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%Promedio y STD de todo el espectro normalizado
% 


Cub_meanNS = mean(NsCub,2);
Cub_STDNS = std(NsCub,0,2);



subplot(2,2,4)

plot(wavelengths,Light_meanNS, '-','linewidth',3);
xlim([450 1000]);
xlabel('Wavelength (nm)');
ylabel('Absorbancia (OD)');
%title('Microchannel - 0^{\circ}');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;


title('Average of Cuvette Normalized absorbance')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% %%Grafica los picos calculados


% spects = (1:length(files));
% subplot(2,2,2)
% plot(spects,peak_spectrum, 'ro')
% xlabel('Measurement N^{\circ}');
% ylabel('Peak (nm)');
% %title('Microchannel - 0^{\circ}');
% 
% set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
% grid on;
% 
% 
% 
% subplot(2,2,4)
% plot(spects,peak_spectrum_Ns, 'ro')
% xlabel('Measurement N^{\circ}');
% ylabel('Normalized Peak (nm)');
% %title('Microchannel - 0^{\circ}');
% set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
% grid on;
% 
% 
% 
% sgtitle('Direct light','FontSize',16)




%% Agua




files = dir('C:\Matias\PosDoc\Publicaciones\Reproducciones\Measurements\SPR-setup\SpectroThorlabs-SLS201\03112021\Absorbancia\3)Agua\*.txt');



%prealoca los vectores mean y std para cada angulo
intensities1= zeros(2048,5); %para el Prom y STD de todo el espectro
peak_spectrum_Ns = zeros(1,length(files));
peak_spectrum = zeros(1,length(files));
NsWater =  zeros(2048,5);


figure(3)
clf


for i = 1:length(files)
a=importdata([files(i).folder '\' files(i).name]);
wavelengths = a(:,1);
absorbancia = a(:,2);
%m = 1:length(wavelengths);


subplot(2,2,1)

plot(wavelengths, absorbancia, '-','linewidth',3);
hold on
xlim([450 1000]);
xlabel('Wavelength (nm)');
ylabel('Absorbancia (OD)');
title('Water');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;




% %%Calculo del pico de cada medicion
% n=10       ;   % polinomion grade
% [p,~,mu] = polyfit(wavelengths,absorbancia,n); %fiteo con centering and scaling to improve the numerical properties.
% fitted = polyval(p, wavelengths,[],mu);   % Create polynom
% plot(wavelengths, fitted, '-.');
% [Max,LOCS] = max(fitted);              % Find index of peak in fitted function
% peak_spectrum(i) = wavelengths(LOCS);



%Normalizacion de los espectros

Ns = normalize (absorbancia,'range');

%%%Normalizacion para hacer promedio
NsWater(:,i) = Ns(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(2,2,3)
plot(wavelengths,Ns,'linewidth',3)
hold on
xlim([450 1000]);
ylim([0 1.1])
xlabel('Wavelength (nm)');
ylabel('Normalized Absorbance');
title('Normalized Water Absorbance');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;







%Calculo de los picos sobre la normalizacion

% n=15       ;   % polinomion grade
% [pNs,~,muNs] = polyfit(wavelengths,Ns,n); %fiteo con centering and scaling to improve the numerical properties.
% fittedNs = polyval(pNs, wavelengths,[],muNs);   % Create polynom
% plot(wavelengths, fittedNs, '-.');
% [MaxNs,LOCSNs] = max(fittedNs);              % Find index of peak in fitted function
% peak_spectrum_Ns(i) = wavelengths(LOCSNs);


%Para prom y STD de todo el espectro
intensities1(:,i) = a(:,2);




end



%%Promedio y STD de todo el espectro
% 


Water_mean = mean(intensities1,2);
Water_STD = std(intensities1,0,2);



subplot(2,2,2)

plot(wavelengths,Water_mean, '-','linewidth',3);
hold on
xlim([450 1000]);
xlabel('Wavelength (nm)');
ylabel('Absorbancia (OD)');
%title('Microchannel - 0^{\circ}');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;


title('Average of Water absorbance')



Wat_meanNS = mean(NsWater,2);
Wat_STDNS = std(NsWater,0,2);



subplot(2,2,4)

plot(wavelengths,Wat_meanNS, '-','linewidth',3);
xlim([450 1000]);
xlabel('Wavelength (nm)');
ylabel('Normalized Absorbance');
%title('Microchannel - 0^{\circ}');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;


title('Average of Water Normalized absorbance')






% %%Grafica los picos calculados


% spects = (1:length(files));
% subplot(2,2,2)
% plot(spects,peak_spectrum, 'ro')
% xlabel('Measurement N^{\circ}');
% ylabel('Peak (nm)');
% %title('Microchannel - 0^{\circ}');
% 
% set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
% grid on;
% 
% 
% 
% subplot(2,2,4)
% plot(spects,peak_spectrum_Ns, 'ro')
% xlabel('Measurement N^{\circ}');
% ylabel('Normalized Peak (nm)');
% %title('Microchannel - 0^{\circ}');
% set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
% grid on;
% 
% 
% 
% sgtitle('Direct light','FontSize',16)




%% AuNPs a concentracion original

files = dir('C:\Matias\PosDoc\Publicaciones\Reproducciones\Measurements\SPR-setup\SpectroThorlabs-SLS201\03112021\Absorbancia\4)AuNPs\*.txt');



%prealoca los vectores mean y std para cada angulo
intensities1= zeros(2048,5); %para el Prom y STD de todo el espectro
peak_spectrum_Ns = zeros(1,length(files));
peak_spectrum = zeros(1,length(files));
NsAuNPs = zeros(2048,5);


figure(4)
clf


for i = 1:length(files)
a=importdata([files(i).folder '\' files(i).name]);
wavelengths = a(:,1);
absorbancia = a(:,2);
%m = 1:length(wavelengths);


subplot(2,2,1)

plot(wavelengths, absorbancia, '-','linewidth',3);
hold on
xlim([450 1000]);
xlabel('Wavelength (nm)');
ylabel('Absorbancia (OD)');
title('AuNPS');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;




% %%Calculo del pico de cada medicion
% n=10       ;   % polinomion grade
% [p,~,mu] = polyfit(wavelengths,absorbancia,n); %fiteo con centering and scaling to improve the numerical properties.
% fitted = polyval(p, wavelengths,[],mu);   % Create polynom
% plot(wavelengths, fitted, '-.');
% [Max,LOCS] = max(fitted);              % Find index of peak in fitted function
% peak_spectrum(i) = wavelengths(LOCS);



%Normalizacion de los espectros

Ns = normalize (absorbancia,'range');

%%%Normalizacion para hacer promedio
NsAuNPs(:,i) = Ns(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(2,2,3)
plot(wavelengths,Ns,'linewidth',3)
hold on
xlim([450 1000]);
ylim([0 1.1])
xlabel('Wavelength (nm)');
ylabel('Normalized Absorbance');
title('Normalized AuNPs Absorbance');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;




%Calculo de los picos sobre la normalizacion

% n=15       ;   % polinomion grade
% [pNs,~,muNs] = polyfit(wavelengths,Ns,n); %fiteo con centering and scaling to improve the numerical properties.
% fittedNs = polyval(pNs, wavelengths,[],muNs);   % Create polynom
% plot(wavelengths, fittedNs, '-.');
% [MaxNs,LOCSNs] = max(fittedNs);              % Find index of peak in fitted function
% peak_spectrum_Ns(i) = wavelengths(LOCSNs);


%Para prom y STD de todo el espectro
intensities1(:,i) = a(:,2);




end



%%Promedio y STD de todo el espectro
% 


AuNPs_mean = mean(intensities1,2);
AuNPs_STD = std(intensities1,0,2);



subplot(2,2,2)

plot(wavelengths,AuNPs_mean, '-','linewidth',3);
hold on
xlim([450 1000]);
xlabel('Wavelength (nm)');
ylabel('Absorbancia (OD)');
%title('Microchannel - 0^{\circ}');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;


title('Average of AuNPs absorbance')



AuNPs_meanNS = mean(NsAuNPs,2);
AuNPs_STDNS = std(NsAuNPs,0,2);



subplot(2,2,4)

plot(wavelengths,AuNPs_meanNS, '-','linewidth',3);
xlim([450 1000]);
xlabel('Wavelength (nm)');
ylabel('Normalized Absorbance');
%title('Microchannel - 0^{\circ}');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;


title('Average of AuNPs Normalized absorbance')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% %%Grafica los picos calculados


% spects = (1:length(files));
% subplot(2,2,2)
% plot(spects,peak_spectrum, 'ro')
% xlabel('Measurement N^{\circ}');
% ylabel('Peak (nm)');
% %title('Microchannel - 0^{\circ}');
% 
% set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
% grid on;
% 
% 
% 
% subplot(2,2,4)
% plot(spects,peak_spectrum_Ns, 'ro')
% xlabel('Measurement N^{\circ}');
% ylabel('Normalized Peak (nm)');
% %title('Microchannel - 0^{\circ}');
% set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
% grid on;
% 
% 
% 
% sgtitle('Direct light','FontSize',16)

