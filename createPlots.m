function [final_data,Max_all,Locs_all,Max_Ns,Locs_Ns] =  createPlots(x)
  files = dir(x); %% change the date folder for processing all experiments



k=1;
for i = 3:length(files) 
Names{k} = files(i).name;

k = k+1;
end %loading the names of each experiment


k=1;
for i=3:length(Names)
m{k} = [files(i).folder '\'   files(i).name];

k = k+1;
end

%%


% 

%prealoca los vectores mean y std para cada angulo
intensities1= zeros(2048,5); %para el Prom y STD de todo el espectro
% peak_spectrum_Ns = zeros(1,length(files));
% peak_spectrum = zeros(1,length(files));

% figure(i)
% clf

for j = 1:length(m)
files_1exp = [m{j} '\' '*.txt'];
files = dir(files_1exp);

for i = 1:length(files)
    figure(j)
    
clf


a=importdata([files(i).folder '\' files(i).name]);
wavelengths = a(700:900,1);
intensities = a(700:900,2);



intensities2(:,i) = a(700:900,2);


subplot(2,2,1)



plot(wavelengths, intensities2, '-','linewidth',3);
hold on
xlim([600 700]);
%ylim([2.8E4 5E4])
xlabel('Wavelength (nm)');
ylabel('Intensity (counts)');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;


n=10       ;   % polinomion grade
[p,~,mu] = polyfit(wavelengths,intensities,n); %fiteo con centering and scaling to improve the numerical properties.
fitted = polyval(p, wavelengths,[],mu);   % Create polynom
plot(wavelengths, fitted, '-.');

% hold off

[Max,LOCS] = max(fitted);              % Find index of peak in fitted function
peak_spectrum(i) = wavelengths(LOCS);



Max_all{1,j} = char(Names(j));
Max_all{2,j} = Max;


Locs_all{1,j} = char(Names(j));
Locs_all{2,j} = wavelengths(LOCS);




%Normalizacion de los espectros

Ns = normalize (intensities,'range');



%Calculo de los picos sobre la normalizacion

n=15       ;   % polinomion grade
[pNs,~,muNs] = polyfit(wavelengths,Ns,n); %fiteo con centering and scaling to improve the numerical properties.
fittedNs = polyval(pNs, wavelengths,[],muNs);   % Create polynom
plot(wavelengths, fittedNs, '-.');
% hold off


[MaxNs,LOCSNs] = max(fittedNs);              % Find index of peak in fitted function
peak_spectrum_Ns(i) = wavelengths(LOCSNs);


Max_Ns{1,j} = char(Names(j));
Max_Ns{2,j} = MaxNs;


Locs_Ns{1,j} = char(Names(j));
Locs_Ns{2,j} = wavelengths(LOCS);





%Para prom y STD de todo el espectro
intensities1(:,i) = a(:,2);









end





for l = 1:length(files)

Ns2(:,l) = normalize (intensities2(:,l),'range');



end


subplot(2,2,3)
plot(wavelengths,Ns2,'linewidth',3)
hold on
xlim([500 700]);
ylim([0 1.1])
xlabel('Wavelength (nm)');
ylabel('Normalized Intensity');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;








spects = (1:length(files));
subplot(2,2,2)
plot(spects,peak_spectrum, 'ro')
% hold on
xlabel('Measurement N^{\circ}');
ylabel('Peak (nm)');
set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;
% hold off


subplot(2,2,4)
plot(spects,peak_spectrum_Ns, 'ro')
% hold on
xlabel('Measurement N^{\circ}');
ylabel('Normalized Peak (nm)');

set(gca,'XMinorTick','on','YMinorTick','off','fontsize',14);
grid on;



sgtitle(char(Names(j)),'FontSize',16)




%%%Guarda los promedios y desviacion de los espectros%%%%

final_data{1,j} = char(Names(j));
final_data{2,j} = mean(intensities1,2);
final_data{3,j} = std(intensities1,0,2);




end




end %Function's END