close all;
clc;

% Pomocna skripta za odredjivanje broja bita kvantizatora kojim se
% kvantizuju odbirci RRC filtra sistema kako bi zadovoljio potrebne
% specifikacije projekta
%  - tako da slabljenje u nepropusnom opsegu bude bolje od 30dB

% Parametri rrc filtra
rrc_alpha = 0.22;                   
rrc_M = 8;                          
rrc_k = 10;                         
hrrc = rrc(rrc_k,rrc_M,rrc_alpha);  % izracunavanje imp. odziva RRC filtra

% Parametri kvantizatora
adc_N = 4;
adc_Vfs = 1;
hrrc_k = adc(hrrc,adc_N,adc_Vfs);   % kvantizovanje odbiraka RRC filtra

figure(1);                          % crtanje impulsnog odziva
% plot(1:1:length(hrrc),hrrc);
% hold on;
plot(1:1:length(hrrc_k),hrrc_k);    % kvantizovan imp. odziv
title('Impulsni odziv kvantizovanog RRC filtra');
xlim([0,81]);
ylabel('h(t)');
% legend('Pre kvantizacije','Posle kvantizacije');

[HRRC, w] = freqz(hrrc,1,2048);     % frekv. odziv RRC filtra
[HRRC_k] = freqz(hrrc_k,1,2048);    % frekv. odziv RRC filtra

figure(2);                          % crtanje amplitudske karakteristike
% plot(w./(2*pi) , 20*log10(abs(HRRC)./max(abs(HRRC))));
% hold on;                            % kvantizovan
plot(w./(2*pi) , 20*log10(abs(HRRC_k)./max(abs(HRRC_k))), ...
    'LineWidth',1);
% plot([0,0.5],[-30,-30],'r --','LineWidth',1);
title('Amplitudska karakteristika kvantizovanog RRC filtra')
xlabel('F');
ylabel('20log(|H|/Hmax) [dB] ');
% legend('Pre kvantizacije','Posle kvantizacije');

