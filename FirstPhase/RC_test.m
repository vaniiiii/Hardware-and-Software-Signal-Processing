close all;
clc

alpha = 0.22;                   % uslov projekta
M = 4;                          % uslov projekta
k = 10;                         % duzina impulsnog odziva u broju simbola

hrrc = rrc(k,M,alpha);          % izracunavanje imp. odziva RRC filtra
hrc=conv(hrrc,hrrc);  

[HRC, w] = freqz(hrc,1,2048);% izracunavanje frekv. odziva RRC filtra

figure(2);                      % crtanje amplitudske karakteristike
plot(w./(2*pi) , 20*log10(abs(HRC)./max(abs(HRC))));
hold on;                        % crtanje granica po uslova projekta
title('Amplitudska karakteristika RC filtra')
xlabel('F');
ylabel('20log(|H|/Hmax) [dB] ');