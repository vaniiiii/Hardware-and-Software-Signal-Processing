close all;
clc;
% ************************************************************************
% [1]
% Odrediti red RRC filtra N = kM za vrednost parametra alpha = 0.22 
% tako da je amplituda najveceg bocnog luka u nepropusnom opsegu 
% F > (1 + alpha) / (2M), gde je faktor interpolacije M = 8, 
% manja od 30 dB i EVM < 1 %.

                                % Inicijalizacija parametara RRC filtra
alpha = 0.22;                   % uslov projekta
M = 8;                          % uslov projekta
k = 10;                         % duzina impulsnog odziva u broju simbola

hrrc = rrc(k,M,alpha);          % izracunavanje imp. odziva RRC filtra

figure(1);                      % crtanje impulsnog odziva
plot(1:1:length(hrrc),hrrc);
title('Impulsni odziv RRC filtra');
xlim([0,81]);
%xlabel('')
ylabel('h(t)');

[HRRC, w] = freqz(hrrc,1,2048);% izracunavanje frekv. odziva RRC filtra

figure(2);                      % crtanje amplitudske karakteristike
plot(w./(2*pi) , 20*log10(abs(HRRC)./max(abs(HRRC))));
hold on;                        % crtanje granica po uslova projekta
plot([(1+alpha)/(2*M),(1+alpha)/(2*M)],[-60,0],'g');
text((1+alpha)/(2*M),-20,'\leftarrow (1+\alpha)/(2*M)')
plot([0,0.5],[-30,-30],'r','LineWidth',1);
title('Amplitudska karakteristika RRC filtra')
xlabel('F');
ylabel('20log(|H|/Hmax) [dB] ');


% ************************************************************************
% [2] 
% Odrediti odnos vrsne i srednje snage za 4-QAM i 16-QAM koji su uobliceni
% projektovanim RRC filtrom na uzorku od 1024 slucajnih simbola.

qam4_map = [-1-1j, -1+1j, 1-1j, 1+1j];          % 4-QAM mapa
syms_qam4_tx = randi([0 3],1,1024);             % generisanje sluc. simbola
syms_qam4_map_tx = qam4_map(syms_qam4_tx + 1);  % mapiranje simbola
syms_qam4_up = expander(syms_qam4_map_tx,M);    % povecanje uc. odabiranja
qam4_tx_out = filter(hrrc,1,syms_qam4_up);      % filtriranje simbola

rrc_delay = k*M/2;                  % kasnjenje RRC filtra
nplot = 1:1:512;                   % broj odbiraka za crtanje
nplot_len = length(nplot);
figure(3);                          
subplot(2,1,1);                     % crtanje I komponente signala QAM4 TX
plot(nplot,real(syms_qam4_up(1:1:nplot_len)),'o');
title ('I komponenta QAM4 predajnika')
ylabel('Amplituda');
xlim([0,nplot_len]);
ylim([-1.2,1.2]);
hold on;                            
plot(real(qam4_tx_out(rrc_delay:nplot_len + rrc_delay)),'r');

subplot(2,1,2);                     % crtanje Q komponente signala QAM4 TX
plot(nplot,imag(syms_qam4_up(1:nplot_len)),'o');
title ('Q komponenta QAM4 predajnika');
ylabel('Amplituda');
xlim([0,nplot_len]);
ylim([-1.2,1.2]);
hold on;                            
plot(imag(qam4_tx_out(rrc_delay:nplot_len + rrc_delay)),'r');

% Racunanje "Peak to Average Ratio" - PAR
p_avg_qam4 = mean(abs(qam4_tx_out).^2);         % Srednja snaga
p_peak_qam4 = max(abs(qam4_tx_out).^2);         % Vrsna snaga
p_par_qam4 = 10*log10(p_peak_qam4/p_avg_qam4);  % PAR
fprintf('\nPAR QAM4 = %.4f dB\n\n', p_par_qam4);      % Ispisivanje rez.


% Crtanje spektra izlaznog signala QAM4 predajnika
figure(4);

nspec_len = length(qam4_tx_out);        % 
nspec = ((1:nspec_len) - nspec_len./2) / nspec_len;

win = transpose(hanning(nspec_len));    % prozorska funkcija
plot(nspec,20*log10(fftshift(abs(fft(qam4_tx_out.*win)))));
title('Spektar signala na izlazu QAM4 predajnika')
ylabel('|H(F)| [dB]');
xlabel('F');


% "Povezivanje" predajnika i prijemnika
% - dodatno filtriranje signala istim RRC filtrom
% - dodaje se dvostruko kasnjenje pri crtanju zbog dvostrukog filtriranja

qam4_rx_out = filter(hrrc,1,qam4_tx_out); %signal na izlazu prijemnika

figure(5);                          
subplot(2,1,1);                     % crtanje I komponente signala QAM4 RX
plot(nplot,real(syms_qam4_up(1:1:nplot_len)),'o');
title ('I komponenta QAM4 prijemnika')
ylabel('Amplituda');
xlim([0,nplot_len]);
ylim([-2,2]);
hold on;                            
plot(real(qam4_rx_out(2*rrc_delay:nplot_len + 2*rrc_delay)),'r');

subplot(2,1,2);                     % crtanje Q komponente signala QAM4 RX
plot(nplot,imag(syms_qam4_up(1:nplot_len)),'o');
title ('Q komponenta QAM4 prijemnika');
ylabel('Amplituda');
xlim([0,nplot_len]);
ylim([-2,2]);
hold on;                            
plot(imag(qam4_rx_out(2*rrc_delay:nplot_len + 2*rrc_delay)),'r');


% Crtanje spektra izlaznog signala QAM4 prijemnika
figure(6);

plot(nspec,20*log10(fftshift(abs(fft(qam4_rx_out.*win)))));
title('Spektar signala na izlazu QAM4 prijemnika')
ylabel('|H(F)| [dB]');
xlabel('F');

% ************************************************************************
% [3]
% Nacrtati konstelaciju signala i izracunati EVM~1/sqrt(SNR) za 4-QAM koji 
% su uobliceni projektovanim RRC filtrom na uzorku od 1024 sluc. simbola.

ofs = 1;
% izdvajanje I i Q komponenti iz signala na izlazu prijemnika
qam4_rx_out_I = real(qam4_rx_out(2*rrc_delay+ofs:M:length(qam4_rx_out)));
qam4_rx_out_Q = imag(qam4_rx_out(2*rrc_delay+ofs:M:length(qam4_rx_out)));

figure(7);                              % crtanje konstelacije signala
plot(qam4_rx_out_I,qam4_rx_out_Q,'o');
hold on;
plot([-1.33, 1.33], [0, 0],'r--');      % crtanje granica simbola
plot([0, 0],[-1.33,1.33],'r--');
title('Konstelacija signala QAM4 prijemnika bez suma');
xlabel('I');
ylabel('Q');

% gornja granica za niz simbola - "iseceni" pri filtriranju
cut = length(syms_qam4_map_tx) - 2*rrc_delay/8;

% racunanje EVM-a za QAM4
qam4_evm = sqrt(sum((qam4_rx_out_I - real(syms_qam4_map_tx(1:cut))).^2 +... 
    (qam4_rx_out_Q - imag(syms_qam4_map_tx(1:cut))).^2)  / ... 
    sum(abs(qam4_rx_out_I + 1j*qam4_rx_out_Q ).^2)) * 100;

%ispisivanje rezultata u konzoli
fprintf('EVM QAM4 (bez suma) = %.4f \n', qam4_evm);


% ************************************************************************
% [4]
% Ponoviti prethodnu tacku u prisustvu aditivnog belog Gausovog šuma (AWGN) 
% koji se dodaje izmedju predajnika i prijemnika, odnosno dva RRC filtra, 
% za SNR = {30, 20, 10} dB.

SNR_s = [30,20,10];

% alokacija nizova
qam4_tx_out_noise = zeros(length(SNR_s), length(qam4_tx_out));
qam4_rx_out_noise = zeros(length(SNR_s), length(qam4_tx_out));
qam4_rx_out_I_noise = zeros(length(SNR_s), ...
                (length(qam4_tx_out)- 2*rrc_delay)/8);
qam4_rx_out_Q_noise = zeros(length(SNR_s), ...
                (length(qam4_tx_out)- 2*rrc_delay)/8);
qam4_evm_noise = zeros(length(SNR_s));
for i = 1:length(SNR_s);
    % dodavanje suma sa i-tim SNRom
    qam4_tx_out_noise(i,:) = awgn(qam4_tx_out, SNR_s(i));
    
    % filtriranje signala sa dodatim sumom
    qam4_rx_out_noise(i,:) = filter(hrrc,1,qam4_tx_out_noise(i,:));
    
    %izdvajanje komponenti I i Q iz "primljenog" signala
    qam4_rx_out_I_noise(i,:) = real(qam4_rx_out_noise(i,2*rrc_delay+ofs:... 
        M:length(qam4_rx_out_noise(i,:))));
    qam4_rx_out_Q_noise(i,:) = imag(qam4_rx_out_noise(i,2*rrc_delay+ofs:... 
        M:length(qam4_rx_out_noise(i,:))));

    figure;                             % crtanje konstelacije
    plot(qam4_rx_out_I_noise(i,:), qam4_rx_out_Q_noise(i,:), 'o');
    title(sprintf('Konstelacija signala QAM4 prijemnika, SNR = %d',... 
        SNR_s(i)));
    hold on
    plot([-1.33, 1.33], [0, 0],'r--');  % crtanje granica simbola
    plot([0, 0],[-1.33,1.33],'r--');
    xlabel('I');
    ylabel('Q');
    
    % Izracunavanje EVM
    qam4_evm_noise(i) = sqrt(sum((qam4_rx_out_I_noise(i,:) - ... 
        real(syms_qam4_map_tx(1:cut))).^2 + (qam4_rx_out_Q_noise(i,:)...
        - imag(syms_qam4_map_tx(1:cut))).^2) / ... 
        sum(abs(qam4_rx_out_I_noise(i,:) + 1j * ... 
        qam4_rx_out_Q_noise(i,:)).^2)) * 100;
    
    %ispisivanje rezultata u konzoli
    fprintf('EVM QAM4 (SNR: %d dB) = %.4f \n', SNR_s(i),qam4_evm_noise(i));
end



% ************************************************************************
% [5]
% Odrediti minimalnu rezoluciju A/D konvertora za koju se obezbedjuje 
% minimalni SNR od 15 dB za 4-QAM, odnosno 18 dB za 16-QAM, 
% ukoliko je u susednom kanalu prisutan bloker istog tipa modulacije 
% i snage 30 dB vece od željenog signala. 
% Pretpostaviti da je maksimalna trenutna vrednost signala jednaka 
% naponu pune skale A/D konvertora, ucestanost odabiranja 
% fs = 8 * 3.84 MHz. Izracunati ENBW iz odbiraka RRC filtra.

qam4_snr = 15;                              % minimalni SNR za QAM4
bloker = 30;                                % snaga blokera (iznad zeljene)
% Racunanje snage kvantizacionog suma i konverzija iz dB
p_n = 10^((0 - p_par_qam4 - bloker - qam4_snr)/10);
% Ekvivalenti propusni opseg suma
enbw = sum(abs(hrrc).^2) / abs(sum(hrrc)).^2;
step = sqrt(12 * p_n / enbw);               % korak kvantizacije ADCa
adc_vfs = 1;                                % napon pune skale ADCa
adc_n = ceil (log2(adc_vfs / step));        % broj bita
fprintf('\nADC_N QAM4 = %d bita\n', adc_n); % Ispisivanje rez.


