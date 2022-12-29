close all;
clear;
clc;
% ************************************************************************
% [1]
% Implementirati predajnik, simulirati rad predajnika na du≈æini od 1024 
% simbola i nacrtati spektar izlaznog signala.

qam16_map = [-3-3j, -3-1j, -3+3j, -3+1j, -1-3j, -1-1j, -1+3j, -1+1j, ...
    3-3j, 3-1j, 3+3j, 3+1j, 1-3j, 1-1j, 1+3j, 1+1j,]/3;
qam16_syms_tx = randi([0 15],1,1024);               % generisanje simbola
qam16_syms_map_tx = qam16_map(qam16_syms_tx + 1);   % mapiranje simbola

qam16_tx_I = real(qam16_syms_map_tx);           % Izdvajanje I komponente
qam16_tx_Q = imag(qam16_syms_map_tx);           % Izdvajanje Q komponente

% Inicijalizacija parametara RRC filtra
rrc_alpha = 0.22;
rrc_M = 8;
rrc_k = 10;
rrc_delay = rrc_k * rrc_M / 2;

hrrc = rrc(rrc_k,rrc_M,rrc_alpha);  % izracunavanje imp. odziva RRC filtra
hrrc_poly = make_poly(hrrc,rrc_M);  % polifazno dekomponovan filtar

% Filtriranje polifaznom dekompozicijom filtra
qam16_tx_I_filt = filter_poly(qam16_tx_I,hrrc_poly,rrc_M);
qam16_tx_Q_filt = filter_poly(qam16_tx_Q,hrrc_poly,rrc_M);
% Formiranje signala na izlazu predajnika
qam16_tx = qam16_tx_I_filt+ 1j * qam16_tx_Q_filt; 

win=transpose(hann(length(qam16_tx)));      % prozorska funkcija

nspec_len = length(qam16_tx);
nspec = ((1:nspec_len) - nspec_len./2) / nspec_len;
% Iscrtavanje spektra signala na izlazu predajnika
plot(nspec, 20*log10(fftshift(abs(fft(qam16_tx.*win)))));
title('Spektar signala na izlazu 16-QAM predajnika ');
ylabel('|16-QAM-TX(F)|');
xlabel('F');

% ************************************************************************
% [2]
% Napraviti model A/D konvertora sa brojem bita odredjenim u prvoj fazi

adc_Vfs = 2;            % napon pune skale A/D konvertora
adc_N = 7;          

qam16_rx_I = adc(qam16_tx_I_filt, adc_N, adc_Vfs); % kvantizacija I komp.
qam16_rx_Q = adc(qam16_tx_Q_filt, adc_N, adc_Vfs); % kvantizacija Q komp.

qam16_rx_kvant = qam16_rx_I + 1j * qam16_rx_Q;
hold on;
plot(nspec, 20*log10(fftshift(abs(fft(qam16_rx_kvant.*win)))));
legend('Pre kvantizacije', 'Posle kvantizacije');

% ************************************************************************
% [3]
% Implementirati numericki kontrolisani oscilator i kompleksni mesac 
% sirine reci od 14 bita koriscenjem aritmetike u fiksnom zarezu 
tic;
% Definisanje broja bita za aritmetiku sa fiksnim zarezom
total_length = 19;
frac_length = 17;
% 14 bita po uslovu projektnog zadatka, +3 bita
% da greöka usled konacne duûine reci bude manja od LSB
% dodavanjem zaötitnih bitova (nguard < log2n)
num_of_bits = struct('SIGNAL_BITLENGTH', total_length,...
    'SIGNAL_FRACLENGTH', frac_length);

% Racunanje ukupnog multiplikativnog faktora za zadati broj mikrorotacija
K = 1;
for i = 0:(num_of_bits.SIGNAL_FRACLENGTH-1);
        K = K*(1/sqrt(1 + 2^(-2*(i))));
end

% Skaliranje komponenti izracunatim faktorom
qam16_rx_I_k = K * qam16_rx_I;
qam16_rx_Q_k = K * qam16_rx_Q;

% Frekvencijska rezolucija NCO-a i frekvencija signala
fnco = 10.6667;  
fs = 3.84 * 10^6;  
% NCO kontrolna rec
W = 2 ^ num_of_bits.SIGNAL_FRACLENGTH * (fnco/fs); 
% Pocetne vrednosti za fazni akumulator
A = 0;
B = W;
% Nizovi u koje se upisuju izlazne vrednosti
qam16_rx_I_c =  zeros(1,length(qam16_rx_I_k));
qam16_rx_Q_c =  zeros(1,length(qam16_rx_Q_k));

% Broj iteracija jednak je broju simbola i za svaku se izröava CORDIC
% algoritam(f-ja CORDIC)
for i = 1:length(qam16_rx_I_k)
    [qam16_rx_I_c(i), qam16_rx_Q_c(i)] = ...
        CORDIC(qam16_rx_I_k(i), qam16_rx_Q_k(i), A+B, num_of_bits); 
    A = A+B;          
end
toc;
% ************************************************************************
% [4]
% Implementirati filtar za uoblicavanje signala i 
% programabilni pomerac faze.

% Moguce je iskoristiti impulsni odziv RRC filtra predajnika za
% uoblicavanje signala
qam16_rx_I_filt = filter(hrrc,1,qam16_rx_I_c);
qam16_rx_Q_filt = filter(hrrc,1,qam16_rx_Q_c);

% Racunanje impulsnog odziva RC filtra koji se koristi u pomeracu_faze
rc_alpha = 0.22;                   
rc_M = 4;                          
rc_k = 10;                         
hrrc_rc = rrc(rc_k,rc_M,rc_alpha);  % izracunavanje imp. odziva RRC filtra
hrc = conv(hrrc_rc,hrrc_rc);        % formiranje imp. odziva RC filtra

% Parametri pomeraca faze
zk = 6;
zl = 0;
% Signali na izlazu pomeraca faze
qam16_rx_I_ps = pom_faze( qam16_rx_I_filt, hrc, rc_M, zk, zl);
qam16_rx_Q_ps = pom_faze( qam16_rx_Q_filt, hrc, rc_M, zk, zl);

% Signal na ulazu decimatora
qam16_rx_ps = qam16_rx_I_ps + 1j * qam16_rx_Q_ps;  


% Decimacija signala biranjem svakog M-tog odbirka
qam16_rx_out = qam16_rx_ps( 3*rrc_delay + 1 :rrc_M: length(qam16_rx_ps));  
% Izdvajanje I i Q komponenti
qam16_rx_out_I = real(qam16_rx_out);
qam16_rx_out_Q = imag(qam16_rx_out);


%Crtanje konstelacionog dijagrama
figure();
plot(qam16_rx_out_I, qam16_rx_out_Q,'o');
hold on;
% crtanje granice
for k = 0:2
    plot([-1.33, 1.33], [-0.667+k*0.667,-0.667+k*0.667],'r--')
    plot([-0.667+k*0.667,-0.667+k*0.667],[-1.33,1.33],'r--')
end
title('Konstelacioni dijagram 16-QAM');
xlabel('I');
ylabel('Q');
