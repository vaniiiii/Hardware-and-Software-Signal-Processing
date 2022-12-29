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
for zl = 1:4
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
figure(4);
hold on;
plot(qam16_rx_out_I, qam16_rx_out_Q,'o');
end

legend('zl = 1','zl = 2','zl = 3','zl = 4','Location', 'northeast')
% crtanje granice
for k = 0:2
    plot([-1.33, 1.33], [-0.667+k*0.667,-0.667+k*0.667],'r--')
    plot([-0.667+k*0.667,-0.667+k*0.667],[-1.33,1.33],'r--')
end
title('Konstelacioni dijagram 16-QAM');
xlabel('I');
ylabel('Q');