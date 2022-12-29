% Moguce je iskoristiti impulsni odziv RRC filtra predajnika za
% uoblicavanje signala
qam4_rx_I_filt = filter(hrrc,1,qam4_rx_I_c);
qam4_rx_Q_filt = filter(hrrc,1,qam4_rx_Q_c);

% Racunanje impulsnog odziva RC filtra koji se koristi u pomeracu_faze
rc_alpha = 0.22;                   
rc_M = 4;                          
rc_k = 10;                         
hrrc_rc = rrc(rc_k,rc_M,rc_alpha);  % izracunavanje imp. odziva RRC filtra
hrc = conv(hrrc_rc,hrrc_rc);        % formiranje imp. odziva RC filtra

% Parametri pomeraca faze
zk = 6;
zl = 4;

for zl = 0
    % Signali na izlazu pomeraca faze
    qam4_rx_I_ps = pom_faze( qam4_rx_I_filt, hrc, rc_M, zk, zl);
    qam4_rx_Q_ps = pom_faze( qam4_rx_Q_filt, hrc, rc_M, zk, zl);

    % Signal na ulazu decimatora
    qam4_rx_ps = qam4_rx_I_ps + 1j * qam4_rx_Q_ps;  


    % Decimacija signala biranjem svakog M-tog odbirka
    qam4_rx_out = qam4_rx_ps( 3 * rrc_delay + 1 : rrc_M : length(qam4_rx_ps));  
    % Izdvajanje I i Q komponenti
    qam4_rx_out_I = real(qam4_rx_out);
    qam4_rx_out_Q = imag(qam4_rx_out);


    %Crtanje konstelacionog dijagrama
    figure(3);
    plot(qam4_rx_out_I, qam4_rx_out_Q,'o');
    hold on;
end

plot([-1.33, 1.33], [0,0],'r--');
plot([0,0],[-1.33,1.33],'r--');
legend('zl = 1','zl = 2','zl = 3','zl = 4','Location', 'east')
title('Konstelacioni dijagram 4-QAM za razlicite vrednosti zl');
xlabel('I');
ylabel('Q');