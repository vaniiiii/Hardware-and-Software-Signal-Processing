function [xnew,ynew] = CORDIC(x0,y0,z0,num_of_bits) 
% CORDIC - Funkcija racuna CORDIC algoritam za ulazne vrednosti simbola 
% i ugla na izlazu faznog akumulatora

% x0,y0,z0 - ulazi
% xnew,ynew - izlazi
% num_of_bits - broja bita za aritmetiku sa fiksnim zarezom
% mod: Rotacioni

FixedPointAttributes = fimath ( 'RoundingMethod', 'Round',...
    'OverflowAction', 'Saturate', 'ProductMode', 'SpecifyPrecision', ...
    'ProductWordLength', num_of_bits.SIGNAL_BITLENGTH, ...
    'ProductFractionLength', num_of_bits.SIGNAL_FRACLENGTH, ... 
    'SumMode', 'SpecifyPrecision', 'SumWordLength', ...
    num_of_bits.SIGNAL_BITLENGTH,'SumFractionLength', ...
    num_of_bits.SIGNAL_FRACLENGTH ) ;

% Normalizacija ugla z0:

z0 = 2*pi*z0/(2^num_of_bits.SIGNAL_FRACLENGTH);

% opseg [0,2pi]
if(z0>=2*pi)
    z0 = z0-2*pi;
end    

% Prerotacija vektora - U rotacionom modu je potrebno obezbediti z<z_max 
% Pod pretpostavkom da je  2pi?theta?0 
if ((z0>pi/2) && (z<=3*pi/2))
    z = z0 - pi;
    x = -x0;
    y = -y0;


else if (z0>3*pi/2)
    z = z0-3*pi/2;
    x = y0;
    y = -x0;

else   
    z = z0;   
    x = x0;
    y = y0;

    end
end

% Broj iteracija CORDIC-a = Duzina razlomljenog dela reci
for i=0:num_of_bits.SIGNAL_FRACLENGTH-1

    % sigma u rotacionom modu
    sigma = -1;
    if (z >= 0)  %proveri ovo, da l' je nebitno koja se vrednost uzme ako je 0
        sigma = 1;
    end
    % Zaokruzivanje vrednosti, aritmetika sa fiksnim zarezom
    x = fi (x , true , num_of_bits.SIGNAL_BITLENGTH, ...
        num_of_bits.SIGNAL_FRACLENGTH, FixedPointAttributes);
    y = fi (y , true , num_of_bits.SIGNAL_BITLENGTH, ...
        num_of_bits.SIGNAL_FRACLENGTH, FixedPointAttributes);
    xnew = fi (zeros(1, length(x)) , true , num_of_bits.SIGNAL_BITLENGTH, ...
        num_of_bits.SIGNAL_FRACLENGTH, FixedPointAttributes);
    ynew = fi (zeros(1, length(y)) , true , num_of_bits.SIGNAL_BITLENGTH, ...
        num_of_bits.SIGNAL_FRACLENGTH, FixedPointAttributes);
    z = fi (z , true , num_of_bits.SIGNAL_BITLENGTH, ...
        num_of_bits.SIGNAL_FRACLENGTH, FixedPointAttributes);

    % U rotacinom modu izlazne vrednosti se racunaju kao:
    xnew = x - sigma * y * 2^(-i);
    ynew = y + sigma * x * 2^(-i);
    z = z - sigma * atan(2^(-i));

    % Azuriranje vrednosti za sledecu iteraciju
    x = xnew;
    y = ynew;
end


end
