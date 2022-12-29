function [ y ] = adc( x, N,Vfs )
%ADC Funkcija simulira rad A/D konvertora prema ulaznim parametrima
%   y - kvantizovan signal
%   x - ulazni signal
%   N - rezolucija konvertora
%   Vfs - napon pune skale - pretpostavlja se centriranje oko 0

step = Vfs / 2^(N);                             % jedan kvant
min_val = -Vfs/2;                               % granice kvantizatora
max_val = Vfs/2;
ofs = step * 0.5;                               % ofset

x(x >= max_val) = max_val;                      % zasicenje
x(x <= min_val) = min_val;                      

y = (round((x - min_val)./ step)) * step;          %kvantizacija
y = y + min_val;
y (y > max_val - ofs) = max_val - step;


end
