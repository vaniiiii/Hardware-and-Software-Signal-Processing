function [ xI ] = expander( x,I )
%EXPANDER - Funkcija povecava ucestanost odabiranja ulaznog signala I puta
%   Signal xI se sastoji od odbiraka ulaznog signala x izmedju kojih je
%   ubaceno I nula.
%   x - odbirci ulaznog signala
%   I - faktor povecanja ucestanosti odabiranja

len = length(x)*I;
xI = zeros(1,len);
xI(1:I:len) = x;

end

