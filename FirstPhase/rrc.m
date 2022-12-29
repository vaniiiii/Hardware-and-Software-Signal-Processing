function [ h ] = rrc( k, M, alpha )
%RRC -  Funkcija vraca impulsni odziv filtra "koren podignutog kosinusa" 
%       h prema zadatim parametrima k, m, alpha gde su:
%       k - duzina impulsnog odziva u broj simbola
%       M - faktor interpolacije
%       alpha - parametar definise oblik/sirinu filtra

nsym = k/2;
n = -M*nsym:1:M*nsym;
h = (sin(pi*n./M*(1-alpha))+4*alpha*n./M.*cos(pi*n./M*(1+alpha))) ./ ...
    (pi*n.*(1-(4*alpha*n./M).^2));
h(nsym*M+1) = (1+alpha*(4/pi-1))/M;
h = h./sqrt(sum(h.^2));

end
