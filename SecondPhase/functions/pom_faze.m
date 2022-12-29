function [ y ] = pom_faze( x, hrc, rc_M, zk, zl )
%POM_FAZE Funckija modeluje programabilni pomerac faze prijemnika
%   Ulazni signal se prvenstveno zakasni za onoliko odbiraka koliko je
%   definisano parametrom zk, a onda se u izlazni niz smestaju ti odbirci
%   filtrirani polifaznom komponentom RC filtra koja je odredjena
%   parametrom zl
%   x       -   ulazni signal
%   hrc     -   imp. odziv RC filtra
%   rc_M    -   faktor interpolacije RC filtra
%   zk      -   kasnjenje
%   zl      -   izbor polifazne komponente
%   y       -   izlazni signal

x = [ zeros(1,zk) x(1:(length(x)-zk)) ];  % kasnjenje ulaznog signala

hrc_poly = make_poly(hrc,rc_M);     % polifazna dekompozicija RC filtra

% Filtriranje ulaznog signala odredjenom polifaznom komponentom i 
% smestanje rezultata u izlazni niz
y = filter( hrc_poly(mod((4-zl),4)+1,:), 1, x );

end

