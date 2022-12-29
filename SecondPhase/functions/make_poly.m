function y = make_poly(x,M)
%MAKE_POLY Ova funkcija pravi polifaznu dekompoziciju filtra
%   x - Koeficijenti ulaznog filtra
%   M - broj zeljenih polifaznih komponenti

n = ceil(length(x)/M);      % najmanji broj polifaznih dekompozicija
tmp = zeros(1,n*M);           
tmp(1:1:length(x)) = x;     % niz koji sadrzi x i popunjen je nulama 

y = zeros(M,n);             % alokacija izlaznog niza

%U i-tu polifaznu dekompozicija smestamo svaki M-ti koeficijent ulaznog
% filtra pocev od i-tog
for i = 1:1:M
    y(i,:) = tmp(i:M:length(tmp));
end

end

