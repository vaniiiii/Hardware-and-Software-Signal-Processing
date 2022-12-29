function [ y ] = filter_poly( sig,FIR_poly,M )
%FILTER_POLY Filtriranje polifaznom dekompozicijom filtra
%   sig - ulazni signal
%   FIR_poly - polifazna dekompozicija FIR filtra
%   M - faktor interpolacije

sig_filt = zeros(M,length(sig));    % alociranje niza za filtrirane odbirke

% Filtriranje ulaznog signala svakom polifaznom komponentom
for i = 1:1:M
    sig_filt(i,:) = filter(FIR_poly(i,:),1,sig);
end

y=zeros(1,M*length(sig_filt));      % alociranje izlaznog niza

% Preuredjivanje redosleda odbiraka tako da se u prvih M pozicija stave
% prvi odbirci svih polifazno filtriranih komponenti, u sledecih M pozicija
% drugi odbirci svih polifazno filtriranih komponenti itd...
for k=1:length(sig_filt)
    for i=1:M
        y( (k-1) * M + i ) = sig_filt(i,k);
    end
end

end

