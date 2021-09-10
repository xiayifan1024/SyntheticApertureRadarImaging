function fs=fty(s)
%对矩阵行做fft
fs = fftshift(fft(fftshift(s.'))).';
%fs=fftshift(fft(fftshift(s,2),[],2),2);