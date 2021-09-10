function fs = ftx(s)
%FTX 此处显示有关此函数的摘要
%   此处显示详细说明
fs = fftshift(fft(fftshift(s)));
%fs=fftshift(fft(fftshift(s,1)),1);
end

