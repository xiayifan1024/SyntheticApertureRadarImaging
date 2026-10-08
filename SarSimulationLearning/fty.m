function out = fty(s)
% 中心化FFT，显式指定第 2 维，支持奇数长度。
out = fftshift(fft(ifftshift(s,2),[],2),2);
end
