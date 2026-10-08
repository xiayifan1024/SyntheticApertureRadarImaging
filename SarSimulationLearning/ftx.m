function out = ftx(s)
% 中心化FFT，显式指定第 1 维，支持奇数长度。
out = fftshift(fft(ifftshift(s,1),[],1),1);
end
