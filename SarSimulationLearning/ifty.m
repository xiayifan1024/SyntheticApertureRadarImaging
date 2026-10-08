function out = ifty(s)
% 中心化IFFT，显式指定第 2 维，支持奇数长度。
out = fftshift(ifft(ifftshift(s,2),[],2),2);
end
