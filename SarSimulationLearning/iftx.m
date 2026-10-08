function out = iftx(s)
% 中心化IFFT，显式指定第 1 维，支持奇数长度。
out = fftshift(ifft(ifftshift(s,1),[],1),1);
end
