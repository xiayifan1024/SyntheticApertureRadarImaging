function s = iftx(fs)
%IFTX 此处显示有关此函数的摘要
%   此处显示详细说明
s = ifftshift(ifft(ifftshift(fs)));
%s = ifftshift(ifft(ifftshift(fs,1)),1);
end

