function s=ifty(fs)
%IFTY 此处显示有关此函数的摘要
%   对矩阵行做ifft
s=ifftshift(ifft(ifftshift(fs.'))).';
%s=ifftshift(ifft(ifftshift(fs,2),[],2),2);
end

