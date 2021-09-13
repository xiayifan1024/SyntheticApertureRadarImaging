%-------------omega-K算法---------------
% 目前只适用于正侧视
%-------------二维FFT-------------------
%S_f = fft2(s);
S_f = ftx(fty(s));
%-------------参考函数相乘---------------
theta_ref = 4*pi*R_nc/c*sqrt((f0+ones(Naz,1)*fr).^2-c^2*(fa'*ones(1,Nrg)).^2/4/Vr^2)+pi*ones(Naz,1)*fr.^2/Kr;
Href = exp(1i*theta_ref);
S1 = S_f.*Href;
S1=S1.*exp(-1j*4*pi*R_nc/c*(ones(Naz,1)*fr+f0));  % 由公式5.11,将频谱移动到中心
% 思考:为什么wk算法需要这样
% %--------------精确实现算法--------------
% %----------------变量代换----------------
% % 将原来的距离项频率fr替换成新的距离项频率frNew
% % 在fa-frNew中的频谱图像发生了形变
% % 此时再将发生形变的频谱再放入fa-fr中
% 
% fr_ = sqrt((f0+ones(Naz,1)*fr).^2-c^2*(fa'*ones(1,Nrg)).^2/4/Vr^2)-f0;     % 公式8.5
% % plot(fr,fr_)                % 新距离向频率fr_与原距离向频率之间的关系
% % 新的fr_与fr呈现线性关系,此外fr_这里是一个矩阵的形式,实际上这个矩阵是fr_与fa的耦合
% frmax_ = max(max(fr_));
% frmin_ = min(min(fr_));
% frNew = linspace(frmin_,frmax_,Nrg);    % 重建fr_序列,本质上也是重采样的过程
% frOrigin = sqrt((f0+ones(Naz,1)*frNew).^2+c^2*(fa'*ones(1,Nrg)).^2/4/Vr^2)-f0; 
% % 逆用公式8.5,获得fr_对应在老fr域中的序列;
% % frOrigin矩阵的值是不同(fa,fr_)所对应的原本的fr值
% 
% % frOrigin = sqrt((f0+ones(Naz,1)*fr).^2+c^2*(fa'*ones(1,Nrg)).^2/4/Vr^2)-f0; %
% % 直接公式8.5逆变量代换好像也行
% 
% %-----------------Stolt插值---------------
% % 将(fa,fr)坐标的值付给(fa,frOrigin)
% % 注意此处有两种意义上的矩阵
% % 一个是信号矩阵S,一个是频率矩阵frOrigin,两者大小都是1024*512,意义不同
% % 有点类似于RD算法中的sinc插值
% 
% 
% 
% 
% S2 = zeros(Naz,Nrg);
% 
% delta = frOrigin - ones(Naz,1)*fr;      % 频率差
% dFr = Fr/(Nrg-1);                       % 单位距离项频率差
% deltaPix = delta/dFr;                   % 像素差矩阵
% 
% P = 8;
% for a = 1:Naz
%     for r = P/2:Nrg                                       %防止左越界
%         for i = -P/2+1:P/2
%             if r+deltaPix(a,r)+i > Nrg
%                 S2(a,r) = S2(a,r)+S1(a,Nrg)*sinc(deltaPix(a,r)-i);
%             elseif r+deltaPix(a,r)+i <= 0
%                 S2(a,r) = S2(a,r)+S1(a,1)*sinc(deltaPix(a,r)-i);
%             else
%                 S2(a,r) = S2(a,r)+S1(a,r+fix(deltaPix(a,r))+i)*sinc(deltaPix(a,r)-i);
%             end
%         end
%     end
% end
% 
% %-------------二维FFT逆变换---------------
% S = iftx(ifty(S2));
%--------------近似实现算法---------------
%---------------距离向IFFT----------------
Sr = ifty(S1);
%------------补余方位匹配滤波相乘----------
D_nr = sqrt(1-c^2*fa.^2/4/Vr^2/f0^2)'*ones(1,Nrg);   
theta_ = 4*pi*(ones(Naz,1)*R0-R_nc)*f0.*D_nr/c;                    % 式8.35
Src = Sr.*exp(1i*theta_);
%---------------方位向IFFT----------------
S = iftx(Src);
%-----------------------------------------
imagesc(abs(S))
