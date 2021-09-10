%-------------RD算法-------------
% 适用于正侧视或小斜视角
%-------------距离压缩-----------
% 距离压缩 方式一 直接生成
% Hr = exp(1j*pi/Kr*fr.^2);
% s_mf = ifty(fty(s).*Hr);

% 距离压缩 方式二 复制脉冲补零进行DFT,结果取复共轭
Hr = conj(fty(exp(1i*pi*Kr*(tr-2*R_nc/c).^2).*(abs(tr-2*R_nc/c)<=(Tr/2))));
s_mf = ifty(fty(s).*Hr);

% % 距离压缩 方式三 时间反褶取共轭，计算补零DFT
% Hr = fty(conj(exp(1i*pi*Kr*fliplr(tr-2*R_nc/c).^2).*(abs(tr-2*R_nc/c)<=(Tr/2))));
% s_mf = ifty(fty(s).*Hr);


figure(2)
imagesc(abs(s_mf))
%colormap(gray)
%------------距离徙动矫正-------------
s_afft = ftx(s_mf);                 % 方位fft
RCM_matrix = zeros(Naz,Nrg);
Drcm = lambda^2*(ones(Naz,1)*R0).*(fa'*ones(1,Nrg)).^2/8/Vr^2;      
% 距离徙动矫正需要偏移的差
RCM = Drcm/dr;                          % 距离徙动矫正需要偏移像素数
% % 由于距离徙动,本来在(X,Y)位置的落入了(X,Y+Drcm)位置,因而要将(X,Y+Drcm)的插值给(X,Y)
% % 此时X不变,Y改变,看做X(即Naz)个一维插值
% % 由于是网格化,所以坐标用网格位置表示,插值也需要取整插值
% % 可以看做重新重采样插值的过程
% %----------最近邻法-----------------
% for m = 1:Naz
%     for n = 1:Nrg 
%         dec = RCM(m,n)-fix(RCM(m,n)); % 截断小数
%         if RCM(m,n)+ n > Nrg 
%             %判断Y+Drcm相加后是否超出距离Ymax
%             RCM_matrix(m,n) = s_afft(m,Nrg);
%         else
%             if dec >= 0.5
%                 RCM_matrix(m,n) = s_afft(m,n+fix(RCM(m,n))+1);
%                 %这里RCM>0
%             else
%                 RCM_matrix(m,n) = s_afft(m,n+fix(RCM(m,n)));
%             end
%         end
%     end
% end
%------------sinc插值函数法--------------
% 插值需要求出非整数点(X,Y+Drcm)的值,先对原信号进行插值
% 求出非整数点后的值,在进行插值
P = 8;
for a = 1:Naz
    for r = P/2:Nrg                                       %防止左越界
        for i = -P/2+1:P/2
            if r+RCM(a,r)+i > Nrg
                RCM_matrix(a,r) = RCM_matrix(a,r)+s_afft(a,Nrg)*sinc(RCM(a,r)-i);
            else
                RCM_matrix(a,r) = RCM_matrix(a,r)+s_afft(a,r+fix(RCM(a,r))+i)*sinc(RCM(a,r)-i);
            end
        end
    end
end
              
% % ---------------------------------------
s_rmc = iftx(RCM_matrix);
figure(3)
imagesc(abs(s_rmc))
%------------方位向压缩---------------
Ha = exp(-1j*pi/Ka*fa.^2);
s_rmf = iftx(ftx(s_rmc).*Ha.');
figure(4)
imagesc(abs(s_rmf))
