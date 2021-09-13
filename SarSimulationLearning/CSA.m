%-------------CS算法---------------
% 目前只适用于正侧视
%------------方位向FFT-------------
Srd = ftx(s);
%-----------建立变标方程------------
R_0=ones(Naz,1)*R0;                                                           % R0矩阵化
D_nr = sqrt(1-c^2*fa.^2/4/Vr^2/f0^2)'*ones(1,Nrg);                            % 徙动参数
Km = Kr./(1-Kr*c*R_0.*(fa'*ones(1,Nrg)).^2/2/Vr^2/f0^3./D_nr.^3);             % 调整后的距离向调频率
fa_ref = fa(1,Naz/2+1);                                                       % 参考点(中心点)方位向频率
Vrref = Vr;                                                                   % 假设参考点速度不变
D_nref_rref = sqrt(1-c^2*fa_ref^2/4/Vrref^2/f0^2);
D_n_rref = sqrt(1-c^2*fa.^2/4/Vrref^2/f0^2)'*ones(1,Nrg);
ta_ = ta'*ones(1,Nrg)- 2*R_nc/c./D_n_rref;                                    % 时间坐标参考点变换
Ssc = exp(1j*pi*Km.*(D_nref_rref./D_n_rref-1).*(ta_).^2);                     % 建立变标方程
S1 = Ssc.*Srd;
% %------------距离项FFT-----------------
S_d = fty(Srd);
% %------------距离压缩------------------
I1 = exp(1i*pi*D_nr./Km./D_nref_rref.*(ones(Naz,1)*fr).^2);                   % 公式7.32第2指数项,注意正负取反
I2 = exp(1i*4*pi/c*(1./D_n_rref-1./D_nref_rref)*R_nc.*(ones(Naz,1)*fr));      % 公式7.32第4指数项,同上
S2_d = S_d.*I1.*I2;
% %------------距离项IFFT----------------
S2 = ifty(S2_d);
% %------------方位压缩及相位矫正---------
I3 = exp(1i*4*pi*R_0*f0.*D_nr/c);                                             % 公式7.34第2指数项,同上                                      
I4 = exp(-1i*4*pi*Km/c^2.*(1-D_n_rref./D_nref_rref).*(R_0./D_nr-R_nc./D_nr)); % 公式7.34第3指数项,同上
S3_rd = S2.*I3.*I4;
% %------------方位项IFFT----------------
S3 = ftx(S3_rd);                                                             % 为什么这里直接就用ifft就行?

figure(5)
imagesc(abs(S3))

%--------------点目标分析------------------
