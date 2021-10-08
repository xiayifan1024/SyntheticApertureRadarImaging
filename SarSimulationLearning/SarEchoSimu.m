% SAR点目标仿真
% 公式全部来自Cumming书
%----------------------------
% 后期可以添加sinc包络
%---------------------------
% 参考
% https://gitee.com/zhaofei2048/sar-algorithm 代码最规范,美观
% https://github.com/xiayifan1024/SAR-Synthetic-Aperture-Radar/ RD,CS完整
% https://wenku.baidu.com/view/2ed300dbcc1755270622080f?fr=uc RD,CS简单易懂
%---------------------------
clc;close all;clear;
%---------------------------
% 参数见表5.2
R_nc = 20e3;            % 景中心斜距
% H = 5e3;                % 高度
Tr = 25e-6;             % 发射脉冲时宽
Kr = 0.25e12;             % 距离调频率,正扫频
% Kr = -0.25e12;           % 负扫频
% Bw = 100e6             % 信号带宽,由Tr*Kr直接推
Fr = 7.5e6;             % 距离采样率
% X = 10e3               % 斜距条带宽度 
Vr = 150;               % 雷达有效速度
f0 = 5.3e9;             % 雷达工作频率
% lambda = 0.032;         % 波长(波长和载频可以相互推)
% Ka = 131;               % 方位向调频率
% La = 1                 % 真实孔径长度
% Ls = 850               % 合成孔径长度
% Ta = 3.4               % 目标照射时间
df_dop = 80;           % 多普勒带宽
Fa = 104;               % 方位采样率(PRF)
theta_rc = 0*pi/180;        % 波束斜视角,弧度,目前斜视角不等于零存在问题  

Naz = 256;          	% 距离线数（即数据矩阵，行数）
Nrg = 256;              % 距离线采样点数（即数据矩阵，列数）
c = 3e8;                % 光速
% 
% f_nc = 2055;            % 多普勒中心频率
% eta_c = -51.7;          % 波束中心偏移时间
% 
% alpha_os_r = 1.2;  % 距离项过采样系数
% alpha_os_a = 1.3;  % 方位向过采样系数
%----------------------------
PRF = Fa;
PRT = 1/PRF;
lambda = c/f0;          %波长(波长和载频可以相互推)
Ka = 2*Vr^2*(cos(theta_rc))^2/lambda/R_nc;  %方位向调频率
f_nc = 2*Vr*sin(theta_rc)/lambda;   % 多普勒中心频率，根据公式（4.34）计算

dr=c/2/Fr;              %系统的距离向分辨率
da=Vr*PRT;              %系统的方位向分辨率

Ymin = R_nc*cos(theta_rc)-dr*Nrg/2;   %距离项最近点(由坐标确定采样点数和由采样点数确定坐标)
Ymax = R_nc*cos(theta_rc)+dr*Nrg/2;   %距离项最远点
Xmin = 0;                                    %方位向最小位置
Xmax = da*Naz;                               %方位向最大位置
R0 = linspace(Ymin,Ymax,Nrg)*cos(theta_rc);              %R0

La = 0.886*2*Vr*cos(theta_rc)/df_dop;                  %真实孔径 4.36
Ls = 0.886*lambda*R_nc/La;               %合成孔径
%----------------------------
% ta = linspace(Xmin/Vr,Xmax/Vr,Naz);   %方位向时间域(慢时间域)
% % ta = (0:Naz-1)/Fa; 
% % ta = linspace((Xmin-Xmax)/Vr,(Xmax-Xmin)/Vr,Naz);
ta = (0:Naz-1)/Fa - R_nc*sin(theta_rc)/Vr;   % 方位时间轴
% % 注意要 - R_nc*sin(theta_rc)/Vr,这是景中心点对应的相对波束中心穿越时刻
fa = ((0:(Naz-1))-Naz/2)/Naz*Fa;  %方位向频率域(此处提前转置之后可以不需要再考虑构建矩阵的问题了)


%tr = linspace(2*Ymin/c,2*Ymax/c,Nrg); %距离项时间域(快时间域)
tr = 2*R_nc/c+(-Nrg/2:(Nrg/2-1))/Fr;
fr = f_nc+((0:(Nrg-1))-Nrg/2)/Nrg*Fr;              %距离项频率域
%-----------------------------
Target = [(Xmax-Xmin)/2,R_nc*cos(theta_rc),1
%     3*(Xmax-Xmin)/8,R_nc+dr*Nrg/4,1
%     3*(Xmax-Xmin)/8,R_nc-dr*Nrg/4,1
%     5*(Xmax-Xmin)/8,R_nc+dr*Nrg/4,1
%     5*(Xmax-Xmin)/8,R_nc-dr*Nrg/4,1
];                                   %目标(X方位向,Y距离项,反射系数)
%-----------------------------
%构造数据
s = zeros(Naz,Nrg);
[Ntarget,temp] = size(Target);%获取目标数
for k = 1:1:Ntarget
    dX = Vr*ta-Target(k,1);%相当于将坐标轴零点移动到目标点
    Rt = sqrt(dX.^2+Target(k,2)^2);%斜距(行向量)
    dt = 2*Rt/c;
    %发送与接收回波延时(行向量),dt包含了慢时间域ta的信息
    %通过dt将慢时间域ta与快时间域tr联系
    tm = ones(Naz,1)*tr - dt'*ones(1,Nrg);
    %建立时间矩阵time matrix,慢时间向量*快时间向量,行向量为快时间轴,列向量为慢时间轴
    %大小(Naz,Nrg),"耦合"
    phase = pi*Kr*tm.^2 - (4*pi/lambda)*(Rt'*ones(1,Nrg));%相位
    
    flag1 = (abs(ta*Vr-(Target(k,1)-Target(k,2)*tan(theta_rc)))<=Ls/2);%包络,目标束缚,排除超过合成孔径长度,方位向约束
    flag2 = (abs(tm)<=Tr/2);%包络,排除时间矩阵中超过脉宽一半的,距离向约束
    s = s+Target(k,3)*exp(1j*phase).*(flag1'*ones(1,Nrg)).*flag2;
end

%-------------------------------
% % 幅度图,没有像图5.16渐变是由于未添加包络
% figure(1)
% imagesc(255-abs(s))
% colormap(gray)
% % 相位图,正负扫频图像不同
% figure(2)
% imagesc(255-angle(s))
% colormap(gray)
% % 二维频谱幅度
% figure(3)
% imagesc(255-abs(fft2(s)))
% colormap(gray)
% % 二维频谱相位
% figure(4)
% imagesc(255-angle(fft2(s)))
% colormap(gray)

% RDA
% CSA
% wkA
