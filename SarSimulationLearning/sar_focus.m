function [image, axesOut] = sar_focus(s,tr,ta,f0,Kr,Vr,c,theta_rc,method)
%SAR_FOCUS 对匀速直线、固定斜视的条带 SAR 回波进行聚焦。
% s: 方位行、距离列；tr/ta: 实际快/慢时间 (s)；theta_rc: rad。
% 输出距离为波束中心斜距 R0/cos(theta_rc)，方位为目标沿轨坐标。
% RD/CS 使用二阶距离频率展开；WK 使用逐多普勒中心化的 Stolt 插值。
% 原始回波必须覆盖目标照射时间和脉冲，去中心后的频谱不能混叠。
[Naz,Nrg] = size(s);
tr = tr(:).';
ta = ta(:).';
assert(Naz>=2 && Nrg>=2 && length(ta)==Naz && length(tr)==Nrg, ...
    '回波尺寸必须与时间轴一致，且每个方向至少有两个采样点。');
assert(c>0 && Vr>0 && f0>0 && Kr~=0 && abs(theta_rc)<pi/2, ...
    '载频、速度和光速必须为正，调频率非零，斜视角绝对值小于 90 度。');
dtRange = tr(2)-tr(1);
dtAz = ta(2)-ta(1);
assert(dtRange>0 && dtAz>0, '时间轴必须递增。');
assert(max(abs(diff(tr)-dtRange))<1e-6*dtRange && ...
    max(abs(diff(ta)-dtAz))<1e-6*dtAz, '时间轴必须均匀采样。');
Fr = 1/dtRange;
Fa = 1/dtAz;
lambda = c/f0;
Dref = cos(theta_rc);
f_nc = 2*Vr*sin(theta_rc)/lambda;
faBase = (-floor(Naz/2):ceil(Naz/2)-1).'*Fa/Naz;
fa = faBase+f_nc;                 % 物理多普勒频率，不能用折叠后的中心频率
fr = (-floor(Nrg/2):ceil(Nrg/2)-1)*Fr/Nrg; % 距离基带不加 f_nc
radicand = 1-(lambda*fa/(2*Vr)).^2;
assert(all(radicand>0), '方位频率超出传播波数范围。');
D = sqrt(radicand);
Rnc = c*tr(floor(Nrg/2)+1)/2;
Rref = Rnc*Dref;                  % 参考最近斜距
range = c*tr/2;                   % 波束中心斜距网格
R0 = range*Dref;
azShift = Rnc*sin(theta_rc)/Vr;    % 把零多普勒坐标平移回采集窗口
axesOut.range = range;
axesOut.azimuth = Vr*(ta+azShift);

% 先去多普勒中心，再 FFT；非零斜视下 f_nc 可以远大于 PRF。
Srd = ftx(bsxfun(@times,s,exp(-2j*pi*f_nc*ta.')));
Ha = exp(4j*pi/lambda*(D*R0));
shift = exp(2j*pi*faBase*azShift);
invKm = 1/Kr-c*Rref*fa.^2./(2*Vr^2*f0^3*D.^3);
assert(all(abs(invKm)>eps/abs(Kr)), '有效距离调频率出现奇点。');
Km = 1./invKm;

switch upper(method)
    case 'RD'
        % 同时完成距离压缩和二次距离压缩 (SRC)。
        Hrange = exp(1j*pi*(invKm*fr.^2));
        compressed = ifty(fty(Srd).*Hrange);
        % 精确徙动量 R0/D-R0/Dref；不再用小多普勒二次近似。
        query = (bsxfun(@rdivide,R0,D)-range(1))/(c/(2*Fr))+1;
        corrected = sar_sinc_interp(compressed,query);
        image = iftx(bsxfun(@times,corrected.*Ha,shift));
    case 'CS'
        % 变标必须使用快时间 tr，且距离 FFT 的输入为变标后的 S1。
        timeOffset = bsxfun(@minus,tr,2*Rref./(c*D));
        Hscale = exp(1j*pi*bsxfun(@times,Km.*(Dref./D-1),timeOffset.^2));
        S1 = Srd.*Hscale;
        Hrange = exp(1j*pi*((D.*invKm/Dref)*fr.^2) + ...
            4j*pi*Rref/c*((1./D-1/Dref)*fr));
        S2 = ifty(fty(S1).*Hrange);
        % 变标引入的剩余相位含距离差的平方。
        residualRange = bsxfun(@rdivide,R0-Rref,D);
        Hresidual = exp(-4j*pi/c^2* ...
            bsxfun(@times,Km.*(1-D/Dref),residualRange.^2));
        image = iftx(bsxfun(@times,S2.*Ha.*Hresidual,shift));
    case 'WK'
        root = sqrt(bsxfun(@minus,(f0+fr).^2,(c*fa/(2*Vr)).^2));
        assert(isreal(root), '距离/方位频率超出传播波数范围。');
        % 去掉快时间 FFT 原点，再用参考目标抵消快速振荡相位。
        S1 = bsxfun(@times,fty(Srd),exp(-2j*pi*fr*(2*Rnc/c)));
        Href = exp(bsxfun(@plus,4j*pi*Rref/c*root,1j*pi*fr.^2/Kr));
        S1 = S1.*Href;
        % 每个多普勒处取 k_r=f0*D+fr_new/Dref，避免斜视下固定
        % Stolt 网格平移数十 MHz 后完全超出原始距离频谱。
        newRoot = bsxfun(@plus,f0*D,fr/Dref);
        sourceFrequency = sqrt(bsxfun(@plus,newRoot.^2, ...
            (c*fa/(2*Vr)).^2))-f0;
        S2 = sar_sinc_interp(S1,(sourceFrequency-fr(1))/(Fr/Nrg)+1);
        Srange = ifty(S2);
        Hresidual = exp(4j*pi*f0/c*(D*(R0-Rref)));
        image = iftx(bsxfun(@times,Srange.*Hresidual,shift));
    otherwise
        error('未知算法 %s；应为 RD、CS 或 WK。',method);
end
end
