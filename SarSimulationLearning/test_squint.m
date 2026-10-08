function results = test_squint()
%TEST_SQUINT 正/负斜视、偏心目标、奇数/非方阵 FFT 回归，无需绘图。
% 在本目录运行 results = test_squint;
rng(0);
for shape = [256,257;320,319]
    data = randn(shape.')+1j*randn(shape.');
    restored = iftx(ifty(ftx(fty(data))));
    assert(norm(restored-data,'fro')/norm(data,'fro')<1e-12, ...
        '二维 FFT 往返失败。');
end
% 插值应保持整数位置，并正确处理两端越界和分数偏移。
signal = ones(2,64);
query = repmat([-1,1,20.25,40.75,64,66],2,1);
resampled = sar_sinc_interp(signal,query);
assert(max(max(abs(resampled(:,[2,3,4,5])-1)))<1e-12);
assert(all(all(resampled(:,[1,6])==0)));

c = 3e8; Vr = 150; f0 = 5.3e9; lambda = c/f0;
Rnc = 20e3; Fr = 7.5e6; Fa = 104; Tr = 25e-6;
dr = c/(2*Fr); da = Vr/Fa;
methods = {'RD','CS','WK'};
angles = [0,5,-5,15,-15,30,-30];
results = struct('angle',{},'method',{},'shape',{},'target',{}, ...
    'positionErrorPixels',{},'width3dB',{},'peak',{});
for dimensions = [256,257;256,320]
    Naz = dimensions(1); Nrg = dimensions(2);
    for angle = angles
        theta = angle*pi/180;
        ta = (0:Naz-1)/Fa-Rnc*sin(theta)/Vr;
        tr = 2*Rnc/c+(-floor(Nrg/2):ceil(Nrg/2)-1)/Fr;
        Ls = lambda*Rnc*80/(2*Vr*cos(theta));
        targets = [floor(Naz/2)*da,Rnc*cos(theta); ...
            (floor(Naz/2)-32)*da,(Rnc-4*dr)*cos(theta); ...
            (floor(Naz/2)+32)*da,(Rnc+4*dr)*cos(theta)];
        for Kr = [0.25e12,-0.25e12]
            for k = 1:size(targets,1)
                x = targets(k,1); y = targets(k,2);
                % 精确几何距离生成回波；不使用成像算法的近似公式。
                Rt = sqrt((Vr*ta-x).^2+y^2);
                tm = bsxfun(@minus,tr,2*Rt.'/c);
                beam = abs(Vr*ta-(x-y*tan(theta)))<=Ls/2;
                phase = bsxfun(@minus,pi*Kr*tm.^2,4*pi*Rt.'/lambda);
                echo = bsxfun(@times,exp(1j*phase).*(abs(tm)<=Tr/2),beam.');
                for m = 1:length(methods)
                    [image,coords] = sar_focus(echo,tr,ta,f0,Kr,Vr,c,theta,methods{m});
                    assert(all(isfinite(image(:))), '成像结果含 NaN/Inf。');
                    magnitude = abs(image);
                    [peak,index] = max(magnitude(:));
                    [a,r] = ind2sub(size(image),index);
                    errorPixels = [abs(coords.azimuth(a)-x)/da, ...
                        abs(coords.range(r)-y/cos(theta))/dr];
                    width = [sum(magnitude(:,r)>=peak/sqrt(2)), ...
                        sum(magnitude(a,:)>=peak/sqrt(2))];
                    assert(peak>0 && all(errorPixels<=1), ...
                        '%s 在 %g 度斜视下位置错误。',methods{m},angle);
                    assert(all(width<=2), ...
                        '%s 在 %g 度斜视下未聚焦。',methods{m},angle);
                    item.angle = angle;
                    item.method = methods{m};
                    item.shape = [Naz,Nrg];
                    item.target = [x,y,Kr];
                    item.positionErrorPixels = errorPixels;
                    item.width3dB = width;
                    item.peak = peak;
                    results(end+1) = item; %#ok<AGROW>
                end
            end
        end
    end
end
fprintf('通过 %d 个聚焦用例，以及 FFT 和插值检查。\n',length(results));
end
