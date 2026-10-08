function out = sar_sinc_interp(data, query, taps)
%SAR_SINC_INTERP 逐行重采样；query 为从 1 开始的浮点采样位置。
% Hann 加窗 sinc，分数偏移相对于 floor(query)，越界零填充。
if nargin < 3
    taps = 16;
end
assert(mod(taps,2)==0 && taps>=2, '插值长度必须为正偶数。');
[naz,nrg] = size(data);
assert(size(query,1)==naz, '查询矩阵的行数必须与数据相同。');
base = floor(query);
out = complex(zeros(size(query)));
weightSum = zeros(size(query));
rows = repmat((1:naz).',1,size(query,2));
for k = -taps/2+1:taps/2
    index = base+k;
    delta = query-index;
    weight = ones(size(delta));
    nonzero = delta~=0;
    weight(nonzero) = sin(pi*delta(nonzero))./(pi*delta(nonzero));
    weight = weight.*(0.5+0.5*cos(2*pi*delta/taps));
    valid = index>=1 & index<=nrg;
    index = min(max(index,1),nrg);
    values = data(sub2ind([naz,nrg],rows,index));
    out = out+values.*weight.*valid;
    weightSum = weightSum+weight;
end
out = out./weightSum;
out(query<1 | query>nrg) = 0;
end
