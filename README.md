# SyntheticApertureRadarImaging

3 algorithms of matlab codes —— RangeDoppler,ChirpScaling,Omega-k

RD, CS, wK MATLAB point-target simulations for broadside and fixed-squint stripmap SAR.

Some of the code needs to be optimized

Comments are Chinese, and some comments need to be optimized




RD、CS、wK 三种 MATLAB 点目标成像算法，已修正固定斜视下的频率轴、几何坐标与聚焦处理。仓库内的入口用于仿真；下方真实数据链接需要另行读取数据和配置参数。

部分代码需要优化

注释为中文,部分注释需要优化

------------------------------------------------------------------------------------------------------------

成像SAR数据,该数据是经过运动补偿后的数据,原数据也有,但是太大没有上传.(懒得上传,笑)

本数据是国科大课程《SAR信号处理与运动补偿》数据,仅供学习参考.

链接: https://pan.baidu.com/s/1Lmh1Y-atdWo701DvFllnXw 提取码: p4n2 复制这段内容后打开百度网盘手机App，操作更方便哦

------------------------------------------------------------------------------------------------------------

参考:
Ian G. Cumming, Frank H.W 著《合成孔径雷达成像算法与实现》(主要参考)

Achim Hein 著 Processing of SAR Data——Fundamentals,Signal Processing,Interferomentry

https://gitee.com/zhaofei2048/sar-algorithm 代码规范,美观(人家不给看了)

https://github.com/denkywu/SAR-Synthetic-Aperture-Radar RD,CS完整

https://wenku.baidu.com/view/2ed300dbcc1755270622080f?fr=uc RD,CS简单易懂



|            |The schedule/进度表|           |
|-------------| :-----------: |------------|
| algorithms  |   Simulation  |   Imaging  |
|RangeDoppler |        O      |            |
|ChirpScaling |        O      |            |
|   omega-K   |        O       |            |

----------------------------------------------------------------------------------------
咕了，搞SAR目标检测去了。

----------------------------------------------------------------------------------------
2026年10月08日

后记

在 codex 辅助下更新了代码中的错误，也算是结了一个这么多年的遗憾

## 2026-10-08 斜视成像修正记录

### 运行方式

在 MATLAB 中进入 `SarSimulationLearning`，运行 `SarEchoSimu`，会生成三种算法的成像结果。修改脚本的 `theta_rc`（弧度），例如 `15*pi/180` 或 `-30*pi/180`，即可模拟斜视。运行 `results = test_squint` 可检查 0°、±5°、±15°、±30°、偏心目标、正负距离调频率和奇数/非方阵数据，共 252 个聚焦用例。

### 修改内容

正侧视时 `f_nc = 0`、`cos(theta_rc) = 1`，部分频率轴和坐标错误被掩盖；引入斜视后会出现位置偏移、距离徙动校正不足或失焦。

| 文件 | 修改内容 |
| --- | --- |
| `SarSimulationLearning/SarEchoSimu.m` | 距离频率 `fr` 保持基带，方位物理频率使用 `fa_base + f_nc`；按快时间采样网格计算最近斜距，消除重复乘 `cos(theta_rc)` 和 `linspace` 步长不一致；区分采样间隔与分辨率；默认运行三种成像算法。 |
| `SarSimulationLearning/sar_focus.m`（新增） | 集中实现三种算法的公共坐标与频率处理：方位 FFT 前去多普勒中心，滤波使用物理多普勒频率，统一输出坐标并补偿方位平移。 |
| `SarSimulationLearning/RDA.m` | 调用修正后的 RD：加入二次距离压缩（SRC），使用精确徙动量 `R0/D(fa) - R0/Dref`，方位滤波随距离变化。 |
| `SarSimulationLearning/CSA.m` | 调用修正后的 CS：变标使用快时间 `tr`，距离 FFT 使用变标结果 `S1`，剩余相位包含距离差的平方，最后使用方位逆 FFT。 |
| `SarSimulationLearning/wkA.m` | 调用修正后的 wK：启用逐多普勒中心化的 Stolt 重采样，并补偿剩余方位相位，避免固定 Stolt 网格在斜视下移出原始距离频谱。 |
| `SarSimulationLearning/sar_sinc_interp.m`（新增） | 提供 16 点 Hann 加窗 sinc 插值；以整数位置计算正确的分数偏移，两端越界零填充。 |
| `SarSimulationLearning/ftx.m`、`iftx.m`、`fty.m`、`ifty.m` | 显式指定变换维度，使用匹配的中心化 FFT/IFFT，支持奇数长度和非方阵。 |
| `SarSimulationLearning/test_squint.m`（新增） | 无需绘图的回归检查：精确几何回波、聚焦位置、−3 dB 宽度、有限值、FFT 往返和插值边界。 |

### 坐标与适用范围

三种算法统一输出波束中心斜距和目标沿轨坐标。`Target(:,2)` 仍表示**最近斜距**，因此斜视时应把波束中心斜距乘以 `cos(theta_rc)` 后填写。距离采样间隔 `c/(2*Fr)` 与距离分辨率 `c/(2*abs(Kr)*Tr)` 不同。

本实现采用匀速直线、固定斜视、窄带条带模型；RD/CS 的 SRC 使用参考距离和二阶距离频率展开。大带宽、大范围或更高斜视角仍需评估高阶耦合。采样必须覆盖回波脉冲与照射时间，去中心后的多普勒带宽必须小于 PRF；截断或混叠不能通过成像滤波恢复。

公式背景：[Davidson、Cumming、Ito，A Chirp Scaling Approach for Processing Squint Mode SAR Data (1996)](https://sar.ece.ubc.ca/papers/Davidson%2BCumming%2BIto_1996.pdf)。


### 验证记录

使用 Python 3.13 / NumPy 2.2.4 对与 MATLAB 实现等价的数值公式进行复核：

- 斜视角：0°、±5°、±15°、±30°。
- 数据尺寸：256×256、257×320；每组包含中心目标及两个偏心目标。
- 距离调频率：`+0.25e12` 和 `-0.25e12` Hz/s。
- 三种算法共 252 个聚焦用例通过；离散峰值位置误差为 0 像素，距离和方位的 −3 dB 宽度均为 1 个采样像素。
- FFT 往返和插值边界检查通过。

以上结果是 NumPy 等价公式复核，当前环境未安装 MATLAB/Octave，尚未直接执行 MATLAB 回归脚本。可在 MATLAB 中进入 `SarSimulationLearning`，运行 `results = test_squint` 进行原生验证。验证结论仅覆盖上述参数与目标配置。
