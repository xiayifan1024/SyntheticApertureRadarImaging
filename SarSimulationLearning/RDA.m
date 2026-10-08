% Range Doppler：支持固定斜视、匀速直线条带 SAR。
[s_rmf,imageAxes] = sar_focus(s,tr,ta,f0,Kr,Vr,c,theta_rc,'RD');
figure(4)
imagesc(imageAxes.range,imageAxes.azimuth,abs(s_rmf))
xlabel('波束中心斜距 (m)'); ylabel('沿轨方位 (m)'); title('Range Doppler');
axis xy
