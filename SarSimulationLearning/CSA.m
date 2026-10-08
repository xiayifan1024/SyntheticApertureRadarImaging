% Chirp Scaling：支持固定斜视、匀速直线条带 SAR。
[S3,imageAxes] = sar_focus(s,tr,ta,f0,Kr,Vr,c,theta_rc,'CS');
figure(5)
imagesc(imageAxes.range,imageAxes.azimuth,abs(S3))
xlabel('波束中心斜距 (m)'); ylabel('沿轨方位 (m)'); title('Chirp Scaling');
axis xy
