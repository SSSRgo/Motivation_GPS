function [fitresult, gof,f] = createFit2(trayontime, PeakInternalTrayon)
%CREATEFIT2(TRAYONTIME,PEAKINTERNALTRAYON)
%  创建一个拟合。
%
%  要进行 '无标题拟合 1' 拟合的数据:
%      X 输入: trayontime
%      Y 输出: PeakInternalTrayon
%  输出:
%      fitresult: 表示拟合的拟合对象。
%      gof: 带有拟合优度信息的结构体。
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 19-Aug-2023 19:22:20 自动生成


%% 拟合: '无标题拟合 1'。
[xData, yData] = prepareCurveData( trayontime, PeakInternalTrayon );

% 设置 fittype 和选项。
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.000224470822980021;

% 对数据进行模型拟合。
[fitresult, gof] = fit( xData, yData, ft, opts );

% 绘制数据拟合图。
f=figure( 'Name', 'Fitting Curve' );
h = plot( fitresult, xData, yData );
legend( h, 'PeakInternalTrayon vs. trayontime', 'Fitting Curve', 'Location', 'NorthEast', 'Interpreter', 'none' );
% 为坐标区加标签
xlabel( 'trayontime', 'Interpreter', 'none' );
% ylabel( 'PeakInternalTrayon', 'Interpreter', 'none' );
grid on


