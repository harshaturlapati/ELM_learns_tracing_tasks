function plot_curve(C,curvature_threshold)
my2Dplot = @(X,msg) plot(X(1,:), X(2,:), msg) ;
mynorm = @(V) sqrt(V(1,:).^2 + V(2,:).^2);
figure % illustrating the figure of 8 and the areas with high curvature
%subplot(1,2,1)
    my2Dplot(C(pi/4 * (0:7)), 'ob'), hold on,
    my2Dplot(C([0:0.01:2*pi]), 'r'), grid on, axis equal, % axis(.15 *[-1 1 -1 1])
    xlabel('x-axis [m]'), ylabel('y-axis [m]')
    p = -2*pi:0.001:2*pi;
    [tmp1, tmp2, C_curv] = f_powerlaw_Zago_v4 (C(p), p, 0.001, nan);

    range_p = logical((p>0)'.*(C_curv>curvature_threshold ));
    my2Dplot(C(p(range_p)),'.b')
    title('High curvature regions highlghted')
%subplot(1,2,2)
figure
    plot(p(p>0), C_curv(p>0), 'r'), grid on, hold on, xlim([0,2.01*pi])
    plot(p(range_p), C_curv(range_p), '.b')
    xlabel('p [rad]'), ylabel('curvature [1/m]')
    set(gca, 'XTick', pi/4 * (0:8), ...
        'XTickLabel', {'0', '', '\pi/2', '', '\pi', '', '3\pi/2', '', '2\pi'})
    title('Curvature vs parameter p')
    
end