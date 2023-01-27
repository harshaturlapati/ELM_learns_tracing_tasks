function set_training_settings
clc; clear; close all;
rng('default');
%% Graphing options
graph_opt.drawVector = @(c,f,color) quiver(c(1),c(2),f(1),f(2),'Color',color);
graph_opt.my2DscatterXY = @(V,args) scatter(V(1,:),V(2,:),args);

%% Defining curvature threshold
curvature_threshold = 1;                %% select minorant for curvature

%% Defining training data length
half_training_length = 15; % seconds, i.e., 2*half_training_length number of seconds of training data will be chosen from the center of the overall dataset
no_samples = 20; % To randomly pick 20 initial pstar and Pt samples to perform generalisation from

%% Defining Zago filtering parameters
sampling_time = 0.01;   %% sec
cutoff_frequency = nan;   %% Hz
beta = 0.05;
%% Defining training curve
my2Dplot = @(X,msg) plot(X(1,:), X(2,:), msg) ;
mynorm = @(V) sqrt(V(1,:).^2 + V(2,:).^2);
syms p real;
origin = [0; 0];
l_x = 0.1*0.5;
l_y = 0.1;
C = @(p)[l_x*sin(2*p); l_y*sin(p)]; % curve definition
simplify(C(p))
figure % illustrating the figure of 8 and the areas with high curvature
subplot(1,2,1)
    my2Dplot(C(pi/4 * (0:7)), 'ob'), hold on,
    my2Dplot(C([0:0.01:2*pi]), 'r'), grid on, axis equal, % axis(.15 *[-1 1 -1 1])
    xlabel('x-axis [m]'), ylabel('y-axis [m]')
    p = -pi:0.001:2*pi;
    [tmp1, tmp2, C_curv] = f_powerlaw_Zago_v4 (C(p), p, 0.001, cutoff_frequency);

    range_p = logical((p>0)' .* (C_curv>curvature_threshold ));
    my2Dplot(C(p(range_p)),'.b')
    title('High curvature regions highlghted')
subplot(1,2,2)
    plot(p(p>0), C_curv(p>0), 'r'), grid on, hold on, xlim([0,2.01*pi])
    plot(p(range_p), C_curv(range_p), '.b')
    xlabel('p [rad]'), ylabel('curvature [1/m]')
    set(gca, 'XTick', pi/4 * (0:8), ...
        'XTickLabel', {'0', '', '\pi/2', '', '\pi', '', '3\pi/2', '', '2\pi'})
    title('Curvature vs parameter p')

    
p_range = [0:0.1:2*pi];
curve_pts = C(p_range);

C_min = min(curve_pts,[],2);
C_max = max(curve_pts,[],2);

x_size = C_max(1) - C_min(1);
y_size = C_max(2) - C_min(2);

percent = 1/100;
margin = 40*percent;  % 20% margin
axis_x_min = C_min(1) - margin*x_size;
axis_y_min = C_min(2) - margin*y_size;
axis_x_max = C_max(1) + margin*x_size;
axis_y_max = C_max(2) + margin*y_size;


save('training_settings.mat');
end