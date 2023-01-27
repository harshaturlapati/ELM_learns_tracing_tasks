
function [time, speed, curvature, A, C,Cd, Cdd] = f_powerlaw_Zago_v4(Pt,time,sampling_time, cutoff_frequency)

% % pathname = [pwd '/datafiles/'];     %% DC: changed slashes
% % filename = 'harsha.fig8.hman.mat';
% % work = load( strcat(pathname,filename) );
% % 
% % SET_Pt = work.data(6:7,:); % Positions
% % SET_Vt = work.data(10:11,:); % Velocities
% % time =  work.data(1,:)'; % Time
% % sampling_time=0.01
% % Pt = SET_Pt';

if size(Pt,1)>size(Pt,2)
    Pt = Pt';
end
switch nargin
    case 1
        time=ones(1, length(Pt));
        sampling_time = 1;
        cutoff_frequency = 0.1/sampling_time;    %% 
    case 2
        time = time(:)';
        sampling_time = min(diff(time));
        cutoff_frequency = 5;    %% Hz as in Zago
    case 3
        time = time(:)';
        cutoff_frequency = 5;    %% Hz as in Zago
end

%% Uniform time re-sampling
time_uniform = time(1):sampling_time:time(end);
Pt = (interp1(time,Pt',time_uniform,'spline'))' ;
time = time_uniform;


%% Butterworth filtering - 5 Hz Cutoff, second-order, zero-phase lag
if not( isnan(cutoff_frequency) )        %% skip filtering if cutoff_freq = NaN
    fc = cutoff_frequency;
    fs = 1/sampling_time;

    [a,b] = butter(2,2*pi*fc, 's');
    SYS = ss(tf(a,b));
    Out1 = lsim(SYS,Pt(1,:),time, [0 Pt(1,1)]/SYS.C(2));
    Out2 = lsim(SYS,Pt(2,:),time, [0 Pt(2,1)]/SYS.C(2));
        % note: lsim accepts initial conditions only if SYS is a
        % state-space model and 
        % the initial conditions is actually the Y0 = C*X0 variable
    Pt = [Out1, Out2]';
end


%% from (filtered) time-samples to splines 

X = spline(time,Pt(1,:));
Y = spline(time,Pt(2,:));

%% compute (spline) higher derivatives

Xd = fnder(X,1);
Yd = fnder(Y,1);

Xdd = fnder(X,2);
Ydd = fnder(Y,2);

%% resample 
X = ppval(X,time');
Y = ppval(Y,time');
Xd = ppval(Xd,time');
Yd = ppval(Yd,time');
Xdd = ppval(Xdd,time');
Ydd = ppval(Ydd,time');

%% compute powerlaw variabels as in Zago et al. (2019)
D = abs( Xd .* Ydd - Xdd .* Yd );
speed = sqrt(Xd.^2 + Yd.^2);
curvature =  D ./ speed.^3;             %% Zago eq.(5)
A = D.^(1/3) .* curvature.^(2/3);       %% Zago eq.(6)

C = [X'; Y'];
Cd = [Xd'; Yd'];
Cdd = [Xdd'; Ydd'];
return
