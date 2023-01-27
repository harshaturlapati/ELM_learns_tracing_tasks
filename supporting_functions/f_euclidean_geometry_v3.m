%% Setup Geometry
    % % REDUCE CODE
%     off lower; on expanddf; 
%     depend s,t; depend C, s;
%     eq1:= C_p = df(C,t,1);
%     eq2:= C_pp = df(C,t,2); 
%     eq3:= C_ppp = df(C,t,3);
%     SOL:= solve({eq1, eq2, eq3}, {df(C,s), df(C,s,2), df(C,s,3)});
%     SOL where{df(~x,t,2)=>mkid(x,!_pp), df(~x,t)=>mkid(x,!_p)};

function GEO = f_euclidean_geometry_v2(C)

syms p real
Pt = sym('Pt', [2 1],'real');

C_p = diff(C(p));
C_pp = diff(C_p);
C_ppp = diff(C_pp);

s_p = sqrt(dot(C_p, C_p));   %% Euclidean velocity
s_pp = diff(s_p);
s_ppp = diff(s_pp);

C_s   = C_p/s_p;
C_ss  = (- C_p*s_pp + C_pp*s_p)/s_p^3;
C_sss = (3*C_p*s_pp^2 - 3*C_pp*s_p*s_pp + C_ppp*s_p^2)/s_p^5;

mu = sqrt(dot(C_ss,C_ss));   %% Euclidean Curvature CHECK SIGN
dist2 = (Pt-C)' * (Pt-C) ;    %% Euclidean squared Distance
grad = -C_p' * (Pt-C) ;       %% Euclidean Distance gradient

JJ = [0 -1; 1 0];       
% technically , the normalization factor is not needed since |C_s|=1
%% Put MF inside GEO - done in v2
MF = [C_s, JJ*C_s] / sqrt(C_s' * C_s);  
GEO.MF = matlabFunction(MF, 'Vars',{p});
GEO.C = C;
GEO.C_s = matlabFunction(C_s, 'Vars',{p});
GEO.C_ss = matlabFunction(C_ss, 'Vars',{p});
GEO.mu = matlabFunction(mu, 'Vars',{p});
GEO.dist2 = matlabFunction(dist2, 'Vars',{p, Pt});
GEO.grad = matlabFunction(grad, 'Vars',{p, Pt});



range = [0 2*pi];

% figure
% subplot(4,2,[2,4,6,8])
%     fplot(@(p) [1 0]*GEO.C(p), @(p) [0 1]*GEO.C(p), range), hold on
%     fplot(@(p) [1 0]*GEO.C(p), @(p) [0 1]*GEO.C(p), range/8, 'LineWidth',2)
% %     plot([1 0]*GEO.C(0), [0 1]*GEO.C(0), 'ro')
% %     plot([1 0]*GEO.C(pi/4), [0 1]*GEO.C(pi/4), 'rx')
%     legend('0<p<2 \pi', '0< p< \pi/4')
%     grid on, axis equal, xlabel('x [m]'), ylabel('y [m]')
% % fplot3(@(p) [1 0]*GEO.C(p), @(p) [0 1]*GEO.C(p), p, [0 2*pi], 'r')
% subplot(4,2,1)
%     fplot(@(p) [1 0]*GEO.C(p), range), hold on
%     fplot(@(p) [0 1]*GEO.C(p), range), grid on, ylabel('C')
% subplot(4,2,3)
%     fplot(@(p) [1 0]*GEO.C_s(p), range), hold on
%     fplot(@(p) [0 1]*GEO.C_s(p), range), grid on, ylabel('C_s')
% subplot(4,2,5)
%     fplot(@(p) [1 0]*GEO.C_ss(p), range), hold on
%     fplot(@(p) [0 1]*GEO.C_ss(p), range), grid on, ylabel('C_{ss}')
% subplot(4,2,7)
%     fplot(GEO.mu, range), grid on, ylabel('\mu')
%     xlabel('curve parameter')

end