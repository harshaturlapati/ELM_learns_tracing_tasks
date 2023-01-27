function [farthest_point, C_pstar_farthest] = find_extrema(SET_Pt,GEO)
%% Reparamaterisation
SET_pstar = []; SET_err = [];

Pt = SET_Pt(:,1);
pstar = fminbnd(@(p) GEO.dist2(p,Pt), 0, 2*pi);   %% NOTE [0 pi] is the upper lobe
 
for Pt = SET_Pt    %% looping through colums
    pstar = fminbnd(@(p) GEO.dist2(p,Pt), pstar-pi/4, pstar+pi/4);
    %pstar-pi/2, pstar+pi/2);   %% NOTE [0 pi] is the upper lobe
    SET_pstar = [SET_pstar, pstar];
    SET_err = [SET_err, sqrt( GEO.dist2(pstar, Pt))];
end

farthest_point_idx = find(SET_err==max(SET_err));

farthest_point = SET_Pt(:,farthest_point_idx);
C_pstar_farthest = GEO.C(SET_pstar(farthest_point_idx));

end