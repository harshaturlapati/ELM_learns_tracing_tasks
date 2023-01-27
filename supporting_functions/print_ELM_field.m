function [h_fig1_elm_gen_fig8,h_elm_field] = print_ELM_field(GEO_gen,C_gen,h_fig1_elm_gen_fig8,Z_BUILDER,ELM,mean_Z,std_Z,mean_Y,std_Y,Xt,sampling_time,beta,my2Dplot)

p_range = [0:0.1:2*pi];
            curve_pts = GEO_gen.C(p_range);
            
            C_min = min(curve_pts,[],2);
            C_max = max(curve_pts,[],2);
            
            x_size = C_max(1) - C_min(1);
            y_size = C_max(2) - C_min(2);
            
            percent = 1/100;
            margin = 50*percent;  % 20% margin
            axis_x_min = C_min(1) - margin*x_size;
            axis_y_min = C_min(2) - margin*y_size;
            axis_x_max = C_max(1) + margin*x_size;
            axis_y_max = C_max(2) + margin*y_size;
            
            set(h_fig1_elm_gen_fig8, 'Position', get(0, 'Screensize')); 
            size_fig = 2;
            field_origin = C_gen(0);
            myXgrid = axis_x_min: 0.02 :axis_x_max;
            myYgrid = axis_y_min: 0.02 :axis_y_max;
            [x_coord,y_coord] = meshgrid(myXgrid, myYgrid);
            Pt_grid = [x_coord(:) y_coord(:)]';
            axis equal, grid on, hold on
            p_0 = 0;
            
               for Pt = Pt_grid
                pstar = fminbnd (@(p) GEO_gen.dist2(p,Pt), 0 , 2*pi);
                Zt = Z_BUILDER(Pt, mod(pstar,2*pi), GEO_gen);
                Zt_norm = (Zt - mean_Z)./std_Z;
                Xt = ELM.Dyn (Xt, Zt_norm);
                Yt_norm = ELM.Eval(Xt) ;
                Yt = (Yt_norm.*std_Y)+mean_Y;
                Ptplusone = Pt + GEO_gen.MF(mod(pstar,2*pi))*sampling_time*Yt + beta*(GEO_gen.C(pstar)-Pt);
                %+ beta*(GEO_gen.C(pstar)-Pt);
                h_elm_field = my2Dplot([Pt, Ptplusone], 'b'); hold on;
                my2Dplot([Pt, Pt], '.b'); hold on; %% this plots a thick mark at the base of vector 
            end
end
        