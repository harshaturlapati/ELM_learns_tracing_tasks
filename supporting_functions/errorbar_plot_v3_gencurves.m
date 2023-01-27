function h_ppl = errorbar_plot_v3_gencurves(total_subjects,pathname,real_K,real_V,real_K_CI,real_V_CI,ELM_K,ELM_V,ELM_K_CI,ELM_V_CI,curve_idx)

h_ppl = figure;
%real_K = [11.2513 11.2823 14.17637 10.9074 15.6335];
%real_V = [0.1446 0.2214 0.2169 0.2059 0.0771];
real_V_neg = real_V_CI{curve_idx};
real_V_pos = real_V_neg;
real_K_neg = real_K_CI{curve_idx};
real_K_pos = real_K_neg;
e_real = errorbar(real_K{curve_idx},real_V{curve_idx},real_V_neg,real_V_pos,real_K_neg,real_K_pos,'o'); hold on;
e_real.Marker = 'square';
e_real.MarkerSize = 10;
e_real.Color = 'red';
e_real.CapSize = 15;
e_real.MarkerEdgeColor = 'red';
e_real.MarkerFaceColor = 'red';



%ELM_K = [21.2717 13.8305 16.8864 20.4570 8.4609];
%ELM_V = [0.1296 0.1920 0.1964 0.1712 0.1060];
ELM_V_neg = ELM_V_CI{curve_idx};
ELM_V_pos = ELM_V_neg;
ELM_K_neg = ELM_K_CI{curve_idx};
ELM_K_pos = ELM_K_neg;
e_elm = errorbar(ELM_K{curve_idx},ELM_V{curve_idx},ELM_V_neg,ELM_V_pos,ELM_K_neg,ELM_K_pos,'o'); hold on;
e_elm.Marker = 'diamond';
e_elm.MarkerSize = 10;
e_elm.Color = 'blue';
e_elm.CapSize = 15;
e_elm.MarkerEdgeColor = 'blue';
e_elm.MarkerFaceColor = 'blue';

subject_order = [1 10 2 3 4 5 6 7 8 9];
for subject_idx = 1 : length(real_V{curve_idx})
    %text(real_K(subject_idx)+0.5*real_K_neg(subject_idx),real_V(subject_idx)+0.5*real_V_neg(subject_idx),sprintf('$S_%d^{real}$',subject_idx),'Interpreter','Latex');
    %text(ELM_K(subject_idx)+0.5*ELM_K_neg(subject_idx),ELM_V(subject_idx)+0.5*ELM_V_neg(subject_idx),sprintf('$S_%d^{ELM}$',subject_idx),'Interpreter','Latex');
    %text_h = text(0.5*(real_K(subject_idx)+ELM_K(subject_idx)),0.5*(real_V(subject_idx)+ELM_V(subject_idx)),sprintf('$S_%d$',subject_idx),'Interpreter','Latex','fontsize',12);
    if(~mod(subject_idx,2))

    text_h = text(real_K{curve_idx}(subject_idx)+2.5,real_V{curve_idx}(subject_idx),sprintf('$S_{%d}$',subject_order(subject_idx)),'Interpreter','Latex','fontsize',20);
    else
    text_h = text(ELM_K{curve_idx}(subject_idx)+2.5,ELM_V{curve_idx}(subject_idx),sprintf('$S_{%d}$',subject_order(subject_idx)),'Interpreter','Latex','fontsize',20);
    end
    plot_h = plot([real_K{curve_idx}(subject_idx) ELM_K{curve_idx}(subject_idx)],[real_V{curve_idx}(subject_idx) ELM_V{curve_idx}(subject_idx)],'Color',rand(1,3));
end
    legend([e_real, e_elm],{'Real','ELM'})
    ylabel('speed (m/s)'),
    xlabel('Curvature (1/m)')
    title(sprintf('Power - law comparison for curve %d',curve_idx))
    
    
    table_allS = [];
    for subject_idx = 1 : length(real_V{curve_idx})
        % || v_c^real v_CI^real k_c^real k_CI^real | v_c^ELM v_CI^ELM k_c^ELM k_CI^ELM ||
        subject_chars = [real_V{curve_idx}(subject_idx) real_V_CI{curve_idx}(subject_idx) real_K{curve_idx}(subject_idx) real_K_CI{curve_idx}(subject_idx) ELM_V{curve_idx}(subject_idx) ELM_V_CI{curve_idx}(subject_idx) ELM_K{curve_idx}(subject_idx) ELM_K_CI{curve_idx}(subject_idx)];
        table_allS = [table_allS; subject_chars];
    end
    digits(2);
    latex_table = latex(sym(vpa(table_allS)));
    allS_table_filename = strcat(pathname,'all_S.table.mat');
    save(allS_table_filename,'latex_table');
    
end