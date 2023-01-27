function h_ppl = errorbar_plot_v2(total_subjects,pathname,real_K,real_V,real_K_CI,real_V_CI,ELM_K,ELM_V,ELM_K_CI,ELM_V_CI)

h_ppl = figure;
real_V_neg = real_V_CI;
real_V_pos = real_V_neg;
real_K_neg = real_K_CI;
real_K_pos = real_K_neg;
e_real = errorbar(real_K,real_V,real_V_neg,real_V_pos,real_K_neg,real_K_pos,'o'); hold on;
e_real.Marker = 'square';
e_real.MarkerSize = 10;
e_real.Color = 'red';
e_real.CapSize = 15;
e_real.MarkerEdgeColor = 'red';
e_real.MarkerFaceColor = 'red';


ELM_V_neg = 0*ELM_V_CI;
ELM_V_pos = ELM_V_neg;
ELM_K_neg = 0*ELM_K_CI;
ELM_K_pos = ELM_K_neg;
e_elm = plot(ELM_K,ELM_V,'d','color','b','markerfacecolor','b','markersize',10);
% e_elm = errorbar(ELM_K,ELM_V,ELM_V_neg,ELM_V_pos,ELM_K_neg,ELM_K_pos,'o'); hold on;
% e_elm.Marker = 'diamond';
% e_elm.MarkerSize = 10;
% e_elm.Color = 'blue';
% e_elm.CapSize = 15;
% e_elm.MarkerEdgeColor = 'blue';
% e_elm.MarkerFaceColor = 'blue';

subject_order = [1 10 2 3 4 5 6 7 8 9];
for subject_idx = 1 : length(real_V)
    if(~mod(subject_idx,2))

    text_h = text(real_K(subject_idx)*1.05,real_V(subject_idx)*1.05,sprintf('$S_{%d}$',subject_order(subject_idx)),'Interpreter','Latex','fontsize',20);
    else
    text_h = text(ELM_K(subject_idx)*1.05,ELM_V(subject_idx)*1.05,sprintf('$S_{%d}$',subject_order(subject_idx)),'Interpreter','Latex','fontsize',20);
    end
    plot_h = plot([real_K(subject_idx) ELM_K(subject_idx)],[real_V(subject_idx) ELM_V(subject_idx)],'Color',rand(1,3));
end
    legend([e_real, e_elm],{'Real','ELM'})
    ylabel('speed (m/s)'),
    xlabel('Curvature (1/m)')
    
    table_allS = [];
    for subject_idx = 1 : length(real_V)
        subject_chars = [real_V(subject_idx) real_V_CI(subject_idx) real_K(subject_idx) real_K_CI(subject_idx) ELM_V(subject_idx) ELM_V_CI(subject_idx) ELM_K(subject_idx) ELM_K_CI(subject_idx)];
        table_allS = [table_allS; subject_chars];
    end
    digits(2);
    latex_table = latex(sym(vpa(table_allS)));
    allS_table_filename = strcat(pathname,'all_S.table.mat');
    save(allS_table_filename,'latex_table');
    
end