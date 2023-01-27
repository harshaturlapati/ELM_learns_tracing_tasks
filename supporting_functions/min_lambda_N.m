function [er_os_total] = min_lambda_N(arg_goal, ELM, N_train, Z_train_norm, Y_train, pstar_all, Pt_all, GEO, Z_BUILDER, sampling_time)
% Xt contains the last value of the washout set
    c = arg_goal(1);
    l = arg_goal(2);
    ELM.Dyn = @(Xt, Zt)radbas(ELM.A * Xt + c * ELM.C * Zt + ELM.zeta);
    j = 1;
    Xt = zeros(ELM.Nx, 1);
    X_train = zeros(ELM.Nx, N_train);
    for Zt = Z_train_norm
        Xt = ELM.Dyn(Xt, Zt);
        X_train(:, j) = Xt;
        j = j + 1;
    end


    Nx = ELM.Nx;
    if N_train > 8000
        hhh = 100;
        coeff = 1/2;
        inds = (10:hhh:floor((N_train-100)*1/4));
        len_int = length(inds);
    else
        %hhh = 10;
        coeff = 3/4;
        len_int = min([floor((N_train-100)*1/4),400]);
        rng('default');
        inds = randperm(floor((N_train-100)*1/4), len_int);

    end


    er_is_temp=[];
    er_os_temp = [];
    for k = 1:len_int
        N_cur = inds(k):inds(k)+floor((N_train-100)*coeff)-50;
        W4 = (X_train(:,N_cur)*X_train(:,N_cur)'+ l * eye(Nx))\(X_train(:,N_cur) * Y_train(:,N_cur)');
        b = mean(Y_train(:,N_cur),2) - W4' * mean(X_train(:,N_cur), 2);
        Y_train_hat = W4' * X_train(:,N_cur) + b;
        %X_t_cur_prev = X_train(:,N_cur);
        Y_t_cur(:,1) = Y_train_hat(:, end);
        pstar = pstar_all(N_cur(end));
        %Pt = res(ll).Y_t_cur(:,1);
        Pt = Pt_all(:, N_cur(end));
        % 10-steps ahead forecasting in autonomous regime
        for h = 1:100
            pstar = fminbnd (@(p) GEO.dist2(p,Pt), pstar-pi/8 , pstar+pi/8);
            Zt = Z_BUILDER(Pt, mod(pstar,2*pi), GEO);
            X_t_cur(:,h) = radbas(c* ELM.C * Zt + ELM.zeta );
            Y_t_cur(:,h) = W4' * X_t_cur(:,h) + b;
            %X_t_cur_prev = X_t_cur(:,h);
            Pt = Pt + GEO.MF(mod(pstar,2*pi))*sampling_time*Y_t_cur(:,h);
        end
        er_is_temp(k) = mean((mean((Y_train_hat-Y_train(:,N_cur)).^2,2)));
        Y_test_cur = Y_train(:,N_cur(end)+1:N_cur(end)+100);
        Y_test_hat = Y_t_cur;
        er_os_temp(k) = mean((mean((Y_test_hat-Y_test_cur).^2,2)));
    end
    er_is_total = mean(er_is_temp,2);
    er_os_total = mean(er_os_temp,2);


end

