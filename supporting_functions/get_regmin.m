function ridgeLambda = get_regmin(X_train, Z_wash)
[N, T_ridge_regression] = size(X_train);
 x1 = X_train;
 y1 = Z_wash;
 x = X_train - repmat(mean(X_train,2), 1, T_ridge_regression);
 y = Z_wash - mean(Z_wash,2);
 x = x./repmat(std(x1, [], 2), 1, T_ridge_regression);
 [U, s, V] = csvd(x');
 disp('size of U is ')
 size(U)
 disp('size of s is ') 
 size(s)
 disp('size of V is ') 
 size(V)
 [reg_min, G, reg_param] = gcvHansen(U, s, y', 'Tikh'); % performs cross validation, while respecting the temporal properties
 display(reg_min);
 
ridgeLambda = reg_min;
end