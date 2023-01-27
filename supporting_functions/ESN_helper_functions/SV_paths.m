function [S_terminal,S_next,Sigma_next,F,Paths,Sigma2,Sigma2p,Sigma2k,Sigma2u]...
    =SV_paths(parameters,So,T,N,Tnext,method,bo,epsilon,epsilonw)
% 
% Creates paths out of the GJR model with NIG innovations.
% Requires the use of randraw and Werner's NIG package
% Inputs: 
% So initial price, T time to maturity, N number of paths.
% kernel: measure with respect to which the simulations are carried out. 
% The options are:
%     'historical'
%     'mcmm_ems': mean correcting martingale measure and the resulting paths are
%     corrected using Duan's Empirical Martingale Simulation technique.
% sigma2o (optional) variance at T=0, epsilono (optional) innovation at T=0. 
% If sigma2o and
% epsilono are not specified the following default values are used:
% epsilono=0 and sigma2o the marginal variance of the model
% 
% Outputs: 
% S_terminal: vector containing the terminal values (t=T) of the N paths.
% S_next: vector containing the values at t=1 of the N paths.
% Sigma_next: vector containing the values of the volatility sigma at t=1 of the N paths.
% F: vector containing the values of the martingale correction factor 
% (same notation as in the paper) corresponding to the N paths.
% Paths: Optional, NxT matrix contaning all the paths.
% Sigma2: Optional, NxT matrix contaning the conditional variances of all the paths.
% Epsilon: Optional, NxT matrix contaning the residuals of all the paths.


r=cell2mat(parameters(1));
gamma=cell2mat(parameters(2));
phi=cell2mat(parameters(3));
sigmaw=cell2mat(parameters(4));
nu=0;
alphasv=gamma/(2*(1-phi));
muxi=-.63518;
H=pi^2/8;
sigmaeta2=sigmaw^2/4;
% tic

if nargin<7
%   Default value for bo
    bo=gamma/(1-phi)*ones(3,1);
    epsilon=randn(N,T);
    epsilonw=randn(N,T);
elseif nargin<8
    epsilon=randn(N,T);
    epsilonw=randn(N,T);
elseif isempty(bo)==1
    bo=gamma/(1-phi)*ones(3,1);
end 
w=sigmaw*epsilonw;
S_terminal=zeros(N,1);
S_next=zeros(N,1);
Sigma_next=zeros(N,3);
Paths=zeros(N,T);
Sigma2p=zeros(N,T);
Sigma2u=zeros(N,T);
Sigma2=zeros(N,T);
Sigma2k=zeros(N,T);
F=zeros(N,2);
if strcmp(method,'h-likelihood')==1
    for i=1:N 
        b=bo(1);
        bu=bo(2);
        s=So;
        np=1;
        nminp=1;
        for j=1:Tnext    
            bp=gamma+phi*bu;
            sigmap=exp(.5*bp);
            b=gamma+phi*b+w(i,j);
            sigma=exp(.5*b);
            ys=sigma*epsilon(i,j);
            bu=fzero(@(bunext) vola_update_zero(bunext,bu,ys,sigmaw^2,phi,gamma),bp);
            sigmau=sqrt(exp(bu));
            s=s*exp(r+ys);
%             np=np*n_mcmm(nu,epsilon(i,j),sigma,sigmap);
%             nminp=nminp*n_minimal(r,epsilon(i,j),sigma,sigmap);
        end
        Sigma_next(i,1)=sigma;
        Sigma_next(i,2)=sigmap;
        Sigma_next(i,3)=sigmau;
        S_next(i,1)=s;
        for j=Tnext+1:T    
            bp=gamma+phi*bu;
            sigmap=exp(.5*bp);
            b=gamma+phi*b+w(i,j);
            sigma=exp(.5*b);
            ys=sigma*epsilon(i,j);
            bu=fzero(@(bunext) vola_update_zero(bunext,bu,ys,sigmaw^2,phi,gamma),bp);
            sigmau=sqrt(exp(bu));
            s=s*exp(r+ys);
%             np=np*n_mcmm(nu,epsilon(i,j),sigma,sigmap);
%             nminp=nminp*n_minimal(r,epsilon(i,j),sigma,sigmap);
        end
       S_terminal(i)=s;
%        F(i,:)=[np nminp];
    end
elseif strcmp(method,'kalman')==1
    for i=1:N 
        b=bo(1);
        s=So;
        Fk=H+(sigmaeta2/(1-phi^2));
        nk=1;
        nmink=1;
        for j=1:Tnext    
            b=gamma+phi*b+w(i,j);
            sigma=exp(.5*b);
            ys=sigma*epsilon(i,j);
            s=s*exp(r+ys);
            l=log(abs(ys));
            if j==1
                lp=.5*bo(1)+muxi;
            else
                lp=(1-phi)*(alphasv+muxi)+phi*l-(phi*H*(l-lp)/Fk);
            end
            Fk=H+(phi^2)*(Fk-H)*H*Fk^(-1)+sigmaeta2;
            sigmak=exp(lp-muxi);
%             nk=nk*n_mcmm(nu,epsilon(i,j),sigma,sigmak);
%             nmink=nmink*n_minimal(r,epsilon(i,j),sigma,sigmak);
        end
        Sigma_next(i,1)=sigma;
        Sigma_next(i,2)=sigmak;
        S_next(i,1)=s;
        for j=Tnext+1:T    
            b=gamma+phi*b+w(i,j);
            sigma=exp(.5*b);
            ys=sigma*epsilon(i,j);
            s=s*exp(r+ys);
            l=log(abs(ys));
            lp=(1-phi)*(alphasv+muxi)+phi*l-(phi*H*(l-lp)/Fk);
            Fk=H+(phi^2)*(Fk-H)*H*Fk^(-1)+sigmaeta2;
            sigmak=exp(lp-muxi);
%             nk=nk*n_mcmm(nu,epsilon(i,j),sigma,sigmak);
%             nmink=nmink*n_minimal(r,epsilon(i,j),sigma,sigmak);
        end
       S_terminal(i)=s;
%        F(i,:)=[nk nmink];
    end
elseif strcmp(method,'all')==1  
    Sigma_next=zeros(N,4);
    F=zeros(N,4);
    w=sigmaw*epsilonw;
    for i=1:N 
        b=bo(1);
        bu=bo(2);
        s=So;
        Fk=H+(sigmaeta2/(1-phi^2));
        np=1;
        nk=1;
        nminp=1;
        nmink=1;
        for j=1:Tnext    
            bp=gamma+phi*bu;
            sigmap=exp(.5*bp);
            b=gamma+phi*b+w(i,j);
            sigma=exp(.5*b);
            ys=sigma*epsilon(i,j);
            bu=fzero(@(bunext) vola_update_zero(bunext,bu,ys,sigmaw^2,phi,gamma),bp);
            sigmau=sqrt(exp(bu));
            s=s*exp(r+ys);
            l=log(abs(ys));
            if j==1
                lp=.5*bo(3)+muxi;
            else
                lp=(1-phi)*(alphasv+muxi)+phi*l-(phi*H*(l-lp)/Fk);
            end
            Fk=H+(phi^2)*(Fk-H)*H*Fk^(-1)+sigmaeta2;
            sigmak=exp(lp-muxi);
%             np=np*n_mcmm(nu,epsilon(i,j),sigma,sigmap);
%             nk=nk*n_mcmm(nu,epsilon(i,j),sigma,sigmak);
%             nminp=nminp*n_minimal(r,epsilon(i,j),sigma,sigmap);
%             nmink=nmink*n_minimal(r,epsilon(i,j),sigma,sigmak);
            Paths(i,j)=s;
            Sigma2p(i,j)=sigmap^2;
            Sigma2u(i,j)=sigmau^2;
            Sigma2(i,j)=sigma^2;
            Sigma2k(i,j)=sigmak^2;
        end
        Sigma_next(i,1)=sigma;
        Sigma_next(i,2)=sigmap;
        Sigma_next(i,3)=sigmak;
        Sigma_next(i,4)=sigmau;
        S_next(i,1)=s;
        for j=Tnext+1:T    
            bp=gamma+phi*bu;
            sigmap=exp(.5*bp);
            b=gamma+phi*b+w(i,j);
            sigma=exp(.5*b);
            ys=sigma*epsilon(i,j);
            bu=fzero(@(bunext) vola_update_zero(bunext,bu,ys,sigmaw^2,phi,gamma),bp);
            sigmau=sqrt(exp(bu));
            s=s*exp(r+ys);
            l=log(abs(ys));
            lp=(1-phi)*(alphasv+muxi)+phi*l-(phi*H*(l-lp)/Fk);
            Fk=H+(phi^2)*(Fk-H)*H*Fk^(-1)+sigmaeta2;
            sigmak=exp(lp-muxi);
%             np=np*n_mcmm(nu,epsilon(i,j),sigma,sigmap);
%             nk=nk*n_mcmm(nu,epsilon(i,j),sigma,sigmak);
%             nminp=nminp*n_minimal(r,epsilon(i,j),sigma,sigmap);
%             nmink=nmink*n_minimal(r,epsilon(i,j),sigma,sigmak);
            Paths(i,j)=s;
            Sigma2p(i,j)=sigmap^2;
            Sigma2u(i,j)=sigmau^2;
            Sigma2(i,j)=sigma^2;
            Sigma2k(i,j)=sigmak^2;
        end
        S_terminal(i,1)=s;
%         F(i,:)=[np nk nminp nmink];
    end
end

function value=vola_update_zero(b,bprevious,y,sigmaw2,phi,gamma)
value=-y^2*exp(-b)+1+(2/sigmaw2)*(b-gamma-phi*bprevious);
