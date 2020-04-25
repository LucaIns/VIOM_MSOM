% saveas(gca, 'newMMwDer3.eps','epsc');

clc, clear, clf, close all
rng(10) % seed for general simulation 
% NOTE:
% for the single simulation example: 
    % uncomment 2nd line frac_cont
    % n=1000 
    % monit=1 (and stop at the first iteration)

%% shortcuts

% if monit=1 each iteration is monitored (and we add a stop at each iteration)
% if monit=2 results over all iterations are monitored (without stopping)
% if monit=3 final results function of n are showed
% [TBA: if monit=4 final results function of contamination mechanism are showed]
monit = 3;
    % monit = 1;
% progress bar (it might become annoying)
progr = 1;

% adopted estimator
% please respect the order as the full vector est_nam_f (see draft, described also later)
% [TBA: general ordering]
% est_nam_f = {'opt', 'OLS', 'LTS', 'MM', 'FSR', 'FSRwj', 'FSRws'};
% est_nam = est_nam_f;
est_nam = {'opt', 'OLS', 'LMS', 'MM', 'FSR', 'FSRws'}; % 'FSRwj'
lest = length(est_nam);

%% Simulation design

% num. iterations (it has to be >1)
iter_tot = 500;
% range of samples sizes
n_tot = round(linspace(100, 1000, 19));
    % n_tot = 1000;
% signal-to-noise ratio
SNR_tot = 3; %[3, 4, 6]; % (SNR=5 is a fairly high Rsquared, approx 80%)
% total fraction of contamination
frac_cont = [0.25, 0, 0.5];
%     frac_cont = 0.5;
% fraction of contamination for MSOM (for VIOM it is the complement to cont_frac)
fracMSOM_vec = frac_cont/2;

% shifting parameters
shiftX = 3; 
shiftY= -3;
% inflation parameters (isometric) 
inflV = 10; % strictly > 1
                
% intercept (T/F)
intercept = 1;
% relevant predictors (different from 0)
s = 1 + intercept;
% irrelevant predictors
irr_p = 0;  
% total predictors
p = s + irr_p;

% parameters
beta_fixed = 2;
beta_true = ones(s, 1) * beta_fixed;
beta_true = [beta_true; zeros(irr_p, 1)];

% uncontaminated X (multivariate normal)
mu_X_u = zeros(p, 1);
sigma_X_u = 1;
rho_X_u = eye(p) * sigma_X_u;

%% initialization

% beta estimates
sol = nan(lest, p, iter_tot);
% sigma estimates (using WLS and non)
sols = nan(lest, 1, iter_tot);
% time
tim = nan(lest, iter_tot);
% time where each n is averaged over iterr
timn = nan(lest, length(n_tot));

% beta MSE=variance+bias^2 accross iterr iterations and for different sample size
solVB = nan(length(n_tot), lest, p);
solBB = nan(length(n_tot), lest, p);
% sigma MSE
solVS = nan(length(n_tot), lest, 1);
solBS = nan(length(n_tot), lest, 1);
            
% progress bar
if progr == 1
    hwait = waitbar(0,'Please wait...');
end

%% beginning

for SNR_i = 1:length(SNR_tot)
    
    for cont_i = 1:length(frac_cont)
        
        for n_i = 1:length(n_tot)
            
            % obs
            n = round(n_tot(n_i));
            
            % weights estimates
            solw = nan(lest, n, iter_tot);
            
            for iterr_i = 1:iter_tot
                
                % Data generation
                
                % normal X
                X_u = mvnrnd(mu_X_u, rho_X_u, n);
                % add intercept
                if intercept == 1
                    X_u(:, 1) = ones(n, 1);
                end
                % errors' std dev (to obtain the pre-specified SNR)
                sigma_e_u = sqrt(var(X_u * beta_true) / SNR_tot(SNR_i));
                % normal errors
                e_u = randn(n, 1) .* sigma_e_u;
                                
                % fraction of obs arising from a MSOM
                cont_fracMSOM_i = fracMSOM_vec(cont_i);
                % number of units arising from a MSOM
                n_MSOM = round(n * cont_fracMSOM_i);
                % indexes MSOM
                ind_MSOM = randperm(n, n_MSOM);
                % MSOM (also on X)
                X_c = X_u;
                X_c(ind_MSOM, intercept+1:end) = X_u(ind_MSOM, intercept+1:end) + ...
                    shiftX; % .* rand(round(n*cont_fracMSOM_i), 1);
                e_c = e_u;
                e_c(ind_MSOM) = e_u(ind_MSOM) + shiftY; % .* rand(round(n*cont_fracMSOM_i), 1);
                
                % fraction arising from a VIOM
                frac_VIOM = frac_cont(cont_i) - cont_fracMSOM_i;                
                % number of units arising from a VIOM
                n_VIOM = round(n * frac_VIOM);
                % indexes VIOM
                ind_VIOM = datasample(setdiff(1:n, ind_MSOM), n_VIOM, 'replace', false);
                % ind_VIOM = randperm(n, n_VIOM);
                %VIOM
                % inflVV = 1 + inflV; % rand(round(n*frac_VIOM), 1) .* inflV; % different var infl
                inflVV = inflV;
                e_c(ind_VIOM) = e_u(ind_VIOM) .* sqrt(inflVV);
                % contaminated model
                y_c = X_u * beta_true + e_c;
                
                % total contaminated units
                ind_cont = unique([ind_MSOM, ind_VIOM]);
                % clean obs indexes
                ind_keep = setdiff(1:n, ind_cont);
                
                %% Fitting
                               
                % use contaminated data
                X = X_c;
                y = y_c;
                % estimator counter
                est_i = 1;
                
                % oracle
                if any(strcmp('opt', est_nam))
                    ttic = tic;
                    w_opt = ones(n, 1);
                    w_opt(ind_VIOM) = 1 ./ inflVV; 
                    w_opt(ind_MSOM) = 0; % MSOM after VIOM (so in case of overlaps we set w=0)
                    W = w_opt .* eye(n);
                    solOpt = (X'*W*X)\X'*W*y;
                    sol(est_i, :, iterr_i) = solOpt';
                    sols(est_i, 1, iterr_i) = 1/(n-p) * (w_opt' * (y - X * solOpt).^2) / (sum(w_opt)/n);
                    solw(est_i, :, iterr_i) = w_opt;
                    tim(est_i, iterr_i) = toc(ttic);
                    est_i = est_i + 1;
                end
                
                
                % OLS
                % [TBA: it might be computed from FSR, but we need to time it]
                if any(strcmp('OLS', est_nam))
                    ttic = tic;
                    [solOLS, ~, ~, ~, statOLS] = regress(y, X);
                    sol(est_i, :, iterr_i) = solOLS';
                    solw(est_i, :, iterr_i) = ones(n,1);
                    sols(est_i, 1, iterr_i) =  1/(n-p) * sum((y - X * solOLS).^2);
                    tim(est_i, iterr_i) = toc(ttic);
                    est_i = est_i + 1;
                end
                
                % LXS
                if any(ismember({'LTS', 'LMS', 'MM', 'FSR'}, est_nam))
                    if any(strcmp(est_nam, 'LMS' ))
                        LXSk = 1;
                    else
                        LXSk = 2;
                    end
                    ttic = tic;
                    solLXS = LXS(y, X, 'lms', LXSk, 'h', floor(0.5*(n+p+1)), 'intercept', 0, ...
                        'msg', 0, 'plots', monit==1); 
                    if ismember('LTS', est_nam) || ismember('LMS', est_nam)
                        sol(est_i, :, iterr_i) = solLXS.beta';
                        solw(est_i, :, iterr_i) = solLXS.weights;
                        sols(est_i, 1, iterr_i) = 1/(n-p) * (solLXS.weights' * (y - X * solLXS.beta).^2) / (sum(solLXS.weights)/n);
                        %solLXS.scale.^2;
                        timLTS = toc(ttic);
                        tim(est_i, iterr_i) = timLTS;
                        est_i = est_i + 1;
                    end
                end
                
                % MM
                if any(strcmp('MM', est_nam))
                    ttic = tic;
                    InitialEst.beta=solLXS.beta;
                    InitialEst.scale=solLXS.scale;
                    [solMM]=MMreg(y,X,'intercept', 0, 'InitialEst', InitialEst, ...
                        'rhofunc', 'bisquare', 'eff', 0.85, 'plots', monit==1);% 
                    solMMest = solMM.beta;
                    sol(est_i, :, iterr_i) = solMMest';
                    solw(est_i, :, iterr_i) = solMM.weights;
                    sols(est_i, 1, iterr_i) = 1/(n-p) * (solMM.weights' * (y - X * solMM.beta).^2) / (sum(solMM.weights)/n);
                    % solMM.auxscale^2;
                    tim(est_i, iterr_i) = toc(ttic) + timLTS; 
                    est_i = est_i + 1;
                end
                
                % FSR
                if any(ismember({'FSR', 'FSRw', 'FSRwj'}, est_nam))
                    ttic = tic;
                    FSRout = FSR(y, X, 'lms', solLXS.bs, 'intercept', 0, ...
                        'init', floor(n/2)-1, 'msg', 0, 'plots', monit==1); 
                    if any(ismember('FSR', est_nam))
                        sol(est_i, :, iterr_i) = FSRout.beta';
                        tmp = ones(n, 1);
                        if ~isnan(FSRout.outliers)
                            tmp(FSRout.outliers) = 0;
                        end
                        solw(est_i,:, iterr_i) = tmp; 
                        sols(est_i, 1, iterr_i) =  1/(n-p) * (tmp' * (y - X * FSRout.beta).^2) / (sum(tmp)/n);
                        % FSRout.scale^2;
                        tim(est_i, iterr_i) = toc(ttic) + timLTS;
                        est_i = est_i + 1;
                    end              
                end
                
                % FSRwj & FSRws
                % we get two signals from FSR, and respectively estimate the weights
                % jointly and singularly
                if any(ismember({'FSRwj', 'FSRws'}, est_nam)) 
                    ttic = tic;
                    % 'weakened' FSR to detect a double signal
                    FSRoutW = FSR(y, X, 'lms', solLXS.bs, 'intercept', 0, ...
                        'init', floor(n/2), 'msg', 0, 'plots', monit==1, 'weak',1);
                    %FSRoutW = FSRout; % use FSR as usual
                    % stronger signal
                    trim_FSR =  FSRoutW.outliers';
                    if any(isnan(FSRoutW.outliers))
                        trim_FSR =  [];
                    end
                    % weaker signal
                    down_FSR = FSRoutW.VIOMout';
                    tempt = toc(ttic);
                    
                    % FSRwj: joint weights
                    if ismember('FSRwj', est_nam)
                        ttic = tic;
                        solFSRwj =  VIOM(y, X, down_FSR, 'trim', trim_FSR, ...
                            'mult', 1, 'intercept', 0);
                        sol(est_i, :, iterr_i) = solFSRwj.beta';
                        solw(est_i, :, iterr_i) = solFSRwj.w;
                        sols(est_i, 1, iterr_i) =  1/(n-p) * (solFSRwj.w' * (y - X * solFSRwj.beta).^2) / (sum(solFSRwj.w)/n);
                        % FSRoutW.scale^2;
                        tim(est_i, iterr_i) = toc(ttic) + tempt + timLTS;
                        est_i = est_i + 1;
                    end

                    % FSRws: single independent weights
                    if ismember('FSRws', est_nam)
                        ttic = tic;
                        solFSRws =  VIOM(y, X, down_FSR, 'trim', trim_FSR, ...
                            'intercept', 0);
                        sol(est_i, :, iterr_i) = solFSRws.beta';
                        solw(est_i, :, iterr_i) = solFSRws.w;
                        sols(est_i, 1, iterr_i) =  1/(n-p) * (solFSRws.w' * (y - X * solFSRws.beta).^2) / (sum(solFSRws.w)/n);
                        % FSRoutW.scale^2;
                        tim(est_i, iterr_i) = toc(ttic) + tempt + timLTS;
                        % est_i = est_i + 1;
                    end
                end
                            
                %% Monitoring each iteration
                if monit==1

                    % plot data and FSRws solution
                    pl1_dat

                    % plot diagnostics (MM weights monitoring, their
                    % derivative, FSR weights monitorin)
                    pl1_diag

                    %% Console messages
                    pl1_msg
                    
                end
                           
                % update progress bar
                if progr == 1
                    waitbar((iterr_i + iter_tot * (n_i-1) +  iter_tot * n_i * (SNR_i-1) ...
                        + iter_tot * n_i * SNR_i * (cont_i-1)) / ...
                        (iter_tot * length(n_tot) * length(SNR_tot) * length(frac_cont)));
                end
                
                % print iter
                if rem(iterr_i, 100) == 0 % every 100 iteration it appears by default
                    % [TBA: put it as a parameter at the begin]
                    fprintf('iteration: %d/%d \n', iterr_i, iter_tot);
                    fprintf('SNR: %d/%d  \n', SNR_i, length(SNR_tot));
                    fprintf('cont: %d/%d  \n', cont_i, length(frac_cont));
                    fprintf('n: %d/%d  \n', n_i, length(n_tot));
                end
                
                if monit == 1
                    cascade
                    % put a breakpoint here and uncomment when monit == 1
                    close all
                end
                
                % uncomment to monitor specific interations
                % if any(abs(solFSRws.beta) < 1) 
                % 	fprintf('bad fit'); 
                % end
                
            end % end iter - iterating over sample size
            
            %% store for each sample size
            
            % beta MSE decomp in var and bias^2
            vvB = mean((sol - mean(sol, 3)).^2, 3);
            bbB = (mean(sol, 3) - beta_true').^2;
            % store each coeff separately
            solVB(n_i, :, :) = vvB;
            solBB(n_i, :, :) = bbB;
                        
            % sigma MSE decomp in var and bias^2
            vvS = mean((sols - mean(sols, 3)).^2, 3);
            bbS = (mean(sols, 3) - sigma_e_u^2').^2;
            % store each coeff separately
            solVS(n_i, :, :) = vvS;
            solBS(n_i, :, :) = bbS;
            
            % time over iterr
            % averaged tiem over iterr for each n
            timn(:, n_i) = mean(tim, 2);        

            
            %% Monitoring iterations for each n
            if (monit==2 || monit==1)
                
                % plot bias (beta and weights [TBA_sigma]) and time
                pl2_bias
                
                % put a breakpoint here when monit == 2
                cascade
                % close all
            end
            
        end  % end sample size - iterating over contamination mechanism
        
        
        %% Monitoring each contamination mechanism
        if monit == 3
            
            % plot MSE results and time
            pl3_resN
            % put a breakpoint here when monit == 3
            cascade 
            
        end
        
    end  % end contamination mechanism - iterating over contamination fraction
    
        %if monit == 4
           % [TBA] 
           %pl4_resCont
           %cascade
        %end
end       
 


if progr == 1 % (it also closes figures)
    % delete progress bar
    close(hwait)
end

