
clc, clear, clf, close all

% choose dataset 
datch = 3;
% intercept
intercept = 1;

% adopted estimators (here the order does not matter)
% est_nam_f = {OLS, 'LTS', 'MM', 'FSR', 'FSRwj', 'FSRws', 'OLS'};
est_nam = {'OLS', 'MM', 'FSRws'}; % 'LMS', FSR
len_est = length(est_nam);

%% load data

if datch == 1
    %% Cook et al. (1982) data
    X = [116, 132, 104, 139, 114, 129, 720, 174, 312, 338, 465]';
    y = [105, 120, 85, 121, 115, 127, 630, 155, 248.4, 310, 443]';
    tit = sprintf('Data used by Cook et al. (1982)');
    xla = 'Destructive thickness measure';
    yla = 'Non-destructive thickness measure';
elseif datch == 2
    %% Gumedze (2018) data
    data = importdata('seedData.csv');
    y = data.data(:,2);
    X = data.data(:,3);
    tit = sprintf('Data used by Gumedze (2018)');
    xla = 'Seed length';
    yla = 'Seed weight';
elseif datch == 3
    %% loyalty data
    data = importdata('loyalty.mat');
    y = data.data(:,end);
    X = data.data(:,1);
    tit = sprintf('Loyalty data');
    xla = 'Number of visits';
    yla = 'Amount spent (in €)';
    intercept = 0;
elseif datch == 4
    %% bank data
    data = importdata('BalanceSheets.mat');
    y = data.data(:,end);
    X = data.data(:,1);
    % X = data.data(:,1:end-1);
    tit = sprintf('Balance sheets data');
    xla = 'Labour Share';
    yla = 'Return on Sales';
end 

n = size(X, 1);
if intercept == 1
    X = [ones(n, 1), X];
end
p = size(X, 2);

% collect estimation results
betaHat = nan(len_est, p);
h = nan(len_est, 1);

%% estimators

% iterate accross estimators 
% respecting the pre-specified order
for i = 1:len_est
    
    est_nam_i = est_nam(i);
    
    % OLS
    if strcmp('OLS', est_nam_i)
        betaHat(i,:) = regress(y, X);
    end

    % LTS
    if any(ismember({'LTS', 'LMS', 'FSR', 'FSRws', 'FSRwj', 'MM'}, est_nam_i))
        % LTS
        solLXS = LXS(y, X, 'lms', 1, 'h', floor(0.5*(n+p+1)), 'intercept', 0, ...
            'msg', 0, 'plots', 1, 'conflev', 0.9); % 'rew', 1,
        betaHat(i,:) = solLXS.beta';
    end

    % MM
    if strcmp('MM', est_nam_i)
        InitialEst.beta=solLXS.beta;
        InitialEst.scale=solLXS.scale;
        [solMM]=MMreg(y,X,'intercept', 0, 'InitialEst', InitialEst, ...
            'rhofunc', 'bisquare', 'eff', 0.85, 'plots', 1);
        betaHat(i,:) = solMM.beta;  
    end

    % FSR
    if strcmp('FSR', est_nam_i)
        FSRout = FSR(y, X, 'lms', solLXS.bs, 'intercept', 0, ...
            'init', floor((n+p)/2), 'msg', 0, 'plots', 1); % , 'bonflev', 0.95<
        betaHat(i,:) = FSRout.beta';
    end
    
    % FSRwj & FSRws
    % we get two signals from FSR, and respectively estimate the weights
    % jointly and singularly
    if any(ismember({'FSRwj', 'FSRws'}, est_nam_i)) 
        % 'weakeed' FSR to detect a double signal
        FSRoutw = FSR(y, X, 'lms', solLXS.bs, 'intercept', 0, ...
            'init', floor(n/2)-1, 'msg', 0, 'plots', 1, 'weak',1);
        % stronger signal
         trim_FSR =  FSRoutw.outliers;
        if any(isnan(FSRoutw.outliers))
            trim_FSR =  [];
        end
        % weaker signal
        down_FSR = FSRoutw.VIOMout; % setdiff(1:n, [FSRoutw.ListCl'; trim_FSR])';  
        if any(isnan(down_FSR))
            down_FSR = solLXS.outliers;
            warning('FSR has found no outliers; LMS is used to detect VIOM');
        end

        % FSRwj: joint weights
        if ismember('FSRwj', est_nam_i)
            solFSRwj =  VIOM(y, X, down_FSR, 'trim', trim_FSR, ...
                'mult', 1, 'intercept', 0);
            betaHat(i,:) = solFSRwj.beta';
        end

        % FSRws: single independent weights
        if ismember('FSRws', est_nam_i)
            ttic = tic;
            solFSRws =  VIOM(y, X, down_FSR, 'trim', trim_FSR, ...
                'intercept', 0);
            betaHat(i,:) = solFSRws.beta';
        end
    end    
end

%% diagnostics

cleans = setdiff(1:n, [down_FSR, trim_FSR]);

if any(strcmp('MM', est_nam))
    % MM
    InitialEst.beta=solLXS.beta;
    InitialEst.scale=solLXS.scale;
    k = 100;
    eff = linspace(0.5, 0.999, k);
    w = nan(n, k);
    for j = 1:length(eff)
        eff_i = eff(j);
        [solMMdiag]=MMreg(y,X,'intercept', 0, 'rhofunc', 'bisquare', 'InitialEst', InitialEst, ...
            'eff', eff_i, 'plots', 0, 'msg', 0);
        w(:,j) = solMMdiag.weights;
    end
    %w = wMM;
    figure % try with resfwplot to get an interactive plot
    h1=plot(eff, w(cleans,:), 'b-', 'DisplayName', 'clean', 'LineWidth', 1);
    if ~isnan(trim_FSR)
        hold on
        h2=plot(eff, w(trim_FSR,:), 'r-.', 'DisplayName', 'MSOM', 'LineWidth', 1);
    end
    hold on
    h3=plot(eff, w(down_FSR,:), 'g--', 'DisplayName', 'VIOM', 'LineWidth', 1);
    if ~isnan(trim_FSR)
        clickableMultiLegend([h1(1), h2(1), h3(1)], 'clean', 'MSOM', 'VIOM', 'Location', 'northwest');
    %     clb = legend([h1(1), h2(1), h3(1)], 'clean', 'MSOM', 'VIOM', 'Location', 'northwest');
    %     set(clb,'FontSize',12);
    else
        clickableMultiLegend([h1(1), h3(1)], 'clean', 'VIOM', 'Location', 'northwest');
    end
    %title(sprintf('Weights monitoring as a function of efficiency in MM-estimator'), ...
    %    'FontSize', 12);
    title('')
    xlabel('Efficiency', 'FontSize', 15); 
    ylabel('MM weights derivatives estimates', 'FontSize', 15); 
    axis manual;
end

%% FSR for weights

[sol_FS] = FSReda(y, X, solLXS.bs, 'intercept', 0, 'wREML', 1); 
sol_FS.RES = sol_FS.w;
datatooltip.DisplayStyle = 'datatip'; 
fground.fthresh = 50; %1.5; %2.5;
resfwdplot(sol_FS,'datatooltip',datatooltip,'tag','FSwdiag33', 'fground', fground);
    % databrush.persist = 'on';
    % databrush.RemoveLabels='on';
    % databrush.RemoveLabels='off';
    % databrush.Label='on'; 
    % databrush.selectionmode='rect';
    % resfwdplot(sol_FS, 'databrush', databrush); 
%title(sprintf('weights monitoring using FSR'), 'FontSize', 12);
    % title('');
xlabel('iteration b') %, 'FontSize', 25); 
ylabel('FSRws weights estimates') %, 'FontSize', 25);
if datch ~= 1
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',25)
    set(gca,'YTickLabel',get(gca,'YTickLabel'),'fontsize',25)
    xlabel('iteration b', 'FontSize', 25); 
    ylabel('FSRws weights estimates', 'FontSize', 25);
else
    xlabel('iteration b'); 
    ylabel('FSRws weights estimates');
end
% text(6,0.5, num2str(1), 'fontsize',12);

%% MM derivative

if any(strcmp('MM', est_nam))
    h = 1/k;       % step size
    wMM = [diff(w, 1, 2), zeros(n, 1)]/h;   % 1st derivative
    wMM(abs(wMM)>7) = NaN;
    % % wMM = w; % plot simple weights
    sol_FS.RES = wMM; 
        % databrush.persist = 'on';
        % databrush.RemoveLabels='on';
        % databrush.RemoveLabels='off';
        % databrush.Label='on'; 
        % databrush.selectionmode='rect';
        % resfwdplot(sol_FS, 'databrush', databrush); 
    resfwdplot(sol_FS, 'datatooltip', 'datatip'); 
    % title(sprintf('Weights derivative monitoring'), 'FontSize', 10);
    xlabel('Efficiency'); 
    ylabel('MM weights estimates derivatives'); 
end
%% data

figure
hold on
plot(X(:, intercept+1), y, 'b.', 'MarkerSize', 15);
if ~isnan(trim_FSR)
    hold on
    plot(X(trim_FSR, intercept+1), y(trim_FSR), 'r.', 'MarkerSize', 15);
%     text(X(trim_FSR, intercept+1),y(trim_FSR), num2str(trim_FSR), 'fontsize',10);
end
hold on
plot(X(down_FSR, intercept+1), y(down_FSR), 'g.', 'MarkerSize', 15);
% text(X(down_FSR, intercept+1), y(down_FSR), num2str(down_FSR), 'fontsize',10);

    % txun = [5, 6, 7, 9, 11]';
    % text(X(txun, intercept+1), y(txun)+25, num2str(txun), 'fontsize',15);

lst = {'--','-', '-.', ':', '--','-','-.',':',};

col = 'bkgmycbkgmyc';
for i = 1:len_est
    hold on
    if intercept == 1
        a = betaHat(i, 1);
        b = betaHat(i, 2);
    else
        a = 0;
        b = betaHat(i);
    end 
    h(i) = refline(b,a);
    set(h(i),'color', col(i),'LineWidth', 1.2, 'LineStyle', lst{i});
end
clb = clickableMultiLegend(h, est_nam, 'Location', 'northeast');
set(clb,'FontSize',12);
    % title(tit, 'FontSize', 10);
xlabel(xla); 
ylabel(yla); 
% grid
box
cascade
