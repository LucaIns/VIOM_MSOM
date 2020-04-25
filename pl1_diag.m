                    %% diagnostic: monitoring weights vs efficiency for MM
                    
                    if monit==1 

                        % monitoring efficiency
                        InitialEst.beta=solLXS.beta;
                        InitialEst.scale=solLXS.scale;
                        k = 1000;
                        eff = linspace(0.5, 0.99, k);
                        w = nan(n, k);
                        bMM = nan(p, k);
                        for j = 1:length(eff)
                            eff_i = eff(j);
                            [solMMest]=MMreg(y,X,'InitialEst',InitialEst, 'intercept', 0, 'rhofunc', 'bisquare', ...
                                'eff', eff_i, 'plots', 0);
                            w(:,j) = solMMest.weights;
                            bMM(:,j) = solMMest.beta;
                        end
                        
                        % shortcut to plot derivatives
                        % w = wMM;
                        figure 
                        h1=plot(eff, w(~ismember(1:n,ind_cont),:), 'b-', 'DisplayName', 'clean', 'LineWidth', 1);
                        hold on
                        if ~isempty(ind_VIOM)
                            h2=plot(eff, w(ind_VIOM, :), 'g--', 'DisplayName', 'VIOM', 'LineWidth', 1);
                            hold on
                        end
                        if ~isempty(ind_MSOM)
                            h3=plot(eff, w(ind_MSOM,:), 'r--', 'DisplayName', 'MSOM', 'LineWidth', 1);
                            hold on
                        end
                        if ~isempty(ind_VIOM) && ~isempty(ind_MSOM)
                            % clickableMultiLegend([h1(1), h2(1), h3(1)], 'clean', 'VIOM', 'MSOM', 'Location', 'southeast');
                            % axis manual
                            clb = legend([h1(1), h2(1), h3(1)], 'clean', 'VIOM', 'MSOM', 'Location', 'northwest');
                            set(clb,'FontSize',8);
                        else 
                            % legend()
                        end
%                         title(sprintf('monitoring weights deriv. as a function of efficiency in MM-estimator \n (SNR: %d, MSOM: %.2f, VIOM: %.2f)', ... 
%                             SNR_tot(SNR_i), cont_fracMSOM_i, frac_cont(cont_i)-cont_fracMSOM_i), ...
%                             'FontSize', 12);
                        xlabel('Efficiency', 'FontSize', 20); 
                        ylabel('MM weights estimates', 'FontSize', 20);
                     
                        % beta
                        % figure 
                        % h1=plot(eff, bMM, 'b-');
                        % title(sprintf('beta monitoring as a function of efficiency in MM-estimator (SNR: %d, MSOM: %.2f, VIOM: %.2f)', ... 
                        %     SNR_tot(SNR_i), cont_fracMSOM_i, cont_frac(cont_i)-cont_fracMSOM_i), ...
                        %     'FontSize', 20);
                        % xlabel('Efficiency', 'fontweight', 'bold', 'FontSize', 18); 
                        % ylabel('beta estimates', 'fontweight', 'bold', ...
                        %     'FontSize', 18);

                        %% Diagnostics: Weights monitoring using FSR

                        % we added weight estimation in FSReda (only few lines of code)
                        % replace the FSDA function and check the differences
                        wwlab = 1;
                        [sol_FS] = FSReda(y, X, solLXS.bs, 'intercept', 0, 'wREML', wwlab);
                        sol_FS.RES = sol_FS.w;
                            % databrush.persist = 'on';
                            % databrush.RemoveLabels='on';
                            % databrush.RemoveLabels='off';
                            % databrush.Label='on'; 
                            % databrush.selectionmode='rect';
                            % resfwdplot(sol_FS, 'databrush', databrush); % 'datatooltip', 'datatip', 
                        datatooltip.DisplayStyle = 'datatip';
                        fground.fthresh = 10;
                        resfwdplot(sol_FS,'datatooltip',datatooltip,'fground',fground, 'tag','myplot');
%                         title(sprintf('Monitoring using FSR (SNR: %d, MSOM: %.2f, VIOM: %.2f)', ... 
%                             SNR_tot(SNR_i), cont_fracMSOM_i, frac_cont(cont_i)-cont_fracMSOM_i), ...
%                             'FontSize', 20);
                        xlabel('iteration b', 'FontSize', 30); 
                        if wwlab == 1
                            ylabel('FSRws weights estimates', 'FontSize', 30);
                        else
                            ylabel('Scaled residuals', 'FontSize', 30);
                        end
                        xtn = get(gca,'XTickLabel');
                        set(gca,'XTickLabel', xtn,'fontsize',27)
                        ytn = get(gca,'YTickLabel');
                        set(gca,'yTickLabel', ytn,'fontsize',27)
                        
                        
                    %% MM DERIVATIVE
                                         
                        % interactive
                        h = 1/k;       % step size
                        wMM = [diff(w, 1, 2), zeros(n, 1)]/h;   % 1st derivative
                        wMM(abs(wMM)>2) = NaN;
                        w = wMM;
                        figure 
                        h1=plot(eff, w(~ismember(1:n,ind_cont),:), 'b-', 'DisplayName', 'clean', 'LineWidth', 1);
                        hold on
                        if ~isempty(ind_VIOM)
                            h2=plot(eff, w(ind_VIOM, :), 'g--', 'DisplayName', 'VIOM', 'LineWidth', 1);
                            hold on
                        end
                        if ~isempty(ind_MSOM)
                            h3=plot(eff, w(ind_MSOM,:), 'r--', 'DisplayName', 'MSOM', 'LineWidth', 1);
                            hold on
                        end
                        if ~isempty(ind_VIOM) && ~isempty(ind_MSOM)
                            % clickableMultiLegend([h1(1), h2(1), h3(1)], 'clean', 'VIOM', 'MSOM', 'Location', 'southeast');
                            % axis manual
                            clb = legend([h1(1), h2(1), h3(1)], 'clean', 'VIOM', 'MSOM', 'Location', 'northwest');
                            set(clb,'FontSize',8);
                        else 
                            % legend()
                        end
%                         title(sprintf('monitoring weights deriv. as a function of efficiency in MM-estimator \n (SNR: %d, MSOM: %.2f, VIOM: %.2f)', ... 
%                             SNR_tot(SNR_i), cont_fracMSOM_i, frac_cont(cont_i)-cont_fracMSOM_i), ...
%                             'FontSize', 12);
                        xlabel('Efficiency', 'FontSize', 20); 
                        ylabel('MM weights derivatives estimates', 'FontSize', 20);
                        
                        % using FSReda
    %                         wMM(abs(wMM)>2) = NaN;
    %                         sol_FS.RES = wMM'; 
    %                             % databrush.persist = 'off';
    %                             % databrush.RemoveLabels='on';
    %                             % databrush.RemoveLabels='off';
    %                             % databrush.Label='on'; 
    %                             % databrush.selectionmode='rect';
    %                             % resfwdplot(sol_FS, 'databrush', databrush); % 'datatooltip', 'datatip', 
    %                         resfwdplot(sol_FS, 'datatooltip', 'datatip'); 
    %                         %title(sprintf('Weights derivative monitoring (SNR: %d, MSOM: %.2f, VIOM: %.2f)', ... 
    %                          %   SNR_tot(SNR_i), cont_fracMSOM_i, frac_cont(cont_i)-cont_fracMSOM_i), ...
    %                          %   'FontSize', 10);
    %                         xlabel('Efficiency', 'FontSize', 20); 
    %                         ylabel('MM weights estimates derivatives', 'FontSize', 20);
    %                         xticks('manual')
    %                         xticklabels(round(linspace(0.5, 1, 11), 2));
    %                         % ylim([-1, 10])
    
                    


                    end