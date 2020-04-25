            %% Betas: Variance bias decomposition 
            solVBp = mean(solVB, 3);
            solBBp = mean(solBB, 3);
            figure()
            lestR = length(est_nam);
            for es_i = 1:lestR
                subplot(2, round((lestR)/2), es_i);
                h1 = plot(n_tot,solVBp(:, es_i),  'g-', ...
                    'DisplayName', 'Var', 'LineWidth', 3);
                hold on
                h2 = plot(n_tot,solBBp(:, es_i), 'b-', ...
                    'DisplayName', 'Sq. bias', 'LineWidth', 3);
                hold on
                h3 = plot(n_tot,solVBp(:, es_i) + solBBp(:, es_i), 'r--', ...
                    'DisplayName', 'MSE', 'LineWidth', 3);
                title(sprintf('%s', est_nam{es_i}), ...
                    'FontSize', 8, 'fontweight', 'bold');
                xlabel('n', 'FontSize', 5);
                ylabel('MSE', 'FontSize', 5);
                hnd = gca;
                hnd.YAxis.FontSize = 10;
                hnd.XAxis.FontSize= 10;
                limax = max(solVBp+solBBp);
                if es_i > 3 || es_i == 1 % || es_i == 3
                    ylim([0,max(limax([1,4:lest]))]);
                else
                    ylim([0,max(limax(es_i), max(limax([1,4:lest])))]);
                end
                xlim([n_tot(1), n_tot(end)]);
                % axis manual
                grid();
            end
            % for MATLAB2019b
            sgtitle(sprintf('MSE beta \n (SNR: %d, MSOM: %.3f, VIOM: %.3f)', ... 
               SNR_tot(SNR_i), cont_fracMSOM_i, frac_cont(cont_i)-cont_fracMSOM_i));
            % for earlier than MATLAB2019b
            %suptitle(sprintf('MSE beta \n (SNR: %d, MSOM: %.3f, VIOM: %.3f)', ... 
            %   SNR_tot(SNR_i), cont_fracMSOM_i, frac_cont(cont_i)-cont_fracMSOM_i));
            %
            % clb = clickableMultiLegend(gca, 'Var', 'Sq. bias', 'MSE', 'Location', 'northwest');
            clb = legend(gca, 'Var', 'Sq. bias', 'MSE', 'Location', 'northeast');
            set(clb,'FontSize',8);
            axis manual



            %% sigma (WLS): Variance bias decomposition 
            solVSp = solVS(:,:,1);
            solBSp = solBS(:,:,1);
            figure()
            lestR = length(est_nam);
            for es_i = 1:lestR
                subplot(2, round((lestR)/2), es_i);
                h1 = plot(n_tot,solVSp(:, es_i),  'g-', ...
                    'DisplayName', 'Var', 'LineWidth', 3);
                hold on
                h2 = plot(n_tot,solBSp(:, es_i), 'b-', ...
                    'DisplayName', 'Sq. bias', 'LineWidth', 3);
                hold on
                h3 = plot(n_tot,solVSp(:, es_i) + solBSp(:, es_i), 'r--', ...
                    'DisplayName', 'MSE', 'LineWidth', 3);
                title(sprintf('%s', est_nam{es_i}), ...
                   'FontSize', 8, 'fontweight', 'bold');
                xlabel('n', 'FontSize', 5);
                ylabel('MSE', 'FontSize', 5);
                hnd = gca;
                hnd.YAxis.FontSize = 10;
                hnd.XAxis.FontSize= 10;
                limax = max(solVSp+solBSp);
                if es_i > 3 || es_i == 1
                    ylim([0,max(limax([1,4:lest]))]);
                else
                    ylim([0,max(limax(es_i), max(limax([1,4:lest])))]);
                end
                xlim([n_tot(1), n_tot(end)]);
                % axis manual
                grid();
            end
%             suptitle(sprintf('MSE sigma with WLS res \n (SNR: %d, MSOM: %.3f, VIOM: %.3f)', ... 
%                 SNR_tot(SNR_i), cont_fracMSOM_i, frac_cont(cont_i)-cont_fracMSOM_i));
            sgtitle(sprintf('MSE sigma with WLS res \n (SNR: %d, MSOM: %.3f, VIOM: %.3f)', ... 
                SNR_tot(SNR_i), cont_fracMSOM_i, frac_cont(cont_i)-cont_fracMSOM_i));
            % clb = clickableMultiLegend(gca, 'Var', 'Sq. bias', 'MSE', 'Location', 'northwest');
            clb = legend(gca, 'Var', 'Sq. bias', 'MSE', 'Location', 'northeast');
            set(clb,'FontSize',8);
            axis manual


            %% time vs n

            col = 'bkgmycrbkgmycbkgmyc';
            figure
            for t_i = 1:lest
                col_ch = col(t_i);
                h1 = plot(n_tot, timn(t_i, :), col_ch, ...
                    'DisplayName', string(est_nam(t_i)), 'LineWidth', 2);
                hold on
            end
            hold on
            % clickableMultiLegend(gca, est_nam, 'Location', 'northwest')
            clb = legend(gca, est_nam, 'Location', 'northwest');
            set(clb,'FontSize',10);
            title(sprintf('Computing time \n (SNR: %d, MSOM: %.3f, VIOM: %.3f)', ... 
                SNR_tot(SNR_i), cont_fracMSOM_i, frac_cont(cont_i)-cont_fracMSOM_i)); 
            xlabel('n', 'FontSize', 10); 
            ylabel('time (in sec.)', 'FontSize', 10);
            % axis manual
            grid();
