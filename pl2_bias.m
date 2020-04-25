                %% plot beta bias over iterations
                solBiasB = (sol - beta_true');
                % recombine solution
                % (iteration * estimator * coeff)
                ssolB = permute(solBiasB, [3 1 2]);
                % col = 'bkgmycbkgmycbkgmyc';
                for j = 1:p
                    figure()
                    violin(ssolB(:,:,j),'xlabel', est_nam, 'facecolor', parula(lest), ...
                        'edgecolor','k', 'mc','k-.','medc','r-.');
                    ylabel('Bias','FontSize',14)
                    hold on
                    plot(get(gca,'xlim'),[0 0], 'g--', 'DisplayName', 'Baseline');
                    title(sprintf('Bias in Beta %d for n=%d and %d iterations', j-1, n, iter_tot));
                    %clickableMultiLegend(L)
                    
                    % NOTE: compare it with a boxplot
%                         figure
%                         boxplot(ssol(:,:,j), est_nam)
%                         title(sprintf('Beta %d', j-1));
%                         hold on
%                         xlim=get(gca,'xlim');
%                         plot(xlim,[0 0], '--')
                end
                

                
                %% plot time over iterations
                
%                 col = 'bkgmycrbkgmycbkgmyc';
%                 figure
%                 for t_i = 1:lest
%                     col_ch = col(t_i);
%                     h1 = plot(1:iterr, tim(t_i, :), col_ch, ...
%                         'DisplayName', string(est_nam(t_i)), 'LineWidth', 2);
%                     hold on
%                 end
%                 hold on
%                 clickableMultiLegend(gca, est_nam, 'Location', 'northwest')
%                 title(sprintf('Computing time for n=%d', n), ...
%                     'Interpreter', 'latex');
%                 xlabel('iteration', 'fontweight', 'bold', 'FontSize', 15, 'Interpreter', 'latex');
%                 ylabel('time (in sec.)', 'fontweight', 'bold', ...
%                     'FontSize', 15, 'Interpreter', 'latex');
%                 % axis manual
