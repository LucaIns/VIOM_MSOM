                    % plot data
                    if monit==1 && p<=2  %  (monit==1 || iterr_i ==iterr)

                        % Real model
                        figure
                        h1 = plot(X(ind_keep,intercept+1), y(ind_keep), 'b.', 'DisplayName', 'clean');
                        % VIOM
                        hold on
                        h2 = plot(X(ind_VIOM,intercept+1), y(ind_VIOM), 'g.', 'DisplayName', 'VIOM');
                        if n < 200
                            text(X(ind_VIOM,intercept+1),y(ind_VIOM),num2str(((ind_VIOM))'));
                            % highlight selected units
                            % msel = [11, 24, 43, 50];
                            % text(X(msel,intercept+1),y(msel),num2str(((msel))'));
                        end
                        % MSOM
                        hold on
                        h3 = plot(X(ind_MSOM,intercept+1), y(ind_MSOM), 'r.', 'DisplayName', 'MSOM');
                        if n < 200
                            text(X(ind_MSOM,intercept+1),y(ind_MSOM),num2str(((ind_MSOM))'));
                        end
                        ylim([min(y)-1, max(y)+1]);
                        xlim([min(X(:,2))-1, max(X(:,2))+1]);
                        % estimates
                        a = sol(:, 1, iterr_i);
                        b = sol(:, 2, iterr_i);
                        h=refline(b, a);
%                         h.LineWidth =  cell2mat(get(h, 'LineWidth')) .* 2;
                        % legend
                        hold on
                        hleg = legend('Location', 'northwest');
                        hleg.String(end-length(est_nam)+1:end) = est_nam;
                        if monit == 11 % ~isempty(ind_VIOM) && ~isempty(ind_MSOM)
                            clickableMultiLegend([h1, h2, h3], 'clean', 'VIOM', 'MSOM', 'Location', 'northwest');
                            axis manual
                        end
                        title('Real model')


                        % FSRws solution
                        ind_MSOMfsrws = trim_FSR';
                        ind_VIOMfsrws = down_FSR';
                        ind_keepfsrws = find(~ismember(1:n, down_FSR));
                        figure
                        h1 = plot(X(ind_keepfsrws,intercept+1), y(ind_keepfsrws), 'b.', 'DisplayName', 'clean');
                        % VIOM
                        hold on
                        h2 = plot(X(ind_VIOMfsrws,intercept+1), y(ind_VIOMfsrws), 'g.', 'DisplayName', 'VIOM');
                        if n < 200
                            text(X(ind_VIOMfsrws,intercept+1),y(ind_VIOMfsrws),num2str(((ind_VIOMfsrws))'));
                        end
                        % MSOM
                        hold on
                        h3 = plot(X(ind_MSOMfsrws,intercept+1), y(ind_MSOMfsrws), 'r.', 'DisplayName', 'MSOM');
                        if n < 200
                            text(X(ind_MSOMfsrws,intercept+1),y(ind_MSOMfsrws),num2str(((ind_MSOMfsrws))'));
                        end
                        hold on
                        if ~isempty(ind_VIOMfsrws) && ~isempty(ind_MSOMfsrws)
                            clickableMultiLegend([h1, h2, h3], 'clean', 'VIOM', 'MSOM', 'Location', 'northwest');
                            axis manual
                        else %FIX: or add all cases with clickablelegend
                            legend()
                        end
                        title('FSRws solution')
                        ylim([min(y)-1, max(y)+1]);
                        xlim([min(X(:,2))-1, max(X(:,2))+1]);

                    end