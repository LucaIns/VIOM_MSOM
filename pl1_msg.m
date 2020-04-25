

                        % oracle
                        disp('--------------------------')
                        disp('oracle estimates: ')
                        disp(solOpt')
                        disp('--------------------------')

                        % LTS
                        disp('--------------------------')
                        disp('estimates LTS: ')
                        disp(solLXS.beta')
                        disp('--------------------------')

                        % MM
                        disp('--------------------------')
                        disp('estimates MM: ')
                        disp(solMM.beta')
                        disp('--------------------------')

                        % FSR
                        disp('--------------------------')
                        disp('estimates FSR: ')
                        disp(FSRout.beta')
                        disp('--------------------------')

                        % FSRwj
                        if ismember('FSRwj', est_nam)
                            disp('estimates FSRwj: ')
                            disp(solFSRwj.beta')
                            disp('--------------------------')
                        end

                        % FSRws
                        if ismember('FSRws', est_nam)
                            disp('estimates FSRws: ')
                            disp(solFSRws.beta')
                            disp('--------------------------')
                        end