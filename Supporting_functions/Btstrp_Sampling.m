<<<<<<< HEAD
function [Btstrp_NucVol, Btstrp_Acc,...
          Btstrp_Vis_Nuc, Btstrp_Vis_Cyto,...
          Btstrp_Vis_PeriCore,Btstrp_Vis_Domain,...
          Best_HyperPara] = Btstrp_Sampling(btstrp_index,num_per_loop,nWorkers,maxRetries,...
                                            Store_track, Store_mask, Store_fits,...
                                            bootstrp_list, morph_by_mov, num_pix_y,...
                                            locFitCI_UB, pixel_size, cell_avg_width,...
                                            StepSize, GapSize, weight_ratio)
    % ------------------------------------------------------------------
    % Turn worker-abort warnings to errors 
    % ------------------------------------------------------------------
    abortWarn1 = 'parallel:client:WorkerAborted';
    abortWarn2 = 'parallel:client:SpmdSessionGone';
    
    oldWarn1   = warning('query', abortWarn1);   % remember previous modes
    oldWarn2   = warning('query', abortWarn2);
    
    warning('error', abortWarn1);                % make them error
    warning('error', abortWarn2);
    % ------------------------------------------------------------------
    
    for block = 1:numel(btstrp_index)
        ind_s   = btstrp_index(block);
        ind_e   = ind_s + num_per_loop - 1;
        blk_len = ind_e - ind_s + 1;
    
        retries = 0;
        while retries < maxRetries
            try
                % ---------- Pool Initiation----------------------------
                if exist('pool','var') && ~isempty(pool)
                    evalc('delete(pool)');
                end
                [~, pool] = evalc( ...
                sprintf('parpool(''local'',%d,''IdleTimeout'',Inf);', nWorkers));
    
                % ---------- Broadcast for pool -------------------
                STC = parallel.pool.Constant(Store_track);   
                SMC = parallel.pool.Constant(Store_mask);
                SFC = parallel.pool.Constant(Store_fits);
    
                fprintf('\nBootstrap %d–%d  (attempt %d)\n',ind_s,ind_e,retries+1);
                tic
                [Btstrp_NucVol1, Btstrp_Acc1, Btstrp_Vis_Nuc1, ...
                 Btstrp_Vis_Cyto1, Btstrp_Vis_PeriCore1, Btstrp_Vis_Domain1, ...
                 Best_HyperPara1] = ...
                     Btstrp_PhyProp( bootstrp_list(ind_s:ind_e), ...
                         STC.Value, SMC.Value, SFC.Value, ...
                         morph_by_mov, num_pix_y, locFitCI_UB, pixel_size, ...
                         cell_avg_width,'all','SSIM',50, StepSize, ...
                         GapSize,{[0 1 0],[1 1 0]},0.05,{},weight_ratio );
                toc
    
                % ---------- Store Results -----------------------------
                Btstrp_NucVol(ind_s:ind_e)       = num2cell(Btstrp_NucVol1,2);
                Btstrp_Acc(ind_s:ind_e)          = num2cell(Btstrp_Acc1,2);
                Btstrp_Vis_Nuc(ind_s:ind_e)      = num2cell(Btstrp_Vis_Nuc1,2);
                Btstrp_Vis_Cyto(ind_s:ind_e)     = num2cell(Btstrp_Vis_Cyto1,2);
                Best_HyperPara(ind_s:ind_e)      = num2cell(Best_HyperPara1,2);
                Btstrp_Vis_PeriCore(ind_s:ind_e) = ...
                    [num2cell(Btstrp_Vis_PeriCore1,2) ; ...
                     repmat({[]},blk_len-size(Btstrp_Vis_PeriCore1,1),1)];
                Btstrp_Vis_Domain(ind_s:ind_e)   = ...
                    [num2cell(Btstrp_Vis_Domain1,2)   ; ...
                     repmat({[]},blk_len-size(Btstrp_Vis_Domain1,1),1)];
    
                evalc('delete(pool)');            % Delete pool after batch
                break                             
    
            catch ME
                warning('Block %d (%d–%d) failed: %s', ...
                         block,ind_s,ind_e,ME.message);
                evalc('delete(gcp(''nocreate''))');   
                retries = retries + 1;
                pause(2)                              
                if retries == maxRetries
                    error('Aborting: block %d failed %d times.',block,maxRetries);
                end
                % numWorkers = max(numWorkers-2,4);      
            end
        end
    end
    
    %% Flatten cells
    Btstrp_NucVol       = vertcat(Btstrp_NucVol{:});
    Btstrp_Acc          = vertcat(Btstrp_Acc{:});
    Btstrp_Vis_Nuc      = vertcat(Btstrp_Vis_Nuc{:});
    Btstrp_Vis_Cyto     = vertcat(Btstrp_Vis_Cyto{:});
    Btstrp_Vis_PeriCore = vertcat(Btstrp_Vis_PeriCore{:});
    Btstrp_Vis_Domain   = vertcat(Btstrp_Vis_Domain{:});
    Best_HyperPara      = vertcat(Best_HyperPara{:});
    
    clear ind_s ind_e blk_len btstrp_index num_per_loop
    clear Btstrp_*1
    
    % ------------------------------------------------------------------
    % Restore original warning modes 
    % ------------------------------------------------------------------
    warning(oldWarn1.state, abortWarn1);
    warning(oldWarn2.state, abortWarn2);
=======
function [Btstrp_NucVol, Btstrp_Acc,...
          Btstrp_Vis_Nuc, Btstrp_Vis_Cyto,...
          Btstrp_Vis_PeriCore,Btstrp_Vis_Domain,...
          Best_HyperPara] = Btstrp_Sampling(btstrp_index,num_per_loop,nWorkers,maxRetries,...
                                            Store_track, Store_mask, Store_fits,...
                                            bootstrp_list, morph_by_mov, num_pix_y,...
                                            locFitCI_UB, pixel_size, cell_avg_width,...
                                            StepSize, GapSize, weight_ratio)
    % ------------------------------------------------------------------
    % Turn worker-abort warnings to errors 
    % ------------------------------------------------------------------
    abortWarn1 = 'parallel:client:WorkerAborted';
    abortWarn2 = 'parallel:client:SpmdSessionGone';
    
    oldWarn1   = warning('query', abortWarn1);   % remember previous modes
    oldWarn2   = warning('query', abortWarn2);
    
    warning('error', abortWarn1);                % make them error
    warning('error', abortWarn2);
    % ------------------------------------------------------------------
    
    for block = 1:numel(btstrp_index)
        ind_s   = btstrp_index(block);
        ind_e   = ind_s + num_per_loop - 1;
        blk_len = ind_e - ind_s + 1;
    
        retries = 0;
        while retries < maxRetries
            try
                % ---------- Pool Initiation----------------------------
                if exist('pool','var') && ~isempty(pool)
                    evalc('delete(pool)');
                end
                [~, pool] = evalc( ...
                sprintf('parpool(''local'',%d,''IdleTimeout'',Inf);', nWorkers));
    
                % ---------- Broadcast for pool -------------------
                STC = parallel.pool.Constant(Store_track);   
                SMC = parallel.pool.Constant(Store_mask);
                SFC = parallel.pool.Constant(Store_fits);
    
                fprintf('\nBootstrap %d–%d  (attempt %d)\n',ind_s,ind_e,retries+1);
                tic
                [Btstrp_NucVol1, Btstrp_Acc1, Btstrp_Vis_Nuc1, ...
                 Btstrp_Vis_Cyto1, Btstrp_Vis_PeriCore1, Btstrp_Vis_Domain1, ...
                 Best_HyperPara1] = ...
                     Btstrp_PhyProp( bootstrp_list(ind_s:ind_e), ...
                         STC.Value, SMC.Value, SFC.Value, ...
                         morph_by_mov, num_pix_y, locFitCI_UB, pixel_size, ...
                         cell_avg_width,'all','SSIM',50, StepSize, ...
                         GapSize,{[0 1 0],[1 1 0]},0.05,{},weight_ratio );
                toc
    
                % ---------- Store Results -----------------------------
                Btstrp_NucVol(ind_s:ind_e)       = num2cell(Btstrp_NucVol1,2);
                Btstrp_Acc(ind_s:ind_e)          = num2cell(Btstrp_Acc1,2);
                Btstrp_Vis_Nuc(ind_s:ind_e)      = num2cell(Btstrp_Vis_Nuc1,2);
                Btstrp_Vis_Cyto(ind_s:ind_e)     = num2cell(Btstrp_Vis_Cyto1,2);
                Best_HyperPara(ind_s:ind_e)      = num2cell(Best_HyperPara1,2);
                Btstrp_Vis_PeriCore(ind_s:ind_e) = ...
                    [num2cell(Btstrp_Vis_PeriCore1,2) ; ...
                     repmat({[]},blk_len-size(Btstrp_Vis_PeriCore1,1),1)];
                Btstrp_Vis_Domain(ind_s:ind_e)   = ...
                    [num2cell(Btstrp_Vis_Domain1,2)   ; ...
                     repmat({[]},blk_len-size(Btstrp_Vis_Domain1,1),1)];
    
                evalc('delete(pool)');            % Delete pool after batch
                break                             
    
            catch ME
                warning('Block %d (%d–%d) failed: %s', ...
                         block,ind_s,ind_e,ME.message);
                evalc('delete(gcp(''nocreate''))');   
                retries = retries + 1;
                pause(2)                              
                if retries == maxRetries
                    error('Aborting: block %d failed %d times.',block,maxRetries);
                end
                % numWorkers = max(numWorkers-2,4);      
            end
        end
    end
    
    %% Flatten cells
    Btstrp_NucVol       = vertcat(Btstrp_NucVol{:});
    Btstrp_Acc          = vertcat(Btstrp_Acc{:});
    Btstrp_Vis_Nuc      = vertcat(Btstrp_Vis_Nuc{:});
    Btstrp_Vis_Cyto     = vertcat(Btstrp_Vis_Cyto{:});
    Btstrp_Vis_PeriCore = vertcat(Btstrp_Vis_PeriCore{:});
    Btstrp_Vis_Domain   = vertcat(Btstrp_Vis_Domain{:});
    Best_HyperPara      = vertcat(Best_HyperPara{:});
    
    clear ind_s ind_e blk_len btstrp_index num_per_loop
    clear Btstrp_*1
    
    % ------------------------------------------------------------------
    % Restore original warning modes 
    % ------------------------------------------------------------------
    warning(oldWarn1.state, abortWarn1);
    warning(oldWarn2.state, abortWarn2);
>>>>>>> 5e4cf7b (Initial commit)
end