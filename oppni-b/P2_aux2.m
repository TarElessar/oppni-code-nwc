function P2_aux2( broadcast )
% check existence of subject specific struct file
    if exist(fullfile( outpath,broadcast.subj,'InputStruct_ssa.mat'),'file')  
        InputStruct_ssa = load( fullfile( outpath,broadcast.subj,'InputStruct_ssa.mat') );
    else
        error('cannot find Input struct file for subject: %s \n');
    end

    % quick formatting stuff, again assuRImes that directory structure was already constructed in "P0" pipeline step 
    opath0 = fullfile(outpath,InputStruct_ssa.PREFIX,'rawdata');
    %
    opath1a = fullfile( outpath, InputStruct_ssa.PREFIX,'anat_proc');
    opath2a = fullfile( opath1a,sprintf('subpipe_%03u',    broadcast.PipeStruct.pipe_idx.Warp));
    opath3a = fullfile( opath2a,sprintf('seg_subpipe_%03u',broadcast.PipeStruct.pipe_idx.Seg ));
    %
    opath1p = fullfile(outpath,InputStruct_ssa.PREFIX,'phys_proc');
    %
    opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
    opath2f = fullfile(opath1f,sprintf('subpipe_%03u',broadcast.PipeStruct.pipe_idx.P1));
    opath3f = fullfile(opath2f,sprintf('seg_subpipe_%03u',broadcast.PipeStruct.pipe_idx.Seg ));
    %
    opath4f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',['pipe_',broadcast.PipeStruct.PNAME{1}]); 

    fprintf('\n===> subj %s (%u), func. processing (part-2)...\n',broadcast.subj,broadcast.index),
    %% =======================================================================
    %%      FUNC Processing, Block-2 ...
    %% =======================================================================



    %%%%%========== compatibility adjustments for BLOCK2 ... TASKREG only
    missing_physio=0;
    for nr=1:InputStruct_ssa.N_func
        if isempty(InputStruct_ssa.frun(nr).PHYSIO_filename)
            missing_physio=missing_physio+1;
        end
    end
    if analysis_struct.uses_taskfile==0  && ~strcmpi(broadcast.PipeStruct.TASKREG{1},'OFF')
        warning('your analysis model does not specify a task design. Turning TASKREG off!')
        broadcast.PipeStruct.TASKREG{1}={'OFF'};
    end
    %%%%%========== compatibility adjustments, done

    anymiss=0;
    for nr = 1:InputStruct_ssa.N_func
        if ~exist( [opath4f,'/func',num2str(nr),'_fullproc.mat'],'file')
            anymiss=1;
        end
    end

    % appending info about run-length (for regressor construction and interp-contrast-list)
    for nr = 1:InputStruct_ssa.N_func
        %
        opath0   = fullfile(outpath,InputStruct_ssa.PREFIX,'rawdata');
        unix(sprintf('gunzip %s/func%u.nii.gz',opath0,nr));
        hdr      = load_nii_hdr(sprintf('%s/func%u.nii',opath0,nr)); %***%
        unix(sprintf('gzip %s/func%u.nii',opath0,nr));
        InputStruct_ssa.frun(nr).Nt_raw   = hdr.dime.dim(5);
        InputStruct_ssa.frun(nr).Nt_adj   = hdr.dime.dim(5) - InputStruct_ssa.frun(nr).DROP_first - InputStruct_ssa.frun(nr).DROP_last;
    end

    % extract contrasts specified for analysis --> appended into InputStruct, delete unformatted split-info field 
    InputStruct_ssa = interpret_contrast_list( InputStruct_ssa, analysis_struct, ParamStruct_aug.CONTRAST); % generate contrast list for each subject and run    
    %InputStruct_ssa = rmfield(InputStruct_ssa,'task_unformat');
    % extract seed contrasts SAA
    InputStruct_ssa = interpret_seed_list( InputStruct_ssa, analysis_struct, ParamStruct_aug.CONTRAST); % generate contrast list for each subject and run 
    %InputStruct_ssa = rmfield(InputStruct_ssa,'seed_unformat');

    if anymiss==0 && maskisnew==0 % we can skip all of this, if data are processed and masks havent been re-constructed
        disp('Files found, no update to mask(s). Skipping last proc. stage!')
    else
        %--> previously specified opaths over here?
        
        for nr=1:InputStruct_ssa.N_func
    
            if ~exist( [opath4f,'/func',num2str(nr),'_fullproc.mat'],'file') || maskisnew==1 %only do runs with missing outputs / or if mask was just updated
        
                % loading data into mats:
                VSstr = sprintf('%s/postwarp/func%u_warped_smo.nii.gz',opath2f,nr); %***%
                VS = load_untouch_niiz(VSstr);
                volmatS = nifti_to_mat(VS,MB); 
            
                %% first part of regression proc: constructing regressors from DETREG,GSREG,MOTREG,ROIREG,TASKREG 
            
                % >>> Temporal Detrending
                Step = 'DETREG';
                if strcmpi(broadcast.PipeStruct.(Step){1},'OFF')
                    xdet = ones( InputStruct_ssa.frun(nr).Nt_adj, 1 );
                    statd = [];
                else
                    % get function handle for analysis model of interest
                    currPath=pwd;                               % get current path
                    cd(pipeline_struct.(Step).filepath);               % jump to module directory
                    pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
                    cd(currPath);                               % jump back to current path
                    % execute step:
                    [xdet,statd] = pfun( InputStruct_ssa.frun(nr).Nt_adj, InputStruct_ssa.TR_MSEC, broadcast.PipeStruct.(Step)(2:end) );  
                end

                % >>> Global Signal Regression
                Step = 'GSREG';
                if strcmpi(broadcast.PipeStruct.(Step){1},'OFF')
                    xglb = [];
                    statg = [];
                else
                    % get function handle for analysis model of interest
                    currPath=pwd;                               % get current path
                    cd(pipeline_struct.(Step).filepath);               % jump to module directory
                    pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
                    cd(currPath);                               % jump back to current path
                    % execute step:
                    [xglb,statg] = pfun( volmatS, broadcast.PipeStruct.(Step)(2:end) );  
                end

                % >>> Motion Parameter Regression
                Step = 'MOTREG';
                if strcmpi(broadcast.PipeStruct.(Step){1},'OFF')
                    xmot = [];
                    statm = [];
                else
                    % get function handle for analysis model of interest
                    currPath=pwd;                               % get current path
                    cd(pipeline_struct.(Step).filepath);               % jump to module directory
                    pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
                    cd(currPath);                               % jump back to current path
                    % execute step:
                    [xmot,statm] = pfun( sprintf('%s/warp/func%u_mpe',opath2f,nr), broadcast.PipeStruct.(Step)(2:end) );  
                end

                % >>> Noise ROI Regression
                Step = 'ROIREG';
                if strcmpi(broadcast.PipeStruct.(Step){1},'OFF')
                    xroi = [];
                    statr = [];
                else
                    maskpath  = [outpath,'/_group_level/masks/pipe_',broadcast.PipeStruct.PNAME{1}];
                    parcpath  = [outpath,'/_group_level/parcellations/pipe_',broadcast.PipeStruct.PNAME{1}];
                    % get function handle for analysis model of interest
                    currPath=pwd;                               % get current path
                    cd(pipeline_struct.(Step).filepath);               % jump to module directory
                    pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
                    cd(currPath);                               % jump back to current path
                    % execute step:
                    [xroi,statr] = pfun( sprintf('%s/postwarp/func%u_warped.nii.gz',opath2f,nr), {maskpath,parcpath}, broadcast.PipeStruct.(Step)(2:end) );  
                end
                
                % >>> Regression of Task Design
                Step = 'TASKREG';
                if strcmpi(broadcast.PipeStruct.(Step){1},'OFF')
                    Xsignal = [];
                elseif strcmpi(broadcast.PipeStruct.TASKREG{1},'OP1')
                    Xsignal = InputStruct_ssa.task_contrast{nr}.design_mat;
                    Xsignal = bsxfun(@rdivide, Xsignal, sqrt(sum(Xsignal.^2)));
                else
                    error('unknown taskreg option!');
                end

            
                %% second part of regression proc: concat vectors and apply to data matrix 
        
                % --> some regressional stats
                ztmp = zscore(volmatS')';
                if nr==1 && ~isempty(xdet) && size(xdet,2)>1
                    xtmp = zscore(xdet(:,2:end));
                    RegCorr2.det(:,:,broadcast.index) = (ztmp*xtmp)./(InputStruct_ssa.frun(nr).Nt_adj-1);
                    RegStat.det(broadcast.index,:) = statd;
                end
                if nr==1 && ~isempty(xglb)
                    xtmp = zscore(xglb);
                    RegCorr2.glb(:,:,broadcast.index) = (ztmp*xtmp)./(InputStruct_ssa.frun(nr).Nt_adj-1);
                    RegStat.glb(broadcast.index,:) = statg;
                end
                if nr==1 && ~isempty(xmot)
                    xtmp = zscore(xmot);
                    RegCorr2.mot(:,:,broadcast.index) = (ztmp*xtmp)./(InputStruct_ssa.frun(nr).Nt_adj-1);
                    RegStat.mot(broadcast.index,:) = statm;
                end
                if nr==1 && ~isempty(xroi)
                    xtmp = zscore(xroi);
                    RegCorr2.roi(:,:,broadcast.index) = (ztmp*xtmp)./(InputStruct_ssa.frun(nr).Nt_adj-1);
                    RegStat.roi(broadcast.index,:) = statr;
                end
                % --> some regressional stats
        
                Xnoise = [xdet(:,2:end),xglb,xmot,xroi]; % full noise matrix, normed to unit length
                Xnoise = bsxfun(@rdivide,Xnoise,sqrt(sum(Xnoise.^2)));
                [ output ] = GLM_model_fmri( volmatS, [0], [Xnoise], [Xsignal], 1, 1 );
                volmatF = output.vol_denoi;
            
                % >>> Low-pass Filtering
                Step = 'LOPASS';
                if strcmpi(broadcast.PipeStruct.(Step){1},'OFF')
                    disp('skipping lopass');
                else
                    % get function handle for analysis model of interest
                    currPath=pwd;                               % get current path
                    cd(pipeline_struct.(Step).filepath);               % jump to module directory
                    pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
                    cd(currPath);                               % jump back to current path
                    % execute step:
                    volmatF = pfun( volmatF, InputStruct_ssa.TR_MSEC, broadcast.PipeStruct.(Step)(2:end) );  
                end
                
                % >>> Component-based filtering
                if strcmpi(broadcast.PipeStruct.COMPFILT{1},'OFF')
                    disp('skipping compfilt');
                else
                    error('no compfilt options enabled yet!!')
                end
            
                %% Export fully processed data as .mat file and full .nii file
                clear TMPVOL;
                for(p=1:size(volmatF,2) )
                    tmp=double(MB.img);
                    tmp(tmp>0)= volmatF(:,p);
                    TMPVOL(:,:,:,p) = tmp;
                end
                nii=VS;
                nii.img = TMPVOL;
                nii.hdr.dime.datatype = 16;
                nii.hdr.hist = VS.hdr.hist;
                nii.hdr.dime.dim(5) = size(volmatF,2);
                save_untouch_niiz(nii,[opath4f,'/func',num2str(nr),'_fullproc.nii.gz']); 
                save([opath4f,'/func',num2str(nr),'_fullproc.mat'],'volmatF');
            end
        end
    end


    if strcmpi( ParamStruct_aug.ANALYSIS, 'NONE')
        disp('no further analysis, just doing one last quick check...')
        % extract and put the runs into cell array
        clear volcel;
        for nr=1:InputStruct_ssa.N_func
           x=load([opath4f,'/func',num2str(nr),'_fullproc.mat']);
           % quick check in case mask somehow doesnt match matfile
           if size(x.volmatF,1) ~= sum(maskS(:))
                error('fully processed matfile in %s does not match functional mask! Try deleting group level folders and rerunning P2!',opath4f)
           end
        end
        disp('ok!');

    else
        disp('now doing subject-level analysis...')

        % construct relevant path, load runfiles
        opath5f = fullfile(opath4f, ParamStruct_aug.Variable_ID);
        mkdir_r(opath5f);

        if ~isempty(ParamStruct_aug.GMMASK_ANL)
            if strcmpi(ParamStruct_aug.GMMASK_ANL,'AUTO')
                gmmask_anl_file = [outpath,'/_group_level/masks/pipe_',broadcast.PipeStruct.PNAME{1},'/func_GM_mask_grp.nii'];
            else
                if exist(ParamStruct_aug.GMMASK_ANL,'file')
                    gmmask_anl_file = ParamStruct_aug.GMMASK_ANL;
                else
                    error('analysis grey matter mask %s cannot be located?',ParamStruct_aug.GMMASK_ANL)
                end
            end
            if  contains(gmmask_anl_file,'.nii')
                unix(sprintf('cp %s %s/gmmask_anl_pre.nii',gmmask_anl_file,opath5f));
                unix(sprintf('gzip %s/gmmask_anl_pre.nii',opath5f))
            elseif contains(gmmask_anl_file,'.nii.gz')
                unix(sprintf('cp %s %s/gmmask_anl_pre.nii',gmmask_anl_file,opath5f));
            else
                error('non-nifti format of analysis grey matter mask')
            end
            VGA = load_untouch_niiz(gmmask_anl_file);
            VGA.img = VGA.img .* maskS;
            save_untouch_nii(VGA,sprintf('%s/gmmask_anl.nii',opath5f));
            unix(sprintf('rm %s/gmmask_anl_pre.nii',opath5f));
            vga_vec = nifti_to_mat(VGA,MB);
            kepix_vga = find(vga_vec>0);
            maskSsub = VGA.img;
        else
            kepix_vga = 1:sum(maskS(:));
            maskSsub = maskS;
        end
        
        % extract and put the runs into cell array
        clear volcel;
        for nr=1:InputStruct_ssa.N_func
           x=load([opath4f,'/func',num2str(nr),'_fullproc.mat']);

           % quick check in case mask somehow doesnt match matfile
           if size(x.volmatF,1) ~= sum(maskS(:))
                error('fully processed matfile in %s does not match functional mask! Try deleting group level folders and rerunning P2!',opath4f)
           end

           volcel{nr} = x.volmatF(kepix_vga,:);
        end

        % if seed-based analysis, make sure to produce appropriately resampled copies.....
        if analysis_struct.uses_roifile>0
            seedpath = fullfile(outpath,InputStruct_ssa.PREFIX,'func_seeds',sprintf('subpipe_%03u',broadcast.PipeStruct.pipe_idx.P1));
            for i=1:numel(InputStruct_ssa.seed_info.contrast)
                if ~exist(sprintf('%s/%s.resam.nii.gz',seedpath,InputStruct_ssa.seed_info.contrast(i).label),'file')

                    % copy over base mask and gzip if unzipped
                    if contains(InputStruct_ssa.seed_info.contrast(i).location,'.nii.gz')
                        unix(sprintf('cp %s %s/%s.nii.gz',InputStruct_ssa.seed_info.contrast(i).location, seedpath, InputStruct_ssa.seed_info.contrast(i).label))
                    elseif contains(InputStruct_ssa.seed_info.contrast(i).location,'.nii')
                        unix(sprintf('cp %s %s/%s.nii',InputStruct_ssa.seed_info.contrast(i).location, seedpath, InputStruct_ssa.seed_info.contrast(i).label));
                        unix(sprintf('gzip %s/%s.nii', seedpath, InputStruct_ssa.seed_info.contrast(i).label));
                    end
                    unix(sprintf('3dresample -master %s/postwarp/func1_warped_mask_clean.nii.gz -input %s/%s.nii.gz -prefix %s/%s.resam.nii.gz',opath2f, seedpath,InputStruct_ssa.seed_info.contrast(i).label, seedpath,InputStruct_ssa.seed_info.contrast(i).label));
                else
                    disp('seedmap found! skipping to next');
                end
                % loading data into mats:
                VSstr = sprintf('%s/%s.resam.nii.gz',seedpath, InputStruct_ssa.seed_info.contrast(i).label); %***%
                VS = load_untouch_niiz(VSstr);
                stmp = nifti_to_mat(VS,MB); 
                if size(stmp,2)>1
                    error('seed nifti files must be 3D (single-volume)!')
                end
                InputStruct_ssa.seed_info.seedmat(:,i) = stmp(kepix_vga); 
            end
        end

        % get function handle for analysis model of interest
        currPath=pwd;                               % get current path
        cd(analysis_struct.filepath);               % jump to module directory
        pfun= str2func(analysis_struct.model_name); % get function handle
        cd(currPath);                               % jump back to current path
        
        % call function
        ParamStruct_aug.TR_MSEC = InputStruct_ssa.TR_MSEC;
        ParamStruct_aug.mask    = maskSsub;
        %
        if analysis_struct.uses_taskfile>0
            out_analysis = pfun( volcel, InputStruct_ssa.task_info, ParamStruct_aug );  
        elseif analysis_struct.uses_roifile>0
            out_analysis = pfun( volcel, InputStruct_ssa.seed_info, ParamStruct_aug );  
        else
            out_analysis = pfun( volcel, ParamStruct_aug );  
        end
        % store submask:
        out_analysis.submask = maskSsub;
        %
        save([opath5f,'/out_analysis.mat'],'out_analysis');

        fi = fieldnames( out_analysis.image );

        for i=1:numel(fi)
            clear TMPVOL;
            for p=1:size(out_analysis.image.(fi{i}),2)
                tmp=double(maskS);
                tmp(tmp>0)= out_analysis.image.(fi{i})(:,p);
                TMPVOL(:,:,:,p) = tmp;
            end
            nii=MB;
            nii.img = TMPVOL;
            nii.hdr.dime.datatype = 16;
            nii.hdr.hist = MB.hdr.hist;
            nii.hdr.dime.dim(5) = size(out_analysis.image.(fi{i}),2);
            save_untouch_niiz(nii,[opath5f,'/',fi{i},'.nii.gz']); 
        end
    end
return;
