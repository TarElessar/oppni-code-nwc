function P2_aux1( broadcast )

    outpath = broadcast.outpath;
    pipeline_struct = broadcast.pipeline_struct;

    % check existence of subject specific struct file
    if exist(fullfile( outpath,broadcast.subj,'InputStruct_ssa.mat'),'file')  
        S = load( fullfile( outpath,broadcast.subj,'InputStruct_ssa.mat') );
        InputStruct_ssa = S.InputStruct_ssa;
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

    fprintf('\n===> subj %s (%u), anat. processing...\n',broadcast.subj,broadcast.index),
    %% =======================================================================
    %%      ANAT Processing ...
    %% =======================================================================
    
    nr=1; %-fixing to run=1 for now

    % tag for rawdata in case z-clipping was performed
    if (ischar(InputStruct_ssa.arun(nr).ZCLIP_thr) && strcmpi(InputStruct_ssa.arun(nr).ZCLIP_thr,'AUTO')) || (isnumeric(InputStruct_ssa.arun(nr).ZCLIP_thr) && isfinite(InputStruct_ssa.arun(nr).ZCLIP_thr))
        zclip_tag = '_zclip';
    else
        zclip_tag = '';
    end

    % fixing orientation stuff -- adjusting for obliquity, switching to standard mni-compatible orientation 
    if ~exist(sprintf('%s/anat%u_2std.nii.gz',opath1a,nr),'file')
        unix(sprintf('3dWarp -oblique2card -prefix %s/anat%u_deob.nii.gz -cubic %s/anat%u%s.nii.gz', opath1a,nr, opath0,nr, zclip_tag)); %-wsinc5
        unix(sprintf('fslreorient2std %s/anat%u_deob.nii.gz %s/anat%u_2std.nii.gz', opath1a,nr, opath1a,nr));        
        if exist(sprintf('%s/anat%u_2std.nii.gz',opath1a,nr),'file')
            unix(sprintf('rm %s/anat%u_deob.nii.gz', opath1a,nr)); % not a really useful intermediate - delete it
        else
            error('failed to create deob/reoriented anatomical file')
        end
    end

    % >>> Spatial Masking
    Step = 'AMASK';
    if strcmpi(broadcast.PipeStruct.(Step){1},'OFF')
        error('cannot turn this step off!');
    else
        % get function handle for analysis model of interest
        currPath=pwd;                               % get current path
        cd(pipeline_struct.(Step).filepath);               % jump to module directory
        pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
        cd(currPath);                               % jump back to current path
        % execute step:
        pfun( sprintf('%s/anat%u_2std.nii.gz', opath1a,nr), ParamStruct_aug.TEMPLATE_loc, opath2a, broadcast.PipeStruct.(Step)(2:end) );  
    end

    % >>> Spatial Warping
    Step = 'AWARP';
    if strcmpi(broadcast.PipeStruct.(Step){1},'OFF')
        error('cannot turn this step off!');
    else
        % get function handle for analysis model of interest
        currPath=pwd;                               % get current path
        cd(pipeline_struct.(Step).filepath);               % jump to module directory
        pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
        cd(currPath);                               % jump back to current path
        % execute step:
        pfun( sprintf('%s/anat%u_2std.nii.gz', opath1a,nr), sprintf('%s/anatBrainMask.nii.gz',opath2a), ParamStruct_aug.TEMPLATE_loc, opath2a, broadcast.PipeStruct.(Step)(2:end) );  
    end

    % >>> Anatomic segmentation
    Step = 'ASEG';
    if strcmpi(broadcast.PipeStruct.(Step){1},'OFF')
        error('cannot turn this step off!');
    else
        % get function handle for analysis model of interest
        currPath=pwd;                               % get current path
        cd(pipeline_struct.(Step).filepath);               % jump to module directory
        pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
        cd(currPath);                               % jump back to current path
        % execute step:
        pfun( sprintf('%s/anat_procss.nii.gz',opath2a), sprintf('%s/anatBrainMask.nii.gz',opath2a), opath3a, broadcast.PipeStruct.(Step)(2:end) );  
    end
    
    % extra step: warping the anatomical segmentations into template space
    tisslist = {'CSF','GM','WM'}; % list tissues in increasing order of T1 intensity
    for i=1:3
        if ~exist(sprintf('%s/anat_seg_%s_warped.nii.gz',opath3a,tisslist{i}),'file')
            if contains(broadcast.PipeStruct.AWARP{1},'AF') %afni-styles warp
                unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -source %s/anat_seg_%s.nii.gz -nwarp "%s/anatQQ_WARP.nii.gz %s/anatQQ.aff12.1D" -prefix %s/anat_seg_%s_warped.nii.gz',...
                    opath2a,opath3a,tisslist{i},opath2a,opath2a,opath3a,tisslist{i}));
            elseif contains(broadcast.PipeStruct.AWARP{1},'AN') %ants-styles warp
                unix(sprintf('antsApplyTransforms -d 3 -i %s/anat_seg_%s.nii.gz -r %s/anat_warped.nii.gz -n linear -t %s/anatQQ_WARP.nii.gz -t %s/anatQQ_GenericAffine.mat -o %s/anat_seg_%s_warped.nii.gz',...
                    opath3a,tisslist{i}, opath2a, opath2a,opath2a, opath3a,tisslist{i}));
            else
                error('unrecognized warping style!')
            end
        else
            disp('skipping tissue seg warping...')
        end
    end

    fprintf('\n===> subj %s (%u/%u), physio. processing...\n',broadcast.subj,broadcast.index,numel(subject_list)),
    %% =======================================================================
    %%      PHYSIO Processing ...
    %% =======================================================================

    fprintf('\n===> phys-proc. now on subj %u/%u: %s...\n',ns,numel(subject_list),broadcast.subj),
    
    disp('nothin for physio-proc so far. put in soon!')

    fprintf('\n===> subj %s (%u/%u), func. processing (part-1)...\n',broadcast.subj,broadcast.index,numel(subject_list)),
    %% =======================================================================
    %%      FUNC Processing, Block-1 ...
    %% =======================================================================

    %%%%%========== compatibility adjustments for BLOCK1 ... RICOR only
    missing_physio=0;
    for nr=1:InputStruct_ssa.N_func
        if isempty(InputStruct_ssa.frun(nr).PHYSIO_filename)
            missing_physio=missing_physio+1;
        end
    end
    if missing_physio>0 && ~strcmpi(broadcast.PipeStruct.RICOR{1},'OFF')
        warning('subject %s has %u/%u runs with missing physio. Turning RICOR off!',InputStruct_ssa.PREFIX, missing_physio,InputStruct_ssa.N_func)
        broadcast.PipeStruct.RICOR{1}='OFF';
    end
    %%%%%========== compatibility adjustments, done

    % --> NB: block-1 processing gets split into prewarp/warp/postwarp subfolders for ease of debugging
    
    % clear for variable run lengths between subj.
    clear prefix_set Funcfile_set base_set prefix_set_wrp Funcfile_set_wrp;

    for nr=1:InputStruct_ssa.N_func

        % tag for rawdata in case scan dropping was performed
        if InputStruct_ssa.frun(nr).DROP_first>0 || InputStruct_ssa.frun(nr).DROP_last>0
            drop_tag = '_drop';
        else
            drop_tag = '';
        end

        % * the steps below (INIMOT, DESPIKE, RICOR, TSHIFT) do slice-specific processing before 
        % we do any deobliqueing/alignment which destroys slice-based information
        
        % INIMOT is constructing some preliminary estimates of displacement -- finds the minimum-displacement volume for motion correction later 
        if strcmpi(ParamStruct_aug.INIMOT{1},'OP1')
            motref_0rel = inimot_OP1( sprintf('%s/func%u%s.nii.gz',opath0,nr,drop_tag), sprintf('func%u',nr), sprintf('%s/init_mot_estim',opath1f) );
        else
            error('unrecognized initial motion estimator?!')
        end

        % >>> Removing "Spikes" in fMRI data
        Step = 'DESPIKE';
        if strcmpi(broadcast.PipeStruct.(Step){1},'OFF')
            unix(sprintf('cp %s/func%u%s.nii.gz %s/prewarp/func%u_despike.nii.gz',opath0,nr,drop_tag, opath2f,nr));
        else
            % get function handle for analysis model of interest
            currPath=pwd;                               % get current path
            cd(pipeline_struct.(Step).filepath);               % jump to module directory
            pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
            cd(currPath);                               % jump back to current path
            % execute step:
            pfun( sprintf('%s/func%u%s.nii.gz',opath0,nr,drop_tag), sprintf('func%u',nr), sprintf('%s/prewarp',opath2f), motref_0rel, broadcast.PipeStruct.(Step)(2:end) );  
        end

        % >> Pipeline Step #5: "RICOR"
        if strcmpi(broadcast.PipeStruct.RICOR{1},'AF1')
                    %-- c. slice-based physio correction
                    %
                    % get the stripped-down physio data matrix
                    fid = fopen( sprintf('%s.slibase.1D',ostr4) );
                    tline = fgetl(fid);
                    kq=0; ise=0;
                    while ischar(tline) 
                        if    ( contains(tline,'# >') ) ise=1;
                        elseif( contains(tline,'# <') ) ise=2;
                        elseif( ise==1 ) % only if currently flagged on
                            kq=kq+1;
                            ricormat(kq,:) = str2num( tline );
                        end
                        tline = fgetl(fid);
                    end
                    fclose(fid);
                    % take the matrix of slicewise physio covariates and regress slice-by-slice from the despiked data
                    ricor_regress( sprintf('%s/prewarp/func%u_despike.nii.gz',opath2f,nr), ricormat, sprintf('%s/prewarp/func%u_ricor.nii.gz',opath2f,nr) );
        elseif strcmpi(broadcast.PipeStruct.RICOR{1},'OFF')
            unix(sprintf('cp %s/prewarp/func%u_despike.nii.gz %s/prewarp/func%u_ricor.nii.gz',opath2f,nr, opath2f,nr));
        else
            error('unrecognized ricorring?!')
        end

        % >>> Slice-Timing Correction
        Step = 'TSHIFT';
        if strcmpi(broadcast.PipeStruct.(Step){1},'OFF')
            unix(sprintf('cp %s/prewarp/func%u_ricor.nii.gz %s/prewarp/func%u_tshift.nii.gz',opath2f,nr, opath2f,nr));
        else
            % get function handle for analysis model of interest
            currPath=pwd;                               % get current path
            cd(pipeline_struct.(Step).filepath);               % jump to module directory
            pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
            cd(currPath);                               % jump back to current path
            % execute step:
            pfun( sprintf('%s/prewarp/func%u_ricor.nii.gz',opath2f,nr), sprintf('func%u',nr), sprintf('%s/prewarp',opath2f), InputStruct_ssa.TPATTERN, broadcast.PipeStruct.(Step)(2:end) );  
        end

        % fixing orientation stuff -- adjusting for obliquity, switching to standard mni-compatible orientation 
        if ~exist(sprintf('%s/prewarp/func%u_2std.nii.gz',opath2f,nr),'file')
            unix(sprintf('3dWarp -oblique2card -prefix %s/prewarp/func%u_deob.nii.gz -wsinc5 %s/prewarp/func%u_tshift.nii.gz' , opath2f,nr, opath2f,nr));
            unix(sprintf('fslreorient2std %s/prewarp/func%u_deob.nii.gz %s/prewarp/func%u_2std.nii.gz',opath2f,nr,opath2f,nr));
            if exist(sprintf('%s/prewarp/func%u_2std.nii.gz',opath2f,nr),'file')
                unix(sprintf('rm %s/prewarp/func%u_deob.nii.gz', opath2f,nr)); % not really useful intermediate - delete it
            else
                error('failed to create deob/reoriented functional file')
            end
        else
            disp('skipping alignment prep...')
        end

        % store fields for FWARP step later on
        prefix_set{nr} = sprintf('func%u',nr);
        Funcfile_set{nr} = sprintf('%s/prewarp/func%u_2std.nii.gz',opath2f,nr);
        base_set(nr) = motref_0rel;
        % and for SMOTTINGT
        prefix_set_wrp{nr} = sprintf('func%u_warped',nr);
        Funcfile_set_wrp{nr} = sprintf('%s/postwarp/func%u_warped.nii.gz',opath2f,nr);
    end

    % pre-specifying some input/output directories
    odir1 = sprintf('%s/warp',opath2f);
    odir2 = sprintf('%s/postwarp',opath2f);
    Anatloc = opath2a;

    % >>> Functional Warping
    Step = 'FWARP';
    if strcmpi(broadcast.PipeStruct.(Step){1},'OFF')
        error('cannot turn this step off!');
    else
        % get function handle for analysis model of interest
        currPath=pwd;                               % get current path
        cd(pipeline_struct.(Step).filepath);               % jump to module directory
        pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
        cd(currPath);                               % jump back to current path
        % execute step:
        pfun( Funcfile_set, prefix_set, odir1, odir2, base_set, Anatloc, broadcast.PipeStruct.(Step)(2:end) );  
    end

    % >>> Spatial Smoothing
    Step = 'SMOOTH';
    if strcmpi(broadcast.PipeStruct.(Step){1},'OFF')
        error('cannot turn this step off!');
    else
        % get function handle for analysis model of interest
        currPath=pwd;                               % get current path
        cd(pipeline_struct.(Step).filepath);               % jump to module directory
        pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
        cd(currPath);                               % jump back to current path
        % execute step:
        for ni=1:numel(Funcfile_set_wrp)
            pfun( Funcfile_set_wrp{ni}, prefix_set_wrp{ni}, odir2, broadcast.PipeStruct.(Step)(2:end) );  
        end
    end

    % extra step: make mask, mean, sd maps --> run-1 mean and sd used for creating group-level maps. For runs >1, just kept for qc purposes
    for nr=1:InputStruct_ssa.N_func
        if ~exist(sprintf('%s/postwarp/func%u_warped_mask.nii.gz',opath2f,nr),'file') || ...
           ~exist(sprintf('%s/postwarp/func%u_warped_tav.nii.gz',opath2f,nr),'file')  || ...
           ~exist(sprintf('%s/postwarp/func%u_warped_tsd.nii.gz',opath2f,nr),'file') 

            unix(sprintf('3dAutomask -prefix %s/postwarp/func%u_warped_mask.nii.gz %s/postwarp/func%u_warped.nii.gz',opath2f,nr,opath2f,nr));
            unix(sprintf('3dTstat -mean  -prefix %s/postwarp/func%u_warped_tav.nii.gz %s/postwarp/func%u_warped.nii.gz',opath2f,nr,opath2f,nr));
            unix(sprintf('3dTstat -stdev -prefix %s/postwarp/func%u_warped_tsd.nii.gz %s/postwarp/func%u_warped.nii.gz',opath2f,nr,opath2f,nr));
        else
            disp('func mask done - skipping ahead.')
        end
    end
    % extra step: tidied up functional mask using anatomical data - for later group mask construction
    if ~exist( sprintf('%s/postwarp/func%u_warped_mask_clean.nii.gz',opath2f,nr),'file')
        nr=1; % run-1 only
        if ~exist(sprintf('%s/warp/anat_warped_rs.nii.gz',opath2f),'file') || ...
           ~exist(sprintf('%s/warp/anat_warped_rs_mask.nii.gz',opath2f),'file') || ...
           ~exist(sprintf('%s/postwarp/func%u_warped_mask_clean.nii.gz',opath2f,nr),'file')

            unix(sprintf('3dresample -master %s/postwarp/func%u_warped_mask.nii.gz -input %s/anat_warped.nii.gz -prefix %s/warp/anat_warped_rs.nii.gz',opath2f,nr,opath2a,opath2f));
            unix(sprintf('3dmask_tool -dilate_input 5 -5 -fill_holes -input %s/warp/anat_warped_rs.nii.gz -prefix %s/warp/anat_warped_rs_mask.nii.gz',opath2f,opath2f))
            unix(sprintf('3dmask_tool -input %s/postwarp/func%u_warped_mask.nii.gz %s/warp/anat_warped_rs_mask.nii.gz -inter -prefix %s/postwarp/func%u_warped_mask_clean.nii.gz',opath2f,nr,opath2f,opath2f,nr))
        else
            disp('clean func mask done - skipping ahead.')
        end
    else
        disp('skipping newspace masking...')
    end
    % extra step: resampling tissue segmentations into functional space
    tisslist = {'CSF','GM','WM'}; % tissues in increasing order of T1 intensity
    for i=1:3
        if ~exist( sprintf('%s/anat_seg_%s_resam.nii.gz',opath3f,tisslist{i}),'file')
            unix(sprintf('3dresample -master %s/postwarp/func1_warped_mask_clean.nii.gz -input %s/anat_seg_%s_warped.nii.gz -prefix %s/anat_seg_%s_resam.nii.gz',...
                opath2f,opath3a,tisslist{i},opath3f,tisslist{i}));
        else
            disp('skipping tissue seg warping...')
        end
    end

% %     % SCRATCH -- extra step: effex maps
% %     mkdir(sprintf('%s/effex',opath2f));
% %     if ns==20
% %         return;
% %     end
% %     %
% %     V1 = load_untouch_niiz(sprintf('%s/func%u%s.nii',opath0,nr,drop_tag));
% %     V2 = load_untouch_niiz(sprintf('%s/prewarp/func%u_despike.nii',opath2f,nr));
% %     map = kurtosis(double(V2.img),0,4)-kurtosis(double(V1.img),0,4);
% %     map(~isfinite(map))=eps;
% %     mosaic_viewer( map, 3, [], [], 'jet', 1 )
return;
