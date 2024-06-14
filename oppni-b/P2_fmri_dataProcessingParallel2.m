function P2_fmri_dataProcessing( inputfile, pipelinefile, paramlist, outpath, big_skip, mask_subj_idxes, subj_subset, num_cores )
%
% =========================================================================
% P2_FMRI_DATAPROCESSING: this script should be run after the "P0" and "P1"
% steps of the bold pipeline. It does the actual data processing!
% =========================================================================
%
% Syntax:
%
%     P2_dataProcessing( inputfile, pipelinefile, paramlist, outpath, big_skip, mask_subj_idxes )
%
% Input:
%      inputfile : string, giving name of input textfile listing data to process
%      pipelinefiles : string, giving name of pipeline textfile specifying which pipeline steps to apply to the data
%      paramlist : string, giving name of parameter file specifying special arguments to use in pipelines
%      outpath : string, specifying destination directory for processed outputs
%      big_skip : completely bypasses anatomical processing checks / most of functional processing
%                 0=do not skip, 1=do skip. ONLY USE IF YOU KNOW WHAT YOU ARE DOING!
%      mask_subj_idxes : 
%     
%          if mask_subj_idxes=[], it will take all subjects in your input file
%          if mask_subj_idxes=numeric vector, this is will pull corresponding subjects from lines of the input file
%          if mask_subj_idxes=matfile, it will load the file and search for a "mask_subj_idxes" numeric vector, treated as above
%                                      if no such field, will instead seek "importpath" structure
%                                      organized as:
%                                         importpath.brain_maps
%                                         importpath.masks
%                                         importpath.parcellations
%                                      *this will set to IMPORT_MASKFILES --> copies over contents of brain_maps/masks/parcellations
%

% declaring path
CODE_PATH = fileparts(which('P0_fmri_populateDirectories.m'));
if CODE_PATH(end)~='/'
    CODE_PATH = [CODE_PATH '/'];
end

if nargin<5
    big_skip=0;
end
if nargin<6
    mask_subj_idxes=[];
end
if nargin<7
    subj_subset=0;
end
if nargin<8
    num_cores=0;
else
    parpool(num_cores);
end

% initializing structure and running checkes
[subject_list, InputStruct_aug, PipeStruct_aug, ParamStruct_aug] = P0_fmri_populateDirectories( inputfile, pipelinefile, paramlist, outpath );

% now augmenting outpath... do this after P0!
outpath = fullfile(outpath,'fmri_proc'); % subdir should be fmri_proc
if ~exist(outpath,'dir') error('fmri proc directory does not exist!'); end
e=dir([outpath,'*']); % dir to get absolute
outpath = fullfile(e.folder,e.name); % convert to absolute
% check for missing files from previous step too
File_Existence_Checker_fmri(InputStruct_aug,outpath,1); 

% check validity of analysis model etc. --> carry forward information about the analysis too!
analysis_struct = check_fmri_analysis_model( ParamStruct_aug.ANALYSIS );
% also, modify "number of components" field if more than one contrast is specified
if( ~isempty(strfind(ParamStruct_aug.CONTRAST,',' )) )
    analysis_struct.num_comp = 'multi_component';
end
% % % check validity of pipelines --> WILL SOON DO MORE??
pipeline_struct = check_fmri_processing_model( PipeStruct_aug );

% Now configuring param defaults if unspecified
%
if ~isfield(ParamStruct_aug,'INIMOT')
    ParamStruct_aug.INIMOT={'OP1',[]};
end
if ~isfield(ParamStruct_aug,'ROIMASK')
    ParamStruct_aug.ROIMASK={'OP1',[]};
end
if ~isfield(ParamStruct_aug,'TRUNC_ANL')
    ParamStruct_aug.TRUNC_ANL=[];
end
if ~isfield(ParamStruct_aug,'GMMASK_ANL')
    ParamStruct_aug.GMMASK_ANL=[];
end


% list of subjects for constructing group masks
IMPORT_MASKFILES=0;
if isempty(mask_subj_idxes)
    disp('using all subj in current input file for mask construction!')
    mask_subj_idxes = 1:numel(subject_list);
elseif ischar(mask_subj_idxes) 
    
    if ~exist(mask_subj_idxes,'file')
        error('mask id file not found!')
    end
    
    x = load(mask_subj_idxes);
    if isfield(x,'mask_subj_idxes')
        disp('loading list of subject rows in input file for mask construction!');
        if isempty(x.mask_subj_idxes)
            error('list of subject rows for mask-making is empty')
        elseif numel(x.mask_subj_idxes)<10
            warning('not a lot of subjects for group mask-making ... might be unstable')
        end
        mask_subj_idxes = x.mask_subj_idxes; clear x;
    elseif isfield(x,'importpath')
        IMPORT_MASKFILES=1; % case where we import files
        mask_subj_paths.brain_maps    = x.importpath.brain_maps;
        mask_subj_paths.masks         = x.importpath.masks;
        mask_subj_paths.parcellations = x.importpath.parcellations;
        clear mask_subj_idxes x;
    else
        error('mask_subj_idxes matfile contains information in unrecognized format...')
    end
elseif ~isnumeric(mask_subj_idxes)
    error('unrecognized mask id format?')
else
    disp('using numeric list of subj values for mask construction!')
end

if IMPORT_MASKFILES==0
    % store information about masking sublist...
    subject_list_formask = subject_list(mask_subj_idxes);
    % consistency checking
    if exist([outpath,'/_group_level/pipe_',PipeStruct_aug.PNAME{1},'_mask_subj_idxes.mat'],'file')
        x=load([outpath,'/_group_level/pipe_',PipeStruct_aug.PNAME{1},'_mask_subj_idxes.mat']);
        if     ~isempty( setdiff(subject_list_formask,x.subject_list_formask) ) 
            error('custom subject list for masking :: subjects in new list not present in old! delete group level folders if you want to update!')
        elseif ~isempty( setdiff(x.subject_list_formask,subject_list_formask) )
            error('custom subject list for masking :: subjects not in new list that are present in old! delete group level folders if you want to update!')
        else
            disp('custom subject list for masking :: list is consistent with old one ... continuing without modification!')
        end        
    end
    % re-saving, in case indexes need updating ( they may change, as long as subject prefix list doesnt )
    save([outpath,'/_group_level/pipe_',PipeStruct_aug.PNAME{1},'_mask_subj_idxes.mat'],'mask_subj_idxes','subject_list_formask');
else
    % store information about importing paths + consistency checking...
    if exist([outpath,'/_group_level/pipe_',PipeStruct_aug.PNAME{1},'_mask_subj_paths.mat'],'file')
        x=load([outpath,'/_group_level/pipe_',PipeStruct_aug.PNAME{1},'_mask_subj_paths.mat']);
        mismatchpth=0; mistag=[];
        if ~strcmp(mask_subj_paths.brain_maps,x.mask_subj_paths.brain_maps)
            mismatchpth=mismatchpth+1;
            mistag = [mistag, ', brain_maps'];
        elseif ~strcmp(mask_subj_paths.masks,x.mask_subj_paths.masks)
            mismatchpth=mismatchpth+1;
            mistag = [mistag, ', brain_maps'];
        elseif ~strcmp(mask_subj_paths.parcellations,x.mask_subj_paths.parcellations)
            mismatchpth=mismatchpth+1;
            mistag = [mistag, ', brain_maps'];
        end
        if mismatchpth>0
            error('custom imported files for masking :: %s path(s) differ from old! delete group level folders if you want to update!', mistag(3:end));
        else
            disp('custom imported files for masking :: paths are consistent with old ones ... continuing without modification!')
        end
    else
        save([outpath,'/_group_level/pipe_',PipeStruct_aug.PNAME{1},'_mask_subj_paths.mat'],'mask_subj_paths');
    end
    % --> import everything immediately!
    unix(sprintf('cp %s/*.nii* %s/_group_level/brain_maps/pipe_%s',mask_subj_paths.brain_maps,outpath,PipeStruct_aug.PNAME{1}))
    unix(sprintf('cp %s/*.nii* %s/_group_level/masks/pipe_%s',mask_subj_paths.masks,outpath,PipeStruct_aug.PNAME{1}))
    unix(sprintf('cp %s/*.nii* %s/_group_level/parcellations/pipe_%s',mask_subj_paths.parcellations,outpath,PipeStruct_aug.PNAME{1}))
end
% .to reset this part, need to delete group_level folder + strip out P2 processed data! 

% subject listing
if isempty(subj_subset) || subj_subset==0 || subj_subset<0
    subj_list_for_proc = 1:numel(subject_list);
else
    subj_list_for_proc = subj_subset;
end

parfor (ns=1:subj_list_for_proc, num_cores) % step through anat-proc, func-proc (block-1)
    broadcast_aux1 = struct;
    broadcast_aux1.index = ns;
    broadcast_aux1.subj = subject_list{ns};
    broadcast_aux1.PipeStruct = PipeStruct_aug;
    P2_aux1(broadcast_aux1);
end

%% Checkpoint

% check for all completed base-proccing of participants before masking can
% start (in case only a subset were run):
nrmax = 0;
for ns=1:numel(subject_list)
    % check existence of subject specific struct file
    if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')  
        load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end
    nrmax = max([nrmax InputStruct_ssa.N_func]);
end
complet_mat = NaN*ones(numel(subject_list),nrmax);

for ns=1:numel(subject_list)
    
    % check existence of subject specific struct file
    if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')  
        load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end

    % quick formatting stuff, again assuRImes that directory structure was already constructed in "P0" pipeline step 
    opath0 = fullfile(outpath,InputStruct_ssa.PREFIX,'rawdata');
    %
    opath1a = fullfile( outpath, InputStruct_ssa.PREFIX,'anat_proc');
    opath2a = fullfile( opath1a,sprintf('subpipe_%03u',    PipeStruct_aug.pipe_idx.Warp));
    opath3a = fullfile( opath2a,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
    %
    opath1p = fullfile(outpath,InputStruct_ssa.PREFIX,'phys_proc');
    %
    opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
    opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
    opath3f = fullfile(opath2f,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
    %
    opath4f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',['pipe_',PipeStruct_aug.PNAME{1}]); 

    for nr=1:InputStruct_ssa.N_func
        complet_mat(ns,nr)=2; %--fully processed already
        if ~exist( [opath4f,'/func',num2str(nr),'_fullproc.mat'],'file') %only do runs with missing outputs
            complet_mat(ns,nr)=1; %--processed up to masking
            if ~exist(sprintf('%s/postwarp/func%u_warped_smo.nii.gz',opath2f,nr),'file')
                complet_mat(ns,nr)=0; %--not processed
            end
        end
    end
end
if sum(complet_mat(:)==0)>0
    disp('not all subjects processed enough to mask. Halting for now!');
    disp('Unfinished:')
    ix = find(sum(complet_mat==0,2)>0);
    for i=1:numel(ix)
        ix2 = find( complet_mat(ix(i),:)==0);
        for j=1:numel(ix2)
            flagd = complet_mat(ix(i),ix2(j));
            fprintf('%s -- run %u.\n',subject_list{ix(i)},ix2(j));
        end
    end
    return;
end

%% INTERMEZZO-ANATOMICAL: getting group-level maps

if IMPORT_MASKFILES==0
    
    disp('now constructing Anatomical group level brain maps n masks');
    
    % constructing a group-level consensus mask, along with a probabilistic map of brain voxels
    if ~exist( [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/anat_brain_mask_grp.nii'], 'file') || ... % Brain mask
       ~exist( [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_pBRAIN_grp.nii'], 'file')
    
        for ni=1:numel(mask_subj_idxes)
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')  
                load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
            opath1a = fullfile( outpath, InputStruct_ssa.PREFIX,'anat_proc');
            opath2a = fullfile( opath1a,sprintf('subpipe_%03u',    PipeStruct_aug.pipe_idx.Warp));
            Vw = load_untouch_niiz(sprintf('%s/anatBrainMask_warped.nii.gz',opath2a)); %***%
            if ni==1
                volref = double(Vw.img);
            else
                volref = volref + double(Vw.img);
            end
        end
        volref = volref./numel(mask_subj_idxes);
        Vw.img = double(volref>0.50); % included in majority of individuals
        save_untouch_niiz(Vw,[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/anat_brain_mask_grp.nii']);
        Vw2=Vw;
        Vw2.img = volref; % prob
        Vw2.hdr.dime.bitpix=32;
        Vw2.hdr.dime.datatype=16;
        save_untouch_niiz(Vw2,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_pBRAIN_grp.nii']);
        clear Vw Vw2;
    end
    
    % constructing a group-level mean anatomical image
    if ~exist( [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_brain_grp.nii'], 'file') % Mean image
    
        for ni=1:numel(mask_subj_idxes)
        
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')  
                load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
            opath1a = fullfile( outpath, InputStruct_ssa.PREFIX,'anat_proc');
            opath2a = fullfile( opath1a,sprintf('subpipe_%03u',    PipeStruct_aug.pipe_idx.Warp));
            Vw = load_untouch_niiz(sprintf('%s/anat_warped.nii.gz',opath2a)); %***%
            if ni==1
                volref = double(Vw.img);
            else
                volref = volref + double(Vw.img);
            end
        end
        volref = volref./numel(mask_subj_idxes);
        Vw.img = volref;
        save_untouch_niiz(Vw,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_brain_grp.nii']);
        clear Vw;
    end

    % constructing group-level probabilistic tissue maps of CSF, GM, WM
    tisslist = {'CSF','GM','WM'};
    for i=1:3
        if ~exist( [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_p',tisslist{i},'_grp.nii'], 'file') % Mean image
        
            for ni=1:numel(mask_subj_idxes)
            
                % check existence of subject specific struct file
                if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')  
                    load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
                else
                    error('cannot find Input struct file for subject: %s \n');
                end
                opath1a = fullfile( outpath, InputStruct_ssa.PREFIX,'anat_proc');
                opath2a = fullfile( opath1a,sprintf('subpipe_%03u',    PipeStruct_aug.pipe_idx.Warp));
                opath3a = fullfile( opath2a,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
                Vw = load_untouch_niiz(sprintf('%s/anat_seg_%s_warped.nii.gz',opath3a,tisslist{i})); %***%
                if ni==1
                    volref = double(Vw.img);
                else
                    volref = volref + double(Vw.img);
                end
            end
            volref = volref./numel(mask_subj_idxes);
            Vw.img = volref;
            save_untouch_niiz(Vw,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_p',tisslist{i},'_grp.nii']);
            clear Vw;
        end
    end
      
    % rough "ventricular/sulcal map" for alignment checking
    if ~exist([outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/anat_sulcal_mask_grp.nii'],'file')
        Va  = load_untouch_niiz([outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_brain_grp.nii']);
        Vb  = load_untouch_niiz([outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/anat_brain_mask_grp.nii']);
        vsm = double(Va.img);
        vsi = 1 - (vsm - min(vsm(:)))./(max(vsm(:))-min(vsm(:)));
        vsi = vsi .* double( Vb.img );
        V   = Va;
        V.img = double( vsi > 0.5);% prctile(vsi(vsi>0),75));
        save_untouch_niiz(V,[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/anat_sulcal_mask_grp.nii']);
        clear Va Vb V;
    end

else
    disp('using premade/imported Anatomical group level brain maps n masks');

    % pull represntative first file
    if exist(fullfile( outpath,subject_list{1},'InputStruct_ssa.mat'),'file')  
        load( fullfile( outpath,subject_list{1},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end
    opath1a = fullfile( outpath, InputStruct_ssa.PREFIX,'anat_proc');
    opath2a = fullfile( opath1a,sprintf('subpipe_%03u',    PipeStruct_aug.pipe_idx.Warp));
    Vo = load_untouch_niiz(sprintf('%s/anat_warped.nii.gz',opath2a)); %***%
    
    % here is our checklist: subdirs and files
    chklist_a = {'masks',              'brain_maps',     'brain_maps',    'brain_maps',   'brain_maps',  'brain_maps',  'masks'};
    chklist_b = {'anat_brain_mask_grp','anat_pBRAIN_grp','anat_brain_grp','anat_pCSF_grp','anat_pGM_grp','anat_pWM_grp','anat_sulcal_mask_grp'};
      
    for k=1:numel(chklist_a)
        fileimp = [outpath,'/_group_level/',chklist_a{k},'/pipe_',PipeStruct_aug.PNAME{1},'/',chklist_b{k},'.nii'];
        % quick check to make sure everything is there
        % & quick check for compatibility with anaomical processed datas
        if exist( fileimp,'file' )
            Vi = load_untouch_niiz(fileimp);
        elseif exist( [fileimp,'.gz'],'file' )
            Vi = load_untouch_niiz([fileimp,'.gz']);
        else
            error('failed to find imported anatomical mask file %s!',fileimp);
        end

        if sum( Vo.hdr.dime.pixdim(2:5)==Vi.hdr.dime.pixdim(2:5) ) ~= 4
            error('imported masks do not match on vox-res!');
        end
        if sum( Vo.hdr.dime.dim(2:5)==Vi.hdr.dime.dim(2:5) ) ~= 4
            error('imported masks do not match on matrix-size!');
        end        

        Vo.hdr.dime.pixdim(2:5); %vox res
        Vo.hdr.dime.dim(2:5); %matx size
    end
end

%% INTERMEZZO-FUNCTIONAL: getting group-level maps

maskisnew=0;

if IMPORT_MASKFILES==0

    disp('now constructing Functional group level brain maps n masks')
    
    % creating group-level consensus brain mask for restricting analysis
    if ~exist( [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_brain_mask_grp.nii'], 'file') % Brain mask
        for ni=1:numel(mask_subj_idxes)
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
                load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
            opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
            opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
            Vw=load_untouch_niiz(sprintf('%s/postwarp/func1_warped_mask_clean.nii.gz',opath2f)); %***%
            if ni==1
                volref = double(Vw.img);
            else
                volref = volref + double(Vw.img);
            end
        end
        volref = volref./numel(mask_subj_idxes);
        Vw.img = double(volref>0.50); % included in majority of individuals
        save_untouch_niiz(Vw,[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_brain_mask_grp.nii']);
        Vw2=Vw;
        Vw2.img = volref; % prob
        Vw2.hdr.dime.bitpix=32;
        Vw2.hdr.dime.datatype=16;
        save_untouch_niiz(Vw2,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_pBRAIN_grp.nii']);
        clear Vw Vw2;

        maskisnew=1;
    end

    % creating group-level mean average epi image for QC
    if ~exist( [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_tAV_grp.nii'], 'file') % group mean tav epi map
        for ni=1:numel(mask_subj_idxes)
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
                load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
            opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
            opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
            Vw=load_untouch_niiz(sprintf('%s/postwarp/func1_warped_tav.nii.gz',opath2f)); %***%
            if ni==1
                volref = double(Vw.img);
            else
                volref = volref + double(Vw.img);
            end
        end
        volref = volref./numel(mask_subj_idxes);
        Vw.img = volref;
        save_untouch_niiz(Vw,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_tAV_grp.nii']); % group mean sd epi map
        clear Vw;

        maskisnew=1;
    end
    % creating group-level mean temporal SD map (+and variance mask!)
    if ~exist( [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_tSD_grp.nii'], 'file') % t-std mask
        for ni=1:numel(mask_subj_idxes)
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
                load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
            opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
            opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
            Vw=load_untouch_niiz(sprintf('%s/postwarp/func1_warped_tsd.nii.gz',opath2f)); %***%
            if ni==1
                volref = double(Vw.img);
            else
                volref = volref + double(Vw.img);
            end
        end
        volref = volref./numel(mask_subj_idxes);
        Vw.img = volref;
        save_untouch_niiz(Vw,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_tSD_grp.nii']); % group mean sd epi map
        clear Vw;

        maskisnew=1;
    end

    % constructing group-level probabilistic tissue maps of CSF, GM, WM
    tisslist = {'CSF','GM','WM'};
    for i=1:3
        if ~exist( [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_p',tisslist{i},'_grp.nii'], 'file') % Mean image
        
            for ni=1:numel(mask_subj_idxes)
            
                % check existence of subject specific struct file
                if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')  
                    load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
                else
                    error('cannot find Input struct file for subject: %s \n');
                end
                opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
                opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
                opath3f = fullfile(opath2f,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
                Vw=load_untouch_niiz(sprintf('%s/anat_seg_%s_resam.nii.gz',opath3f,tisslist{i})); %***%
                if ni==1
                    volref = double(Vw.img);
                else
                    volref = volref + double(Vw.img);
                end
            end
            volref = volref./numel(mask_subj_idxes);
            Vw.img = volref;
            save_untouch_niiz(Vw,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_p',tisslist{i},'_grp.nii']);
            clear Vw;

            maskisnew=1;
        end
    end

    % ROIMASK is constructing binary functional masks for roi-based regression
    pcsf_file = [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_pCSF_grp.nii'];
    pwm_file  = [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_pWM_grp.nii'];
    pgm_file  = [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_pGM_grp.nii'];
    tsd_file  = [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_tSD_grp.nii'];
    bmsk_file = [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_brain_mask_grp.nii'];
    
    % constructing binary tissue masks -- for ROI regression (And other things)
    if strcmpi(ParamStruct_aug.ROIMASK{1},'OP1')
        roimask_OP1( pcsf_file,pwm_file,pgm_file,tsd_file,bmsk_file, [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1}], ParamStruct_aug.ROIMASK(2:end) );
    else
        error('unrecognized roimask estimator?!')
    end

    % --> FOR ROIREG... special masks/parcellations, depending on what variant you are using!!!
    
    if strcmpi(PipeStruct_aug.ROIREG{1},'OP1') || strcmpi(PipeStruct_aug.ROIREG{1},'OP2') %% creating PCA-based spatial maps
    
        masklist = {'CSF','WM','tSD'};
        parcdir  = [outpath,'/_group_level/parcellations/pipe_',PipeStruct_aug.PNAME{1}];
        parcpref = [parcdir,'/__opptmp_inter_spwt'];
        if ~exist([parcdir,'/Ugrp_WM.nii'], 'file') || ~exist([parcdir,'/Ugrp_CSF.nii'], 'file') || ~exist([parcdir,'/Ugrp_tSD.nii'], 'file')
            %
            for i=1:numel(masklist)
                mkdir_r(parcpref);
                maskname = [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_',masklist{i},'_mask_grp.nii'];
                Mtmp = load_untouch_niiz(maskname);
                ucat=[];
                unix(sprintf('rm %s/volblur.nii',parcpref))
                for ni=1:numel(mask_subj_idxes)
                    [ni numel(mask_subj_idxes)],
                    % check existence of subject specific struct file
                    if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
                        load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
                    else
                        error('cannot find Input struct file for subject: %s \n');
                    end
                    % NB: only uses run-1 per subject
                    opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
                    opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
                    unix(sprintf('3dBlurInMask -input %s/postwarp/func1_warped.nii.gz -prefix %s/volblur.nii -mask %s -fwhm 6',opath2f,parcpref,maskname));
                    Vtmp = load_untouch_niiz(sprintf('%s/volblur.nii',parcpref));
                    volmat = nifti_to_mat(Vtmp,Mtmp); npc_10 = floor(0.1*size(volmat,2));
                    volmat = bsxfun(@rdivide, bsxfun(@minus,volmat,mean(volmat,2)), std(volmat,0,2)+eps);
                    [u,~,~]=svd( volmat,'econ');
                    ucat = [ucat, u(:,1:npc_10)];
                    npc_10_list(ni,1) = npc_10;
                    unix(sprintf('rm %s/volblur.nii',parcpref))
                end
                [Ugrp,~,~] = svd( ucat,'econ');
                Ugrp = Ugrp(:,1:round(median(npc_10_list)));
                save([parcdir,'/Ugrp_',masklist{i},'.mat'],'Ugrp','npc_10_list');
            
                clear TMPVOL;
                for(p=1:size(Ugrp,2) )
                    tmp=double(Mtmp.img);
                    tmp(tmp>0)= Ugrp(:,p);
                    TMPVOL(:,:,:,p) = tmp;
                end
                nii=Vtmp;
                nii.img = TMPVOL;
                nii.hdr.dime.datatype = 16;
                nii.hdr.hist = Vtmp.hdr.hist;
                nii.hdr.dime.dim(5) = size(Ugrp,2);
                save_untouch_niiz(nii,[parcdir,'/Ugrp_',masklist{i},'.nii']); 
            
                unix(sprintf('rm %s/*.nii', parcpref))
                unix(sprintf('rmdir %s', parcpref))
            end
        end
    
    elseif strcmpi(PipeStruct_aug.ROIREG{1},'OP3') %% resampling pre-made parcellations
        
        roipath = [CODE_PATH,'/reference_maps/roimask']; %% HARD-CODED for now
        e=dir([roipath,'/*.nii']);
        for i=1:numel(e)
            [~,epr,esf]=fileparts(e(i).name);
            if ~strcmp(esf,'.nii')
                error('roimask file does not seem to be in nifti format');
            else
                prefixe = ['roimask_resam_',epr];
            end
            if ~exist(sprintf('%s/%s.nii',[outpath,'/_group_level/parcellations/pipe_',PipeStruct_aug.PNAME{1}],prefixe),'file')
                unix(sprintf('3dresample -master %s -input %s -prefix %s/%s.nii -rmode NN',[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_brain_mask_grp.nii'],fullfile(e(i).folder,e(i).name), [outpath,'/_group_level/parcellations/pipe_',PipeStruct_aug.PNAME{1}],prefixe  ));
            end
        end
    elseif strcmpi(PipeStruct_aug.ROIREG{1},'OFF')
        disp('no ROIREG step - no parcellations imported!')
    else
        error('unrecognized ROIREG setting-- not sure how to configure tissue templates!')
    end

else

    disp('using premade/imported Functional group level brain maps n masks');

    % try to pull a represntative first file >> for qc checking
    if exist(fullfile( outpath,subject_list{1},'InputStruct_ssa.mat'),'file')  
        load( fullfile( outpath,subject_list{1},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end
    opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
    opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
    Vo=load_untouch_niiz(sprintf('%s/postwarp/func1_warped_tav.nii.gz',opath2f)); %***%
    
    % here is our checklist: subdirs and files
    chklist_a = {'masks','masks','masks','masks','masks','brain_maps','brain_maps'};
    chklist_b = {'func_brain_mask_grp','func_CSF_mask_grp','func_WM_mask_grp','func_tSD_mask_grp','func_GM_mask_grp','func_tAV_grp','func_tSD_grp'};
      
    for k=1:numel(chklist_a)
        fileimp = [outpath,'/_group_level/',chklist_a{k},'/pipe_',PipeStruct_aug.PNAME{1},'/',chklist_b{k},'.nii'];
        % quick check to make sure everything is there
        % & quick check for compatibility with anaomical processed datas
        if exist( fileimp,'file' )
            Vi = load_untouch_niiz(fileimp);
        elseif exist( [fileimp,'.gz'],'file' )
            Vi = load_untouch_niiz([fileimp,'.gz']);
        else
            error('failed to find imported functional mask file %s!',fileimp);
        end

        if sum( Vo.hdr.dime.pixdim(2:5)==Vi.hdr.dime.pixdim(2:5) ) ~= 4
            error('imported masks do not match on vox-res!');
        end
        if sum( Vo.hdr.dime.dim(2:5)==Vi.hdr.dime.dim(2:5) ) ~= 4
            error('imported masks do not match on matrix-size!');
        end        

        Vo.hdr.dime.pixdim(2:5); %vox res
        Vo.hdr.dime.dim(2:5); %matx size
    end

    maskisnew = 1; % for now, default is to redo if you are using an imported mask set
end

% pre-load mask for blck-2
MBstr = ([outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_brain_mask_grp.nii']);
MB = load_untouch_niiz(MBstr);
maskS = double(MB.img);
% declare empty struct --> nuisance regressor correlation maps
RegCorr2.det=[];
RegCorr2.glb=[];
RegCorr2.mot=[];
RegCorr2.roi=[];

parfor (ns=1:subj_list_for_proc, num_cores)  % step through func-proc (block-2)
    broadcast_aux2 = struct;
    broadcast_aux2.index = ns;
    broadcast_aux2.subj = subject_list{ns};
    broadcast_aux2.PipeStruct = PipeStruct_aug;
    P2_aux2(broadcast_aux2);
end

% TODO: add functionality for this back in
%save([outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/regstat.mat'],'RegCorr2');

disp('funxionale block-2 done');
