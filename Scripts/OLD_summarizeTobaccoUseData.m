%% Script to load in and organize substance use data from pHCP

function [options,data] = OLD_summarizeTobaccoUseData(options)

%% options
if ~exist('options','var')
    options = [];
end
if ~isfield(options,'displayFigs_stats')
    options.displayFigs_stats = 0;
end
if ~isfield(options,'mps_plotSpread_file')
    options.mps_plotSpread_file = '/home/shaw-raid1/data/psychophysics/SAiP.git/scripts/plotSpread/';
end
addpath(genpath(options.mps_plotSpread_file))
if ~isfield(options,'phcp_target_file')
    %     options.phcp_target_file = 'C:\Users\kkill\OneDrive\Desktop\GitRepos\SAiP.git\Scripts\PHCP-SubstanceUse_DATA_LABELS_2022-08-18_1012.csv';
    %     options.phcp_target_file = 'E:\GitRepos\SAiP.git\Scripts\PHCP-SubstanceUse_DATA_LABELS_2022-08-18_1012.csv';
    options.phcp_target_file = '/home/shaw-raid1/data/psychophysics/SAiP.git/scripts/PHCP-SubstanceUse_DATA_LABELS_2022-08-18_1012.csv';
end
if ~isfield(options,'syon_target_file')
    %     options.phcp_target_file = 'C:\Users\kkill\OneDrive\Desktop\GitRepos\SAiP.git\Scripts\SYON-SubstanceUse_DATA_LABELS_2023-04-19_1412.csv';
    %     options.phcp_target_file = 'E:\GitRepos\SAiP.git\Scripts\SYON-SubstanceUse_DATA_LABELS_2023-04-19_1412.csv';
    options.syon_target_file = '/home/shaw-raid1/data/psychophysics/SAiP.git/SYON-SubstanceUse_DATA_LABELS_2023-04-19_1412.csv';
end
% Grab behavioral data for pHCP SFM task
if ~isfield(options,'phcp_sfm_file')
    options.phcp_sfm_file = '/home/shaw-raid1/data/psychophysics/SFM.git/';
end
addpath(genpath(options.phcp_sfm_file))
if ~isfield(options,'phcp_sfm_struct')
    phcp_sfm_opt = [];
    phcp_sfm_opt.dateCutoff = 0;
    phcp_sfm_opt.displayFigs = 0;
    phcp_sfm_opt.run_symptom_corrs = 0;
    phcp_sfm_opt.run_demographics = 0;   % Don't want to run the demo portion as it takes forever
    phcp_sfm_opt.run_symptom_corrs = 0;   % Same w/ symptom analysis
    options.phcp_sfm_struct = summarize_SFM_results( phcp_sfm_opt );
end
% % Grab behavioral data for SYON SFM task
% if ~isfield(options,'syon_sfm_file')
%     options.syon_sfm_file = '/home/shaw-raid1/data/psychophysics/SYON.git/Bistable_Tasks/SFM_Task';
% end
% addpath(genpath(options.syon_sfm_file))
% if ~isfield(options,'syon_sfm_struct')
%     sfm_opt = [];
% %     sfm_opt.dateCutoff = 0;
% %     sfm_opt.displayFigs = 0;
% %     sfm_opt.run_symptom_corrs = 0;
%     options.syon_sfm_struct = analyze_SFM_data( sfm_opt );
% end
if ~isfield(options,'mrs_file')
    error('No options.mrs_file provided!')
    % e.g., options.mrs_file = '/home/shaw-raid1/data/MRS/processed_data/20220830_phcp_OCC_193subj_H2O_scaled.csv';
    % e.g., options.mrs_file = '/home/shaw-raid1/data/MRS/processed_data/20220830_phcp_PFC_147subj_H2O_scaled.csv';    
end
addpath(genpath(options.mrs_file))
if strcmp(options.mrs_file,...
        '/home/shaw-raid1/data/MRS/processed_data/20220830_phcp_OCC_193subj_H2O_scaled.csv')
    options.whichVoxel = 'OCC';
elseif strcmp(options.mrs_file,...
        '/home/shaw-raid1/data/MRS/processed_data/20220830_phcp_PFC_147subj_H2O_scaled.csv')
    options.whichVoxel = 'PFC';
end
if ~isfield(options,'mrs_n_col')
    options.mrs_n_col = 504; % LCM default, if using Gosia's notes = 446
end
if ~isfield(options,'mrs_header_lines')
    options.mrs_header_lines = 6; % LCM default, if using Gosia's notes = 6
end
if ~isfield(options,'mrs_struct')
    mrs_opt.target_file = options.mrs_file;
    mrs_opt.n_col = options.mrs_n_col;
    mrs_opt.header_lines = options.mrs_header_lines;
    options.mrs_struct = read_in_LCModel_results(mrs_opt);
end
if ~isfield(options,'which_metab')
    options.which_metab = {'Glu','GABA','Gln'};
    %         options.which_metab = {'Glu','GABA','Gln','GSH','NAA','NAAG'};
    %     options.which_metab = {'Glu','GABA','Gln','GSH','NAA','NAAG','Asp','Asc','MacY','Glc'};
    %     options.which_metab = {'Glu','GABA','Gln','MacY','NAA','NAAG','Glc'};
    warning('options.which_metab not specified, assuming you want to look at only Glu, Gln, and GABA...')
end

options.curDur = pwd;

%% Load in redcap substance abuse data for pHCP
[data,options] = load_pHCP_subAbuseRedcap(options);

% Combine all data across 7T runs


%% Configure the pHCP SFM data
% Grab the switch rates for  each included subject
data.phcp_sfm.dataAve = (squeeze(nanmean(...
    options.phcp_sfm_struct.illusory_task.Hz.all_data_plot(),3)));
data.phcp_sfm.dataAveLog = (squeeze(nanmean(...
    options.phcp_sfm_struct.illusory_task.Hz.all_data(),3)));
data.phcp_sfm.subjNum = [options.phcp_sfm_struct.subj_number{1};options.phcp_sfm_struct.subj_number{2};...
    options.phcp_sfm_struct.subj_number{3}];
% Also grab block data for SFM to do stats
data.phcp_sfm.dataAllLog = options.phcp_sfm_struct.illusory_task.Hz.all_data_plot;
% Remove nans from SFM list, if there are any
data.phcp_sfm.subjNum(isnan(data.phcp_sfm.dataAve(:,1)) & isnan(data.phcp_sfm.dataAve(:,2))) = [];
data.phcp_sfm.dataAve((isnan(data.phcp_sfm.dataAve(:,1)) & isnan(data.phcp_sfm.dataAve(:,2))),:) = [];
data.phcp_sfm.dataAllLog((isnan(data.phcp_sfm.dataAveLog(:,1)) & isnan(data.phcp_sfm.dataAveLog(:,2))),:,:) = [];
data.phcp_sfm.dataAveLog((isnan(data.phcp_sfm.dataAveLog(:,1)) & isnan(data.phcp_sfm.dataAveLog(:,2))),:) = [];
for iI=3:5
    % If looking at MRIA/B visits, grab the SFM from the first column
    % If looking at MRIZ visits, grab the second column
    if iI==3 || iI==4
        sfmIdx = 1;
    elseif iI==5
        sfmIdx = 2;
    end
    
    % Grab the subjNums of subjects who have tobacco and SFM data
    data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum = ...
        intersect(data.(data.runListUniqShort{iI}).subjNum,data.phcp_sfm.subjNum);
    
    % Now grab the switch rates and tobacco use data for each subject included
    data.(data.runListUniqShort{iI}).tobacco.sfm.subUse = data.(data.runListUniqShort{iI}).tobacco.cigAveDay(...
        ismember(data.(data.runListUniqShort{iI}).subjNum,data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum));
    data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate = data.phcp_sfm.dataAve(...
        ismember(data.phcp_sfm.subjNum,data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum),sfmIdx);
    data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateLog = data.phcp_sfm.dataAveLog(...
        ismember(data.phcp_sfm.subjNum,data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum),sfmIdx);
    data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateAllLog = data.phcp_sfm.dataAllLog(...
        ismember(data.phcp_sfm.subjNum,data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum),sfmIdx,:);
    data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum = data.phcp_sfm.subjNum(...
        ismember(data.phcp_sfm.subjNum,data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum));
    
    % Make group index for plotting
    data.(data.runListUniqShort{iI}).tobacco.sfm.groupIdx(...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<2000000) = 1;
    data.(data.runListUniqShort{iI}).tobacco.sfm.groupIdx(...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>=2000000 & ...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<6000000) = 2;
    data.(data.runListUniqShort{iI}).tobacco.sfm.groupIdx(...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>6000000) = 3;
end





%% Combine the SFM data structs across SYON and pHCP





%% Run sfm correlations
data.corrType = 'Spearman';
% % Correlate for B visit, using all subjs (including 0's)
% [data.sfm.corrs.tobacco.mriB_all.r data.sfm.corrs.tobacco.mriB_all.p] = ...
%     corr(data.mri7TB.tobacco.sfm.subUse',data.mri7TB.tobacco.sfm.switchRateLog,...
%     'type',data.corrType);

% Correlate combining all subjects across runs, but not including 0's
data.sfm.corrs.tobacco.allUsers.allDataSubUse = [];
data.sfm.corrs.tobacco.allUsers.allDataSwitchRate = [];
data.sfm.corrs.tobacco.allUsers.allDataSwitchRateLog = [];
data.sfm.corrs.tobacco.allUsers.allDataSwitchRateLogBlock = [];
data.sfm.corrs.tobacco.allUsers.subjNum = [];
for iI=3:5
    data.sfm.corrs.tobacco.allUsers.allDataSubUse = [data.sfm.corrs.tobacco.allUsers.allDataSubUse, ...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse(data.(data.runListUniqShort{iI}).tobacco.sfm.subUse>0)];
    data.sfm.corrs.tobacco.allUsers.allDataSwitchRate = [data.sfm.corrs.tobacco.allUsers.allDataSwitchRate; ...
        data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate(data.(data.runListUniqShort{iI}).tobacco.sfm.subUse>0)];
    data.sfm.corrs.tobacco.allUsers.allDataSwitchRateLog = [data.sfm.corrs.tobacco.allUsers.allDataSwitchRateLog; ...
        data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateLog(data.(data.runListUniqShort{iI}).tobacco.sfm.subUse>0)];
    data.sfm.corrs.tobacco.allUsers.allDataSwitchRateLogBlock = [data.sfm.corrs.tobacco.allUsers.allDataSwitchRateLogBlock; ...
        data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateAllLog(data.(data.runListUniqShort{iI}).tobacco.sfm.subUse>0,:,:)];
    
    % Make subject lists
    data.sfm.corrs.tobacco.allUsers.subjNum = [data.sfm.corrs.tobacco.allUsers.subjNum; ...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum(...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse>0)];
end
% Make group indices for plotting
data.sfm.corrs.tobacco.allUsers.groupIdx(...
    data.sfm.corrs.tobacco.allUsers.subjNum < 2000000) = 1;
data.sfm.corrs.tobacco.allUsers.groupIdx(...
    data.sfm.corrs.tobacco.allUsers.subjNum >= 2000000 &...
    data.sfm.corrs.tobacco.allUsers.subjNum < 6000000) = 2;
data.sfm.corrs.tobacco.allUsers.groupIdx(...
    data.sfm.corrs.tobacco.allUsers.subjNum > 6000000) = 3;

% Do corrs
[data.sfm.corrs.tobacco.allUsers.r data.sfm.corrs.tobacco.allUsers.p] = ...
    corr(data.sfm.corrs.tobacco.allUsers.allDataSubUse',...
    data.sfm.corrs.tobacco.allUsers.allDataSwitchRateLog,...
    'type',data.corrType);

% Do corrs for PwPP seperately
[data.sfm.corrs.tobacco.groupUsers(3).r data.sfm.corrs.tobacco.groupUsers(3).p] = ...
    corr(data.sfm.corrs.tobacco.allUsers.allDataSubUse(data.sfm.corrs.tobacco.allUsers.groupIdx==3)',...
    data.sfm.corrs.tobacco.allUsers.allDataSwitchRateLog(data.sfm.corrs.tobacco.allUsers.groupIdx==3),...
    'type',data.corrType);

%% Compare switch rates in users vs in non users
% Average switch rates for each group and the average, as a function of
% users vs non-users.
for iI = 3:5
    %% Record individual switches for later averaging
    %% Tobacco
    % Controls
    % Non log
    data.sfm.average.tobacco.allSwitchData{iI-2,1,1} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<2000000));   % User
    data.sfm.average.tobacco.allSwitchData{iI-2,2,1} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<2000000));   % Non user
    % Log
    data.sfm.average.tobacco.allSwitchDataLog{iI-2,1,1} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<2000000));   % User
    data.sfm.average.tobacco.allSwitchDataLog{iI-2,2,1} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<2000000));   % Non user
    % Log w/ blocks
    data.sfm.average.tobacco.allSwitchDataBlockLog{iI-2,1,1} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateAllLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<2000000),:,:);   % User
    data.sfm.average.tobacco.allSwitchDataBlockLog{iI-2,2,1} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateAllLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<2000000),:,:);   % Non user
    % Relatives
    % Non log
    data.sfm.average.tobacco.allSwitchData{iI-2,1,2} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>=2000000 & data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<6000000));   % User
    data.sfm.average.tobacco.allSwitchData{iI-2,2,2} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>=2000000 & data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<6000000));   % Non user
    % Log
    data.sfm.average.tobacco.allSwitchDataLog{iI-2,1,2} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>=2000000 & data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<6000000));   % User
    data.sfm.average.tobacco.allSwitchDataLog{iI-2,2,2} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>=2000000 & data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<6000000));   % Non user
    % Log w/ blocks
    data.sfm.average.tobacco.allSwitchDataBlockLog{iI-2,1,2} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateAllLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>=2000000 & data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<6000000),:,:);   % User
    data.sfm.average.tobacco.allSwitchDataBlockLog{iI-2,2,2} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateAllLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>=2000000 & data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<6000000),:,:);   % Non user
    % PwPp
    % Non log
    data.sfm.average.tobacco.allSwitchData{iI-2,1,3} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>6000000));   % User
    data.sfm.average.tobacco.allSwitchData{iI-2,2,3} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>6000000));  % Non user
    % Log
    data.sfm.average.tobacco.allSwitchDataLog{iI-2,1,3} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>6000000));   % User
    data.sfm.average.tobacco.allSwitchDataLog{iI-2,2,3} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>6000000));  % Non user
    % Log w/ blocks
    data.sfm.average.tobacco.allSwitchDataBlockLog{iI-2,1,3} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateAllLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>6000000),:,:);   % User
    data.sfm.average.tobacco.allSwitchDataBlockLog{iI-2,2,3} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateAllLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>6000000),:,:);  % Non user
    % All participants
    % Non log
    data.sfm.average.tobacco.allSwitchData{iI-2,1,4} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0));   % User
    data.sfm.average.tobacco.allSwitchData{iI-2,2,4} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0));  % Non user
    % Log
    data.sfm.average.tobacco.allSwitchDataLog{iI-2,1,4} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0));   % User
    data.sfm.average.tobacco.allSwitchDataLog{iI-2,2,4} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0));  % Non user
    % Log w/ blocks
    data.sfm.average.tobacco.allSwitchDataBlockLog{iI-2,1,4} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateAllLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0),:,:);   % User
    data.sfm.average.tobacco.allSwitchDataBlockLog{iI-2,2,4} = data.(data.runListUniqShort{iI}).tobacco.sfm.switchRateAllLog((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0),:,:);  % Non user
    
    %% Average across group and session, for users and non users
    % Tobacco
    % Controls
    data.sfm.average.tobacco.aveSwitch(iI-2,1,1) = nanmean(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<2000000)));   % User
    data.sfm.average.tobacco.aveSwitch(iI-2,2,1) = nanmean(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<2000000)));   % Non user
    data.sfm.average.tobacco.steSwitch(iI-2,1,1) = nanstd(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<2000000))) /...
        (numel(data.sfm.average.tobacco.allSwitchData{iI-2,1,1})-1);   % User
    data.sfm.average.tobacco.steSwitch(iI-2,2,1) = nanmean(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<2000000))) /...
        (numel(data.sfm.average.tobacco.allSwitchData{iI-2,2,1})-1);   % Non user
    % Relatives
    data.sfm.average.tobacco.aveSwitch(iI-2,1,2) = nanmean(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>=2000000 & data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<6000000)));   % User
    data.sfm.average.tobacco.aveSwitch(iI-2,2,2) = nanmean(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>=2000000 & data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<6000000)));   % Non user
    data.sfm.average.tobacco.steSwitch(iI-2,1,2) = nanstd(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>=2000000 & data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<6000000))) /...
        (numel(data.sfm.average.tobacco.allSwitchData{iI-2,1,2})-1);   % User
    data.sfm.average.tobacco.steSwitch(iI-2,2,2) = nanstd(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>=2000000 & data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum<6000000))) /...
        (numel(data.sfm.average.tobacco.allSwitchData{iI-2,2,2})-1);   % Non user
    % PwPp
    data.sfm.average.tobacco.aveSwitch(iI-2,1,3) = nanmean(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>6000000)));   % User
    data.sfm.average.tobacco.aveSwitch(iI-2,2,3) = nanmean(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>6000000)));   % Non user
    data.sfm.average.tobacco.steSwitch(iI-2,1,3) = nanstd(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>6000000))) /...
        (numel(data.sfm.average.tobacco.allSwitchData{iI-2,1,3})-1);   % User
    data.sfm.average.tobacco.steSwitch(iI-2,2,3) = nanstd(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0) &...
        (data.(data.runListUniqShort{iI}).tobacco.sfm.subjNum>6000000))) /...
        (numel(data.sfm.average.tobacco.allSwitchData{iI-2,2,3})-1);   % Non user
    % All participants
    data.sfm.average.tobacco.aveSwitch(iI-2,1,4) = nanmean(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0)));   % User
    data.sfm.average.tobacco.aveSwitch(iI-2,2,4) = nanmean(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0)));   % Non user
    data.sfm.average.tobacco.steSwitch(iI-2,1,4) = nanstd(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'>0))) /...
        (numel(data.sfm.average.tobacco.allSwitchData{iI-2,1,4})-1);   % User
    data.sfm.average.tobacco.steSwitch(iI-2,2,4) = nanstd(data.(data.runListUniqShort{iI}).tobacco.sfm.switchRate((...
        data.(data.runListUniqShort{iI}).tobacco.sfm.subUse'==0))) /...
        (numel(data.sfm.average.tobacco.allSwitchData{iI-2,2,4})-1);   % Non user
    
end

%% Run stats to look for group diffs
% Run individual T-tests to look at differences between users/non users for
% each group
% Ttest2
for iI=1:4   % For each of the groups and the average
    [~, data.sfm.stats.tobacco.ttest.p(iI) , ~ , ...
        data.sfm.stats.tobacco.ttest.stats{iI}] = ...
        ttest2(data.sfm.average.tobacco.allSwitchDataLog{2,1,iI},...
        data.sfm.average.tobacco.allSwitchDataLog{2,2,iI});
end

% For users/non users seperately, is there a group difference?
for iJ = 1:2
    nest = zeros(3,3);
    nest(1,2) = 1;
    if options.displayFigs_stats
        show_stats_fig = 'on';
    else
        show_stats_fig = 'off';
    end
    
    % Tobacco
    clear all_data all_subj all_group all_block subUseIdx
    % For the anova make arrays to store the correct data
    all_data = cat(1,data.sfm.average.tobacco.allSwitchDataBlockLog{2,iJ,1:3});
    all_subj = repmat([1:size(all_data,1)]',[1 size(all_data,2) size(all_data,3)]);
    if iJ == 1   % Users
        subUseIdx = data.mri7TB.tobacco.sfm.subUse>0;
    elseif iJ == 2   % Non users
        subUseIdx = data.mri7TB.tobacco.sfm.subUse==0;
    end
    all_group = [ones(numel(data.mri7TB.tobacco.sfm.subjNum(data.mri7TB.tobacco.sfm.groupIdx==1 & ...
        subUseIdx)),...
        size(all_data,2), size(all_data,3)) ; ...
        2*ones(numel(data.mri7TB.tobacco.sfm.subjNum(data.mri7TB.tobacco.sfm.groupIdx==2 & ...
        subUseIdx)),...
        size(all_data,2), size(all_data,3)) ; ...
        3*ones(numel(data.mri7TB.tobacco.sfm.subjNum(data.mri7TB.tobacco.sfm.groupIdx==3 & ...
        subUseIdx)),...
        size(all_data,2), size(all_data,3)) ];
    all_block = [];
    
    for iB = 1:size(all_data,3)
        all_block = cat(3, all_block, iB*ones(size(all_data,1), size(all_data,2)));
    end
    
    [~, data.sfm.stats.tobacco.anovan.table{iJ}] = anovan(all_data(:),{all_subj(:),...
        all_group(:),all_block(:)},'random',1,'continuous',[3],...
        'nested',nest,'model','full','varnames',{'subj','group','block'},...
        'display',show_stats_fig);
end
clear all_data all_subj all_group all_block subUseIdx

%% Load in pupil size data


%% Configure the MRS data
% Get mrs file info
if ~isfield(options,'toss_mrs_quality')
    options.toss_mrs_quality = 1; % 0 = no, 1 = yes (new way, using summarize_phcp_data_status
end
if ~isfield(options,'skip_last_n_rows')
    options.skip_last_n_rows = 4;
end

for iFile = 1:numel(options.mrs_struct.row_name)-options.skip_last_n_rows % skip last 2, mean and sd/mean
    name_idx = regexp(options.mrs_struct.row_name{iFile},'P\d\d\d\d\d\d\d');
    options.mrs_struct.subj_number(iFile,1) = str2num(options.mrs_struct.row_name{iFile}...
        (name_idx+1:name_idx+7));
    date_idx = regexp(options.mrs_struct.row_name{iFile},'\d\d\d\d\d\d\d\d');
    options.mrs_struct.date_number(iFile,1) = datenum(options.mrs_struct.row_name{iFile}...
        (date_idx:date_idx+7),'yyyymmdd');
end
clear date_idx name_idx

% toss MRS data for poor quality
if options.toss_mrs_quality
    
    summ_opt = [];
    summ_opt.displayFigs = 0;
    summ_opt.check_MRS_quality = 1;
    summ_opt.toss_date = 20170701;
    summ_opt.OCC_mrs_file = options.mrs_file; % fudge this here to use whichever
    % data file we are using for both OCC and PFC, then figure
    % out which ROI's data quality metric to use later...
    summ_opt.PFC_mrs_file = options.mrs_file;
    summ_opt.skip_last_n_rows = options.skip_last_n_rows;
    summarize_data = summarize_pHCP_data_status( summ_opt );
    
    mrs_qual_metrics = {'lw_H2O','lw_LCM','SNR','excluded'};
    high_is_good = [0 0 1 0];
    
    % Tobacco
    for iJ=3:5
        % Make subject/date list of participants that have sub use data and MRS
        % data.
        % Convert date to datenum format
        for iI=1:numel(data.(data.runListUniqShort{iJ}).date)
            dateHolder{iI} = num2str(data.(data.runListUniqShort{iJ}).date(iI));
        end
        all_subj_date = [data.(data.runListUniqShort{iJ}).subjNum; datenum(dateHolder,'yyyymmdd')']';
        clear dateHolder
        
        all_mrs_data = nan(size(all_subj_date,1), numel(options.which_metab));
        for iS = 1:size(all_subj_date,1)
            mrs_idx = find( options.mrs_struct.subj_number == all_subj_date(iS,1) & ...
                options.mrs_struct.date_number == all_subj_date(iS,2) );
            
            if ~isempty(mrs_idx) % if in MRS structure, else leave as NaN
                for iM = 1:numel(options.which_metab)
                    all_mrs_data(iS,iM) = options.mrs_struct.(options.which_metab{iM}...
                        )(mrs_idx);
                end
            end
        end
        
        toss_qual = [];
        find_out_mrs = [];
        missing_qual = [];
        for iS = 1:size(all_subj_date,1)
            subj_idx = find(summarize_data.subject_numbers == all_subj_date(iS,1) & ...
                abs(summarize_data.psy_date_numbers - all_subj_date(iS,2)) <= 14); % date within 2 weeks
            if ~isempty(subj_idx)
                toss_qual(iS,1) = squeeze( sum( summarize_data.mrs_quality_failed_binary(...
                    subj_idx,1,:),3)) >= 1;
            else
                toss_qual(iS,1) = NaN;
                missing_qual = [missing_qual ; all_subj_date(iS,:)];
            end
        end
        clear subj_idx
        
        %         if ~isempty(missing_qual)
        %             error('Houston, we have a problem...');
        %         end
        
        find_out_mrs = [find_out_mrs ; find(toss_qual == 1)];
        
        if ~isempty(find_out_mrs)
            %         MRS_subj_date(find_out_mrs,:) = [];
            warning(['Tossing ' num2str(numel(find_out_mrs)) ' data sets, for MRS quality outliers..']);
        end
        all_mrs_data(find_out_mrs,:) = NaN;
        
        data.(data.runListUniqShort{iJ}).tobacco.mrs.all_mrs_data = all_mrs_data(~isnan(all_mrs_data(:,1)),:);
        data.(data.runListUniqShort{iJ}).tobacco.mrs.subjNum = all_subj_date(~isnan(all_mrs_data(:,1)),1);
        data.(data.runListUniqShort{iJ}).tobacco.mrs.date = all_subj_date(~isnan(all_mrs_data(:,1)),2);
        
        clear all_subj_date all_mrs_data mrs_idx missing_qual find_out_mrs toss_qual
    end
    
end

% Find sub use data for subset of MRS subjects
for iJ=3:5
    data.(data.runListUniqShort{iJ}).tobacco.mrs.subUse = data.(data.runListUniqShort{iJ}).tobacco.cigAveDay(...
        ismember(data.(data.runListUniqShort{iJ}).subjNum',...
        data.(data.runListUniqShort{iJ}).tobacco.mrs.subjNum));
end


%% Run MRS correlations
% Correlate combining all subjects across runs, but not including 0's
% Make new lists including all subjects/data from each session
data.mrs.corrs.tobacco.allUsers.allDataSubUse = [];
data.mrs.corrs.tobacco.allUsers.allDataMRS = [];
data.mrs.corrs.tobacco.allUsers.subjNum = [];
for iI=3:5
    data.mrs.corrs.tobacco.allUsers.allDataSubUse = [data.mrs.corrs.tobacco.allUsers.allDataSubUse, ...
        data.(data.runListUniqShort{iI}).tobacco.mrs.subUse(data.(data.runListUniqShort{iI}).tobacco.mrs.subUse>0)];
    data.mrs.corrs.tobacco.allUsers.allDataMRS = cat(1,data.mrs.corrs.tobacco.allUsers.allDataMRS, ...
        data.(data.runListUniqShort{iI}).tobacco.mrs.all_mrs_data((data.(data.runListUniqShort{iI}).tobacco.mrs.subUse>0)',:));
    
    % Make subject lists
    data.mrs.corrs.tobacco.allUsers.subjNum = [data.mrs.corrs.tobacco.allUsers.subjNum; ...
        data.(data.runListUniqShort{iI}).tobacco.mrs.subjNum(...
        data.(data.runListUniqShort{iI}).tobacco.mrs.subUse>0)];
end
% Make group indices for plotting
data.mrs.corrs.tobacco.allUsers.groupIdx(...
    data.mrs.corrs.tobacco.allUsers.subjNum < 2000000) = 1;
data.mrs.corrs.tobacco.allUsers.groupIdx(...
    data.mrs.corrs.tobacco.allUsers.subjNum >= 2000000 &...
    data.mrs.corrs.tobacco.allUsers.subjNum < 6000000) = 2;
data.mrs.corrs.tobacco.allUsers.groupIdx(...
    data.mrs.corrs.tobacco.allUsers.subjNum > 6000000) = 3;

% Glu/GABA/Gln
for iI = 1:length(options.which_metab)
    % Tobacco
    [data.mrs.corrs.tobacco.allUsers.r(iI) data.mrs.corrs.tobacco.allUsers.p(iI)] = ...
        corr(data.mrs.corrs.tobacco.allUsers.allDataSubUse',...
        data.mrs.corrs.tobacco.allUsers.allDataMRS(:,iI),...
        'type',data.corrType);
    
    for iG=1:3   % For the 3 groups seperately
        [data.mrs.corrs.tobacco.groupUsers(iG).r(iI) ...
            data.mrs.corrs.tobacco.groupUsers(iG).p(iI)] = ...
            corr(data.mrs.corrs.tobacco.allUsers.allDataSubUse(...
            data.mrs.corrs.tobacco.allUsers.groupIdx==iG)',...
            data.mrs.corrs.tobacco.allUsers.allDataMRS(...
            data.mrs.corrs.tobacco.allUsers.groupIdx==iG,iI),...
            'type',data.corrType);
    end
end


%% Load in and correlate clinical tests from redcap w/ use data
% Can try looking at PID5 and SGI first (possibly SANS and SAPS and MARS)
% options.symptom_list = {'BACS Composite Z Scores','BPRS Total Score','spq total'};
% options.symptom_short = {'BACS','BPRS','SPQ'};
options.symptom_list = {'BACS Composite Z Scores','BPRS Total Score','BPRS Pos','BPRS Neg','BPRS Disorg','BPRS Depress','BPRS Mania','spq total'};
options.symptom_short = {'BACS','BPRS','BPRS_P','BPRS_N','BPRS_DIS','BPRS_DEP','BPRS_M','SPQ'};
% {'bprs_pos','bprs_negative','bprs_disorg','bprs_dep','bprs_mania'}
for iI = 3:5
    sym_opt.subj_number = data.(data.runListUniqShort{iI}).subjNum';
    for iJ=1:numel(data.(data.runListUniqShort{iI}).date)   % Convert to datenum format
        dateHolder{iJ} = num2str(data.(data.runListUniqShort{iI}).date(iJ));
    end
    sym_opt.date_number = datenum(dateHolder,'yyyymmdd');
    clear dateHolder
    
    for iS = 1:numel(options.symptom_list)
        sym_opt.symptom_measure = options.symptom_list{iS};
        % Call MPS get_phcp_symptoms function
        cd /home/shaw-raid1/matlab_tools/mpsCode/
        sym_out = get_phcp_symptoms(sym_opt);
        cd(options.curDur);
        data.(data.runListUniqShort{iI}).symptom.(options.symptom_short{iS}).sympDataFull = sym_out.psy_list;
    end
    clear sym_opt
end

% Combine the symptom data across runs for correlations
for iJ=1:numel(options.symptom_short)
    data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.allDataSubUse = [];
    data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.allDataSymptom = [];
    data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.subjNum = [];
    for iI=3:5
        % Find nans in the symp data to remove subjects
        clear nanHolderTob 
        nanHolderTob = ~isnan(data.(data.runListUniqShort{iI}).symptom.(options.symptom_short{iJ}).sympDataFull);
        
        % Grab the data to correlate for all users
        data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.allDataSubUse = [data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.allDataSubUse; ...
            data.(data.runListUniqShort{iI}).tobacco.cigAveDay(nanHolderTob' & data.(data.runListUniqShort{iI}).tobacco.cigAveDay>0)'];
        data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.allDataSymptom = [data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.allDataSymptom; ...
            data.(data.runListUniqShort{iI}).symptom.(options.symptom_short{iJ}).sympDataFull(nanHolderTob' & data.(data.runListUniqShort{iI}).tobacco.cigAveDay>0)];

        % Make subject lists
        data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.subjNum = [data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.subjNum; ...
            data.(data.runListUniqShort{iI}).subjNum(...
            nanHolderTob' & data.(data.runListUniqShort{iI}).tobacco.cigAveDay>0)'];
    end
    % Make group indices for plotting
    data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.groupIdx(...
        data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.subjNum < 2000000) = 1;
    data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.groupIdx(...
        data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.subjNum >= 2000000 &...
        data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.subjNum < 6000000) = 2;
    data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.groupIdx(...
        data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.subjNum > 6000000) = 3;
end


% Correlate symptom data with substance use
for iJ=1:numel(options.symptom_short)
    % Tobacco
    [data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.r ...
        data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.p] = ...
        corr(data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.allDataSubUse,...
        data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.allDataSymptom,...
        'type',data.corrType);
    
    for iI=1:3   % For the 3 groups seperately
        [data.symptom.(options.symptom_short{iJ}).corrs.tobacco.groupUsers(iI).r ...
            data.symptom.(options.symptom_short{iJ}).corrs.tobacco.groupUsers(iI).p] = ...
            corr(data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.allDataSubUse(...
            data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.groupIdx==iI),...
            data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.allDataSymptom(...
            data.symptom.(options.symptom_short{iJ}).corrs.tobacco.allUsers.groupIdx==iI),...
            'type',data.corrType);
    end
end


%% Plotting variables
% Define group colors
options.group_def_colors{1} = [0 1 0];
options.group_def_colors{2} = [0 0 1];
options.group_def_colors{3} = [1 0 0];

% %% Plot n's of users vs non users from each group within the last 7 days
% % of clinical visit. Plot number of people who haven't smoked, have smoked,
% % and didn't answer (NaNs).
%
% % Tobacco
% figure()
% figSize.blockAve.baseSize = get(0,'Screensize');   % Base size in pixels
% figSize.blockAve.aspectRatio = [10.9849 9.2814];   % Aspect ratio
% figSize.blockAve.figSize = [0 0 ...
%     figSize.blockAve.baseSize(4)*...
%     (figSize.blockAve.aspectRatio(1)/figSize.blockAve.aspectRatio(2))...
%     figSize.blockAve.baseSize(4)];   % Size/postion of fig
% hb1 = bar(1,[numel(data.Clinical.tobacco.useLast7Days(data.groupTypeUniq==1 & isnan(data.Clinical.tobacco.useLast7Days))),...
%     numel(data.Clinical.tobacco.useLast7Days(data.groupTypeUniq==1 & data.Clinical.tobacco.useLast7Days==0)),...
%     numel(data.Clinical.tobacco.useLast7Days(data.groupTypeUniq==1 & data.Clinical.tobacco.useLast7Days==1))],'FaceColor','flat');
% hold on
% hb2 = bar(2,[numel(data.Clinical.tobacco.useLast7Days(data.groupTypeUniq==2 & isnan(data.Clinical.tobacco.useLast7Days))),...
%     numel(data.Clinical.tobacco.useLast7Days(data.groupTypeUniq==2 & data.Clinical.tobacco.useLast7Days==0)),...
%     numel(data.Clinical.tobacco.useLast7Days(data.groupTypeUniq==2 & data.Clinical.tobacco.useLast7Days==1))],'FaceColor','flat');
% hb3 = bar(3,[numel(data.Clinical.tobacco.useLast7Days(data.groupTypeUniq==3 & isnan(data.Clinical.tobacco.useLast7Days))),...
%     numel(data.Clinical.tobacco.useLast7Days(data.groupTypeUniq==3 & data.Clinical.tobacco.useLast7Days==0)),...
%     numel(data.Clinical.tobacco.useLast7Days(data.groupTypeUniq==3 & data.Clinical.tobacco.useLast7Days==1))],'FaceColor','flat');
% hb4 = bar(4,[NaN, NaN, NaN],'FaceColor','flat');   % Plot nans to create values needed for legend
% for iI = 1:3
%     hb1(iI).CData(:) = [0 1*((1/3)*(4-iI)) 0];
%     hb2(iI).CData = [0 0 1*((1/3)*(4-iI))];
%     hb3(iI).CData = [1*((1/3)*(4-iI)) 0 0;];
%     hb4(iI).CData = [1*((1/3)*(4-iI)) 1*((1/3)*(4-iI)) 1*((1/3)*(4-iI))];
% end
% legend(hb4,'No data','Non smoker','Smoker')
% xlim([0.5 3.5])
% title('Participants who smoked in the last 7 days of clinical visit.')
% xticks(1:3)
% xticklabels({sprintf('%s%d%s','Controls (n=',sum(data.groupTypeUniq==1),')'),...
%     sprintf('%s%d%s','Relatives (n=',sum(data.groupTypeUniq==2),')'),...
%     sprintf('%s%d%s','Probands (n=',sum(data.groupTypeUniq==3),')')})
%
% %% Of the people who do report drinking, plot the average weekly use.
% lineOffset = .3;   % How much to offset the positioning of each line
%
% % First average across days for each group
% % Tobacco
% use7DayAveTob(1,:) = [nanmean(data.Clinical.tobacco.any(data.groupTypeUniq==1 & data.Clinical.tobacco.useLast7Days==1,:),1) ...
%     nanmean(data.Clinical.tobacco.aveDailyUse(data.groupTypeUniq==1 & data.Clinical.tobacco.useLast7Days==1))];
% use7DayAveTob(2,:) = [nanmean(data.Clinical.tobacco.any(data.groupTypeUniq==2 & data.Clinical.tobacco.useLast7Days==1,:),1) ...
%     nanmean(data.Clinical.tobacco.aveDailyUse(data.groupTypeUniq==2 & data.Clinical.tobacco.useLast7Days==1))];
% use7DayAveTob(3,:) = [nanmean(data.Clinical.tobacco.any(data.groupTypeUniq==3 & data.Clinical.tobacco.useLast7Days==1,:),1) ...
%     nanmean(data.Clinical.tobacco.aveDailyUse(data.groupTypeUniq==3 & data.Clinical.tobacco.useLast7Days==1))];
% use7DaySteTob(1,:) = [nanstd(data.Clinical.tobacco.any(data.groupTypeUniq==1 & data.Clinical.tobacco.useLast7Days==1,:),1) ./...
%     sqrt(sum(data.groupTypeUniq==1 & data.Clinical.tobacco.useLast7Days==1)) ...
%     nanstd(data.Clinical.tobacco.aveDailyUse(data.groupTypeUniq==1 & data.Clinical.tobacco.useLast7Days==1)) ./...
%     sqrt(sum(data.groupTypeUniq==1 & data.Clinical.tobacco.useLast7Days==1))];
% use7DaySteTob(2,:) = [nanstd(data.Clinical.tobacco.any(data.groupTypeUniq==2 & data.Clinical.tobacco.useLast7Days==1,:),1) ./...
%     sqrt(sum(data.groupTypeUniq==2 & data.Clinical.tobacco.useLast7Days==1)) ...
%     nanstd(data.Clinical.tobacco.aveDailyUse(data.groupTypeUniq==2 & data.Clinical.tobacco.useLast7Days==1)) ./...
%     sqrt(sum(data.groupTypeUniq==2 & data.Clinical.tobacco.useLast7Days==1))];
% use7DaySteTob(3,:) = [nanstd(data.Clinical.tobacco.any(data.groupTypeUniq==3 & data.Clinical.tobacco.useLast7Days==1,:),1) ./...
%     sqrt(sum(data.groupTypeUniq==3 & data.Clinical.tobacco.useLast7Days==1)) ...
%     nanstd(data.Clinical.tobacco.aveDailyUse(data.groupTypeUniq==3 & data.Clinical.tobacco.useLast7Days==1)) ./...
%     sqrt(sum(data.groupTypeUniq==3 & data.Clinical.tobacco.useLast7Days==1))];
%
% % Plot tobacco
% figure()
% xAxis = [1:2:size(use7DayAveTob,2)*2]';
% hb1 = bar([xAxis-lineOffset],use7DayAveTob(1,:)',.1,'FaceColor','flat');
% hold on
% hb2 = bar([xAxis],use7DayAveTob(2,:)',.1,'FaceColor','flat');
% hb3 = bar([xAxis+lineOffset],use7DayAveTob(3,:)',.1,'FaceColor','flat');
%
% % Plot bee swarm
% for iI=1:size(use7DayAveTob,2)-1
%     bee_bin_width = .1;
%     bee_spread_width = .5;
%     beePlot = plotSpread({data.Clinical.tobacco.any(data.groupTypeUniq==1 & data.Clinical.tobacco.useLast7Days==1,iI),...
%         data.Clinical.tobacco.any(data.groupTypeUniq==2 & data.Clinical.tobacco.useLast7Days==1,iI),...
%         data.Clinical.tobacco.any(data.groupTypeUniq==3 & data.Clinical.tobacco.useLast7Days==1,iI)},...
%         'binWidth', bee_bin_width,...
%         'distributionColors', {options.group_def_colors{1},...
%         options.group_def_colors{2},options.group_def_colors{3}},...
%         'xValues', [xAxis(iI)-lineOffset, xAxis(iI), xAxis(iI)+lineOffset],...
%         'spreadWidth', bee_spread_width);
%     set(beePlot{1},'MarkerSize',10)
%     hold on
% end
%
% % Plot errorbars
% errorbar(xAxis-lineOffset,use7DayAveTob(1,:)',use7DaySteTob(1,:)','.k');
% errorbar(xAxis,use7DayAveTob(2,:)',use7DaySteTob(2,:)','.k');
% errorbar(xAxis+lineOffset,use7DayAveTob(3,:)',use7DaySteTob(3,:)','.k');
% hold on
%
% hb1(1).CData = repmat(options.group_def_colors{1},[8 1]); hb1(1).FaceAlpha = .5;
% hb2(1).CData = repmat(options.group_def_colors{2},[8 1]); hb3(1).FaceAlpha = .5;
% hb3(1).CData = repmat(options.group_def_colors{3},[8 1]); hb2(1).FaceAlpha = .5;
% hold off
% xlim([0.5 16])
% title('Participants average daily tobacco use.')
% xticks(1:2:8*2)
% xticklabels({'Mon','Tue','Wed','Thu','Fri','Sat','Sun','Ave'})
% legend([hb1 hb2 hb3],{sprintf('%s%d%s','Controls (n=',sum(data.groupTypeUniq==1 & data.Clinical.tobacco.useLast7Days==1),')'),...
%     sprintf('%s%d%s','Relatives (n=',sum(data.groupTypeUniq==2 & data.Clinical.tobacco.useLast7Days==1),')'),...
%     sprintf('%s%d%s','Probands (n=',sum(data.groupTypeUniq==3 & data.Clinical.tobacco.useLast7Days==1),')')})
%

%% Plot n's of smokers vs non smokers for each of the 7T sessions.
% Plot tobacco
figure()
% Plot 7TA data
subplot(3,1,1)
hb1 = bar([.85 1.15],[sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum<2000000)==0) ...
    sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum<2000000)>0)],'FaceColor','flat');
text([.85 1.15],...
    [sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum<2000000)==0) ...
    sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum<2000000)>0)],...
    [char({'NO=','YES='}') num2str([sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum<2000000)==0) ...
    sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum<2000000)>0)]')],'vert','bottom','horiz','center');   % Plot n's above bars
hold on
hb2 = bar([1.85 2.15],[sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>=2000000 & data.mri7TA.subjNum<6000000)==0) ...
    sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>=2000000 & data.mri7TA.subjNum<6000000)>0)],'FaceColor','flat');
text([1.85 2.15],...
    [sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>=2000000 & data.mri7TA.subjNum<6000000)==0) ...
    sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>=2000000 & data.mri7TA.subjNum<6000000)>0)],...
    [char({'NO=','YES='}') num2str([sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>=2000000 & data.mri7TA.subjNum<6000000)==0) ...
    sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>=2000000 & data.mri7TA.subjNum<6000000)>0)]')],'vert','bottom','horiz','center');
hb3 = bar([2.85 3.15],[sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>6000000)==0) ...
    sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>6000000)>0)],'FaceColor','flat');
text([2.85 3.15],...
    [sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>6000000)==0) ...
    sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>6000000)>0)],...
    [char({'NO=','YES='}') num2str([sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>6000000)==0) ...
    sum(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>6000000)>0)]')],'vert','bottom','horiz','center');

hb1.FaceColor = options.group_def_colors{1};
hb2.FaceColor = options.group_def_colors{2};
hb3.FaceColor = options.group_def_colors{3};

xlim([0.5 3.5])
title('Number of tobacco users per group for all 7TA scans.')
set(gca,'xtick',1:3,'xticklabels',{sprintf('%s%d%s','Controls (n=',numel(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum<2000000)),')'),...
    sprintf('%s%d%s','Relatives (n=',numel(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>=2000000 & data.mri7TA.subjNum<6000000)),')'),...
    sprintf('%s%d%s','Probands (n=',numel(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>6000000)),')')})

% Plot 7TB data
subplot(3,1,2)
hb1 = bar([.85 1.15],[sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum<2000000)==0) ...
    sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum<2000000)>0)],'FaceColor','flat');
text([.85 1.15],...
    [sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum<2000000)==0) ...
    sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum<2000000)>0)],...
    [char({'NO=','YES='}') num2str([sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum<2000000)==0) ...
    sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum<2000000)>0)]')],'vert','bottom','horiz','center');   % Plot n's above bars
hold on
hb2 = bar([1.85 2.15],[sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>=2000000 & data.mri7TB.subjNum<6000000)==0) ...
    sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>=2000000 & data.mri7TB.subjNum<6000000)>0)],'FaceColor','flat');
text([1.85 2.15],...
    [sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>=2000000 & data.mri7TB.subjNum<6000000)==0) ...
    sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>=2000000 & data.mri7TB.subjNum<6000000)>0)],...
    [char({'NO=','YES='}') num2str([sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>=2000000 & data.mri7TB.subjNum<6000000)==0) ...
    sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>=2000000 & data.mri7TB.subjNum<6000000)>0)]')],'vert','bottom','horiz','center');
hb3 = bar([2.85 3.15],[sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>6000000)==0) ...
    sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>6000000)>0)],'FaceColor','flat');
text([2.85 3.15],...
    [sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>6000000)==0) ...
    sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>6000000)>0)],...
    [char({'NO=','YES='}') num2str([sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>6000000)==0) ...
    sum(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>6000000)>0)]')],'vert','bottom','horiz','center');

hb1.FaceColor = options.group_def_colors{1};
hb2.FaceColor = options.group_def_colors{2};
hb3.FaceColor = options.group_def_colors{3};

xlim([0.5 3.5])
title('Number of tobacco users per group for all 7TB scans.')
set(gca,'xtick',1:3,'xticklabels',{sprintf('%s%d%s','Controls (n=',numel(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum<2000000)),')'),...
    sprintf('%s%d%s','Relatives (n=',numel(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>=2000000 & data.mri7TB.subjNum<6000000)),')'),...
    sprintf('%s%d%s','Probands (n=',numel(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>6000000)),')')})


% Plot 7TZ data
subplot(3,1,3)
hb1 = bar([.85 1.15],[sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum<2000000)==0) ...
    sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum<2000000)>0)],'FaceColor','flat');
text([.85 1.15],...
    [sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum<2000000)==0) ...
    sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum<2000000)>0)],...
    [char({'NO=','YES='}') num2str([sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum<2000000)==0) ...
    sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum<2000000)>0)]')],'vert','bottom','horiz','center');   % Plot n's above bars
hold on
hb2 = bar([1.85 2.15],[sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>=2000000 & data.mri7TZ.subjNum<6000000)==0) ...
    sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>=2000000 & data.mri7TZ.subjNum<6000000)>0)],'FaceColor','flat');
text([1.85 2.15],...
    [sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>=2000000 & data.mri7TZ.subjNum<6000000)==0) ...
    sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>=2000000 & data.mri7TZ.subjNum<6000000)>0)],...
    [char({'NO=','YES='}') num2str([sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>=2000000 & data.mri7TZ.subjNum<6000000)==0) ...
    sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>=2000000 & data.mri7TZ.subjNum<6000000)>0)]')],'vert','bottom','horiz','center');
hb3 = bar([2.85 3.15],[sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>6000000)==0) ...
    sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>6000000)>0)],'FaceColor','flat');
text([2.85 3.15],...
    [sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>6000000)==0) ...
    sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>6000000)>0)],...
    [char({'NO=','YES='}') num2str([sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>6000000)==0) ...
    sum(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>6000000)>0)]')],'vert','bottom','horiz','center');

hb1.FaceColor = options.group_def_colors{1};
hb2.FaceColor = options.group_def_colors{2};
hb3.FaceColor = options.group_def_colors{3};

xlim([0.5 3.5])
title('Number of tobacco users per group for all 7TZ scans.')
set(gca,'xtick',1:3,'xticklabels',{sprintf('%s%d%s','Controls (n=',numel(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum<2000000)),')'),...
    sprintf('%s%d%s','Relatives (n=',numel(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>=2000000 & data.mri7TZ.subjNum<6000000)),')'),...
    sprintf('%s%d%s','Probands (n=',numel(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>6000000)),')')})


%% Of the people who do report average daily use in 7T runs, how much does each group use on average.
% Tobacco
figure()
% 7TA
subplot(3,1,1)
hb1 = bar(1,nanmean(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum<2000000 & data.mri7TA.tobacco.cigAveDay>0)),...
    'FaceColor',options.group_def_colors{1});
hold on
hb2 = bar(2,nanmean(data.mri7TA.tobacco.cigAveDay((data.mri7TA.subjNum>=2000000 & data.mri7TA.subjNum<6000000) & data.mri7TA.tobacco.cigAveDay>0)),...
    'FaceColor',options.group_def_colors{2});
hb3 = bar(3,nanmean(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>6000000 & data.mri7TA.tobacco.cigAveDay>0)),...
    'FaceColor',options.group_def_colors{3});
errorbar(1:3,[nanmean(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum<2000000 & data.mri7TA.tobacco.cigAveDay>0)) ...
    nanmean(data.mri7TA.tobacco.cigAveDay((data.mri7TA.subjNum>=2000000 & data.mri7TA.subjNum<6000000) & data.mri7TA.tobacco.cigAveDay>0)) ...
    nanmean(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>6000000 & data.mri7TA.tobacco.cigAveDay>0))],...
    [nanstd(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum<2000000 & data.mri7TA.tobacco.cigAveDay>0)) / ...
    (numel(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum<2000000 & data.mri7TA.tobacco.cigAveDay>0))-1),...
    nanstd(data.mri7TA.tobacco.cigAveDay((data.mri7TA.subjNum>=2000000 & data.mri7TA.subjNum<6000000) & data.mri7TA.tobacco.cigAveDay>0)) / ...
    (numel(data.mri7TA.tobacco.cigAveDay((data.mri7TA.subjNum>=2000000 & data.mri7TA.subjNum<6000000) & data.mri7TA.tobacco.cigAveDay>0))-1),...
    nanstd(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>6000000 & data.mri7TA.tobacco.cigAveDay>0)) / ...
    (numel(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>6000000 & data.mri7TA.tobacco.cigAveDay>0))-1)],...
    '.k')

% Plot bee swarm
bee_bin_width = .1;
bee_spread_width = .5;
beePlot = plotSpread({data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum<2000000 & data.mri7TA.tobacco.cigAveDay>0),...
    data.mri7TA.tobacco.cigAveDay((data.mri7TA.subjNum>=2000000 & data.mri7TA.subjNum<6000000) & data.mri7TA.tobacco.cigAveDay>0),...
    data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>6000000 & data.mri7TA.tobacco.cigAveDay>0)},...
    'binWidth', bee_bin_width,...
    'distributionColors', {[.2 .2 .2],[.2 .2 .2],[.2 .2 .2]},...
    'xValues', [1, 2, 3],...
    'spreadWidth', bee_spread_width);
% set(beePlot{1},'MarkerSize',10)

xlim([0.5 3.5])
title('Average tobacco use per group for all 7TA scans.')
set(gca,'xtick',1:3,'xticklabels',{sprintf('%s%d%s','Controls (n=',...
    numel(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum<2000000 & data.mri7TA.tobacco.cigAveDay>0)),')'),...
    sprintf('%s%d%s','Relatives (n=',...
    numel(data.mri7TA.tobacco.cigAveDay((data.mri7TA.subjNum>=2000000 & data.mri7TA.subjNum<6000000) & data.mri7TA.tobacco.cigAveDay>0)),')'),...
    sprintf('%s%d%s','Probands (n=',...
    numel(data.mri7TA.tobacco.cigAveDay(data.mri7TA.subjNum>6000000 & data.mri7TA.tobacco.cigAveDay>0)),')')})

% 7TB
subplot(3,1,2)
hb1 = bar(1,nanmean(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum<2000000 & data.mri7TB.tobacco.cigAveDay>0)),...
    'FaceColor',options.group_def_colors{1});
hold on
hb2 = bar(2,nanmean(data.mri7TB.tobacco.cigAveDay((data.mri7TB.subjNum>=2000000 & data.mri7TB.subjNum<6000000) & data.mri7TB.tobacco.cigAveDay>0)),...
    'FaceColor',options.group_def_colors{2});
hb3 = bar(3,nanmean(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>6000000 & data.mri7TB.tobacco.cigAveDay>0)),...
    'FaceColor',options.group_def_colors{3});
errorbar(1:3,[nanmean(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum<2000000 & data.mri7TB.tobacco.cigAveDay>0)) ...
    nanmean(data.mri7TB.tobacco.cigAveDay((data.mri7TB.subjNum>=2000000 & data.mri7TB.subjNum<6000000) & data.mri7TB.tobacco.cigAveDay>0)) ...
    nanmean(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>6000000 & data.mri7TB.tobacco.cigAveDay>0))],...
    [nanstd(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum<2000000 & data.mri7TB.tobacco.cigAveDay>0)) / ...
    (numel(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum<2000000 & data.mri7TB.tobacco.cigAveDay>0))-1),...
    nanstd(data.mri7TB.tobacco.cigAveDay((data.mri7TB.subjNum>=2000000 & data.mri7TB.subjNum<6000000) & data.mri7TB.tobacco.cigAveDay>0)) / ...
    (numel(data.mri7TB.tobacco.cigAveDay((data.mri7TB.subjNum>=2000000 & data.mri7TB.subjNum<6000000) & data.mri7TB.tobacco.cigAveDay>0))-1),...
    nanstd(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>6000000 & data.mri7TB.tobacco.cigAveDay>0)) / ...
    (numel(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>6000000 & data.mri7TB.tobacco.cigAveDay>0))-1)],...
    '.k')

% Plot bee swarm
bee_bin_width = .1;
bee_spread_width = .5;
beePlot = plotSpread({data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum<2000000 & data.mri7TB.tobacco.cigAveDay>0),...
    data.mri7TB.tobacco.cigAveDay((data.mri7TB.subjNum>=2000000 & data.mri7TB.subjNum<6000000) & data.mri7TB.tobacco.cigAveDay>0),...
    data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>6000000 & data.mri7TB.tobacco.cigAveDay>0)},...
    'binWidth', bee_bin_width,...
    'distributionColors', {[.2 .2 .2],[.2 .2 .2],[.2 .2 .2]},...
    'xValues', [1, 2, 3],...
    'spreadWidth', bee_spread_width);
% set(beePlot{1},'MarkerSize',10)

xlim([0.5 3.5])
title('Average tobacco use per group for all 7TB scans.')
ylabel('Number of Cigarettes / Day')
set(gca,'xtick',1:3,'xticklabels',{sprintf('%s%d%s','Controls (n=',...
    numel(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum<2000000 & data.mri7TB.tobacco.cigAveDay>0)),')'),...
    sprintf('%s%d%s','Relatives (n=',...
    numel(data.mri7TB.tobacco.cigAveDay((data.mri7TB.subjNum>=2000000 & data.mri7TB.subjNum<6000000) & data.mri7TB.tobacco.cigAveDay>0)),')'),...
    sprintf('%s%d%s','Probands (n=',...
    numel(data.mri7TB.tobacco.cigAveDay(data.mri7TB.subjNum>6000000 & data.mri7TB.tobacco.cigAveDay>0)),')')})

% 7TZ
subplot(3,1,3)
hb1 = bar(1,nanmean(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum<2000000 & data.mri7TZ.tobacco.cigAveDay>0)),...
    'FaceColor',options.group_def_colors{1});
hold on
hb2 = bar(2,nanmean(data.mri7TZ.tobacco.cigAveDay((data.mri7TZ.subjNum>=2000000 & data.mri7TZ.subjNum<6000000) & data.mri7TZ.tobacco.cigAveDay>0)),...
    'FaceColor',options.group_def_colors{2});
hb3 = bar(3,nanmean(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>6000000 & data.mri7TZ.tobacco.cigAveDay>0)),...
    'FaceColor',options.group_def_colors{3});
errorbar(1:3,[nanmean(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum<2000000 & data.mri7TZ.tobacco.cigAveDay>0)) ...
    nanmean(data.mri7TZ.tobacco.cigAveDay((data.mri7TZ.subjNum>=2000000 & data.mri7TZ.subjNum<6000000) & data.mri7TZ.tobacco.cigAveDay>0)) ...
    nanmean(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>6000000 & data.mri7TZ.tobacco.cigAveDay>0))],...
    [nanstd(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum<2000000 & data.mri7TZ.tobacco.cigAveDay>0)) / ...
    (numel(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum<2000000 & data.mri7TZ.tobacco.cigAveDay>0))-1),...
    nanstd(data.mri7TZ.tobacco.cigAveDay((data.mri7TZ.subjNum>=2000000 & data.mri7TZ.subjNum<6000000) & data.mri7TZ.tobacco.cigAveDay>0)) / ...
    (numel(data.mri7TZ.tobacco.cigAveDay((data.mri7TZ.subjNum>=2000000 & data.mri7TZ.subjNum<6000000) & data.mri7TZ.tobacco.cigAveDay>0))-1),...
    nanstd(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>6000000 & data.mri7TZ.tobacco.cigAveDay>0)) / ...
    (numel(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>6000000 & data.mri7TZ.tobacco.cigAveDay>0))-1)],...
    '.k')

% Plot bee swarm
bee_bin_width = .1;
bee_spread_width = .5;
beePlot = plotSpread({data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum<2000000 & data.mri7TZ.tobacco.cigAveDay>0),...
    data.mri7TZ.tobacco.cigAveDay((data.mri7TZ.subjNum>=2000000 & data.mri7TZ.subjNum<6000000) & data.mri7TZ.tobacco.cigAveDay>0),...
    data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>6000000 & data.mri7TZ.tobacco.cigAveDay>0)},...
    'binWidth', bee_bin_width,...
    'distributionColors', {[.2 .2 .2],[.2 .2 .2],[.2 .2 .2]},...
    'xValues', [1, 2, 3],...
    'spreadWidth', bee_spread_width);
% set(beePlot{1},'MarkerSize',10)

xlim([0.5 3.5])
title('Average tobacco use per group for all 7TZ scans.')
set(gca,'xtick',1:3,'xticklabels',{sprintf('%s%d%s','Controls (n=',...
    numel(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum<2000000 & data.mri7TZ.tobacco.cigAveDay>0)),')'),...
    sprintf('%s%d%s','Relatives (n=',...
    numel(data.mri7TZ.tobacco.cigAveDay((data.mri7TZ.subjNum>=2000000 & data.mri7TZ.subjNum<6000000) & data.mri7TZ.tobacco.cigAveDay>0)),')'),...
    sprintf('%s%d%s','Probands (n=',...
    numel(data.mri7TZ.tobacco.cigAveDay(data.mri7TZ.subjNum>6000000 & data.mri7TZ.tobacco.cigAveDay>0)),')')})


%% Plot number n of heavy vs light smokers
% Look up how others define heavy vs light...


%% Plot group diffs for SFM switch rates by group and user/non user
% Tobacco
figure; hold on
% Set font sizes
titleFontSize = 12;
axisTitleFontSize = 12;
axisLabelFontSize = 12;
statsFontSize = 12;
% Set figure size
figSize.sfmGrpDiffs.baseSize = get(0,'Screensize');   % Base size in pixels
figSize.sfmGrpDiffs.aspectRatio = [3 3];   % Aspect ratio
figSize.sfmGrpDiffs.figSize = [0 0 ...
    figSize.sfmGrpDiffs.aspectRatio];   % Size/postion of fig

subplot(1,4,1)
% Controls
hb1 = bar(1:2,[data.sfm.average.tobacco.aveSwitch(2,1,1) data.sfm.average.tobacco.aveSwitch(2,2,1)],'FaceColor','flat');
hold on
errorbar(1:2,[data.sfm.average.tobacco.aveSwitch(2,1,1) data.sfm.average.tobacco.aveSwitch(2,2,1)],...
    [data.sfm.average.tobacco.steSwitch(2,1,1),data.sfm.average.tobacco.steSwitch(2,2,1)],'.k')
hb1.FaceColor = [options.group_def_colors{1}];

% Plot bee swarm
bee_bin_width = .1;
bee_spread_width = .5;
beePlot = plotSpread({data.sfm.average.tobacco.allSwitchData{2,1,1},data.sfm.average.tobacco.allSwitchData{2,2,1}},...
    'binWidth', bee_bin_width,...
    'distributionColors', {[.2 .2 .2],[.2 .2 .2]},...
    'xValues', [1, 2],...
    'spreadWidth', bee_spread_width);
set(beePlot{1},'MarkerSize',10)

% Plot stats
text(1,.45,...
    ['t(' num2str(data.sfm.stats.tobacco.ttest.stats{1}.df) ')=' ...
    num2str(data.sfm.stats.tobacco.ttest.stats{1}.tstat)],'fontsize',statsFontSize)
text(1,.4,...
    ['p=' num2str(data.sfm.stats.tobacco.ttest.p(1))],'fontsize',statsFontSize)

ylim([0 .5])
set(gca,'YScale','log')
set(gca,'ytick',[0.005 0.01 0.025 0.05 0.1 0.25 0.5])
ylabel('Switch Rate (Hz)')
title(sprintf('%s\n%s','Average Switch Rates','Controls (B; Tobacco).'))
set(gca,'xticklabels',({sprintf('%s%d%s','Users (n=',numel(data.sfm.average.tobacco.allSwitchData{2,1,1}),')'),...
    sprintf('%s%d%s','Non users (n=',numel(data.sfm.average.tobacco.allSwitchData{2,2,1}),')')}),...
    'xticklabelrotation',45)

subplot(1,4,2)
% Relatives
hb1 = bar(1:2,[data.sfm.average.tobacco.aveSwitch(2,1,2) data.sfm.average.tobacco.aveSwitch(2,2,2)],'FaceColor','flat');
hold on
errorbar(1:2,[data.sfm.average.tobacco.aveSwitch(2,1,2) data.sfm.average.tobacco.aveSwitch(2,2,2)],...
    [data.sfm.average.tobacco.steSwitch(2,1,2),data.sfm.average.tobacco.steSwitch(2,2,2)],'.k')
hb1.FaceColor = [options.group_def_colors{2}];

% Plot bee swarm
bee_bin_width = .1;
bee_spread_width = .5;
beePlot = plotSpread({data.sfm.average.tobacco.allSwitchData{2,1,2},data.sfm.average.tobacco.allSwitchData{2,2,2}},...
    'binWidth', bee_bin_width,...
    'distributionColors', {[.2 .2 .2],[.2 .2 .2]},...
    'xValues', [1, 2],...
    'spreadWidth', bee_spread_width);
set(beePlot{1},'MarkerSize',10)

% Plot stats
text(1,.45,...
    ['t(' num2str(data.sfm.stats.tobacco.ttest.stats{2}.df) ')=' ...
    num2str(data.sfm.stats.tobacco.ttest.stats{2}.tstat)],'fontsize',statsFontSize)
text(1,.4,...
    ['p=' num2str(data.sfm.stats.tobacco.ttest.p(2))],'fontsize',statsFontSize)

ylim([0 .5])
set(gca,'YScale','log')
set(gca,'ytick',[0.005 0.01 0.025 0.05 0.1 0.25 0.5])
title(sprintf('%s\n%s','Average Switch Rates','Relatives (B; Tobacco)'))
set(gca,'xticklabels',({sprintf('%s%d%s','Users (n=',numel(data.sfm.average.tobacco.allSwitchData{2,1,2}),')'),...
    sprintf('%s%d%s','Non users (n=',numel(data.sfm.average.tobacco.allSwitchData{2,2,2}),')')}),...
    'xticklabelrotation',45)


subplot(1,4,3)
% PwPP
hb1 = bar(1:2,[data.sfm.average.tobacco.aveSwitch(2,1,3) data.sfm.average.tobacco.aveSwitch(2,2,3)],'FaceColor','flat');
hold on
errorbar(1:2,[data.sfm.average.tobacco.aveSwitch(2,1,3) data.sfm.average.tobacco.aveSwitch(2,2,3)],...
    [data.sfm.average.tobacco.steSwitch(2,1,3),data.sfm.average.tobacco.steSwitch(2,2,3)],'.k')
hb1.FaceColor = [options.group_def_colors{3}];

% Plot bee swarm
bee_bin_width = .1;
bee_spread_width = .5;
beePlot = plotSpread({data.sfm.average.tobacco.allSwitchData{2,1,3},data.sfm.average.tobacco.allSwitchData{2,2,3}},...
    'binWidth', bee_bin_width,...
    'distributionColors', {[.2 .2 .2],[.2 .2 .2]},...
    'xValues', [1, 2],...
    'spreadWidth', bee_spread_width);
set(beePlot{1},'MarkerSize',10)

% Plot stats
text(1,.45,...
    ['t(' num2str(data.sfm.stats.tobacco.ttest.stats{3}.df) ')=' ...
    num2str(data.sfm.stats.tobacco.ttest.stats{3}.tstat)],'fontsize',statsFontSize)
text(1,.4,...
    ['p=' num2str(data.sfm.stats.tobacco.ttest.p(3))],'fontsize',statsFontSize)

ylim([0 .5])
set(gca,'YScale','log')
set(gca,'ytick',[0.005 0.01 0.025 0.05 0.1 0.25 0.5])
title(sprintf('%s\n%s','Average Switch Rates','PwPP (B; Tobacco)'))
set(gca,'xticklabels',({sprintf('%s%d%s','Users (n=',numel(data.sfm.average.tobacco.allSwitchData{2,1,3}),')'),...
    sprintf('%s%d%s','Non users (n=',numel(data.sfm.average.tobacco.allSwitchData{2,2,3}),')')}),...
    'xticklabelrotation',45)

subplot(1,4,4)
% All participants
hb1 = bar(1:2,[data.sfm.average.tobacco.aveSwitch(2,1,4) data.sfm.average.tobacco.aveSwitch(2,2,4)],'FaceColor','flat');
hold on
errorbar(1:2,[data.sfm.average.tobacco.aveSwitch(2,1,4) data.sfm.average.tobacco.aveSwitch(2,2,4)],...
    [data.sfm.average.tobacco.steSwitch(2,1,4),data.sfm.average.tobacco.steSwitch(2,2,4)],'.k')
hb1.FaceColor = [.8 .8 .8];

% Plot bee swarm
bee_bin_width = .1;
bee_spread_width = .5;
beePlot = plotSpread({data.sfm.average.tobacco.allSwitchData{2,1,4},data.sfm.average.tobacco.allSwitchData{2,2,4}},...
    'binWidth', bee_bin_width,...
    'distributionColors', {[.2 .2 .2],[.2 .2 .2]},...
    'xValues', [1, 2],...
    'spreadWidth', bee_spread_width);
set(beePlot{1},'MarkerSize',10)

% Plot stats
text(1,.45,...
    ['t(' num2str(data.sfm.stats.tobacco.ttest.stats{4}.df) ')=' ...
    num2str(data.sfm.stats.tobacco.ttest.stats{4}.tstat)],'fontsize',statsFontSize)
text(1,.4,...
    ['p=' num2str(data.sfm.stats.tobacco.ttest.p(4))],'fontsize',statsFontSize)

ylim([0 .5])
set(gca,'YScale','log')
set(gca,'ytick',[0.005 0.01 0.025 0.05 0.1 0.25 0.5])
title(sprintf('%s\n%s','Average Switch Rates','All (B; Tobacco)'))
set(gca,'xticklabels',({sprintf('%s%d%s','Users (n=',numel(data.sfm.average.tobacco.allSwitchData{2,1,4}),')'),...
    sprintf('%s%d%s','Non users (n=',numel(data.sfm.average.tobacco.allSwitchData{2,2,4}),')')}),...
    'xticklabelrotation',45)

%% Plot User and Non user groups seperately for each participant group
% Tobacco
figure; hold on
% Set font sizes
titleFontSize = 12;
axisTitleFontSize = 12;
axisLabelFontSize = 12;
statsFontSize = 12;
% Set figure size
figSize.sfmGrpDiffs.baseSize = get(0,'Screensize');   % Base size in pixels
figSize.sfmGrpDiffs.aspectRatio = [3 3];   % Aspect ratio
figSize.sfmGrpDiffs.figSize = [0 0 ...
    figSize.sfmGrpDiffs.aspectRatio];   % Size/postion of fig

subplot(1,2,1)
% Users
hb1 = bar(1,[data.sfm.average.tobacco.aveSwitch(2,1,1)],'FaceColor','flat');
hold on
hb2 = bar(2,[data.sfm.average.tobacco.aveSwitch(2,1,2)],'FaceColor','flat');
hb3 = bar(3,[data.sfm.average.tobacco.aveSwitch(2,1,3)],'FaceColor','flat');
errorbar(1:3,[data.sfm.average.tobacco.aveSwitch(2,1,1) data.sfm.average.tobacco.aveSwitch(2,1,2) data.sfm.average.tobacco.aveSwitch(2,1,3)],...
    [data.sfm.average.tobacco.steSwitch(2,1,1) data.sfm.average.tobacco.steSwitch(2,1,2) data.sfm.average.tobacco.steSwitch(2,1,3)],'.k')
hb1.FaceColor = [options.group_def_colors{1}];
hb2.FaceColor = [options.group_def_colors{2}];
hb3.FaceColor = [options.group_def_colors{3}];

% Plot bee swarm
bee_bin_width = .1;
bee_spread_width = .5;
beePlot = plotSpread({data.sfm.average.tobacco.allSwitchData{2,1,1},data.sfm.average.tobacco.allSwitchData{2,1,2},data.sfm.average.tobacco.allSwitchData{2,1,3}},...
    'binWidth', bee_bin_width,...
    'distributionColors', {[.2 .2 .2],[.2 .2 .2],[.2 .2 .2]},...
    'xValues', [1, 2, 3],...
    'spreadWidth', bee_spread_width);
set(beePlot{1},'MarkerSize',10)

% Plot stats
text(1,.45,...
    ['f(' num2str(data.sfm.stats.tobacco.anovan.table{1}{2,3}) ',' ...
    num2str(data.sfm.stats.tobacco.anovan.table{1}{3,3}) ')=' ...
    num2str(data.sfm.stats.tobacco.anovan.table{1}{3,6})],'fontsize',statsFontSize)
text(1,.4,...
    ['p=' num2str(data.sfm.stats.tobacco.anovan.table{1}{3,7})],'fontsize',statsFontSize)

ylim([0 .5])
set(gca,'YScale','log')
set(gca,'ytick',[0.005 0.01 0.025 0.05 0.1 0.25 0.5])
ylabel('Switch Rate (Hz)')
title(sprintf('%s\n%s','Average Switch Rates','Users (B; Tobacco).'))
set(gca,'xticklabels',({sprintf('%s%d%s','Controls (n=',numel(data.sfm.average.tobacco.allSwitchData{2,1,1}),')'),...
    sprintf('%s%d%s','Relatives (n=',numel(data.sfm.average.tobacco.allSwitchData{2,1,2}),')'),...
    sprintf('%s%d%s','PwPP (n=',numel(data.sfm.average.tobacco.allSwitchData{2,1,3}),')')}),...
    'xticklabelrotation',45)

subplot(1,2,2)
% Non users
hb1 = bar(1,[data.sfm.average.tobacco.aveSwitch(2,2,1)],'FaceColor','flat');
hold on
hb2 = bar(2,[data.sfm.average.tobacco.aveSwitch(2,2,2)],'FaceColor','flat');
hb3 = bar(3,[data.sfm.average.tobacco.aveSwitch(2,2,3)],'FaceColor','flat');
errorbar(1:3,[data.sfm.average.tobacco.aveSwitch(2,2,1) data.sfm.average.tobacco.aveSwitch(2,2,2) data.sfm.average.tobacco.aveSwitch(2,2,3)],...
    [data.sfm.average.tobacco.steSwitch(2,2,1) data.sfm.average.tobacco.steSwitch(2,2,2) data.sfm.average.tobacco.steSwitch(2,2,3)],'.k')
hb1.FaceColor = [options.group_def_colors{1}];
hb2.FaceColor = [options.group_def_colors{2}];
hb3.FaceColor = [options.group_def_colors{3}];
% Plot bee swarm
bee_bin_width = .1;
bee_spread_width = .5;
beePlot = plotSpread({data.sfm.average.tobacco.allSwitchData{2,2,1},data.sfm.average.tobacco.allSwitchData{2,2,2},data.sfm.average.tobacco.allSwitchData{2,2,3}},...
    'binWidth', bee_bin_width,...
    'distributionColors', {[.2 .2 .2],[.2 .2 .2],[.2 .2 .2]},...
    'xValues', [1, 2, 3],...
    'spreadWidth', bee_spread_width);
set(beePlot{1},'MarkerSize',10)

% Plot stats
text(1,.45,...
    ['f(' num2str(data.sfm.stats.tobacco.anovan.table{2}{2,3}) ',' ...
    num2str(data.sfm.stats.tobacco.anovan.table{2}{3,3}) ')=' ...
    num2str(data.sfm.stats.tobacco.anovan.table{2}{3,6})],'fontsize',statsFontSize)
text(1,.4,...
    ['p=' num2str(data.sfm.stats.tobacco.anovan.table{2}{3,7})],'fontsize',statsFontSize)

ylim([0 .5])
set(gca,'YScale','log')
set(gca,'ytick',[0.005 0.01 0.025 0.05 0.1 0.25 0.5])
ylabel('Switch Rate (Hz)')
title(sprintf('%s\n%s','Average Switch Rates','Non users (B; Tobacco).'))
set(gca,'xticklabels',({sprintf('%s%d%s','Controls (n=',numel(data.sfm.average.tobacco.allSwitchData{2,2,1}),')'),...
    sprintf('%s%d%s','Relatives (n=',numel(data.sfm.average.tobacco.allSwitchData{2,2,2}),')'),...
    sprintf('%s%d%s','PwPP (n=',numel(data.sfm.average.tobacco.allSwitchData{2,2,3}),')')}),...
    'xticklabelrotation',45)


%% Plot correlations between SFM and use data.
% % Tobacco
% % Fit the data
% [poly_fit] = polyfit(data.mri7TB.tobacco.sfm.subUse', ...
%     data.mri7TB.tobacco.sfm.switchRate, 1);
% fit_x = [min(data.mri7TB.tobacco.sfm.subUse') max(data.mri7TB.tobacco.sfm.subUse')];
% fit_y = poly_fit(1).*fit_x + poly_fit(2);
% y_range = [min(data.mri7TB.tobacco.sfm.switchRate) max(data.mri7TB.tobacco.sfm.switchRate)];
% 
% figure; hold on
% % Set font sizes
% titleFontSize = 12;
% axisTitleFontSize = 12;
% axisLabelFontSize = 12;
% statsFontSize = 12;
% % Set figure size
% figSize.sfmTobCorr.baseSize = get(0,'Screensize');   % Base size in pixels
% figSize.sfmTobCorr.aspectRatio = [3 3];   % Aspect ratio
% figSize.sfmTobCorr.figSize = [0 0 ...
%     figSize.sfmTobCorr.aspectRatio];   % Size/postion of fig
% 
% plot(fit_x,fit_y,'k-','linewidth',2)
% hold on
% 
% plot(data.mri7TB.tobacco.sfm.subUse(data.mri7TB.tobacco.sfm.groupIdx==1)', ...
%     data.mri7TB.tobacco.sfm.switchRate(data.mri7TB.tobacco.sfm.groupIdx==1), ...
%     'go','MarkerFaceColor','w','linewidth',2,...
%     'MarkerSize',6)
% plot(data.mri7TB.tobacco.sfm.subUse(data.mri7TB.tobacco.sfm.groupIdx==2)', ...
%     data.mri7TB.tobacco.sfm.switchRate(data.mri7TB.tobacco.sfm.groupIdx==2), ...
%     'bo','MarkerFaceColor','w','linewidth',2,...
%     'MarkerSize',6)
% plot(data.mri7TB.tobacco.sfm.subUse(data.mri7TB.tobacco.sfm.groupIdx==3)', ...
%     data.mri7TB.tobacco.sfm.switchRate(data.mri7TB.tobacco.sfm.groupIdx==3), ...
%     'ro','MarkerFaceColor','w','linewidth',2,...
%     'MarkerSize',6)
% 
% text(fit_x(2),y_range(2),...
%     ['r = ' num2str(data.sfm.corrs.tobacco.mriB_all.r)],'fontsize',statsFontSize)
% text(fit_x(2),y_range(2)-.05,...
%     ['p = ' num2str(data.sfm.corrs.tobacco.mriB_all.p)],'fontsize',statsFontSize)
% 
% set(gca,'xlim',[0 fit_x(2)+1])
% % set(gca,'ylim',[0 y_range(2)+.1])
% ylim([0 .5])
% set(gca,'YScale','log')
% set(gca,'ytick',[0.005 0.01 0.025 0.05 0.1 0.25 0.5])
% set(gca,'xcolor','k','ycolor','k')
% xlabel('Servings of Tobacco / Day','color','k','fontsize',axisTitleFontSize)
% ylabel('Average Switch Rate (Hz)','color','k','fontsize',axisTitleFontSize)
% title(sprintf('%s\n%s','SFM Switch Rate Vs','Tobacco Use'),'fontsize',titleFontSize)
% 
% set(gcf,'Units','inches')
% set(gcf,'Position',figSize.sfmTobCorr.figSize,'color','w')
% 


%% Plot correlations between SFM and use data excluding 0s.
% Tobacco
% Fit the data
[poly_fit] = polyfit(data.sfm.corrs.tobacco.allUsers.allDataSubUse', ...
    data.sfm.corrs.tobacco.allUsers.allDataSwitchRate, 1);
fit_x = [min(data.sfm.corrs.tobacco.allUsers.allDataSubUse') max(data.sfm.corrs.tobacco.allUsers.allDataSubUse')];
fit_y = poly_fit(1).*fit_x + poly_fit(2);
y_range = [min(data.sfm.corrs.tobacco.allUsers.allDataSwitchRate) max(data.sfm.corrs.tobacco.allUsers.allDataSwitchRate)];

% Fit the data for PwPP seperatly
[poly_fit_group] = polyfit(data.sfm.corrs.tobacco.allUsers.allDataSubUse(data.sfm.corrs.tobacco.allUsers.groupIdx==3)', ...
    data.sfm.corrs.tobacco.allUsers.allDataSwitchRate(data.sfm.corrs.tobacco.allUsers.groupIdx==3), 1);
fit_x_group = [min(data.sfm.corrs.tobacco.allUsers.allDataSubUse(data.sfm.corrs.tobacco.allUsers.groupIdx==3)') ...
    max(data.sfm.corrs.tobacco.allUsers.allDataSubUse(data.sfm.corrs.tobacco.allUsers.groupIdx==3)')];
fit_y_group = poly_fit_group(1).*fit_x + poly_fit_group(2);

figure; hold on
% Set font sizes
titleFontSize = 12;
axisTitleFontSize = 12;
axisLabelFontSize = 12;
statsFontSize = 12;
% Set figure size
figSize.sfmTobCorr.baseSize = get(0,'Screensize');   % Base size in pixels
figSize.sfmTobCorr.aspectRatio = [3 3];   % Aspect ratio
figSize.sfmTobCorr.figSize = [0 0 ...
    figSize.sfmTobCorr.aspectRatio];   % Size/postion of fig

plot(fit_x,fit_y,'k-','linewidth',2)
hold on
plot(fit_x_group,fit_y_group,'r-','linewidth',2)

plot(data.sfm.corrs.tobacco.allUsers.allDataSubUse(data.sfm.corrs.tobacco.allUsers.groupIdx==1)', ...
    data.sfm.corrs.tobacco.allUsers.allDataSwitchRate(data.sfm.corrs.tobacco.allUsers.groupIdx==1), ...
    'go','MarkerFaceColor','w','linewidth',2,...
    'MarkerSize',6)
plot(data.sfm.corrs.tobacco.allUsers.allDataSubUse(data.sfm.corrs.tobacco.allUsers.groupIdx==2)', ...
    data.sfm.corrs.tobacco.allUsers.allDataSwitchRate(data.sfm.corrs.tobacco.allUsers.groupIdx==2), ...
    'bo','MarkerFaceColor','w','linewidth',2,...
    'MarkerSize',6)
plot(data.sfm.corrs.tobacco.allUsers.allDataSubUse(data.sfm.corrs.tobacco.allUsers.groupIdx==3)', ...
    data.sfm.corrs.tobacco.allUsers.allDataSwitchRate(data.sfm.corrs.tobacco.allUsers.groupIdx==3), ...
    'ro','MarkerFaceColor','w','linewidth',2,...
    'MarkerSize',6)

text(fit_x(2),y_range(2),...
    ['r = ' num2str(data.sfm.corrs.tobacco.allUsers.r)],'fontsize',statsFontSize)
text(fit_x(2),y_range(2)-.05,...
    ['p = ' num2str(data.sfm.corrs.tobacco.allUsers.p)],'fontsize',statsFontSize)
% Plot stats seperatly for PwPP
text(fit_x(2),y_range(2)-.10,...
    ['r = ' num2str(data.sfm.corrs.tobacco.groupUsers(3).r) ' (PwPP)'],'fontsize',statsFontSize)
text(fit_x(2),y_range(2)-.15,...
    ['p = ' num2str(data.sfm.corrs.tobacco.groupUsers(3).p) ' (PwPP)'],'fontsize',statsFontSize)

set(gca,'xlim',[0 fit_x(2)+1])
% set(gca,'ylim',[0 y_range(2)+.1])
ylim([0 .5])
set(gca,'YScale','log')
set(gca,'ytick',[0.005 0.01 0.025 0.05 0.1 0.25 0.5])
set(gca,'xcolor','k','ycolor','k')
xlabel('Servings of Tobacco / Day','color','k','fontsize',axisTitleFontSize)
ylabel('Average Switch Rate (Hz)','color','k','fontsize',axisTitleFontSize)
title(sprintf('%s\n%s','SFM Switch Rate Vs','Tobacco Use (only users, all scan sessions)'),'fontsize',titleFontSize)

set(gcf,'Units','inches')
set(gcf,'Position',figSize.sfmTobCorr.figSize,'color','w')


%% Plot correlations between MRS and use data excluding 0s.
% Tobacco
figure; hold on
% Set font sizes
titleFontSize = 12;
axisTitleFontSize = 12;
axisLabelFontSize = 12;
statsFontSize = 12;
% Set figure size
figSize.mrsTobCorr.baseSize = get(0,'Screensize');   % Base size in pixels
figSize.mrsTobCorr.aspectRatio = [3 3];   % Aspect ratio
figSize.mrsTobCorr.figSize = [0 0 ...
    figSize.mrsTobCorr.aspectRatio];   % Size/postion of fig

% Glu/GABA/Gln
for iI = 1:length(options.which_metab)
    % Fit the data
    [poly_fit] = polyfit(data.mrs.corrs.tobacco.allUsers.allDataSubUse', ...
        data.mrs.corrs.tobacco.allUsers.allDataMRS(:,iI), 1);
    fit_x = [min(data.mrs.corrs.tobacco.allUsers.allDataSubUse') max(data.mrs.corrs.tobacco.allUsers.allDataSubUse')];
    fit_y = poly_fit(1).*fit_x + poly_fit(2);
    y_range = [min(data.mrs.corrs.tobacco.allUsers.allDataMRS(:,iI)) max(data.mrs.corrs.tobacco.allUsers.allDataMRS(:,iI))];
    
    % Fit the data for patients w/ psychosis seperatly
    [poly_fit_group] = polyfit(data.mrs.corrs.tobacco.allUsers.allDataSubUse(data.mrs.corrs.tobacco.allUsers.groupIdx==3)', ...
        data.mrs.corrs.tobacco.allUsers.allDataMRS(data.mrs.corrs.tobacco.allUsers.groupIdx==3,iI), 1);
    fit_x_group = [min(data.mrs.corrs.tobacco.allUsers.allDataSubUse(data.mrs.corrs.tobacco.allUsers.groupIdx==3)') ...
        max(data.mrs.corrs.tobacco.allUsers.allDataSubUse(data.mrs.corrs.tobacco.allUsers.groupIdx==3)')];
    fit_y_group = poly_fit_group(1).*fit_x_group + poly_fit_group(2);
    
    subplot(1,3,iI)
    plot(fit_x,fit_y,'k-','linewidth',2)
    hold on
    plot(fit_x_group,fit_y_group,'r-','linewidth',2)
    
    plot(data.mrs.corrs.tobacco.allUsers.allDataSubUse(data.mrs.corrs.tobacco.allUsers.groupIdx==1)', ...
        data.mrs.corrs.tobacco.allUsers.allDataMRS((data.mrs.corrs.tobacco.allUsers.groupIdx==1)',iI), ...
        'go','MarkerFaceColor','w','linewidth',2,...
        'MarkerSize',6)
    plot(data.mrs.corrs.tobacco.allUsers.allDataSubUse(data.mrs.corrs.tobacco.allUsers.groupIdx==2)', ...
        data.mrs.corrs.tobacco.allUsers.allDataMRS((data.mrs.corrs.tobacco.allUsers.groupIdx==2)',iI), ...
        'bo','MarkerFaceColor','w','linewidth',2,...
        'MarkerSize',6)
    plot(data.mrs.corrs.tobacco.allUsers.allDataSubUse(data.mrs.corrs.tobacco.allUsers.groupIdx==3)', ...
        data.mrs.corrs.tobacco.allUsers.allDataMRS((data.mrs.corrs.tobacco.allUsers.groupIdx==3)',iI), ...
        'ro','MarkerFaceColor','w','linewidth',2,...
        'MarkerSize',6)
    
    % Plot stats
    text(fit_x(2)-4.5,((y_range(2))*.125),...
        ['r = ' num2str(data.mrs.corrs.tobacco.allUsers.r(iI))],'fontsize',statsFontSize)
    text(fit_x(2)-4.5,((y_range(2))*.125)-(y_range(2))*.03,...
        ['p = ' num2str(data.mrs.corrs.tobacco.allUsers.p(iI))],'fontsize',statsFontSize)
    % Plot Stats for PwPP group seperate
    text(fit_x(2)-4.5,((y_range(2))*.125)-(y_range(2))*.06,...
        ['r = ' num2str(data.mrs.corrs.tobacco.groupUsers(3).r(iI)) ' (PwPP)'],'fontsize',statsFontSize)
    text(fit_x(2)-4.5,((y_range(2))*.125)-(y_range(2))*.09,...
        ['p = ' num2str(data.mrs.corrs.tobacco.groupUsers(3).p(iI)) ' (PwPP)'],'fontsize',statsFontSize)
    
    set(gca,'xlim',[0 fit_x(2)+1])
    set(gca,'ylim',[0 y_range(2)+.1])
    set(gca,'xcolor','k','ycolor','k')
    xlabel('Servings of Tobacco / Day','color','k','fontsize',axisTitleFontSize)
    ylabel(sprintf('%s%s',options.which_metab{iI},' Concentration'),'color','k','fontsize',axisTitleFontSize)
    title(sprintf('%s%s\n%s',options.which_metab{iI},' Vs Tobacco Use','(only users, all scan sessions, ',options.whichVoxel,')'),'fontsize',titleFontSize)
end

set(gcf,'Units','inches')
set(gcf,'Position',figSize.mrsTobCorr.figSize,'color','w')



%% Plot symptom correlations
% Make seperate lists for bprs and other measures
options.symptom_general_list = {options.symptom_list{[1 2 numel(options.symptom_list)]}};
options.symptom_general_short = {options.symptom_short{[1 2 numel(options.symptom_short)]}};
options.symptom_bprs_list = {options.symptom_list{[3:7]}};
options.symptom_bprs_short = {options.symptom_short{[3:7]}};
tobAlcSwitch = {'tobacco'};
tobAlcTitle = {'Tobacco'};
for iJ = 1   % Tob/Alc
    figure; hold on
    % Set font sizes
    titleFontSize = 12;
    axisTitleFontSize = 12;
    axisLabelFontSize = 12;
    statsFontSize = 12;
    % Set figure size
    figSize.symptomTobCorr.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.symptomTobCorr.aspectRatio = [3 3];   % Aspect ratio
    figSize.symptomTobCorr.figSize = [0 0 ...
        figSize.symptomTobCorr.aspectRatio];   % Size/postion of fig
    
    % For SPQ/BPRS total/BACS
    for iI = 1:numel(options.symptom_general_short)
        % Fit the data
        clear poly_fit fit_x fit_y y_range
        [poly_fit] = polyfit(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse, ...
            data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSymptom, 1);
        fit_x = [min(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse') max(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse')];
        fit_y = poly_fit(1).*fit_x + poly_fit(2);
        y_range = [min(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSymptom) max(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSymptom)];
        
        % Fit the data for PwPP alone
        clear poly_fit_group fit_x_group fit_y_group y_range_group
        [poly_fit_group] = polyfit(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse(...
            data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==3), ...
            data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSymptom(...
            data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==3), 1);
        fit_x_group = [min(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse(...
            data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==3)') ...
            max(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse(...
            data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==3)')];
        fit_y_group = poly_fit_group(1).*fit_x_group + poly_fit_group(2);
        
        
        subplot(1,numel(options.symptom_general_list),iI)
        plot(fit_x,fit_y,'k-','linewidth',2)
        hold on
        plot(fit_x_group,fit_y_group,'r-','linewidth',2)
        
        plot(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==1), ...
            data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSymptom((data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==1)'), ...
            'go','MarkerFaceColor','w','linewidth',2,...
            'MarkerSize',6)
        plot(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==2), ...
            data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSymptom((data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==2)'), ...
            'bo','MarkerFaceColor','w','linewidth',2,...
            'MarkerSize',6)
        plot(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==3), ...
            data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSymptom((data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==3)'), ...
            'ro','MarkerFaceColor','w','linewidth',2,...
            'MarkerSize',6)
        
        if iJ==1
            xInd = 5;
        elseif iJ==2
            if strcmp(options.whichAlcMeasure,'aveConsumPerDay')
                xInd = 3;
            elseif strcmp(options.whichAlcMeasure,'aveConsumPerWeek')
                xInd = 25;
            end
        end
        if y_range(1)>=0
            % Plot the stats
            text(fit_x(2)-xInd,((y_range(2))*.125),...
                ['r = ' num2str(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.r)],'fontsize',statsFontSize)
            text(fit_x(2)-xInd,((y_range(2))*.125)-(y_range(2))*.03,...
                ['p = ' num2str(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.p)],'fontsize',statsFontSize)
            % Plot the stats for PwPP
            text(fit_x(2)-xInd,((y_range(2))*.125)-(y_range(2))*.06,...
                ['r = ' num2str(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).groupUsers(3).r) ' (PwPP)'],'fontsize',statsFontSize)
            text(fit_x(2)-xInd,((y_range(2))*.125)-(y_range(2))*.09,...
                ['p = ' num2str(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).groupUsers(3).p) ' (PwPP)'],'fontsize',statsFontSize)
        elseif y_range(1) < 0
            % Plot the stats
            text(fit_x(2)-xInd,((y_range(2)-y_range(1))*.125)+y_range(1),...
                ['r = ' num2str(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.r)],'fontsize',statsFontSize)
            text(fit_x(2)-xInd,(((y_range(2)-y_range(1))*.125)-(y_range(2)-y_range(1))*.03)+y_range(1),...
                ['p = ' num2str(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.p)],'fontsize',statsFontSize)
            % Plot the stats for PwPP
            text(fit_x(2)-xInd,(((y_range(2)-y_range(1))*.125)-(y_range(2)-y_range(1))*.06)+y_range(1),...
                ['r = ' num2str(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).groupUsers(3).r) ' (PwPP)'],'fontsize',statsFontSize)
            text(fit_x(2)-xInd,(((y_range(2)-y_range(1))*.125)-(y_range(2)-y_range(1))*.09)+y_range(1),...
                ['p = ' num2str(data.symptom.(options.symptom_general_short{iI}).corrs.(tobAlcSwitch{iJ}).groupUsers(3).p) ' (PwPP)'],'fontsize',statsFontSize)
        end
        
        set(gca,'xlim',[0 fit_x(2)+1])
        if y_range(1)>=0
            set(gca,'ylim',[0 y_range(2)+.1])
        elseif y_range(1)<0
            set(gca,'ylim',[y_range(1)-.1 y_range(2)+.1])
        end
        set(gca,'xcolor','k','ycolor','k')
        xlabel('Servings of Tobacco / Day','color','k','fontsize',axisTitleFontSize)
        ylabel(sprintf('%s%s',options.symptom_general_list{iI},' '),'color','k','fontsize',axisTitleFontSize)
        title(sprintf('%s%s\n%s%s\n%s',options.symptom_general_list{iI},' Vs',tobAlcTitle{iJ},' Use','(users, all scan sessions)'),'fontsize',titleFontSize)
    end
    
    set(gcf,'Units','inches')
    set(gcf,'Position',figSize.symptomTobCorr.figSize,'color','w')
    
    figure; hold on
    % Set font sizes
    titleFontSize = 12;
    axisTitleFontSize = 12;
    axisLabelFontSize = 12;
    statsFontSize = 12;
    % Set figure size
    figSize.symptomTobCorr.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.symptomTobCorr.aspectRatio = [3 3];   % Aspect ratio
    figSize.symptomTobCorr.figSize = [0 0 ...
        figSize.symptomTobCorr.aspectRatio];   % Size/postion of fig
    
    % For BPRS sub measures
    for iI = 1:numel(options.symptom_bprs_short)
        % Fit the data
        clear poly_fit fit_x fit_y y_range
        [poly_fit] = polyfit(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse, ...
            data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSymptom, 1);
        fit_x = [min(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse') ...
            max(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse')];
        fit_y = poly_fit(1).*fit_x + poly_fit(2);
        y_range = [min(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSymptom) ...
            max(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSymptom)];
        
        % Fit the data for PwPP alone
        clear poly_fit_group fit_x_group fit_y_group y_range_group
        [poly_fit_group] = polyfit(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse(...
            data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==3), ...
            data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSymptom(...
            data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==3), 1);
        fit_x_group = [min(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse(...
            data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==3)') ...
            max(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse(...
            data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==3)')];
        fit_y_group = poly_fit_group(1).*fit_x_group + poly_fit_group(2);
        
        
        subplot(1,numel(options.symptom_bprs_list),iI)
        plot(fit_x,fit_y,'k-','linewidth',2)
        hold on
        plot(fit_x_group,fit_y_group,'r-','linewidth',2)
        
        plot(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==1), ...
            data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSymptom((data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==1)'), ...
            'go','MarkerFaceColor','w','linewidth',2,...
            'MarkerSize',6)
        plot(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==2), ...
            data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSymptom((data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==2)'), ...
            'bo','MarkerFaceColor','w','linewidth',2,...
            'MarkerSize',6)
        plot(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSubUse(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==3), ...
            data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.allDataSymptom((data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.groupIdx==3)'), ...
            'ro','MarkerFaceColor','w','linewidth',2,...
            'MarkerSize',6)
        
        if iJ==1
            xInd = 5;
        elseif iJ==2
            if strcmp(options.whichAlcMeasure,'aveConsumPerDay')
                xInd = 3;
            elseif strcmp(options.whichAlcMeasure,'aveConsumePerWeek')
                xInd = 25;
            end
        end
        if y_range(1)>=0
            % Plot the stats
            text(fit_x(2)-xInd,((y_range(2))*.125),...
                ['r = ' num2str(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.r)],'fontsize',statsFontSize)
            text(fit_x(2)-xInd,((y_range(2))*.125)-(y_range(2))*.03,...
                ['p = ' num2str(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.p)],'fontsize',statsFontSize)
            % Plot the stats for PwPP
            text(fit_x(2)-xInd,((y_range(2))*.125)-(y_range(2))*.06,...
                ['r = ' num2str(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).groupUsers(3).r) ' (PwPP)'],'fontsize',statsFontSize)
            text(fit_x(2)-xInd,((y_range(2))*.125)-(y_range(2))*.09,...
                ['p = ' num2str(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).groupUsers(3).p) ' (PwPP)'],'fontsize',statsFontSize)
        elseif y_range(1) < 0
            % Plot the stats
            text(fit_x(2)-xInd,((y_range(2)-y_range(1))*.125)+y_range(1),...
                ['r = ' num2str(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.r)],'fontsize',statsFontSize)
            text(fit_x(2)-xInd,(((y_range(2)-y_range(1))*.125)-(y_range(2)-y_range(1))*.03)+y_range(1),...
                ['p = ' num2str(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).allUsers.p)],'fontsize',statsFontSize)
            % Plot the stats for PwPP
            text(fit_x(2)-xInd,(((y_range(2)-y_range(1))*.125)-(y_range(2)-y_range(1))*.06)+y_range(1),...
                ['r = ' num2str(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).groupUsers(3).r) ' (PwPP)'],'fontsize',statsFontSize)
            text(fit_x(2)-xInd,(((y_range(2)-y_range(1))*.125)-(y_range(2)-y_range(1))*.09)+y_range(1),...
                ['p = ' num2str(data.symptom.(options.symptom_bprs_short{iI}).corrs.(tobAlcSwitch{iJ}).groupUsers(3).p) ' (PwPP)'],'fontsize',statsFontSize)
        end
        
        set(gca,'xlim',[0 fit_x(2)+1])
        if y_range(1)>=0
            set(gca,'ylim',[0 y_range(2)+.1])
        elseif y_range(1)<0
            set(gca,'ylim',[y_range(1)-.1 y_range(2)+.1])
        end
        set(gca,'xcolor','k','ycolor','k')
        xlabel('Servings of Tobacco / Day','color','k','fontsize',axisTitleFontSize)
        ylabel(sprintf('%s%s',options.symptom_bprs_list{iI},' '),'color','k','fontsize',axisTitleFontSize)
        title(sprintf('%s%s\n%s%s\n%s',options.symptom_bprs_list{iI},' Vs ',tobAlcTitle{iJ},' Use','(users, all sessions)'),'fontsize',titleFontSize)
    end
    
    set(gcf,'Units','inches')
    set(gcf,'Position',figSize.symptomTobCorr.figSize,'color','w')
end


end

% Subject number, recruitment pHCP, consensus diagnosis (consensus by all
% clinical staff)
% 2 fields under subject number
% Primary and secondary conses diagnosis
% Under secondary might be info about alcohol use (?)
% SCID - paper files somewhere in the cold storage
% Info on drug/nicotene/alcohol
% Can ask Caroline, would know whats in there or how many people would have
% informative data there.




