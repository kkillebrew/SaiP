% Load in and sort out the substance use data from the redcap report for
% syon.
%
%

function [data,options] = load_SYON_subAbuseRedcap(options)

%% Load in redcap csv file for substance use:
data.rawCSV = readtable(options.syon_target_file,'FileType','text','Delimiter',...
    'comma','HeaderLines',0,'ReadRowNames',0,'ReadVariableNames',1,...
    'TreatAsEmpty','cormat');

% Find all and unique subjects
% Remove any subjects that don't have the correctly formated subj num 'SXXXXXXX'
validSubjInd = regexp(data.rawCSV.RecordID, '\w\d\d\d\d\d\d\d','once');
validSubjInd = ~cellfun('isempty',validSubjInd);
data.rawCSV = data.rawCSV(validSubjInd,:);
% Index subjects and day type ('Clinical', 'EEG', 'MRI')
data.runListRaw = data.rawCSV.EventName;
data.runListUniq = unique(data.runListRaw);
data.runListUniqShort = {'MRI'};
for iI = 1:length(data.runListUniq)
    data.runTypeRawIdx(strcmp(data.runListRaw,data.runListUniq(iI))) = iI;   % Index run types
end
% We only care about the mri day here
data.mriIdx = zeros([length(data.rawCSV.EventName) 1]);
data.mriIdx(data.runTypeRawIdx == 3) = 1;

data.subjListRaw = cellfun(@(x) x(1:8),data.rawCSV.RecordID,'UniformOutput',0);   % Subj list array including all runs
data.subjListUniq = unique(data.subjListRaw);   % Unique subj list
daysOtW = {'Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday'};

% Index subject group
for iI = 1:length(data.subjListUniq)
    data.subjRawIdx(strcmp(data.subjListRaw,data.subjListUniq(iI))) = iI;   % Index raw subj list
    if str2num(data.subjListRaw{iI}(2:end)) < 2000000
        data.groupTypeRaw(iI) = 1;
    elseif str2num(data.subjListRaw{iI}(2:end)) >= 2000000 && str2num(data.subjListRaw{iI}(2:end)) < 6000000
        data.groupTypeRaw(iI) = 2;
    elseif str2num(data.subjListRaw{iI}(2:end)) >= 6000000
        data.groupTypeRaw(iI) = 3;
    end
end

% Manual exlusion list:
% For some reasone this person has nothing entered in their sub ab redcap (not 0's just empty).
% Removing as the NaNs they create mess up stats.
options.manualExclude = [];

%% Grab all relevant data for each subject
mriCounter = 0;
for iI = 1:size(data.subjListUniq,1)
    if ~ismember(['S' num2str(options.manualExclude)],data.subjListUniq(iI))   % Make sure subj is not in manual exclude list
        %% MRI
        if (sum(data.subjRawIdx==iI & data.runTypeRawIdx==3) >= 1)   % Check that this person has data
            % Find the index value for this participant for this run
            partRunIdx = find(data.subjRawIdx==iI & data.runTypeRawIdx==3);
            % If there are more than 1 of this type of visit, then take the
            % first. Assuming earlier runs contained paperwork at least, w/
            % later runs having MRI.
            if numel(partRunIdx) == 1
            elseif numel(partRunIdx) > 1
                partRunIdx = partRunIdx(1);
            end
            
            % Make sure there is actual data here for this subject
            if ~isempty(data.rawCSV.Date_8{partRunIdx})
                
                mriCounter = mriCounter + 1;
                
                %% Grab nicotine use data from 7TA
                % Subnum and Date of 7TA sessions
                data.(data.runListUniqShort{1}).date(mriCounter) = yyyymmdd(datetime(data.rawCSV.Date_8(partRunIdx)));
                data.(data.runListUniqShort{1}).subjNum(mriCounter) = str2num(data.subjListUniq{iI}(2:end));
                % Total cigarettes/cigars on 7TA day
                data.(data.runListUniqShort{1}).tobacco.cigCurrDay(mriCounter) = ...
                    data.rawCSV.x6_HowManyCigarettesAnd_orE_cigarettesHaveYouHadToday_(partRunIdx);
                % Total pinches on 7TA day
                data.(data.runListUniqShort{1}).tobacco.pinchCurrDay(mriCounter) = ...
                    data.rawCSV.x7_HowManyPinchesOfTobaccoHaveYouHadToday_(partRunIdx);
                % Ave cigarettes/cigars on 7TA day
                data.(data.runListUniqShort{1}).tobacco.cigAveDay(mriCounter) = ...
                    data.rawCSV.x8_HowManyCigarettesAnd_orE_cigarettesDoYouUsuallyHaveByThisTim(partRunIdx);
                % Ave pinches on 7TA day
                data.(data.runListUniqShort{1}).tobacco.pinchAveDay(mriCounter) = ...
                    data.rawCSV.x9_HowManyPinchesOfTobaccoDoYouUsuallyHaveByThisTimeOfTheDayOnA(partRunIdx);
                
                %% Grab alcohol use data from 7TA
                % Total beers/glasses of wine on 7TA day
                data.(data.runListUniqShort{1}).alcohol.beerWineCurrDay(mriCounter) = ...
                    data.rawCSV.x12_HowManyBeersAndGlassesOfWineHaveYouHadToday_(partRunIdx);
                % Total beers/glasses of wine in past 24hrs of 7TA day
                data.(data.runListUniqShort{1}).alcohol.beerWinePast24(mriCounter) = ...
                    data.rawCSV.x14_HowManyBeersAndGlassesOfWineHaveYouHadInTheLast24Hours_(partRunIdx);
                % Total liquor on 7TA day
                data.(data.runListUniqShort{1}).alcohol.beerWineCurrDay(mriCounter) = ...
                    data.rawCSV.x13_HowManyMixedDrinksOrShotsOfLiquorHaveYouHadToday_(partRunIdx);
                % Total liquor in past 24hrs of 7TA day
                data.(data.runListUniqShort{1}).alcohol.beerWinePast24(mriCounter) = ...
                    data.rawCSV.x15_HowManyMixedDrinksOrShotsOfLiquorHaveYouHadInTheLast24Hours(partRunIdx);
                % Last consumption of alcohol
                data.(data.runListUniqShort{1}).alcohol.lastConsumption(mriCounter) = ...
                    yyyymmdd(datetime(data.rawCSV.x16_WhenWasYourLastConsumptionOfAlcohol_mm_dd_yyyy__(partRunIdx)));
                data.(data.runListUniqShort{1}).alcohol.lastConsumptionAmount(mriCounter) = ...
                    data.rawCSV.x17_HowMuchAlcoholDidYouDrink__OfStandardDrinks_DuringYourLastC(partRunIdx);
                % Most consumed in last month
                data.(data.runListUniqShort{1}).alcohol.mostConsumed(mriCounter) = ...
                    data.rawCSV.x18_WhatIsTheMostAlcohol__OfStandardDrinks_YouHaveConsumedInOne(partRunIdx);
                % Average days consumed in a week
                data.(data.runListUniqShort{1}).alcohol.aveDaysConsumed(mriCounter) = ...
                    data.rawCSV.x19_InATypicalWeek_HowManyDaysDoYouHaveSomeKindOfAlcoholToDrink(partRunIdx);
                % 7 Day breakdown
                for iK = 1:7   % For the 7 DOW
                    data.(data.runListUniqShort{1}).alcohol.dailyUse(mriCounter,iK) = ...
                        data.rawCSV.(sprintf('%s%s','x20_',daysOtW{iK}))(partRunIdx);
                end
                % Average consumption / day
                data.(data.runListUniqShort{1}).alcohol.aveConsumPerDay(mriCounter) = ...
                    squeeze(nanmean(data.(data.runListUniqShort{1}).alcohol.dailyUse(mriCounter,:),2));
                % Average consumption / week
                data.(data.runListUniqShort{1}).alcohol.aveConsumPerWeek(mriCounter) = ...
                    nansum(data.(data.runListUniqShort{1}).alcohol.dailyUse(mriCounter,:),2);
            end
        end
    end
end
clear mriCounter



end