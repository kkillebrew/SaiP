% Load in and sort out the substance use data from the redcap report for
% phcp.
%
%

function [data,options] = load_pHCP_subAbuseRedcap(options)

%% Load in redcap csv file for substance use:
data.rawCSV = readtable(options.phcp_target_file,'FileType','text','Delimiter',...
    'comma','HeaderLines',0,'ReadRowNames',0,'ReadVariableNames',1,...
    'TreatAsEmpty','cormat');

% Find all and unique subjects
data.subjListRaw = cellfun(@(x) x(1:8),data.rawCSV.RecordID,'UniformOutput',0);   % Subj list array including all runs
data.subjListUniq = unique(data.subjListRaw);   % Unique subj list
data.runListRaw = data.rawCSV.EventName;
data.runListUniq = unique(data.runListRaw);
data.runListUniqShort = {'mri3TA','mri3TB','mri7TA','mri7TB','mri7TZ','Clinical'};
daysOtW = {'Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday'};

% Index subject group
for iI = 1:length(data.subjListRaw)
    if str2num(data.subjListRaw{iI}(2:end)) < 2000000
        data.groupTypeRaw(iI) = 1;
    elseif str2num(data.subjListRaw{iI}(2:end)) >= 2000000 && str2num(data.subjListRaw{iI}(2:end)) < 6000000
        data.groupTypeRaw(iI) = 2;
    elseif str2num(data.subjListRaw{iI}(2:end)) >= 6000000
        data.groupTypeRaw(iI) = 3;
    end
end

% Index subjects and day type (Clinical', '3T MRI Session 1a', etc.)
% Need in order to seperate data between subjects and runs
for iI = 1:length(data.subjListUniq)
    data.subjRawIdx(strcmp(data.subjListRaw,data.subjListUniq(iI))) = iI;   % Index raw subj list
    if str2num(data.subjListUniq{iI}(2:end)) < 2000000
        data.groupTypeUniq(iI) = 1;   % Index group type
    elseif str2num(data.subjListUniq{iI}(2:end)) >= 2000000 && str2num(data.subjListUniq{iI}(2:end)) < 6000000
        data.groupTypeUniq(iI) = 2;
    elseif str2num(data.subjListUniq{iI}(2:end)) >= 6000000
        data.groupTypeUniq(iI) = 3;
    end
end
for iI = 1:length(data.runListUniq)
    data.runTypeRawIdx(strcmp(data.runListRaw,data.runListUniq(iI))) = iI;   % Index run types
end

% Manual exlusion list:
% For some reasone this person has nothing entered in their sub ab redcap (not 0's just empty).
% Removing as the NaNs they create mess up stats.
options.manualExclude = 1010053;

%% Grab all relevant data for each subject
clinCounter = 0;
aCounter = 0;
bCounter = 0;
zCounter = 0;
for iI = 1:size(data.subjListUniq,1)
    if ~ismember(['P' num2str(options.manualExclude)],data.subjListUniq(iI))
        %% Clinical day
        if (sum(data.subjRawIdx==iI & data.runTypeRawIdx==6) == 1)
            % Find the index value for this participant for this run
            partRunIdx = find(data.subjRawIdx==iI & data.runTypeRawIdx==6);
            % If there are more than 1 of this type of visit, then take the
            % second (or last). Assuming earlier runs were not complete, so
            % grab the most complete run.
            if numel(partRunIdx) == 1
            elseif numel(partRunIdx) > 1
                partRunIdx = partRunIdx(end);
            end
            
            % Make sure there is actual data here for this subject
            if ~isempty(~isempty(data.rawCSV.Date{partRunIdx}))
                
                clinCounter = clinCounter + 1;
                
                %% Grab nicotine use from clinical day
                % Subjnum and Date of clinical
                data.(data.runListUniqShort{6}).date(clinCounter) = yyyymmdd(datetime(data.rawCSV.Date{partRunIdx}));
                data.(data.runListUniqShort{6}).subjNum(clinCounter) = str2num(data.subjListUniq{iI}(2:end));
                % Use in the last 7 days?
                if strcmp(data.rawCSV.HasTheParticipantUsedAnyKindOfTobaccoInThePast7Days_(partRunIdx),'Yes')
                    data.(data.runListUniqShort{6}).tobacco.useLast7Days(clinCounter) = 1;
                elseif strcmp(data.rawCSV.HasTheParticipantUsedAnyKindOfTobaccoInThePast7Days_(partRunIdx),'No')
                    data.(data.runListUniqShort{6}).tobacco.useLast7Days(clinCounter) = 0;
                else
                    data.(data.runListUniqShort{6}).tobacco.useLast7Days(clinCounter) = NaN;
                end
                % Use on clinical day
                data.(data.runListUniqShort{6}).tobacco.useOnClinDay(clinCounter) = ...
                    data.rawCSV.HowManyCigarettes_orOtherTobacco_HaveYouSmokedToday_(partRunIdx);
                % 7 Day breakdown
                for iK = 1:7   % For the 7 DOW
                    data.(data.runListUniqShort{6}).tobacco.cigarettes(clinCounter,iK) = ...
                        data.rawCSV.(sprintf('%s%d%s','TobaccoDay',iK,'Cigarettes'))(partRunIdx);
                    data.(data.runListUniqShort{6}).tobacco.cigars(clinCounter,iK) = ...
                        data.rawCSV.(sprintf('%s%d%s','TobaccoDay',iK,'Cigars'))(partRunIdx);
                    data.(data.runListUniqShort{6}).tobacco.pipe(clinCounter,iK) = ...
                        data.rawCSV.(sprintf('%s%d%s','TobaccoDay',iK,'Pipe'))(partRunIdx);
                    data.(data.runListUniqShort{6}).tobacco.chew(clinCounter,iK) = ...
                        data.rawCSV.(sprintf('%s%d%s','TobaccoDay',iK,'ChewingTobacco'))(partRunIdx);
                    data.(data.runListUniqShort{6}).tobacco.snuff(clinCounter,iK) = ...
                        data.rawCSV.(sprintf('%s%d%s','TobaccoDay',iK,'Snuff'))(partRunIdx);
                    data.(data.runListUniqShort{6}).tobacco.other(clinCounter,iK) = ...
                        data.rawCSV.(sprintf('%s%d%s','TobaccoDay',iK,'Other'))(partRunIdx);
                    data.(data.runListUniqShort{6}).tobacco.any(clinCounter,iK) = nansum([data.rawCSV.(sprintf('%s%d%s','TobaccoDay',iK,'Cigarettes'))(partRunIdx),...
                        data.rawCSV.(sprintf('%s%d%s','TobaccoDay',iK,'Cigars'))(partRunIdx),...
                        data.rawCSV.(sprintf('%s%d%s','TobaccoDay',iK,'Pipe'))(data.subjRawIdx==iI & data.runTypeRawIdx==6),...
                        data.rawCSV.(sprintf('%s%d%s','TobaccoDay',iK,'ChewingTobacco'))(partRunIdx),...
                        data.rawCSV.(sprintf('%s%d%s','TobaccoDay',iK,'Snuff'))(partRunIdx),...
                        data.rawCSV.(sprintf('%s%d%s','TobaccoDay',iK,'Other'))(partRunIdx),]);
                end
                % Average across their 7 days
                data.(data.runListUniqShort{6}).tobacco.aveDailyUse(clinCounter) = nanmean(data.(data.runListUniqShort{6}).tobacco.any(clinCounter,:),2);
                
                %% Grab alcohol use from clinical day
                % Use in the last 7 days?
                if strcmp(data.rawCSV.HasTheParticipantHadAnyAlcoholicDrinksInThePast7Days_(partRunIdx),'Yes')
                    data.(data.runListUniqShort{6}).alcohol.useLast7Days(clinCounter) = 1;
                elseif strcmp(data.rawCSV.HasTheParticipantHadAnyAlcoholicDrinksInThePast7Days_(partRunIdx),'No')
                    data.(data.runListUniqShort{6}).alcohol.useLast7Days(clinCounter) = 0;
                else
                    data.(data.runListUniqShort{6}).alcohol.useLast7Days(clinCounter) = NaN;
                end
                %         % Use on clinical day (no alcohol use day of)
                %         data.(data.runListUniqShort{6}).alcohol.useOnClinDay(clinCounter) = ...
                %             data.rawCSV.HowManyCigarettes_orOtherAlcohol_HaveYouSmokedToday_(partRunIdx);
                % 7 Day breakdown
                for iK = 1:7   % For the 7 DOW
                    data.(data.runListUniqShort{6}).alcohol.beer(clinCounter,iK) = ...
                        data.rawCSV.(sprintf('%s%d%s','AlcoholDay',iK,'Beer_WineCoolers'))(partRunIdx);
                    data.(data.runListUniqShort{6}).alcohol.malt(clinCounter,iK) = ...
                        data.rawCSV.(sprintf('%s%d%s','AlcoholDay',iK,'MaltLiquor'))(partRunIdx);
                    data.(data.runListUniqShort{6}).alcohol.wine(clinCounter,iK) = ...
                        data.rawCSV.(sprintf('%s%d%s','AlcoholDay',iK,'TableWine'))(partRunIdx);
                    data.(data.runListUniqShort{6}).alcohol.liquor(clinCounter,iK) = ...
                        data.rawCSV.(sprintf('%s%d%s','AlcoholDay',iK,'HardLiquor_Spirits'))(partRunIdx);
                    data.(data.runListUniqShort{6}).alcohol.other(clinCounter,iK) = ...
                        data.rawCSV.(sprintf('%s%d%s','AlcoholDay',iK,'Other'))(partRunIdx);
                    data.(data.runListUniqShort{6}).alcohol.any(clinCounter,iK) = nansum([data.rawCSV.(sprintf('%s%d%s','AlcoholDay',iK,'Beer_WineCoolers'))(partRunIdx),...
                        data.rawCSV.(sprintf('%s%d%s','AlcoholDay',iK,'MaltLiquor'))(partRunIdx),...
                        data.rawCSV.(sprintf('%s%d%s','AlcoholDay',iK,'TableWine'))(data.subjRawIdx==iI & data.runTypeRawIdx==6),...
                        data.rawCSV.(sprintf('%s%d%s','AlcoholDay',iK,'HardLiquor_Spirits'))(partRunIdx),...
                        data.rawCSV.(sprintf('%s%d%s','AlcoholDay',iK,'Other'))(partRunIdx)]);
                end
                % Average across their 7 days
                data.(data.runListUniqShort{6}).alcohol.aveDailyUse(clinCounter) = nanmean(data.(data.runListUniqShort{6}).alcohol.any(clinCounter,:),2);
                
                
            end
        end
        
        %% 7TA
        if (sum(data.subjRawIdx==iI & data.runTypeRawIdx==3) >= 1) && ~isempty(data.rawCSV.DateCompleted_mm_dd_yyyy_(iI))
            % Find the index value for this participant for this run
            partRunIdx = find(data.subjRawIdx==iI & data.runTypeRawIdx==3);
            % If there are more than 1 of this type of visit, then take the
            % second (or last). Assuming earlier runs were not complete, so
            % grab the most complete run.
            if numel(partRunIdx) == 1
            elseif numel(partRunIdx) > 1
                partRunIdx = partRunIdx(end);
            end
            
            % Make sure there is actual data here for this subject
            if ~isempty(data.rawCSV.DateCompleted_mm_dd_yyyy_{partRunIdx})
                
                aCounter = aCounter + 1;
                
                %% Grab nicotine use data from 7TA
                % Subnum and Date of 7TA sessions
                data.(data.runListUniqShort{3}).date(aCounter) = yyyymmdd(datetime(data.rawCSV.DateCompleted_mm_dd_yyyy_(partRunIdx)));
                data.(data.runListUniqShort{3}).subjNum(aCounter) = str2num(data.subjListUniq{iI}(2:end));
                % Total cigarettes/cigars on 7TA day
                data.(data.runListUniqShort{3}).tobacco.cigCurrDay(aCounter) = ...
                    data.rawCSV.x6_HowManyCigarettesAndCigarsHaveYouHadToday_(partRunIdx);
                % Total pinches on 7TA day
                data.(data.runListUniqShort{3}).tobacco.pinchCurrDay(aCounter) = ...
                    data.rawCSV.x7_HowManyPinchesOfTobaccoHaveYouHadToday_(partRunIdx);
                % Ave cigarettes/cigars on 7TA day
                data.(data.runListUniqShort{3}).tobacco.cigAveDay(aCounter) = ...
                    data.rawCSV.x8_HowManyCigarettesAndCigarsDoYouUsuallyHaveByThisTimeOfDayOnA(partRunIdx);
                % Ave pinches on 7TA day
                data.(data.runListUniqShort{3}).tobacco.pinchAveDay(aCounter) = ...
                    data.rawCSV.x9_HowManyPinchesOfTobaccoDoYouUsuallyHaveByThisTimeOfTheDayOnA(partRunIdx);
                
                %% Grab alcohol use data from 7TA
                % Total beers/glasses of wine on 7TA day
                data.(data.runListUniqShort{3}).alcohol.beerWineCurrDay(aCounter) = ...
                    data.rawCSV.x12_HowManyBeersAndGlassesOfWineHaveYouHadToday_(partRunIdx);
                % Total beers/glasses of wine in past 24hrs of 7TA day
                data.(data.runListUniqShort{3}).alcohol.beerWinePast24(aCounter) = ...
                    data.rawCSV.x14_HowManyBeersAndGlassesOfWineHaveYouHadInTheLast24Hours_(partRunIdx);
                % Total liquor on 7TA day
                data.(data.runListUniqShort{3}).alcohol.beerWineCurrDay(aCounter) = ...
                    data.rawCSV.x13_HowManyMixedDrinksOrShotsOfLiquorHaveYouHadToday_(partRunIdx);
                % Total liquor in past 24hrs of 7TA day
                data.(data.runListUniqShort{3}).alcohol.beerWinePast24(aCounter) = ...
                    data.rawCSV.x15_HowManyMixedDrinksOrShotsOfLiquorHaveYouHadInTheLast24Hours(partRunIdx);
                % Last consumption of alcohol
                data.(data.runListUniqShort{3}).alcohol.lastConsumption(aCounter) = ...
                    yyyymmdd(datetime(data.rawCSV.x16_WhenWasYourLastConsumptionOfAlcohol_mm_dd_yyyy__(partRunIdx)));
                data.(data.runListUniqShort{3}).alcohol.lastConsumptionAmount(aCounter) = ...
                    data.rawCSV.x17_HowMuchAlcoholDidYouDrink__OfStandardDrinks_DuringYourLastC(partRunIdx);
                % Most consumed in last month
                data.(data.runListUniqShort{3}).alcohol.mostConsumed(aCounter) = ...
                    data.rawCSV.x18_WhatIsTheMostAlcohol__OfStandardDrinks_YouHaveConsumedInOne(partRunIdx);
                % Average days consumed in a week
                data.(data.runListUniqShort{3}).alcohol.aveDaysConsumed(aCounter) = ...
                    data.rawCSV.x19_InATypicalWeek_HowManyDaysDoYouHaveSomeKindOfAlcoholToDrink(partRunIdx);
                % 7 Day breakdown
                for iK = 1:7   % For the 7 DOW
                    data.(data.runListUniqShort{3}).alcohol.dailyUse(aCounter,iK) = ...
                        data.rawCSV.(sprintf('%s%s','x20_',daysOtW{iK}))(partRunIdx);
                end
                % Average consumption / day
                data.(data.runListUniqShort{3}).alcohol.aveConsumPerDay(aCounter) = ...
                    squeeze(nanmean(data.(data.runListUniqShort{3}).alcohol.dailyUse(aCounter,:),2));
                % Average consumption / week
                data.(data.runListUniqShort{3}).alcohol.aveConsumPerWeek(aCounter) = ...
                    nansum(data.(data.runListUniqShort{3}).alcohol.dailyUse(aCounter,:),2);
            end
        end
        
        %% 7TB
        if (sum(data.subjRawIdx==iI & data.runTypeRawIdx==4) >= 1)
            % Find the index value for this participant for this run
            partRunIdx = find(data.subjRawIdx==iI & data.runTypeRawIdx==4);
            % If there are more than 1 of this type of visit, then take the
            % second (or last). Assuming earlier runs were not complete, so
            % grab the most complete run.
            if numel(partRunIdx) == 1
            elseif numel(partRunIdx) > 1
                partRunIdx = partRunIdx(end);
            end
            
            % Make sure there is actual data here for this subject
            if ~isempty(data.rawCSV.DateCompleted_mm_dd_yyyy_{partRunIdx})
                
                bCounter = bCounter + 1;
                
                %% Grab nicotine use data from 7TB
                % Subjnum and Date of 7TB sessions
                data.(data.runListUniqShort{4}).date(bCounter) = yyyymmdd(datetime(data.rawCSV.DateCompleted_mm_dd_yyyy_(partRunIdx)));
                data.(data.runListUniqShort{4}).subjNum(bCounter) = str2num(data.subjListUniq{iI}(2:end));
                % Total cigarettes/cigars on 7TA day
                data.(data.runListUniqShort{4}).tobacco.cigCurrDay(bCounter) = ...
                    data.rawCSV.x6_HowManyCigarettesAndCigarsHaveYouHadToday_(partRunIdx);
                % Total pinches on 7TA day
                data.(data.runListUniqShort{4}).tobacco.pinchCurrDay(bCounter) = ...
                    data.rawCSV.x7_HowManyPinchesOfTobaccoHaveYouHadToday_(partRunIdx);
                % Ave cigarettes/cigars on 7TA day
                data.(data.runListUniqShort{4}).tobacco.cigAveDay(bCounter) = ...
                    data.rawCSV.x8_HowManyCigarettesAndCigarsDoYouUsuallyHaveByThisTimeOfDayOnA(partRunIdx);
                % Ave pinches on 7TA day
                data.(data.runListUniqShort{4}).tobacco.pinchAveDay(bCounter) = ...
                    data.rawCSV.x9_HowManyPinchesOfTobaccoDoYouUsuallyHaveByThisTimeOfTheDayOnA(partRunIdx);
                
                %% Grab alcohol use data from 7TB
                % Total beers/glasses of wine on 7TB day
                data.(data.runListUniqShort{4}).alcohol.beerWineCurrDay(bCounter) = ...
                    data.rawCSV.x12_HowManyBeersAndGlassesOfWineHaveYouHadToday_(partRunIdx);
                % Total beers/glasses of wine in past 24hrs of 7TB day
                data.(data.runListUniqShort{4}).alcohol.beerWinePast24(bCounter) = ...
                    data.rawCSV.x14_HowManyBeersAndGlassesOfWineHaveYouHadInTheLast24Hours_(partRunIdx);
                % Total liquor on 7TB day
                data.(data.runListUniqShort{4}).alcohol.beerWineCurrDay(bCounter) = ...
                    data.rawCSV.x13_HowManyMixedDrinksOrShotsOfLiquorHaveYouHadToday_(partRunIdx);
                % Total liquor in past 24hrs of 7TB day
                data.(data.runListUniqShort{4}).alcohol.beerWinePast24(bCounter) = ...
                    data.rawCSV.x15_HowManyMixedDrinksOrShotsOfLiquorHaveYouHadInTheLast24Hours(partRunIdx);
                % Last consumption of alcohol
                data.(data.runListUniqShort{4}).alcohol.lastConsumption(bCounter) = ...
                    yyyymmdd(datetime(data.rawCSV.x16_WhenWasYourLastConsumptionOfAlcohol_mm_dd_yyyy__(partRunIdx)));
                data.(data.runListUniqShort{4}).alcohol.lastConsumptionAmount(bCounter) = ...
                    data.rawCSV.x17_HowMuchAlcoholDidYouDrink__OfStandardDrinks_DuringYourLastC(partRunIdx);
                % Most consumed in last month
                data.(data.runListUniqShort{4}).alcohol.mostConsumed(bCounter) = ...
                    data.rawCSV.x18_WhatIsTheMostAlcohol__OfStandardDrinks_YouHaveConsumedInOne(partRunIdx);
                % Average days consumed in a week
                data.(data.runListUniqShort{4}).alcohol.aveDaysConsumed(bCounter) = ...
                    data.rawCSV.x19_InATypicalWeek_HowManyDaysDoYouHaveSomeKindOfAlcoholToDrink(partRunIdx);
                % 7 Day breakdown
                for iK = 1:7   % For the 7 DOW
                    data.(data.runListUniqShort{4}).alcohol.dailyUse(bCounter,iK) = ...
                        data.rawCSV.(sprintf('%s%s','x20_',daysOtW{iK}))(partRunIdx);
                end
                % Average consumption / day
                data.(data.runListUniqShort{4}).alcohol.aveConsumPerDay(bCounter) = ...
                    squeeze(nanmean(data.(data.runListUniqShort{4}).alcohol.dailyUse(bCounter,:),2));
                % Average consumption / week
                data.(data.runListUniqShort{4}).alcohol.aveConsumPerWeek(bCounter) = ...
                    nansum(data.(data.runListUniqShort{4}).alcohol.dailyUse(bCounter,:),2);
            end
        end
        
        %% 7TZ
        if (sum(data.subjRawIdx==iI & data.runTypeRawIdx==5) >= 1)
            % Find the index value for this participant for this run
            partRunIdx = find(data.subjRawIdx==iI & data.runTypeRawIdx==5);
            % If there are more than 1 of this type of visit, then take the
            % second (or last). Assuming earlier runs were not complete, so
            % grab the most complete run.
            if numel(partRunIdx) == 1
            elseif numel(partRunIdx) > 1
                partRunIdx = partRunIdx(end);
            end
            
            % Make sure there is actual data here for this subject
            if ~isempty(data.rawCSV.DateCompleted_mm_dd_yyyy_{partRunIdx})
                
                zCounter = zCounter + 1;
                
                %% Grab nicotine use data from 7TZ
                % Date of 7TZ sessions
                data.(data.runListUniqShort{5}).date(zCounter) = yyyymmdd(datetime(data.rawCSV.DateCompleted_mm_dd_yyyy_(partRunIdx)));
                data.(data.runListUniqShort{5}).subjNum(zCounter) = str2num(data.subjListUniq{iI}(2:end));
                % Total cigarettes/cigars on 7TZ day
                data.(data.runListUniqShort{5}).tobacco.cigCurrDay(zCounter) = ...
                    data.rawCSV.x6_HowManyCigarettesAndCigarsHaveYouHadToday_(partRunIdx);
                % Total pinches on 7TZ day
                data.(data.runListUniqShort{5}).tobacco.pinchCurrDay(zCounter) = ...
                    data.rawCSV.x7_HowManyPinchesOfTobaccoHaveYouHadToday_(partRunIdx);
                % Ave cigarettes/cigars on 7TZ day
                data.(data.runListUniqShort{5}).tobacco.cigAveDay(zCounter) = ...
                    data.rawCSV.x8_HowManyCigarettesAndCigarsDoYouUsuallyHaveByThisTimeOfDayOnA(partRunIdx);
                % Ave pinches on 7TZ day
                data.(data.runListUniqShort{5}).tobacco.pinchAveDay(zCounter) = ...
                    data.rawCSV.x9_HowManyPinchesOfTobaccoDoYouUsuallyHaveByThisTimeOfTheDayOnA(partRunIdx);
                
                %% Grab alcohol use data from 7TZ
                % Total beers/glasses of wine on 7TZ day
                data.(data.runListUniqShort{5}).alcohol.beerWineCurrDay(zCounter) = ...
                    data.rawCSV.x12_HowManyBeersAndGlassesOfWineHaveYouHadToday_(partRunIdx);
                % Total beers/glasses of wine in past 24hrs of 7TZ day
                data.(data.runListUniqShort{5}).alcohol.beerWinePast24(zCounter) = ...
                    data.rawCSV.x14_HowManyBeersAndGlassesOfWineHaveYouHadInTheLast24Hours_(partRunIdx);
                % Total liquor on 7TZ day
                data.(data.runListUniqShort{5}).alcohol.beerWineCurrDay(zCounter) = ...
                    data.rawCSV.x13_HowManyMixedDrinksOrShotsOfLiquorHaveYouHadToday_(partRunIdx);
                % Total liquor in past 24hrs of 7TZ day
                data.(data.runListUniqShort{5}).alcohol.beerWinePast24(zCounter) = ...
                    data.rawCSV.x15_HowManyMixedDrinksOrShotsOfLiquorHaveYouHadInTheLast24Hours(partRunIdx);
                % Last consumption of alcohol
                data.(data.runListUniqShort{5}).alcohol.lastConsumption(zCounter) = ...
                    yyyymmdd(datetime(data.rawCSV.x16_WhenWasYourLastConsumptionOfAlcohol_mm_dd_yyyy__(partRunIdx)));
                data.(data.runListUniqShort{5}).alcohol.lastConsumptionAmount(zCounter) = ...
                    data.rawCSV.x17_HowMuchAlcoholDidYouDrink__OfStandardDrinks_DuringYourLastC(partRunIdx);
                % Most consumed in last month
                data.(data.runListUniqShort{5}).alcohol.mostConsumed(zCounter) = ...
                    data.rawCSV.x18_WhatIsTheMostAlcohol__OfStandardDrinks_YouHaveConsumedInOne(partRunIdx);
                % Average days consumed in a week
                data.(data.runListUniqShort{5}).alcohol.aveDaysConsumed(zCounter) = ...
                    data.rawCSV.x19_InATypicalWeek_HowManyDaysDoYouHaveSomeKindOfAlcoholToDrink(partRunIdx);
                % 7 Day breakdown
                for iK = 1:7   % For the 7 DOW
                    data.(data.runListUniqShort{5}).alcohol.dailyUse(zCounter,iK) = ...
                        data.rawCSV.(sprintf('%s%s','x20_',daysOtW{iK}))(partRunIdx);
                end
                data.(data.runListUniqShort{5}).alcohol.aveConsumPerDay(zCounter) = ...
                    squeeze(nanmean(data.(data.runListUniqShort{5}).alcohol.dailyUse(zCounter,:),2));
                % Average consumption / week
                data.(data.runListUniqShort{5}).alcohol.aveConsumPerWeek(zCounter) = ...
                    nansum(data.(data.runListUniqShort{5}).alcohol.dailyUse(zCounter,:),2);
            end
        end
    end
end
clear clinCounter aCounter bCounter zCounter



end