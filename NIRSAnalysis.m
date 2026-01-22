%   Perform NIRS statistical analyses on task or rest (coming soon) data.
%
%   NIRSAnalysis() prompts the user to select a folder containing NIRS
%   data files. It loads all subject/condition data, allows the user to
%   select data type(s) and channels, choose analyses, and then computes
%   various Mass Univariate Permutation-based t-tests.
%
%   NIRSAnalysis(ALLDATA) uses a preloaded ALLDATA struct containing
%   fields:
%       ALLDATATASK - structured task data
%       ALLDATAREST - structured rest data
%       time        - time vector
%
%   The directory MUST be organized as follows:
%
%     Study
%     ├── Group1
%     │   ├── Condition1
%     │   │   ├── data1
%     │   │   └── dataN
%     │   └── Condition2
%     │       ├── data1
%     │       └── dataN
%     ├── Group2
%     │   ├── Condition1
%     │   │   ├── data1
%     │   │   └── dataN
%     │   └── Condition2
%     │       ├── data1
%     │       └── dataN
%     └── Group3
%         ├── Condition1
%         │   ├── data1
%         │   └── dataN
%         └── Condition2
%             ├── data1
%             └── dataN
%
%   When prompted the Study folder should be selected to load data from.
%
%   Inputs:
%       ALLDATA (optional) - Struct containing preloaded task/rest data
%
%   Outputs:
%       Results are saved to files in user-selected directories.
%       The following analysis types are supported:
%           1. Mass Univariate Independent T-test - condition
%           2. Mass Univariate Independent T-test - group
%           3. Mass Univariate Dependent T-test - within groups
%           4. Mass Univariate Dependent T-test - across conditions (no groups)
%           5. Average Mass Univariate Independent T-test - condition
%           6. Average Mass Univariate Independent T-test - group
%           7. Average Mass Univariate Dependent T-test - within groups
%           8. Average Mass Univariate Dependent T-test - across conditions (no groups)
%
%   Notes:
%       - Subtraction of a baseline condition is supported for all analysis types.
%       - Permutation testing (Monte Carlo, 5000 permutations by default) is used
%         to compute p-values for robust statistical inference.
%       - Parallel computing is leveraged if MATLAB Parallel Toolbox is available.
%
%   Example usage:
%       % Run analysis from scratch
%       NIRSAnalysis();
%
%       % Run analysis using preloaded data
%       load('ALLDATA.mat');
%       NIRSAnalysis(ALLDATA);
%
%   Author: Dino Soldic
%   Email: dino.soldic@urjc.es
%   Date: 2025-10-29
%
%   See also plotNIRS, exportNIRS

function NIRSAnalysis(ALLDATA)

    % load data
    if nargin < 1

        % get data path
        rootPath = uigetdir(pwd, "Select folder with NIRS data");
        if rootPath == 0, error("Operation Canceled"), end

        ALLDATATASK = struct();
        ALLDATAREST = struct();

        [ALLDATATASK, ALLDATAREST, time] = loadData(rootPath, ALLDATATASK, ALLDATAREST);

        disp("✅ Done loading data.");
    else
        ALLDATATASK = ALLDATA.ALLDATATASK;
        ALLDATAREST = ALLDATA.ALLDATAREST;
        time = ALLDATA.time;
    end

    % select task
    analDataOpt = questdlg('Do you wish to analyze rest or task data?', 'Data Selection', 'Task', 'Rest', 'Task');
    if strcmp(analDataOpt, 'Task'), AnalData = ALLDATATASK; isTask = true; elseif strcmp(analDataOpt, 'Rest'), AnalData = ALLDATAREST; isTask = false; else, disp('Operation canceled. Shutting down'); return, end

    % save results path
    saveDir = uigetdir(pwd, "Select folder where results will be saved");
    if saveDir == 0, saveDir = pwd; end

    saveDirCSV = [saveDir, '\csv'];
    if ~isfolder(saveDirCSV), mkdir(saveDirCSV); end

    % Extract labels
    groupNames = fieldnames(AnalData);
    condNames = fieldnames(AnalData.(groupNames{1}));
    dataHeaders = split(fieldnames(AnalData.(groupNames{1}).(condNames{1})), "_");

    dataTypes = unique(dataHeaders(:, 1));

    chanLabels = unique(strcat(dataHeaders(:, 2), '-', dataHeaders(:, 3)));
    splitLabels = split(chanLabels, '-');
    num1 = str2double(splitLabels(:, 1));
    num2 = str2double(splitLabels(:, 2));
    [~, sortIdx] = sortrows([num1, num2]);
    chanLabels = chanLabels(sortIdx); % sorted

    % Select data and chan
    [dataOpt, ~] = listdlg('ListString', dataTypes, 'PromptString', 'Select the data type to analyze:', 'SelectionMode', 'multiple');
    if isempty(dataOpt), disp('Operation canceled. Shutting down'); return, end

    [chanOpt, ~] = listdlg('ListString', chanLabels, 'PromptString', 'Select channels to analyze:', 'SelectionMode', 'multiple');
    if isempty(chanOpt), disp('Operation canceled. Shutting down'); return, end

    % Mix types and chans
    selectedDataTypes = dataTypes(dataOpt);
    selectedChannels = chanLabels(chanOpt);
    chanCombos = strings(length(dataOpt) * length(chanOpt), 1);
    idx = 1;

    for d = 1:length(dataOpt)

        for c = 1:length(chanOpt)
            chanCombos(idx) = selectedDataTypes(d) + "_" + replace(selectedChannels(c), "-", "_");
            idx = idx + 1;
        end

    end

    % Select analysis
    if isTask
        analTypes = {"Compare Conditions (between groups)", "Compare Groups (fuse conditions)", "Compare Conditions (within groups)", "Compare Conditions (fuse groups)", ...
                         "Avg. Compare Conditions (between groups)", "Avg. Compare Groups (fuse conditions)", "Avg. Compare Conditions (within groups)", "Avg. Compare Conditions (fuse groups)"};
        [analOpt, ~] = listdlg('ListString', analTypes, 'PromptString', 'Select analyses:', 'SelectionMode', 'multiple');
        if isempty(analOpt), disp('Operation canceled. Shutting down'); return, end

        [condOpt, ~] = listdlg('ListString', condNames, 'PromptString', {'Select condition for Mass', 'Univariate T test:'}, 'SelectionMode', 'multiple');
        if isempty(condOpt), disp('Operation canceled. Shutting down'); return, end
        condOpt = sort(condOpt); % force lowest to highest

        doSubstractCond = questdlg('Do you wish to substract any condition?', 'Condition Substraction', 'Yes', 'No', 'Yes');
        if strcmp(doSubstractCond, 'Yes'), doSubstractCond = true; else, doSubstractCond = false; end

        if doSubstractCond

            while true
                [subCondOpt, ~] = listdlg('ListString', condNames, 'PromptString', {'Select condition to substract', 'from:'}, 'SelectionMode', 'single');
                if ~isempty(subCondOpt), break, end
            end

        end

    end

    % amount of permutation for montecarlo p values (5000 good - Maris & Oostenveld 2007)
    nPerm = 5000;
    p = 0.05;

    try

        if isempty(gcp('nocreate'))
            myCluster = parcluster('local');
            maxWorkers = myCluster.NumWorkers;
            parpool('local', maxWorkers);
        end

        useParallel = true;
        disp("Parallel pool initialized with " + num2str(gcp().NumWorkers) + " workers.");
    catch ME
        warning(ME.identifier, 'Parallel pool could not be started. Running in serial mode.\nReason: %s', ME.message);
        useParallel = false;
    end

    %% run analyses
    % mass uni ind T (muit) - condition
    % Compares condition effects between groups
    if any(analOpt == 1)

        % preallocate
        nCond = length(condOpt);
        nChan = numel(chanCombos);
        nTime = length(time);

        tVals = zeros(nChan, nTime, nCond);
        pVals = zeros(nChan, nTime, nCond);
        hVals = zeros(nChan, nTime, nCond);

        for condIdx = 1:nCond
            cond = condNames{condOpt(condIdx)};

            allG1 = []; allG2 = [];

            for chanIdx = 1:nChan
                chan = chanCombos{chanIdx};

                if doSubstractCond
                    g1 = AnalData.(groupNames{1}).(cond).(chan) - AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan);
                    g2 = AnalData.(groupNames{2}).(cond).(chan) - AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan);
                else
                    g1 = AnalData.(groupNames{1}).(cond).(chan);
                    g2 = AnalData.(groupNames{2}).(cond).(chan);
                end

                allG1(:, :, chanIdx) = g1;
                allG2(:, :, chanIdx) = g2;
            end

            % real t values
            for c = 1:nChan

                for t = 1:nTime
                    [~, ~, ~, stats] = ttest2(allG1(t, :, c), allG2(t, :, c));
                    tVals(c, t, condIdx) = stats.tstat;
                end

            end

            % null t values
            combinedData = cat(2, allG1, allG2); % Combine groups along subject dimension
            nSubj1 = size(allG1, 2);
            maxT_Null = zeros(nPerm, 1);

            if useParallel

                parfor permIdx = 1:nPerm
                    % Shuffle subject labels across the whole dataset
                    permIdxs = randperm(size(combinedData, 2));
                    permG1 = combinedData(:, permIdxs(1:nSubj1), :);
                    permG2 = combinedData(:, permIdxs(nSubj1 + 1:end), :);

                    % Local t-stats for this shuffle
                    shuffTs = zeros(nTime, nChan);

                    for c = 1:nChan

                        for t = 1:nTime
                            % Standard t-test on shuffled subtracted data
                            m1 = mean(permG1(t, :, c)); m2 = mean(permG2(t, :, c));
                            v1 = var(permG1(t, :, c)); v2 = var(permG2(t, :, c));
                            n1 = nSubj1; n2 = size(combinedData, 2) - nSubj1;
                            % Fast t-stat formula for parfor speed
                            t_stat = (m1 - m2) / sqrt(v1 / n1 + v2 / n2);
                            shuffTs(t, c) = t_stat;
                        end

                    end

                    % Store only the highest absolute T from the entire shuffle
                    maxT_Null(permIdx) = max(abs(shuffTs(:)));
                end

            else

                for permIdx = 1:nPerm
                    % Shuffle subject labels across the whole dataset
                    permIdxs = randperm(size(combinedData, 2));
                    permG1 = combinedData(:, permIdxs(1:nSubj1), :);
                    permG2 = combinedData(:, permIdxs(nSubj1 + 1:end), :);

                    % Local t-stats for this shuffle
                    shuffTs = zeros(nTime, nChan);

                    for c = 1:nChan

                        for t = 1:nTime
                            % Standard t-test on shuffled subtracted data
                            m1 = mean(permG1(t, :, c)); m2 = mean(permG2(t, :, c));
                            v1 = var(permG1(t, :, c)); v2 = var(permG2(t, :, c));
                            n1 = nSubj1; n2 = size(combinedData, 2) - nSubj1;
                            % Fast t-stat formula for parfor speed
                            t_stat = (m1 - m2) / sqrt(v1 / n1 + v2 / n2);
                            shuffTs(t, c) = t_stat;
                        end

                    end

                    % Store only the highest absolute T from the entire shuffle
                    maxT_Null(permIdx) = max(abs(shuffTs(:)));
                end

            end

            % compare
            for c = 1:nChan

                for t = 1:nTime
                    pVals(c, t, condIdx) = mean(maxT_Null >= abs(tVals(c, t, condIdx)));
                    hVals(c, t, condIdx) = pVals(c, t, condIdx) < p;
                end

            end

        end

        disp("Mass-Uni Independent T-test for condition " + cond + " done");

        % Make results file and save
        results = struct();
        results.type.analysis = "muit-condition";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.group.code = [];
        results.group.name = [];
        results.condition.code = [];
        results.condition.name = [];
        results.channel.code = [];
        results.channel.name = [];
        results.time = time;

        for condIdx = 1:nCond
            results.condition.code = [results.condition.code; condIdx];
            results.condition.name = [results.condition.name; string(condNames{condOpt(condIdx)})];
        end

        for chanIdx = 1:nChan
            results.channel.code = [results.channel.code; chanIdx];
            results.channel.name = [results.channel.name; string(chanCombos{chanIdx})];
        end

        for grpIdx = 1:numel(groupNames)
            results.group.code = [results.group.code; grpIdx];
            results.group.name = [results.group.name; string(groupNames{grpIdx})];
        end

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_condition_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_condition_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_condition.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_condition.mat"), "results");
            end

        end

    end

    % mass uni ind T (muit) - group
    % Ignores condition separation and checks whether there's an effect group-wise
    if any(analOpt == 2)

        % preallocate
        nCond = numel(condNames(condOpt));
        nChan = numel(chanCombos);
        nTime = length(time);

        allG1 = [];
        allG2 = [];

        for chanIdx = 1:nChan
            chan = chanCombos{chanIdx};
            g1_tmp = []; g2_tmp = [];

            for condIdx = 1:nCond

                if doSubstractCond
                    dat1 = AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan) - ...
                        AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan);
                    dat2 = AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan) - ...
                        AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan);
                else
                    dat1 = AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan);
                    dat2 = AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan);
                end

                g1_tmp = [g1_tmp, dat1]; %#ok<*AGROW>
                g2_tmp = [g2_tmp, dat2];
            end

            allG1(:, :, chanIdx) = g1_tmp;
            allG2(:, :, chanIdx) = g2_tmp;
        end

        % Calculate Real T - Values
        tVals = zeros(nChan, nTime);

        for c = 1:nChan

            for t = 1:nTime
                [~, ~, ~, stats] = ttest2(allG1(t, :, c), allG2(t, :, c));
                tVals(c, t) = stats.tstat;
            end

        end

        % Max - T Permutation
        combinedData = cat(2, allG1, allG2);
        nSubj1 = size(allG1, 2);
        nTotal = size(combinedData, 2);
        maxT_Null = zeros(nPerm, 1);

        if useParallel % run in parallel if toolbox installed for performance

            % compute permuted t tests
            parfor permIdx = 1:nPerm
                permIdxs = randperm(nTotal);
                pG1 = combinedData(:, permIdxs(1:nSubj1), :);
                pG2 = combinedData(:, permIdxs(nSubj1 + 1:end), :);

                % Calculate all T-stats for this shuffle
                m1 = mean(pG1, 2); m2 = mean(pG2, 2);
                v1 = var(pG1, 0, 2); v2 = var(pG2, 0, 2);

                shuffTs = (m1 - m2) ./ sqrt(v1 / nSubj1 + v2 / (nTotal - nSubj1));

                % Store only the absolute maximum T found anywhere in this shuffle
                maxT_Null(permIdx) = max(abs(shuffTs(:)));
            end

        else

            for permIdx = 1:nPerm
                permIdxs = randperm(nTotal);
                pG1 = combinedData(:, permIdxs(1:nSubj1), :);
                pG2 = combinedData(:, permIdxs(nSubj1 + 1:end), :);

                % Calculate all T-stats for this shuffle
                m1 = mean(pG1, 2); m2 = mean(pG2, 2);
                v1 = var(pG1, 0, 2); v2 = var(pG2, 0, 2);

                shuffTs = (m1 - m2) ./ sqrt(v1 / nSubj1 + v2 / (nTotal - nSubj1));

                % Store only the absolute maximum T found anywhere in this shuffle
                maxT_Null(permIdx) = max(abs(shuffTs(:)));
            end

        end

        % Compare
        pVals = zeros(nChan, nTime);

        for c = 1:nChan

            for t = 1:nTime
                pVals(c, t) = mean(maxT_Null >= abs(tVals(c, t)));
            end

        end

        hVals = pVals < p;

        labels = string(strjoin(condNames(condOpt), "_"));

        disp("Mass-Uni Independent T-test for groups done");

        % Make results file and save
        results = struct();
        results.type.analysis = "muit-group";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.group.code = [];
        results.group.name = [];
        results.condition.code = [];
        results.condition.name = [];
        results.channel.code = [];
        results.channel.name = [];
        results.time = time;

        for grpIdx = 1:numel(groupNames)
            results.group.code = [results.group.code; grpIdx];
            results.group.name = [results.group.name; string(groupNames{grpIdx})];
        end

        results.condition.code = 1;
        results.condition.name = labels;

        for chanIdx = 1:nChan
            results.channel.code = [results.channel.code; chanIdx];
            results.channel.name = [results.channel.name; string(chanCombos{chanIdx})];
        end

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_group_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_group_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_group.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_group.mat"), "results");
            end

        end

    end

    % mass uni dep T (mudt) - within groups
    % Compares within each group the effect of condition vs another
    if any(analOpt == 3)

        % preallocate
        nCond = numel(condNames(condOpt));
        if nCond == 2, nCondMax = 1; else, nCondMax = nCond; end
        nGrp = numel(groupNames);
        nChan = numel(chanCombos);
        nTime = length(time);

        tVals = zeros(nChan, nTime, nCondMax, nGrp);
        pVals = zeros(nChan, nTime, nCondMax, nGrp);
        hVals = zeros(nChan, nTime, nCondMax, nGrp);

        labels = strings(nCondMax, 1);

        for grpIdx = 1:nGrp
            group = groupNames{grpIdx};

            for condIdx = 1:nCondMax
                condIdx2 = condIdx + 1;
                if condIdx2 > nCond, condIdx2 = 1; end
                labels(condIdx) = sprintf("%s_%s", condNames{condOpt(condIdx)}, condNames{condOpt(condIdx2)});

                diffData = [];

                for chanIdx = 1:nChan
                    chan = chanCombos{chanIdx};

                    if doSubstractCond
                        g1 = AnalData.(group).(condNames{condOpt(condIdx)}).(chan) - AnalData.(group).(condNames{subCondOpt}).(chan);
                        g2 = AnalData.(group).(condNames{condOpt(condIdx2)}).(chan) - AnalData.(group).(condNames{subCondOpt}).(chan);
                    else
                        g1 = AnalData.(group).(condNames{condOpt(condIdx)}).(chan);
                        g2 = AnalData.(group).(condNames{condOpt(condIdx2)}).(chan);
                    end

                    diffData(:, :, chanIdx) = g1 - g2; % The paired difference
                end

                % Calculate Real T - Stats
                currentTVals = zeros(nTime, nChan);

                for c = 1:nChan

                    for t = 1:nTime
                        [~, ~, ~, stats] = ttest(diffData(t, :, c));
                        currentTVals(t, c) = stats.tstat;
                    end

                end

                % Max - T Permutation (Sign - Flipping)
                nSubj = size(diffData, 2);
                maxT_Null = zeros(nPerm, 1);

                if useParallel

                    parfor permIdx = 1:nPerm
                        % Randomly flip signs for each subject (1 or -1)
                        flipSigns = (rand(1, nSubj) > 0.5) * 2 - 1;
                        % Multiply the same sign-flip across all time and channels for that subject
                        permDiffs = diffData .* reshape(flipSigns, [1, nSubj, 1]);

                        % Fast T-calculation: mean / (std/sqrt(n))
                        m = mean(permDiffs, 2);
                        s = std(permDiffs, 0, 2);
                        shuffTs = m ./ (s / sqrt(nSubj));

                        maxT_Null(permIdx) = max(abs(shuffTs(:)));
                    end

                else

                    for permIdx = 1:nPerm
                        % Randomly flip signs for each subject (1 or -1)
                        flipSigns = (rand(1, nSubj) > 0.5) * 2 - 1;
                        % Multiply the same sign-flip across all time and channels for that subject
                        permDiffs = diffData .* reshape(flipSigns, [1, nSubj, 1]);

                        % Fast T-calculation: mean / (std/sqrt(n))
                        m = mean(permDiffs, 2);
                        s = std(permDiffs, 0, 2);
                        shuffTs = m ./ (s / sqrt(nSubj));

                        maxT_Null(permIdx) = max(abs(shuffTs(:)));
                    end

                end

                % Compare
                for c = 1:nChan

                    for t = 1:nTime
                        tStat = currentTVals(t, c);
                        tVals(c, t, condIdx, grpIdx) = tStat;
                        pVals(c, t, condIdx, grpIdx) = mean(maxT_Null >= abs(tStat));
                    end

                end

                hVals(:, :, condIdx, grpIdx) = pVals(:, :, condIdx, grpIdx) < p;

                disp("Mass-Uni Dependent T-test for " + group + "-" + condNames{condOpt(condIdx)} + "_" + condNames{condOpt(condIdx2)} + " done");
            end

        end

        % Make results file and save
        results = struct();
        results.type.analysis = "mudt-within";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.group.code = [];
        results.group.name = [];
        results.condition.code = [];
        results.condition.name = [];
        results.channel.code = [];
        results.channel.name = [];
        results.time = time;

        for grpIdx = 1:numel(groupNames)
            results.group.code = [results.group.code; grpIdx];
            results.group.name = [results.group.name; string(groupNames{grpIdx})];
        end

        for condIdx = 1:nCondMax
            results.condition.code = [results.condition.code; condIdx];
            results.condition.name = labels;
        end

        for chanIdx = 1:nChan
            results.channel.code = [results.channel.code; chanIdx];
            results.channel.name = [results.channel.name; string(chanCombos{chanIdx})];
        end

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_withinSub_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_withinSub_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_withinSub.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_withinSub.mat"), "results");
            end

        end

    end

    % mass uni dep T (mudt) - no groups
    % Ignores group separation and checks whether there's an effect condition-wise
    if any(analOpt == 4)

        % preallocate
        nCond = numel(condNames(condOpt));
        if nCond == 2, nCondMax = 1; else, nCondMax = nCond; end
        nChan = numel(chanCombos);
        nTime = length(time);

        tVals = zeros(nChan, nTime, nCondMax);
        pVals = zeros(nChan, nTime, nCondMax);
        hVals = zeros(nChan, nTime, nCondMax);

        labels = strings(nCondMax, 1);

        for condIdx = 1:nCondMax
            condIdx2 = condIdx + 1;
            if condIdx2 > nCond, condIdx2 = 1; end
            labels(condIdx) = sprintf("%s_%s", condNames{condOpt(condIdx)}, condNames{condOpt(condIdx2)});

            diffData = [];

            for chanIdx = 1:nChan
                chan = chanCombos{chanIdx};

                % Pool Group 1 and Group 2
                if doSubstractCond
                    g1 = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan) - AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan), ...
                              AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan) - AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan)];
                    g2 = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx2)}).(chan) - AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan), ...
                              AnalData.(groupNames{2}).(condNames{condOpt(condIdx2)}).(chan) - AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan)];
                else
                    g1 = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan), AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan)];
                    g2 = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx2)}).(chan), AnalData.(groupNames{2}).(condNames{condOpt(condIdx2)}).(chan)];
                end

                diffData(:, :, chanIdx) = g1 - g2;
            end

            % Calculate Real T-Stats
            currentTVals = zeros(nTime, nChan);

            for c = 1:nChan

                for t = 1:nTime
                    [~, ~, ~, stats] = ttest(diffData(t, :, c));
                    currentTVals(t, c) = stats.tstat;
                end

            end

            % Max-T Permutation (Sign-Flipping pooled subjects)
            nSubjTotal = size(diffData, 2);
            maxT_Null = zeros(nPerm, 1);

            if useParallel

                parfor permIdx = 1:nPerm
                    % Flip signs for pooled subjects
                    flipSigns = (rand(1, nSubjTotal) > 0.5) * 2 - 1;
                    permDiffs = diffData .* reshape(flipSigns, [1, nSubjTotal, 1]);

                    % Fast T-stat
                    m = mean(permDiffs, 2);
                    s = std(permDiffs, 0, 2);
                    shuffTs = m ./ (s / sqrt(nSubjTotal));

                    maxT_Null(permIdx) = max(abs(shuffTs(:)));
                end

            else

                for permIdx = 1:nPerm
                    % Flip signs for pooled subjects
                    flipSigns = (rand(1, nSubjTotal) > 0.5) * 2 - 1;
                    permDiffs = diffData .* reshape(flipSigns, [1, nSubjTotal, 1]);

                    % Fast T-stat
                    m = mean(permDiffs, 2);
                    s = std(permDiffs, 0, 2);
                    shuffTs = m ./ (s / sqrt(nSubjTotal));

                    maxT_Null(permIdx) = max(abs(shuffTs(:)));
                end

            end

            % Compare Real vs Null
            for c = 1:nChan

                for t = 1:nTime
                    tStat = currentTVals(t, c);
                    tVals(c, t, condIdx) = tStat;
                    pVals(c, t, condIdx) = mean(maxT_Null >= abs(tStat));
                end

            end

            hVals(:, :, condIdx) = pVals(:, :, condIdx) < p;

            disp("Mass-Uni Independent T-test for " + condNames{condOpt(condIdx)} + "_" + condNames{condOpt(condIdx2)} + " done");
        end

        % Make results file and save
        results = struct();
        results.type.analysis = "mudt-nogroup";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.condition.code = [];
        results.condition.name = [];
        results.channel.code = [];
        results.channel.name = [];
        results.time = time;

        for condIdx = 1:nCondMax
            results.condition.code = [results.condition.code; condIdx];
            results.condition.name = labels;
        end

        for chanIdx = 1:nChan
            results.channel.code = [results.channel.code; chanIdx];
            results.channel.name = [results.channel.name; string(chanCombos{chanIdx})];
        end

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_nogroup_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_nogroup_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_nogroup.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_nogroup.mat"), "results");
            end

        end

    end

    %% Averages
    % AVERAGE mass uni ind T (muit) - condition
    % Compares condition effects between groups
    if any(analOpt == 5)

        % preallocate
        nCond = length(condOpt);
        nChan = numel(chanCombos);
        nTime = length(time);

        tVals = zeros(nTime, nCond);
        pVals = zeros(nTime, nCond);
        hVals = zeros(nTime, nCond);

        for condIdx = 1:nCond
            cond = condNames{condOpt(condIdx)};

            % Collect all channels to average them
            g1_allChans = [];
            g2_allChans = [];

            for chanIdx = 1:nChan
                chan = chanCombos{chanIdx};
                g1_allChans(:, :, chanIdx) = AnalData.(groupNames{1}).(cond).(chan);
                g2_allChans(:, :, chanIdx) = AnalData.(groupNames{2}).(cond).(chan);
            end

            % Compute the Average Across Channels
            avgG1 = mean(g1_allChans, 3);
            avgG2 = mean(g2_allChans, 3);

            if doSubstractCond
                g1_sub_all = []; g2_sub_all = [];

                for chanIdx = 1:nChan
                    chan = chanCombos{chanIdx};
                    g1_sub_all(:, :, chanIdx) = AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan);
                    g2_sub_all(:, :, chanIdx) = AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan);
                end

                avgG1 = avgG1 - mean(g1_sub_all, 3);
                avgG2 = avgG2 - mean(g2_sub_all, 3);
            end

            % Real T-stats
            realT = zeros(nTime, 1);

            for t = 1:nTime
                [~, ~, ~, stats] = ttest2(avgG1(t, :), avgG2(t, :));
                realT(t) = stats.tstat;
            end

            % Max-T Permutation
            combinedData = [avgG1, avgG2];
            nSubj1 = size(avgG1, 2);
            nTotal = size(combinedData, 2);
            maxT_Null = zeros(nPerm, 1);

            if useParallel

                parfor permIdx = 1:nPerm
                    permIdxs = randperm(nTotal);
                    pG1 = combinedData(:, permIdxs(1:nSubj1)); %#ok<*PFBNS>
                    pG2 = combinedData(:, permIdxs(nSubj1 + 1:end));

                    % Fast T-stat across all time points for this shuffle
                    m1 = mean(pG1, 2); m2 = mean(pG2, 2);
                    v1 = var(pG1, 0, 2); v2 = var(pG2, 0, 2);

                    shuffTs = (m1 - m2) ./ sqrt(v1 / nSubj1 + v2 / (nTotal - nSubj1));

                    % Max T across the time dimension
                    maxT_Null(permIdx) = max(abs(shuffTs));
                end

            else

                for permIdx = 1:nPerm
                    permIdxs = randperm(nTotal);
                    pG1 = combinedData(:, permIdxs(1:nSubj1));
                    pG2 = combinedData(:, permIdxs(nSubj1 + 1:end));

                    % Fast T-stat across all time points for this shuffle
                    m1 = mean(pG1, 2); m2 = mean(pG2, 2);
                    v1 = var(pG1, 0, 2); v2 = var(pG2, 0, 2);

                    shuffTs = (m1 - m2) ./ sqrt(v1 / nSubj1 + v2 / (nTotal - nSubj1));

                    % Max T across the time dimension
                    maxT_Null(permIdx) = max(abs(shuffTs));
                end

            end

            % Compare
            for t = 1:nTime
                tVals(t, condIdx) = realT(t);
                pVals(t, condIdx) = mean(maxT_Null >= abs(realT(t)));
            end

            hVals(:, condIdx) = pVals(:, condIdx) < p;

            disp("Average Mass-Uni Independent T-test for condition " + cond + " done");
        end

        % Make results file and save
        results = struct();
        results.type.analysis = "muit-condition-avg";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.group.code = [];
        results.group.name = [];
        results.condition.code = [];
        results.condition.name = [];
        results.channel.name = [];
        results.time = time;

        for grpIdx = 1:numel(groupNames)
            results.group.code = [results.group.code; grpIdx];
            results.group.name = [results.group.name; string(groupNames{grpIdx})];
        end

        for condIdx = 1:nCond
            results.condition.code = [results.condition.code; condIdx];
            results.condition.name = [results.condition.name; string(condNames{condOpt(condIdx)})];
        end

        results.channel.name = strjoin(chanCombos, ";");

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_avg_condition_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_condition_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_avg_condition.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_condition.mat"), "results");
            end

        end

    end

    % AVERAGE mass uni ind T (muit) - group
    % Ignores condition separation and checks whether there's an effect group-wise
    if any(analOpt == 6)

        % preallocate
        nCond = numel(condNames(condOpt));
        nChan = numel(chanCombos);
        nTime = length(time);

        g1_allChans = [];
        g2_allChans = [];

        for chanIdx = 1:nChan
            chan = chanCombos{chanIdx};
            g1_fused = [];
            g2_fused = [];

            for condIdx = 1:nCond

                if doSubstractCond
                    dat1 = AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan) - ...
                        AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan);
                    dat2 = AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan) - ...
                        AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan);
                else
                    dat1 = AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan);
                    dat2 = AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan);
                end

                g1_fused = [g1_fused, dat1];
                g2_fused = [g2_fused, dat2];
            end

            g1_allChans(:, :, chanIdx) = g1_fused;
            g2_allChans(:, :, chanIdx) = g2_fused;
        end

        % Average across channels
        avgG1 = mean(g1_allChans, 3);
        avgG2 = mean(g2_allChans, 3);

        % Calculate Real T-Stats
        realT = zeros(nTime, 1);

        for t = 1:nTime
            [~, ~, ~, stats] = ttest2(avgG1(t, :), avgG2(t, :));
            realT(t) = stats.tstat;
        end

        % Permutation
        combinedData = [avgG1, avgG2];
        nSubj1 = size(avgG1, 2);
        nTotal = size(combinedData, 2);
        maxT_Null = zeros(nPerm, 1);

        if useParallel

            parfor permIdx = 1:nPerm
                permIdxs = randperm(nTotal);
                pG1 = combinedData(:, permIdxs(1:nSubj1));
                pG2 = combinedData(:, permIdxs(nSubj1 + 1:end));

                % Fast T-stat calculation for all time points
                m1 = mean(pG1, 2); m2 = mean(pG2, 2);
                v1 = var(pG1, 0, 2); v2 = var(pG2, 0, 2);
                shuffTs = (m1 - m2) ./ sqrt(v1 / nSubj1 + v2 / (nTotal - nSubj1));

                % Capture Max absolute T across the time dimension
                maxT_Null(permIdx) = max(abs(shuffTs));
            end

        else

            for permIdx = 1:nPerm
                permIdxs = randperm(nTotal);
                pG1 = combinedData(:, permIdxs(1:nSubj1));
                pG2 = combinedData(:, permIdxs(nSubj1 + 1:end));
                m1 = mean(pG1, 2); m2 = mean(pG2, 2);
                v1 = var(pG1, 0, 2); v2 = var(pG2, 0, 2);
                shuffTs = (m1 - m2) ./ sqrt(v1 / nSubj1 + v2 / (nTotal - nSubj1));
                maxT_Null(permIdx) = max(abs(shuffTs));
            end

        end

        % Compare
        pVals = mean(maxT_Null >= abs(realT)', 1)';
        tVals = realT;
        hVals = pVals < p;

        labels = string(strjoin(condNames(condOpt), "_"));

        disp("Average Mass-Uni Independent T-test for " + labels + " done");

        % Make results file and save
        results = struct();
        results.type.analysis = "muit-group-avg";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.group.code = [];
        results.group.name = [];
        results.condition.code = [];
        results.condition.name = [];
        results.channel.name = [];
        results.time = time;

        for grpIdx = 1:numel(groupNames)
            results.group.code = [results.group.code; grpIdx];
            results.group.name = [results.group.name; string(groupNames{grpIdx})];
        end

        results.condition.code = 1;
        results.condition.name = labels;

        results.channel.name = strjoin(chanCombos, ";");

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_avg_group_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_group_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_avg_group.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_group.mat"), "results");
            end

        end

    end

    % AVERAGE mass uni dep T (mudt) - within groups
    % Compares within each group the effect of condition vs another
    if any(analOpt == 7)

        % preallocate
        nCond = numel(condNames(condOpt));
        if nCond == 2, nCondMax = 1; else, nCondMax = nCond; end
        nGrp = numel(groupNames);
        nChan = numel(chanCombos);
        nTime = length(time);

        tVals = zeros(nTime, nCondMax, nGrp);
        pVals = zeros(nTime, nCondMax, nGrp);
        hVals = zeros(nTime, nCondMax, nGrp);

        labels = strings(nCondMax, 1);

        for grpIdx = 1:nGrp
            group = groupNames{grpIdx};

            for condIdx = 1:nCondMax
                condIdx2 = condIdx + 1;
                if condIdx2 > nCond, condIdx2 = 1; end

                % Collect and Average Channels
                g1_all = []; g2_all = [];

                for chanIdx = 1:nChan
                    chan = chanCombos{chanIdx};

                    if doSubstractCond
                        g1_all(:, :, chanIdx) = AnalData.(group).(condNames{condOpt(condIdx)}).(chan) - AnalData.(group).(condNames{subCondOpt}).(chan);
                        g2_all(:, :, chanIdx) = AnalData.(group).(condNames{condOpt(condIdx2)}).(chan) - AnalData.(group).(condNames{subCondOpt}).(chan);
                    else
                        g1_all(:, :, chanIdx) = AnalData.(group).(condNames{condOpt(condIdx)}).(chan);
                        g2_all(:, :, chanIdx) = AnalData.(group).(condNames{condOpt(condIdx2)}).(chan);
                    end

                end

                avgG1 = mean(g1_all, 3);
                avgG2 = mean(g2_all, 3);
                labels(condIdx) = sprintf("%s_%s", condNames{condOpt(condIdx)}, condNames{condOpt(condIdx2)});

                % Real Paired T-stats for all time points
                realDiff = avgG1 - avgG2;
                nSubj = size(realDiff, 2);

                mu = mean(realDiff, 2);
                sd = std(realDiff, 0, 2);
                realT = mu ./ (sd ./ sqrt(nSubj));

                % Permutation
                maxT_Null = zeros(nPerm, 1);

                if useParallel

                    parfor permIdx = 1:nPerm
                        % Randomly flip signs of the difference for each subject
                        flipSigns = (rand(1, nSubj) > 0.5) * 2 - 1;
                        permDiff = realDiff .* flipSigns;

                        m_p = mean(permDiff, 2);
                        s_p = std(permDiff, 0, 2);
                        shuffTs = m_p ./ (s_p ./ sqrt(nSubj));

                        maxT_Null(permIdx) = max(abs(shuffTs));
                    end

                else

                    for permIdx = 1:nPerm
                        flipSigns = (rand(1, nSubj) > 0.5) * 2 - 1;
                        permDiff = realDiff .* flipSigns;
                        shuffTs = mean(permDiff, 2) ./ (std(permDiff, 0, 2) ./ sqrt(nSubj));
                        maxT_Null(permIdx) = max(abs(shuffTs));
                    end

                end

                % Compare
                tVals(:, condIdx, grpIdx) = realT;
                pVals(:, condIdx, grpIdx) = mean(maxT_Null >= abs(realT)', 1)';
                hVals(:, condIdx, grpIdx) = pVals(:, condIdx, grpIdx) < p;

                disp("Average Mass-Uni Dependent T-test for " + group + "-" + condNames{condOpt(condIdx)} + "_" + condNames{condOpt(condIdx2)} + " done");
            end

        end

        % Make results file and save
        results = struct();
        results.type.analysis = "mudt-within-avg";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.group.code = [];
        results.group.name = [];
        results.condition.code = [];
        results.condition.name = [];
        results.channel.name = [];
        results.time = time;

        for grpIdx = 1:numel(groupNames)
            results.group.code = [results.group.code; grpIdx];
            results.group.name = [results.group.name; string(groupNames{grpIdx})];
        end

        for condIdx = 1:nCondMax
            results.condition.code = [results.condition.code; condIdx];
            results.condition.name = labels;
        end

        results.channel.name = strjoin(chanCombos, ";");

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_avg_withinSub_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_withinSub_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_avg_withinSub.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_withinSub.mat"), "results");
            end

        end

    end

    % AVERAGE mass uni dep T (mudt) - no groups
    % Ignores group separation and checks whether there's an effect condition-wise
    if any(analOpt == 8)

        % preallocate
        nCond = numel(condNames(condOpt));
        if nCond == 2, nCondMax = 1; else, nCondMax = nCond; end
        nChan = numel(chanCombos);
        nTime = length(time);

        tVals = zeros(nTime, nCondMax);
        pVals = zeros(nTime, nCondMax);
        hVals = zeros(nTime, nCondMax);

        labels = strings(nCondMax, 1);

        for condIdx = 1:nCondMax
            condIdx2 = condIdx + 1;
            if condIdx2 > nCond, condIdx2 = 1; end

            nRows = size(AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chanCombos{1}), 1);
            totalSubj = size(AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chanCombos{1}), 2) + ...
                size(AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chanCombos{1}), 2);

            g1_pooled = zeros(nRows, totalSubj, nChan);
            g2_pooled = zeros(nRows, totalSubj, nChan);

            for chanIdx = 1:nChan
                chan = chanCombos{chanIdx};

                if doSubstractCond
                    dat1 = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan) - AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan), ...
                                AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan) - AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan)];
                    dat2 = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx2)}).(chan) - AnalData.(groupNames{1}).(condNames{subCondOpt}).(chan), ...
                                AnalData.(groupNames{2}).(condNames{condOpt(condIdx2)}).(chan) - AnalData.(groupNames{2}).(condNames{subCondOpt}).(chan)];
                else
                    dat1 = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx)}).(chan), AnalData.(groupNames{2}).(condNames{condOpt(condIdx)}).(chan)];
                    dat2 = [AnalData.(groupNames{1}).(condNames{condOpt(condIdx2)}).(chan), AnalData.(groupNames{2}).(condNames{condOpt(condIdx2)}).(chan)];
                end

                g1_pooled(:, :, chanIdx) = dat1;
                g2_pooled(:, :, chanIdx) = dat2;
            end

            % Average across channels
            avgG1 = mean(g1_pooled, 3);
            avgG2 = mean(g2_pooled, 3);

            % Real T-stats
            realDiff = avgG1 - avgG2;
            nSubj = size(realDiff, 2);

            realT = mean(realDiff, 2) ./ (std(realDiff, 0, 2) ./ sqrt(nSubj));

            % Permutation
            maxT_Null = zeros(nPerm, 1);

            if useParallel

                parfor permIdx = 1:nPerm
                    % Flip signs per subject (same flip across all time points)
                    flipSigns = (rand(1, nSubj) > 0.5) * 2 - 1;
                    permDiff = realDiff .* flipSigns;

                    % Fast T-stat
                    m_p = mean(permDiff, 2);
                    s_p = std(permDiff, 0, 2);
                    shuffTs = m_p ./ (s_p ./ sqrt(nSubj));

                    maxT_Null(permIdx) = max(abs(shuffTs));
                end

            else

                for permIdx = 1:nPerm
                    flipSigns = (rand(1, nSubj) > 0.5) * 2 - 1;
                    permDiff = realDiff .* flipSigns;
                    shuffTs = mean(permDiff, 2) ./ (std(permDiff, 0, 2) ./ sqrt(nSubj));
                    maxT_Null(permIdx) = max(abs(shuffTs));
                end

            end

            % Compare
            tVals(:, condIdx) = realT;
            pVals(:, condIdx) = mean(maxT_Null >= abs(realT)', 1)';
            hVals(:, condIdx) = pVals(:, condIdx) < p;

            labels(condIdx) = sprintf("%s_%s", condNames{condOpt(condIdx)}, condNames{condOpt(condIdx2)});

            disp("Average Mass-Uni Dependent T-test for " + condNames{condOpt(condIdx)} + "_" + condNames{condOpt(condIdx2)} + " done");
        end

        % Make results file and save
        results = struct();
        results.type.analysis = "mudt-nogroup-avg";
        results.type.data = string(analDataOpt);
        results.type.substracted = doSubstractCond;
        if doSubstractCond, results.type.substractedCond = condNames{subCondOpt}; end
        results.stats.t = tVals;
        results.stats.p = pVals;
        results.stats.h = hVals;
        results.condition.code = [];
        results.condition.name = [];
        results.channel.name = [];
        results.time = time;

        for condIdx = 1:nCondMax
            results.condition.code = [results.condition.code; condIdx];
            results.condition.name = labels;
        end

        results.channel.name = strjoin(chanCombos, ";");

        % export & save
        exportNIRS(results, saveDirCSV, time, doSubstractCond);

        if doSubstractCond

            if isTask
                save(fullfile(saveDir, "results_task_avg_nogroup_substracted.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_nogroup_substracted.mat"), "results");
            end

        else

            if isTask
                save(fullfile(saveDir, "results_task_avg_nogroup.mat"), "results");
            else
                save(fullfile(saveDir, "results_rest_avg_nogroup.mat"), "results");
            end

        end

    end

    % save data if no save exists
    if nargin < 1
        ALLDATA.ALLDATATASK = ALLDATATASK;
        ALLDATA.ALLDATAREST = ALLDATAREST;
        ALLDATA.time = time;

        save(fullfile(saveDir, "ALLDATA.mat"), "ALLDATA");
    end

    % Display completion
    fprintf('\n\t-------Analyses completed successfully-------\n');
    fprintf('\n\t\t  /\\_/\\ \t  /\\_/\\ \n\t\t ( o.o )\t ( ^.^ )\n\t\t  > ^ <\t\t  > ^ <\n');

end
