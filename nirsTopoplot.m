%   Author: Dino Soldic
%   Email: dino.soldic@urjc.es
%   Date: 2026-01-12
%
%   See also NIRSAnalysis

function nirsTopoplot(plotData, time)

    if nargin < 1
        [dataFile, dataDir] = uigetfile(".mat", "ALLDATA File", "ALLDATA.mat", "MultiSelect", "off");
        if dataFile == 0, error("Channel coordinates need to be selected"), end

        load(fullfile(dataDir, dataFile), "ALLDATA");
        plotData = ALLDATA.ALLDATATASK;
        time = ALLDATA.time;
    end

    try
        eeglab('nogui');
    catch
        errordlg("EEGLAB toolbox missing or not working", "EEGLAB Error")
        return
    end

    %% load standard channels
    [chansFile, chansDir] = uigetfile("*", "Channel Coordinates", "Standard_Channels.txt", "MultiSelect", "off");
    if chansFile == 0, error("Channel coordinates need to be selected"), end
    chans = readcell(fullfile(chansDir, chansFile));

    shChans = [chans{:, 3}] < 8; % remove short CH
    chans = chans(shChans, :);

    nSens = numel(chans(:, 2));

    labels = cellfun(@(x, y) [num2str(x), '-', num2str(y)], chans(:, 2), chans(:, 3), 'UniformOutput', false);

    %% Adapt nirs to eeg chans for topoplot
    % temp xyz to load
    fidPath = sprintf("%s/chanlocs.xyz", pwd);
    fid = fopen(fidPath, 'w');

    for i = 1:nSens
        fprintf(fid, "%d %.4f %.4f %.4f %s\n", i, chans{i, 4}, chans{i, 5}, chans{i, 6}, labels{i});
    end

    fclose(fid);

    % load and delete temp
    chanlocs = struct();
    chanlocs = pop_chanedit(chanlocs, 'load', fidPath);
    delete(fidPath);

    %% prep data
    % Extract labels
    groupNames = fieldnames(plotData);
    condNames = fieldnames(plotData.(groupNames{1}));
    dataHeaders = split(fieldnames(plotData.(groupNames{1}).(condNames{1})), "_");

    dataTypes = unique(dataHeaders(:, 1));

    chanLabels = unique(strcat(dataHeaders(:, 2), '-', dataHeaders(:, 3)));
    splitLabels = split(chanLabels, '-');
    num1 = str2double(splitLabels(:, 1));
    num2 = str2double(splitLabels(:, 2));
    [~, sortIdx] = sortrows([num1, num2]);
    chanLabels = chanLabels(sortIdx); % sorted

    % Select data and chan
    [dataOpt, ~] = listdlg('ListString', dataTypes, 'PromptString', {'Select the data type for', 'topoplot:'}, 'SelectionMode', 'multiple');
    if isempty(dataOpt), disp('Operation canceled. Shutting down'); return, end

    [chanOpt, ~] = listdlg('ListString', chanLabels, 'PromptString', {'Select channels for', 'topoplot:'}, 'SelectionMode', 'multiple');
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

    % extract chans for topoplot
    chansToPlot = ismember({chanlocs.labels}, selectedChannels);

    % substract?
    nCond = numel(condNames);

    doSubstractCond = questdlg('Do you wish to substract any condition?', 'Condition Substraction', 'Yes', 'No', 'Yes');
    if strcmp(doSubstractCond, 'Yes'), doSubstractCond = true; else, doSubstractCond = false; end

    if doSubstractCond

        while true
            [subCondOpt, ~] = listdlg('ListString', condNames, 'PromptString', {'Select condition to substract', 'from:'}, 'SelectionMode', 'single');
            if ~isempty(subCondOpt), break, end
        end

        subCond = condNames{subCondOpt};
        condNames(subCondOpt) = [];
        nCond = numel(condNames);
    end

    if numel(groupNames) > 1
        doGroupSubstract = questdlg('Do you wish to also plot group differences?', 'Group Substraction', 'Yes', 'No', 'Yes');
        if strcmp(doGroupSubstract, 'Yes'), doGroupSubstract = true; else, doGroupSubstract = false; end
    end

    doGroupMean = questdlg('Do you wish to also plot group means?', 'Group Means', 'Yes', 'No', 'Yes');
    if strcmp(doGroupMean, 'Yes'), doGroupMean = true; else, doGroupMean = false; end

    nGrp = numel(groupNames);

    %% plot
    nRows = ceil(sqrt(nCond));
    nCols = 4;

    for grpIdx = 1:nGrp

        topoFig = figure('Name', sprintf('%s', groupNames{grpIdx}), 'NumberTitle', 'off', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        tiled_topoFig = tiledlayout(topoFig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'loose');

        % determine color limits
        allVals = [];

        for condIdx = 1:nCond

            for chanIdx = 1:numel(chanCombos)
                vals = mean(plotData.(groupNames{grpIdx}).(condNames{condIdx}).(chanCombos{chanIdx}), 2);

                if doSubstractCond
                    vals = vals - mean(plotData.(groupNames{grpIdx}).(subCond).(chanCombos{chanIdx}), 2);
                end

                allVals = [allVals; vals]; %#ok<*AGROW>
            end

        end

        % --- Symmetric + rounded color limits ---
        absMax = max(abs(allVals)); % biggest absolute value
        scale = 10 ^ floor(log10(absMax)); % power-of-ten scaling
        niceMax = floor(absMax / scale) * scale / 2; % round to nearest 0.5 step
        colorLimits = [-niceMax, niceMax];

        % plot each condition
        for condIdx = 1:nCond

            dataToPlot = zeros(numel(chanCombos), length(time));

            for chanIdx = 1:numel(chanCombos)
                vals = mean(plotData.(groupNames{grpIdx}).(condNames{condIdx}).(chanCombos{chanIdx}), 2);

                if doSubstractCond
                    vals = vals - mean(plotData.(groupNames{grpIdx}).(subCond).(chanCombos{chanIdx}), 2);
                end

                dataToPlot(chanIdx, :) = vals;
            end

            nexttile(tiled_topoFig, condIdx);
            topoplot(mean(dataToPlot, 2), chanlocs(chansToPlot), 'emarker', {'.', 'k', 18, 1});
            clim(colorLimits);
            title(condNames{condIdx}, 'FontName', 'Times New Roman', 'FontSize', 24);
        end

        figure(topoFig);

        cb = colorbar('Position', [0.04, 0.6, 0.01, 0.3]);
        cb.FontName = 'Times New Roman';
        cb.FontSize = 24;
        cb.Limits = colorLimits;
        cb.Ticks = [colorLimits(1), 0, colorLimits(2)];

    end

    if doGroupSubstract

        topoFig = figure('Name', 'Group difference', 'NumberTitle', 'off', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        tiled_topoFig = tiledlayout(topoFig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'loose');

        % determine color limits
        allVals = [];

        for condIdx = 1:nCond

            for chanIdx = 1:numel(chanCombos)

                if doSubstractCond
                    vals = (mean(plotData.(groupNames{1}).(condNames{condIdx}).(chanCombos{chanIdx}), 2) ...
                        - mean(plotData.(groupNames{1}).(subCond).(chanCombos{chanIdx}), 2)) ...
                        - (mean(plotData.(groupNames{2}).(condNames{condIdx}).(chanCombos{chanIdx}), 2) ...
                        - mean(plotData.(groupNames{2}).(subCond).(chanCombos{chanIdx}), 2));
                else
                    vals = mean(plotData.(groupNames{1}).(condNames{condIdx}).(chanCombos{chanIdx}), 2) ...
                        - mean(plotData.(groupNames{2}).(condNames{condIdx}).(chanCombos{chanIdx}), 2);
                end

                allVals = [allVals; vals];
            end

        end

        % --- Symmetric + rounded color limits ---
        absMax = max(abs(allVals)); % biggest absolute value
        scale = 10 ^ floor(log10(absMax)); % power-of-ten scaling
        niceMax = floor(absMax / scale) * scale / 2; % round to half
        colorLimits = [-niceMax, niceMax];

        % plot each condition
        for condIdx = 1:nCond

            dataToPlot = zeros(numel(chanCombos), length(time));

            for chanIdx = 1:numel(chanCombos)

                if doSubstractCond
                    vals = (mean(plotData.(groupNames{1}).(condNames{condIdx}).(chanCombos{chanIdx}), 2) ...
                        - mean(plotData.(groupNames{1}).(subCond).(chanCombos{chanIdx}), 2)) ...
                        - (mean(plotData.(groupNames{2}).(condNames{condIdx}).(chanCombos{chanIdx}), 2) ...
                        - mean(plotData.(groupNames{2}).(subCond).(chanCombos{chanIdx}), 2));
                else
                    vals = mean(plotData.(groupNames{1}).(condNames{condIdx}).(chanCombos{chanIdx}), 2) ...
                        - mean(plotData.(groupNames{2}).(condNames{condIdx}).(chanCombos{chanIdx}), 2);
                end

                dataToPlot(chanIdx, :) = vals;
            end

            nexttile(tiled_topoFig, condIdx);
            topoplot(mean(dataToPlot, 2), chanlocs(chansToPlot), 'emarker', {'.', 'k', 18, 1});
            clim(colorLimits);
            title(condNames{condIdx}, 'FontName', 'Times New Roman', 'FontSize', 24);
        end

        figure(topoFig);

        cb = colorbar('Position', [0.04, 0.6, 0.01, 0.3]);
        cb.FontName = 'Times New Roman';
        cb.FontSize = 24;
        cb.Limits = colorLimits;
        cb.Ticks = [colorLimits(1), 0, colorLimits(2)];
    end

    % group mean of conditions
    if doGroupMean

        for grpIdx = 1:nGrp

            topoFig = figure('Name', sprintf('Mean of %s', groupNames{grpIdx}), 'NumberTitle', 'off', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

            dataToPlot = zeros(numel(chanCombos), length(time));
            allVals = [];

            for chanIdx = 1:numel(chanCombos)

                valsAllCond = zeros(nCond, length(time));

                for condIdx = 1:nCond
                    vals = mean(plotData.(groupNames{grpIdx}).(condNames{condIdx}).(chanCombos{chanIdx}), 2);

                    if doSubstractCond
                        vals = vals - mean(plotData.(groupNames{grpIdx}).(subCond).(chanCombos{chanIdx}), 2);
                    end

                    valsAllCond(condIdx, :) = vals;
                end

                % average across conditions
                meanVals = mean(valsAllCond, 1);

                dataToPlot(chanIdx, :) = meanVals;
                allVals = [allVals; meanVals(:)];
            end

            % --- Symmetric + rounded color limits ---
            absMax = max(abs(allVals)); % biggest absolute value
            scale = 10 ^ floor(log10(absMax)); % power-of-ten scaling
            niceMax = floor(absMax / scale) * scale / 2; % round to nearest 0.5 step
            colorLimits = [-niceMax, niceMax];

            % plot
            topoplot(mean(dataToPlot, 2), chanlocs(chansToPlot), 'emarker', {'.', 'k', 18, 1});
            clim(colorLimits);

            figure(topoFig);

            cb = colorbar('Position', [0.04, 0.6, 0.01, 0.3]);
            cb.FontName = 'Times New Roman';
            cb.FontSize = 24;
            cb.Limits = colorLimits;
            cb.Ticks = [colorLimits(1), 0, colorLimits(2)];

        end

    end

end
