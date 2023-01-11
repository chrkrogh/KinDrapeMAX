function StoreResults(Inp,Plt,InpFilePath,x_opt,Dra_S,DT,F,z,OptRes_S,...
    Mold,ObjFun,Opt,OptProb,Set)

% Stop command window logging
diary('off');

if isempty(OptRes_S)
    % Data for display is not available
    return
end

% Promt the user to save figures
fprintf('\nSave results to new folder? \n')
Val = input('Enter name of output folder >> ','s');

if ~isempty(Val)
    
    % Import log as char array
    Log_imported = fileread('./Temporary files/log.txt');
    % Remove bold formatting from tables (trick)
    Log_format = regexprep(Log_imported, '<\x2F?strong>', '');

    % Get string with current date and time
    CurrentDate = char(datetime('now','Format','uuuu-MMMM-dd_HH-mm'));

    % Create output folder
    OutFolder = ['./Output/' CurrentDate '_' Val '/'];
    mkdir(OutFolder);

    % Write new log file to output folder and delete old log file
    LogPath = [OutFolder CurrentDate '_log.txt'];
    writematrix(Log_format,LogPath,'QuoteStrings','none');
    %delete('./Temporary files/log.txt')

    % Save figures as .fig and png (incl. maximization of fig window)
    if isfield(Plt,'ShearFig')
        savefig(Plt.ShearFig,[OutFolder 'ShearFig'],'compact')
        Plt.ShearFig.WindowState = 'max';
        print(Plt.ShearFig,[OutFolder 'ShearFig'],'-dpng','-r300')
        Plt.ShearFig.WindowState = 'normal';
    end

    if isfield(Plt,'FibDevFig')
        savefig(Plt.FibDevFig,[OutFolder 'FibDevFig'],'compact')
        Plt.FibDevFig.WindowState = 'max';
        print(Plt.FibDevFig,[OutFolder 'FibDevFig'],'-dpng','-r300')
        Plt.FibDevFig.WindowState = 'normal';
    end

    if isfield(Plt,'FlatFig')
        savefig(Plt.FlatFig,[OutFolder 'FlatFig'],'compact')
        Plt.FlatFig.WindowState = 'max';
        print(Plt.FlatFig,[OutFolder 'FlatFig'],'-dpng','-r300')
        Plt.FlatFig.WindowState = 'normal';
    end

    if isfield(Plt,'StackFig')
        savefig(Plt.StackFig,[OutFolder 'StackFig'],'compact')
        %Plt.StackFig.WindowState = 'max';
        print(Plt.StackFig,[OutFolder 'StackFig'],'-dpng','-r300')
        %Plt.StackFig.WindowState = 'normal';
    end

    % Copy input file
    copyfile([InpFilePath Inp.FileName],OutFolder)

    % Append a line with x_opt
    FileID = fopen([OutFolder Inp.FileName],'a');

    fprintf(FileID,['\n\n\nInp.x_opt = [' num2str(x_opt) '];']);

    fclose(FileID);

    % Save workspace variables
    % First remove figure handles from Plt struct (they were saved
    % seperately above)
    FigFields = {'ShearFig','FibDevFig','FlatFig','StackFig'};
    for i = 1:length(FigFields)
        if isfield(Plt,FigFields{i})
            Plt = rmfield(Plt,FigFields{i});
        end
    end

    save([OutFolder 'Workspace_vars'],...
        'Inp','Plt','Opt','OptRes_S','Dra_S','Set','x_opt','DT','F','z',...
        'Mold','ObjFun','OptProb');

    fprintf('\nLog, figures and workspace vars. succesfully saved.\n')
else
    % If no saving is requested, delete the temp log file
    %delete('./Temporary files/log.txt')
end

end