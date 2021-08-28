function varargout = ATN_gui(varargin)
% ATN_GUI MATLAB code for ATN_gui.fig
%      ATN_GUI, by itself, creates a new ATN_GUI or raises the existing
%      singleton*.
%
%      H = ATN_GUI returns the handle to a new ATN_GUI or the handle to
%      the existing singleton*.
%
%      ATN_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ATN_GUI.M with the given input arguments.
%
%      ATN_GUI('Property','Value',...) creates a new ATN_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ATN_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ATN_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ATN_gui_OpeningFcn, ...
    'gui_OutputFcn',  @ATN_gui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function ATN_gui_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;
% Save UIControl default values to handles-struct

% handles = addBackgroundImages(handles);
handles = addTextObjects(handles);
handles = groupHandles(handles);
handles = updateEcosystemPopupmenu(handles,[]);
handles.tab_simulate.pushbutton_stop.UserData.simulationStopped = false;

set(handles.tab_guilds.all,'Visible','off')
set(handles.tab_simulate.all,'Visible','on')

handles = saveUIControlValues(handles);

ecosyses = handles.tab_simulate.popupmenu_ecosystem.String;
if ~isempty(ecosyses)
    ecosys = feval(@(h)h.String{h.Value},handles.tab_simulate.popupmenu_ecosystem);
    load(['Data/' ecosys '.mat'])
    handles.Data = Data;
    handles.Data = updateGuildInfo(handles.Data);
    handles.Data = addDerivedQuantities(handles.Data);
    handles.minaxeswidth = 5;   % TODO: gather all gui parameters in one place
    % handles.plots = struct([]);   % Initialize plots-struct
    handles = resetSimulation(handles);
    handles = initializeGuildTab(handles);
    % handles.Game = initializeGameStruct();
    % handles.Game.Settings = setGameSettings();
    handles.isGuildChange = false;
    handles.isSimulationStarted = false;
    handles.isNotSaved = false;   
    
    % Update handles structure
    guidata(hObject, handles)
else
    herror = errordlg('No mat-files found under Data/.');
    uiwait(herror)
    handles.output = [];
    guidata(hObject, handles)
    delete(handles.atn_gui)
end


function varargout = ATN_gui_OutputFcn(~, ~, ~)
% varargout{1} = handles.output;
varargout{1} = [];

% ----------------------------------------------------------------------- %
% ------------------------------ Callbacks ------------------------------ %
% ----------------------------------------------------------------------- %

% ------ MENU
function menu_view_Callback(hObject, ~, handles) %#ok<DEFNU>
if strcmpi(handles.tab_simulate.uipanel_tab_simulate.Visible,'on')
    handles.tab_simulate.visibleStatus = get(handles.tab_simulate.all,'visible');
end
guidata(hObject, handles)
function menu_simulate_tab_Callback(hObject, ~, handles) %#ok<DEFNU>
set(handles.tab_guilds.all,'Visible','off')
if handles.isGuildChange
    handles = resetSimulation(handles);
    for i = 1:length(handles.tab_simulate.all)
        if strcmpi(handles.tab_simulate.all(i).Type,'line')
            if any(handles.tab_simulate.all(i) == handles.plots.BiomassFishGuild)
                set(handles.tab_simulate.all(i),'visible','on')
                handles.tab_simulate.visibleStatus{i} = 'on';
            elseif any(handles.tab_simulate.all(i) == handles.plots.CatchFishGuildYearly)
                set(handles.tab_simulate.all(i),'visible','on')
                handles.tab_simulate.visibleStatus{i} = 'on';
            else
                set(handles.tab_simulate.all(i),'visible','off')
                handles.tab_simulate.visibleStatus{i} = 'off';
            end
        else
            set(handles.tab_simulate.all(i),'visible',handles.tab_simulate.visibleStatus{i})
        end
    end
else
    for i = 1:length(handles.tab_simulate.all)
        set(handles.tab_simulate.all(i),'visible',handles.tab_simulate.visibleStatus{i})
    end
end
str_sel = feval(@(h)h.String{h.Value},handles.tab_simulate.popupmenu_ecosystem);
handles = updateEcosystemPopupmenu(handles,str_sel);
handles = saveParameterValues(handles);
guidata(hObject, handles)
function menu_guilds_tab_Callback(~, ~, handles) %#ok<DEFNU>
set(handles.tab_guilds.all,'Visible','on')
set(handles.tab_simulate.all,'Visible','off')
h = handles.tab_guilds.listbox_guild_list;
h.Value = 1;
h.Callback(h, handles)

% ------ SIMULATE TAB UICONTROLS
function pushbutton_gofishing_Callback(hObject, ~, handles) %#ok<DEFNU>
enableUIControls(hObject, handles, 'off');

handles.Data.F = str2double(handles.uicontrolvalues.edit_fishingmortality);
handles.Data.nGrowthDays = str2double(handles.uicontrolvalues.edit_length_season);
handles.Data.hmax = handles.Data.F/(handles.Data.nGrowthDays+1);
gearStr = handles.tab_simulate.uibuttongroup_gear.SelectedObject.String;
S = selective_fishery(handles.Data, gearStr);
E = S;

if handles.tab_simulate.radiobutton_gear_3.Value == 1 % Fyke net selected
    propRelease = str2double(handles.tab_simulate.edit_releaseproportion.String);
    oldFishInds = find([handles.Data.Guilds(handles.Data.GuildInfo.iAdultFishGuilds).age] == 4);
    E(oldFishInds) = E(oldFishInds)*(1-propRelease);
end

handles.Data.E = E;
handles.Data.P_mat = proba_mature(0);

handles = simulateForward(handles);
enableUIControls(hObject, handles, 'on');
handles.isSimulationStarted = true;
guidata(hObject, handles)
function pushbutton_nofishing_Callback(hObject, ~, handles) %#ok<DEFNU>
enableUIControls(hObject, handles, 'off');
handles.Data.nGrowthDays = str2double(handles.uicontrolvalues.edit_length_season);
handles.Data.E = zeros(handles.Data.GuildInfo.nAdultFishGuilds,1);
handles.Data.P_mat = proba_mature(0);
handles = simulateForward(handles);
enableUIControls(hObject, handles, 'on');
handles.isSimulationStarted = true;
guidata(hObject, handles)
function pushbutton_reset_Callback(hObject, ~, handles) %#ok<DEFNU>
handles = resetSimulation(handles);
enableUIControls(hObject, handles, 'on');
handles.isSimulationStarted = false;
guidata(hObject, handles)
function pushbutton_save_results_Callback(~, ~, handles) %#ok<DEFNU>
[filename, dirname] = uiputfile('*.mat');
Results = handles.Results; %#ok<NASGU>
Game = handles.Game; %#ok<NASGU>
save([dirname filename],'Results','Game')
function pushbutton_rewind_Callback(hObject, ~, handles) %#ok<DEFNU>
if handles.Game.currentyear > 0
    handles.Game.currentyear = handles.Game.currentyear - 1;
    handles.Game.gameover = 0;
    
    handles.Results.B(:,end) = [];
    handles.Results.C(:,end) = [];
    handles.Results.GF(:,end) = [];
    handles.Results.G(:,:,end) = [];
    handles.Results.L(:,:,end) = [];
    
    if size(handles.Results.B,2) > 0
        handles.Simulation.Binit = handles.Results.B(:,end);
    else
        handles.Simulation.Binit = vertcat(handles.Data.Guilds.binit);
    end
    handles = refreshAxes(handles); % Refresh plots at each time step
    
    handles.Results.ARE(:,end) = [];
    handles.Results.allbiomasses(:,end-handles.Data.nGrowthDays:end) = [];
    
    isGameMode = handles.tab_simulate.checkbox_gamemode.Value;
    handles.Game.history(end) = [];
    if isGameMode
        
        handles.Game.balance = handles.Game.history(end).balance;
        
        handles.tab_simulate.text_balance.String = ...
            ['Balance: ' num2str(round(handles.Game.balance))];
    end
    enableUIControls(hObject, handles, 'on');
    guidata(hObject, handles)
end
function pushbutton_stop_Callback(hObject, ~, handles) %#ok<DEFNU>
hObject.UserData.simulationStopped = true;
guidata(hObject, handles)
function pushbutton_planktonfigure_Callback(~, ~, handles) %#ok<DEFNU>

Guilds = handles.Data.Guilds;
gls = {Guilds.label};

iAlg = find(1-cellfun(@isempty,cellfun(@(x)strfind(x,'Alg'),gls,'uniformoutput',false)));
iBac = find(1-cellfun(@isempty,cellfun(@(x)strfind(x,'Bac'),gls,'uniformoutput',false)));
iCil = find(1-cellfun(@isempty,cellfun(@(x)strfind(x,'Cil'),gls,'uniformoutput',false)));
iRot = find(1-(cellfun(@isempty,cellfun(@(x)strfind(x,'Rot'),gls,'uniformoutput',false)) & ...
    cellfun(@isempty,cellfun(@(x)strfind(x,'Asp'),gls,'uniformoutput',false))));
iHCr = find(1-cellfun(@isempty,cellfun(@(x)strfind(x,'Cru'),gls,'uniformoutput',false)));
iCcr = find(1-(cellfun(@isempty,cellfun(@(x)strfind(x,'Cyc'),gls,'uniformoutput',false)) & ...
    cellfun(@isempty,cellfun(@(x)strfind(x,'Lep'),gls,'uniformoutput',false))));

panelwidth = 0.28*ones(1,6);
panelheight = 0.4*ones(1,6);
panelvertoffset = 0.07;
panelhorzoffset = 0.05;
panelbottom0 = 0.08;
panelleft0 = 0.04;
panelbottom = ones(1,6)*panelbottom0 + [1 1 1 0 0 0].*(panelheight+panelvertoffset);
panelleft = ones(1,6)*panelleft0 + [0 1 2 0 1 2].*(panelwidth+panelhorzoffset);

f = figure;
set(f,'paperpositionmode','auto')
set(f,'units','normalized','outerposition',[0.1 0.1 0.8 0.8])

subplot(2,3,1)
B = [[Guilds(iAlg).binit] ; handles.Results.B(iAlg,:)'];
plot(0:size(handles.Results.B,2),B)
set(gca,'tickdir','out')
box off
title('Algae')
set(gca,'position',[panelleft(1) panelbottom(1) panelwidth(1) panelheight(1)])

subplot(2,3,2)
B = [[Guilds(iBac).binit] ; handles.Results.B(iBac,:)'];
plot(0:size(handles.Results.B,2),B)
set(gca,'tickdir','out')
box off
title('Bacteria')
set(gca,'position',[panelleft(2) panelbottom(2) panelwidth(2) panelheight(2)])

subplot(2,3,3)
B = [[Guilds(iCil).binit] ; handles.Results.B(iCil,:)'];
plot(0:size(handles.Results.B,2),B)
set(gca,'tickdir','out')
box off
title('Ciliates')
set(gca,'position',[panelleft(3) panelbottom(3) panelwidth(3) panelheight(3)])

subplot(2,3,4)
B = [[Guilds(iRot).binit] ; handles.Results.B(iRot,:)'];
plot(0:size(handles.Results.B,2),B)
set(gca,'tickdir','out')
box off
title('Rotifiers')
set(gca,'position',[panelleft(4) panelbottom(4) panelwidth(4) panelheight(4)])

subplot(2,3,5)
B = [[Guilds(iHCr).binit] ; handles.Results.B(iHCr,:)'];
plot(0:size(handles.Results.B,2),B)
set(gca,'tickdir','out')
box off
title('Herbivorous crustacean')
set(gca,'position',[panelleft(5) panelbottom(5) panelwidth(5) panelheight(5)])

subplot(2,3,6)
B = [[Guilds(iCcr).binit] ; handles.Results.B(iCcr,:)'];
plot(0:size(handles.Results.B,2),B)
set(gca,'tickdir','out')
box off
title('Carnivorous crustacean')
set(gca,'position',[panelleft(6) panelbottom(6) panelwidth(6) panelheight(6)])
function pushbutton_totalbiomassfigure_Callback(~, ~, handles) %#ok<DEFNU>

GuildInfo = handles.Data.GuildInfo;
Guilds = handles.Data.Guilds;
Iall = 1:GuildInfo.nGuilds;
I = setdiff(Iall,GuildInfo.iDetritusGuilds);

f = figure;
set(f,'paperpositionmode','auto')
set(f,'units','normalized','outerposition',[0.1 0.1 0.8 0.8])

subplot(2,2,1)
plot(1:size(handles.Results.B,2),sum(handles.Results.B(I,:)))
set(gca,'xtick',unique(round(get(gca,'xtick'))))
title('Ecosystem biomass')

subplot(2,2,2)
ecosysindicator = (sum(handles.Results.B(I,:)) / sum([Guilds(I).binit])).^handles.Game.Settings.ECOEXP;
plot(1:size(handles.Results.B,2),ecosysindicator)
hold on
plot(get(gca,'xlim'),repmat(handles.Game.Settings.ECOSYSINDLIMIT,1,2),'--r')
set(gca,'ylim',[0 1])
set(gca,'xtick',unique(round(get(gca,'xtick'))))
title('Ecosystem status indicator')

subplot(2,2,3)
legends = {};
for i = 1:GuildInfo.nFishSpecies
    fishgls = {Guilds(GuildInfo.iFishGuilds).label};
    S.subs = {1,1:3};
    S.type = '()';
    fslabels = unique(cellfun(@(x)subsref(x,S),fishgls,'uniformoutput',false),'stable');
    afishgls = {Guilds(GuildInfo.iAdultFishGuilds).label};
    afishgls = cellfun(@(x)subsref(x,S),afishgls,'uniformoutput',false);
    afsinds = find(strcmpi(fslabels{i},afishgls));
    xdata = 0:size(handles.Results.ARE,2);
    ydata = [0 sum(handles.Results.ARE(afsinds,:))]; %#ok<FNDSB>
    plot(xdata,ydata)
    hold on
    legends{end+1} = fslabels{i}; %#ok<AGROW>
end
hold off
title('Fish egg production')
legend(legends)

function edit_numberofyears_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'int',0,999)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_fishingmortality_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,91)
    set(handles.tab_simulate.slider_F,'Value',str2double(hObject.String))
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_axeswidth_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'int',0,999)
    handles = saveUIControlValues(handles);
end
handles = refreshAxes(handles);
guidata(hObject, handles)
function edit_length_season_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'int',0,365)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_releaseproportion_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,1)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)

function checkbox_biomassplot_Callback(hObject, ~, handles) %#ok<DEFNU>
gls = {handles.Data.Guilds(handles.Data.GuildInfo.iFishGuilds).label};
I = find(strcmpi(gls,hObject.String));
J = find(I==handles.Data.GuildInfo.iAdultFishGuildsInFish);
%J = I - (handles.Data.GuildInfo.nLarvaeFishGuilds+handles.Data.GuildInfo.nJuvenileFishGuilds);
if J < 1
    J = [];
end
switch hObject.Value
    case 0
        set(handles.plots.BiomassFishGuild(I),'Visible','off')
        switch handles.tab_simulate.popupmenu_catchtype.Value
            case 1
                set(handles.plots.CatchFishGuildYearly(J),'Visible','off')
                set(handles.plots.CatchFishGuildCumulative,'Visible','off')
                set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','off')
                set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','off')
            case 2
                set(handles.plots.CatchFishGuildYearly,'Visible','off')
                set(handles.plots.CatchFishGuildCumulative(J),'Visible','off')
                set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','off')
                set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','off')
            case 3
                set(handles.plots.CatchFishGuildYearly,'Visible','off')
                set(handles.plots.CatchFishGuildCumulative,'Visible','off')
                set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','on')
                set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','off')
            case 4
                set(handles.plots.CatchFishGuildYearly,'Visible','off')
                set(handles.plots.CatchFishGuildCumulative,'Visible','off')
                set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','off')
                set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','on')
        end
    case 1
        set(handles.plots.BiomassFishGuild(I),'Visible','on')
        switch handles.tab_simulate.popupmenu_catchtype.Value
            case 1
                set(handles.plots.CatchFishGuildYearly(J),'Visible','on')
                set(handles.plots.CatchFishGuildCumulative,'Visible','off')
                set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','off')
                set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','off')
            case 2
                set(handles.plots.CatchFishGuildYearly,'Visible','off')
                set(handles.plots.CatchFishGuildCumulative(J),'Visible','on')
                set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','off')
                set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','off')
            case 3
                set(handles.plots.CatchFishGuildYearly,'Visible','off')
                set(handles.plots.CatchFishGuildCumulative,'Visible','off')
                set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','on')
                set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','off')
            case 4
                set(handles.plots.CatchFishGuildYearly,'Visible','off')
                set(handles.plots.CatchFishGuildCumulative,'Visible','off')
                set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','off')
                set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','on')
        end
end
guidata(hObject, handles)
function checkbox_gamemode_Callback(~, ~, ~) %#ok<DEFNU>
function checkbox_catchable_Callback(~, ~, ~) %#ok<DEFNU>

function popupmenu_catchtype_Callback(hObject, ~, handles) %#ok<DEFNU>

visPlotInds = find([handles.tab_simulate.uipanel_biomassplots.Children.Value]);
afishgls = {handles.Data.Guilds(handles.Data.GuildInfo.iAdultFishGuilds).label};
visInds = [];
for i = 1:length(visPlotInds)
    child = handles.tab_simulate.uipanel_biomassplots.Children(visPlotInds(i));
    visInds = [visInds find(strcmpi(child.String,afishgls))]; %#ok<AGROW>
end

switch hObject.Value
    case 1
        set(handles.plots.CatchFishGuildYearly(visInds),'Visible','on')
        set(handles.plots.CatchFishGuildCumulative,'Visible','off')
        set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','off')
        set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','off')
    case 2
        set(handles.plots.CatchFishGuildYearly,'Visible','off')
        set(handles.plots.CatchFishGuildCumulative(visInds),'Visible','on')
        set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','off')
        set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','off')
    case 3
        set(handles.plots.CatchFishGuildYearly,'Visible','off')
        set(handles.plots.CatchFishGuildCumulative,'Visible','off')
        set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','on')
        set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','off')
    case 4
        set(handles.plots.CatchFishGuildYearly,'Visible','off')
        set(handles.plots.CatchFishGuildCumulative,'Visible','off')
        set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','off')
        set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','on')
end
function popupmenu_ecosystem_Callback(hObject, ~, handles) %#ok<DEFNU>

if handles.isSimulationStarted
    answer = questdlg('This will reset simulation, and all progress will be lost. Do you want to continue?', ...
        'Confirm simulation reset','Yes','No','No');
end
if ~handles.isSimulationStarted || strcmpi(answer,'yes')
    load(['Data\' hObject.String{hObject.Value} '.mat'])
    handles.Data = Data;
    handles.Data = updateGuildInfo(handles.Data);
    handles.Data = addDerivedQuantities(handles.Data);
    handles.minaxeswidth = 5;   % TODO: gather all gui parameters in one place
    handles = resetSimulation(handles);
    handles = initializeGuildTab(handles);
    handles = saveUIControlValues(handles);
    handles.isGuildChange = false;
    handles.isSimulationStarted = false;
    guidata(hObject, handles)
else
    hObject.Value = handles.uicontrolvalues.(hObject.Tag);
end
function popupmenu_biomasstype_Callback(hObject, ~, handles) %#ok<DEFNU>

visPlotInds = find([handles.tab_simulate.uipanel_biomassplots.Children.Value]);
fishgls = {handles.Data.Guilds(handles.Data.GuildInfo.iFishGuilds).label};
visInds = [];
for i = 1:length(visPlotInds)
    child = handles.tab_simulate.uipanel_biomassplots.Children(visPlotInds(i));
    visInds = [visInds find(strcmpi(child.String,fishgls))]; %#ok<AGROW>
end

switch hObject.Value
    case 1
        set(handles.plots.BiomassFishGuild(visInds),'Visible','on')
        set(handles.plots.BiomassFishGuildPerSpecies,'Visible','off')
    case 2
        set(handles.plots.BiomassFishGuild,'Visible','off')
        set(handles.plots.BiomassFishGuildPerSpecies,'Visible','on')
end
function uibuttongroup_gear_SelectionChangedFcn(hObject, ~, handles) %#ok<DEFNU>
if strcmpi(hObject.String,'Fyke net')
    handles.tab_simulate.edit_releaseproportion.Enable = 'on';
else
    handles.tab_simulate.edit_releaseproportion.Enable = 'off';
end

function slider_F_Callback(hObject, ~, handles) %#ok<DEFNU>
set(handles.tab_simulate.edit_fishingmortality,'String',num2str(hObject.Value))
handles = saveUIControlValues(handles);
guidata(hObject, handles)


% ------ GUILDS TAB UICONTROLS
function pushbutton_import_Callback(hObject, ~, handles) %#ok<DEFNU>
answer1 = 'no';
answer2 = 'yes';
if handles.isNotSaved
    answer1 = questdlg('You have unsaved changes in the current ecosystem definitions. Do you want to save changes first?', ...
        'Confirm save changes','Yes','No','No');
end
if strcmpi(answer1,'yes')
    pushbutton_export_Callback(handles.tab_guilds.pushbutton_export, [], handles)
    answer2 = questdlg('Do you want to proceed to importing ecosystem definitions?', ...
        'Confirm save changes','Yes','No','No');
end
if strcmpi(answer2,'yes')
    if handles.isSimulationStarted
        answer3 = questdlg('This will reset simulation, and all progress will be lost. Do you want to continue?', ...
            'Confirm simulation reset','Yes','No','No');
    end
    if ~handles.isSimulationStarted || strcmpi(answer3,'yes')
        [filename, dirname] = uigetfile('*.mat');
        if ischar(filename) && ischar(dirname)
            load([dirname filename])
            if ~exist('Data','var')
                errordlg('Data file corrupted.')
            else
                handles.Data = Data;
                handles.Data = updateGuildInfo(handles.Data);
                handles.Data = addDerivedQuantities(handles.Data);
                handles = resetSimulation(handles); %%%
                handles = initializeGuildTab(handles);
                
                % Set ecosystem popupmenu selection to imported datafile
                str_sel = filename(1:end-4);
                handles = updateEcosystemPopupmenu(handles,str_sel);
                
                handles.isGuildChange = true;
                handles.isSimulationStarted = false;
                guidata(hObject, handles)
            end
        end
    end
end
function pushbutton_export_Callback(hObject, ~, handles)
handles = saveParameterValues(handles);
[filename, dirname] = uiputfile('*.mat');
if ischar(filename) && ischar(dirname)
    Data = handles.Data;
    Data = rmfield(Data,{'omega','GuildInfo'}); %#ok<NASGU>
    save([dirname filename],'Data')
    
    % Set ecosystem popupmenu selection to exported datafile
    str_sel = filename(1:end-4);
    handles = updateEcosystemPopupmenu(handles,str_sel);
end
handles.isNotSaved = false;
guidata(hObject, handles)
function pushbutton_add_guild_Callback(hObject, ~, handles) %#ok<DEFNU>
if handles.isSimulationStarted
    answer = questdlg('This will reset simulation, and all progress will be lost. Do you want to continue?', ...
        'Confirm simulation reset','Yes','No','No');
end
if ~handles.isSimulationStarted || strcmpi(answer,'yes')
    glabel = handles.tab_guilds.edit_label.String;
    if isempty(glabel)
        herror = errordlg('Guild label cannot be empty.');
        uiwait(herror)
        uicontrol(handles.tab_guilds.edit_label)
        return
    end
    if ismember(glabel,{handles.Data.Guilds.label})
        errordlg('Guild label exists already. Choose another label.')
    else
        gtype = handles.tab_guilds.uibuttongroup_type.SelectedObject.String;
        gname = handles.tab_guilds.edit_name.String;
        if isempty(gname)
            herror = errordlg('Guild name cannot be empty.');
            uiwait(herror)
            uicontrol(handles.tab_guilds.edit_name)
            return
        end
        gbinit = str2double(handles.tab_guilds.edit_initial_biomass.String);
        if isnan(gbinit)
            herror = errordlg('Initial biomass cannot be empty.');
            uiwait(herror)
            uicontrol(handles.tab_guilds.edit_initial_biomass)
            return
        end
        
        gigr = str2double(handles.tab_guilds.edit_intrinsic_growth_rate.String);
        gmbr = str2double(handles.tab_guilds.edit_metabolic_rate.String);
        gavgl = str2double(handles.tab_guilds.edit_average_length.String);
        glw_a = str2double(handles.tab_guilds.edit_lw_a.String);
        glw_b = str2double(handles.tab_guilds.edit_lw_b.String);
        gbx = str2double(handles.tab_guilds.edit_basal_respiration.String);
        gc = str2double(handles.tab_guilds.edit_producer_competition.String);
        gs = str2double(handles.tab_guilds.edit_fraction_of_exudation.String);
        gdiss_rate = str2double(handles.tab_guilds.edit_dissolution_rate.String);
        ghatchery = str2double(handles.tab_guilds.edit_hatchery.String);
        gax = str2double(handles.tab_guilds.edit_activity_respiration_coefficient.String);
        ginvest = str2double(handles.tab_guilds.edit_invest.String);
        gcatchable = handles.tab_guilds.checkbox_catchable.Value;
        
        switch lower(gtype)
            case 'fish'
                if length(glabel)~=4 || isnan(str2double(glabel(4)))
                    herror = errordlg(['Label must be four characters long where the ' ...
                        'last character is a number denoting the age of the guild.']);
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_basal_respiration)
                    return
                end
                gigr = [];
                gc = [];
                gage = str2double(glabel(4));
                gs = [];
                gdiss_rate = [];
                
                if isnan(gbx)
                    herror = errordlg('Basal respiration cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_basal_respiration)
                    return
                end
                if isnan(gavgl)
                    herror = errordlg('Average length of individuals cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_average_length)
                    return
                end
                if isnan(glw_a)
                    herror = errordlg('Length-weight parameter ''a'' cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_lw_a)
                    return
                end
                if isnan(glw_b)
                    herror = errordlg('Length-weight parameter ''b'' cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_lw_b)
                    return
                end
                gmbr = metabolicRate(gavgl,glw_a,glw_b);
                if isnan(ghatchery)
                    hwarn = warndlg('Hatchery was left empty. Using default value of 0.');
                    uiwait(hwarn)
                    handles.tab_guilds.edit_hatchery.String = num2str(0);
                    ghatchery = 0;
                end
                if isnan(gax)
                    herror = errordlg('Activity respiration coefficient cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_activity_respiration_coefficient)
                    return
                end
                if isnan(ginvest)
                    hwarn = warndlg(['Surplus invested to reproduction was left empty.' ...
                        'Using default value of 0.']);
                    uiwait(hwarn)
                    handles.tab_guilds.edit_invest.String = num2str(0);
                    ginvest = 0;
                end
            case 'consumer'
                gage = [];
                gigr = [];
                gc = [];
                gavgl = [];
                glw_a = [];
                glw_b = [];
                gs = [];
                gdiss_rate = [];
                ghatchery = [];
                ginvest = [];
                
                if isnan(gbx)
                    herror = errordlg('Basal respiration cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_basal_respiration)
                    return
                end
                if isnan(gmbr)
                    herror = errordlg('Metabolic rate cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_metabolic_rate)
                    return
                end
                if isnan(gax)
                    herror = errordlg('Activity respiration coefficient cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_activity_respiration_coefficient)
                    return
                end
            case 'producer'
                gage = [];
                gbx = [];
                gmbr = [];
                gavgl = [];
                glw_a = [];
                glw_b = [];
                gdiss_rate = [];
                ghatchery = [];
                gax = [];
                ginvest = [];
                
                if isnan(gigr)
                    herror = errordlg('Intrinsic growth rate cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_intrinsic_growth_rate)
                    return
                end
                if isnan(gc)
                    herror = errordlg('Producer competition cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_producer_competition)
                    return
                end
                if isnan(gs)
                    herror = errordlg('Fraction of exudation cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_fraction_of_exudation)
                    return
                end
            case 'detritus'
                gage = [];
                gigr = [];
                gmbr = [];
                gbx = [];
                gc = [];
                gavgl = [];
                glw_a = [];
                glw_b = [];
                gs = [];
                ghatchery = [];
                gax = [];
                ginvest = [];
                if isnan(gdiss_rate)
                    hwarn = errordlg(['Dissolution rate was left empty. Using default ' ...
                        'value of 0. If this is the pool of undissolved detritus, ' ...
                        'you must change the rate later.']);
                    uiwait(hwarn)
                    handles.tab_guilds.edit_dissolution_rate.String = num2str(0);
                    gdiss_rate = 0;
                end
        end
        
        newguild = struct( ...
            'label',glabel,'name',gname,'type',gtype,'igr',gigr,'mbr',gmbr,'avgl',gavgl, ...
            'lw_a',glw_a,'lw_b',glw_b,'binit',gbinit,'age',gage,'bx',gbx,'c',gc, ...
            's',gs,'diss_rate',gdiss_rate,'hatchery',ghatchery,'ax',gax,'invest',ginvest, ...
            'catchable',gcatchable);
        
        handles.Data.Guilds(end+1) = newguild;
        if isempty(handles.Data.communityMatrix)
            handles.Data.communityMatrix(1,1) = 0;
            handles.Data.B0(1,1) = 0;
            handles.Data.d(1,1) = 0;
            handles.Data.q(1,1) = 0;
            handles.Data.e(1,1) = 0;
            handles.Data.y(1,1) = 0;
        else
            handles.Data.communityMatrix(end+1,:) = 0;
            handles.Data.communityMatrix(:,end+1) = 0;
            handles.Data.B0(end+1,:) = 0;
            handles.Data.B0(:,end+1) = 0;
            handles.Data.d(end+1,:) = 0;
            handles.Data.d(:,end+1) = 0;
            handles.Data.q(end+1,:) = 0;
            handles.Data.q(:,end+1) = 0;
            handles.Data.e(end+1,:) = 0;
            handles.Data.e(:,end+1) = 0;
            handles.Data.y(end+1,:) = 0;
            handles.Data.y(:,end+1) = 0;
        end
        
        value = length(handles.Data.Guilds);
        handles = updateGuildListbox(handles, value);
        handles = updatePreyListbox(handles);
        handles = updateAvailablePreyListbox(handles);
        handles.Data = updateGuildInfo(handles.Data);
        handles.Data = addDerivedQuantities(handles.Data);
        handles = updateEdits(handles);
        handles.isGuildChange = true;
        handles.isNotSaved = true;
        handles.isSimulationStarted = false;
        
        guidata(hObject, handles)
        
    end
    
    
end
function pushbutton_clear_Callback(hObject, ~, handles)
set(findobj(handles.tab_guilds.uipanel_add_guild.Children,'Style','edit'),'String','')
handles.isNotSaved = true;
guidata(hObject, handles)
function pushbutton_delete_Callback(hObject, ~, handles) %#ok<DEFNU>
if handles.isSimulationStarted
    answer = questdlg('This will reset simulation, and all progress will be lost. Do you want to continue?', ...
        'Confirm simulation reset','Yes','No','No');
end
if ~handles.isSimulationStarted || strcmpi(answer,'yes')
    lvalue = handles.tab_guilds.listbox_guild_list.Value;
    lstrings = handles.tab_guilds.listbox_guild_list.String;
    if ~isempty(lstrings)
        slabel = lstrings{lvalue};
        lvalue = max([lvalue-1 1]);
        
        di = strcmp(slabel,{handles.Data.Guilds.label});
        
        handles.Data.Guilds(di) = [];
        handles.Data.communityMatrix(di,:) = [];
        handles.Data.communityMatrix(:,di) = [];
        handles.Data.B0(di,:) = [];
        handles.Data.B0(:,di) = [];
        handles.Data.d(di,:) = [];
        handles.Data.d(:,di) = [];
        handles.Data.q(di,:) = [];
        handles.Data.q(:,di) = [];
        handles.Data.e(di,:) = [];
        handles.Data.e(:,di) = [];
        handles.Data.y(di,:) = [];
        handles.Data.y(:,di) = [];
        
        
        handles.Data = updateGuildInfo(handles.Data);
        handles.Data = addDerivedQuantities(handles.Data);
        handles = updateGuildListbox(handles, lvalue);
        handles = updatePreyListbox(handles);
        handles = updateAvailablePreyListbox(handles);
        pushbutton_clear_Callback(handles.tab_guilds.pushbutton_clear, [], handles);
        handles.isGuildChange = true;
        handles.isNotSaved = true;
        handles.isSimulationStarted = false;
        
        guidata(hObject, handles)
    end
end
function pushbutton_save_changes_Callback(hObject, ~, handles) %#ok<DEFNU>
if handles.isSimulationStarted
    answer = questdlg('This will reset simulation, and all progress will be lost. Do you want to continue?', ...
        'Confirm simulation reset','Yes','No','No');
end
if ~handles.isSimulationStarted || strcmpi(answer,'yes')
    glabel = handles.tab_guilds.edit_label.String;
    if isempty(glabel)
        herror = errordlg('Guild label cannot be empty.');
        uiwait(herror)
        uicontrol(handles.tab_guilds.edit_label)
        return
    end
    i = handles.tab_guilds.listbox_guild_list.Value;
    noti = setdiff(1:handles.Data.GuildInfo.nGuilds,i);
    if ismember(glabel,{handles.Data.Guilds(noti).label})
        errordlg('Guild label exists already. Choose another label.')
    else
        gtype = handles.tab_guilds.uibuttongroup_type.SelectedObject.String;
        gname = handles.tab_guilds.edit_name.String;
        if isempty(gname)
            herror = errordlg('Guild name cannot be empty.');
            uiwait(herror)
            uicontrol(handles.tab_guilds.edit_name)
            return
        end
        gbinit = str2double(handles.tab_guilds.edit_initial_biomass.String);
        if isnan(gbinit)
            herror = errordlg('Initial biomass cannot be empty.');
            uiwait(herror)
            uicontrol(handles.tab_guilds.edit_initial_biomass)
            return
        end
        
        gigr = str2double(handles.tab_guilds.edit_intrinsic_growth_rate.String);
        gmbr = str2double(handles.tab_guilds.edit_metabolic_rate.String);
        gavgl = str2double(handles.tab_guilds.edit_average_length.String);
        glw_a = str2double(handles.tab_guilds.edit_lw_a.String);
        glw_b = str2double(handles.tab_guilds.edit_lw_b.String);
        gbx = str2double(handles.tab_guilds.edit_basal_respiration.String);
        gc = str2double(handles.tab_guilds.edit_producer_competition.String);
        gs = str2double(handles.tab_guilds.edit_fraction_of_exudation.String);
        gdiss_rate = str2double(handles.tab_guilds.edit_dissolution_rate.String);
        ghatchery = str2double(handles.tab_guilds.edit_hatchery.String);
        gax = str2double(handles.tab_guilds.edit_activity_respiration_coefficient.String);
        ginvest = str2double(handles.tab_guilds.edit_invest.String);
        gcatchable = handles.tab_guilds.checkbox_catchable.Value;
        
        switch lower(gtype)
            case 'fish'
                if length(glabel)~=4 || isnan(str2double(glabel(4)))
                    herror = errordlg(['Label must be four characters long where the ' ...
                        'last character is a number denoting the age of the guild.']);
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_basal_respiration)
                    return
                end
                gigr = [];
                gc = [];
                gage = str2double(glabel(4));
                gs = [];
                gdiss_rate = [];
                
                if isnan(gbx)
                    herror = errordlg('Basal respiration cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_basal_respiration)
                    return
                end
                if isnan(gavgl)
                    herror = errordlg('Average length of individuals cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_average_length)
                    return
                end
                if isnan(glw_a)
                    herror = errordlg('Length-weight parameter ''a'' cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_lw_a)
                    return
                end
                if isnan(glw_b)
                    herror = errordlg('Length-weight parameter ''b'' cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_lw_b)
                    return
                end
                gmbr = metabolicRate(gavgl,glw_a,glw_b);
                if isnan(ghatchery)
                    hwarn = warndlg('Hatchery was left empty. Using default value of 0.');
                    uiwait(hwarn)
                    handles.tab_guilds.edit_hatchery.String = num2str(0);
                    ghatchery = 0;
                end
                if isnan(gax)
                    herror = errordlg('Activity respiration coefficient cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_activity_respiration_coefficient)
                    return
                end
                if isnan(ginvest)
                    hwarn = warndlg(['Surplus invested to reproduction was left empty.' ...
                        'Using default value of 0.']);
                    uiwait(hwarn)
                    handles.tab_guilds.edit_invest.String = num2str(0);
                    ginvest = 0;
                end
            case 'consumer'
                gigr = [];
                gc = [];
                gage = [];
                gavgl = [];
                glw_a = [];
                glw_b = [];
                gs = [];
                gdiss_rate = [];
                ghatchery = [];
                ginvest = [];
                
                if isnan(gbx)
                    herror = errordlg('Basal respiration cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_basal_respiration)
                    return
                end
                if isnan(gmbr)
                    herror = errordlg('Metabolic rate cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_metabolic_rate)
                    return
                end
                if isnan(gax)
                    herror = errordlg('Activity respiration coefficient cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_activity_respiration_coefficient)
                    return
                end
            case 'producer'
                gage = [];
                gbx = [];
                gmbr = [];
                gavgl = [];
                glw_a = [];
                glw_b = [];
                gdiss_rate = [];
                ghatchery = [];
                gax = [];
                ginvest = [];
                
                if isnan(gigr)
                    herror = errordlg('Intrinsic growth rate cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_intrinsic_growth_rate)
                    return
                end
                if isnan(gc)
                    herror = errordlg('Producer competition cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_producer_competition)
                    return
                end
                if isnan(gs)
                    herror = errordlg('Fraction of exudation cannot be empty.');
                    uiwait(herror)
                    uicontrol(handles.tab_guilds.edit_fraction_of_exudation)
                    return
                end
            case 'detritus'
                gage = [];
                gigr = [];
                gmbr = [];
                gbx = [];
                gc = [];
                gavgl = [];
                glw_a = [];
                glw_b = [];
                gs = [];
                ghatchery = [];
                gax = [];
                ginvest = [];
                if isnan(gdiss_rate)
                    hwarn = errordlg(['Dissolution rate was left empty. Using default ' ...
                        'value of 0. If this is the pool of undissolved detritus, ' ...
                        'you must change the rate later.']);
                    uiwait(hwarn)
                    handles.tab_guilds.edit_dissolution_rate.String = num2str(0);
                    gdiss_rate = 0;
                end
        end
        
        handles.Data.Guilds(i).name = gname;
        handles.Data.Guilds(i).label = glabel;
        handles.Data.Guilds(i).age = gage;
        handles.Data.Guilds(i).type = gtype;
        handles.Data.Guilds(i).igr = gigr;
        handles.Data.Guilds(i).mbr = gmbr;
        handles.Data.Guilds(i).avgl = gavgl;
        handles.Data.Guilds(i).lw_a = glw_a;
        handles.Data.Guilds(i).lw_b = glw_b;
        handles.Data.Guilds(i).binit = gbinit;
        handles.Data.Guilds(i).bx = gbx;
        handles.Data.Guilds(i).c = gc;
        handles.Data.Guilds(i).s = gs;
        handles.Data.Guilds(i).diss_rate = gdiss_rate;
        handles.Data.Guilds(i).hatchery = ghatchery;
        handles.Data.Guilds(i).ax = gax;
        handles.Data.Guilds(i).invest = ginvest;
        handles.Data.Guilds(i).catchable = gcatchable;
        
        handles = updateGuildListbox(handles, handles.tab_guilds.listbox_guild_list.Value);
        handles.Data = updateGuildInfo(handles.Data);
        handles.Data = addDerivedQuantities(handles.Data);
        handles = updateEdits(handles);
        handles.isGuildChange = true;
        handles.isNotSaved = true;
        handles.isSimulationStarted = false;
        guidata(hObject, handles)
    end
end
function pushbutton_add_food_item_Callback(hObject, ~, handles) %#ok<DEFNU>
if handles.isSimulationStarted
    answer = questdlg('This will reset simulation, and all progress will be lost. Do you want to continue?', ...
        'Confirm simulation reset','Yes','No','No');
end
if ~handles.isSimulationStarted || strcmpi(answer,'yes')
    Guilds = handles.Data.Guilds;
    i = handles.tab_guilds.listbox_guild_list.Value;
    predatorGuild = Guilds(i);
    jj = handles.tab_guilds.listbox_available_prey.Value;
    preyLabel = handles.tab_guilds.listbox_available_prey.String{jj};
    preyGuild = Guilds(strcmpi(preyLabel,{Guilds.label}));
    handles.tab_guilds.listbox_available_prey.String = ...
        setdiff(handles.tab_guilds.listbox_available_prey.String,preyLabel);
    handles.tab_guilds.listbox_available_prey.Value = max([1 jj-1]);
    handles.tab_guilds.listbox_prey.String = ...
        union(handles.tab_guilds.listbox_prey.String,preyLabel);
    handles.tab_guilds.listbox_prey.Value = ...
        find(strcmp(handles.tab_guilds.listbox_prey.String,preyLabel));
    j = find(strcmp({handles.Data.Guilds.label},preyLabel));
    handles.Data.communityMatrix(i,j) = 1;
    
    if strcmpi(predatorGuild.type,'consumer')
        default_B0 = 1500;
        default_d = 0;
        default_q = 1.2;
        default_e = 0.66;
        default_y = 10;
    elseif strcmpi(predatorGuild.type,'fish')
        if predatorGuild.age == 0       % Larvae
            default_B0 = 15000;
            default_d = 0.0001;
            default_q = 1.2;
            default_e = 0.66;
            default_y = 10;
        else                            % Juveniles and adults
            if strcmpi(preyGuild.type,'fish')
                default_B0 = 15000;
                default_d = 0.01;
                default_q = 1.2;
                default_e = 0.66;
                default_y = 10;
            else
                default_B0 = 25000;
                default_d = 0.001;
                default_q = 1.2;
                default_e = 0.66;
                default_y = 10;
            end
        end
    else
        default_B0 = 0;
        default_d = 0;
        default_q = 0;
        default_e = 0;
        default_y = 0;
    end
    handles.Data.B0(i,j) = default_B0;
    handles.Data.d(i,j) = default_d;
    handles.Data.q(i,j) = default_q;
    handles.Data.e(i,j) = default_e;
    handles.Data.y(i,j) = default_y;
    
    handles.tab_guilds.edit_half_saturation_constant.String = num2str(default_B0);
    handles.tab_guilds.edit_feeding_interference_coefficient.String = num2str(default_d);
    handles.tab_guilds.edit_F_q.String = num2str(default_q);
    handles.tab_guilds.edit_assimilation_efficiency.String = num2str(default_e);
    handles.tab_guilds.edit_maximum_consumption_rate.String = num2str(default_y);
    
    handles.Data = updateGuildInfo(handles.Data);
    handles.Data = addDerivedQuantities(handles.Data);
    handles.isGuildChange = true;
    handles.isNotSaved = true;
    handles.isSimulationStarted = false;
    guidata(hObject, handles)
end
function pushbutton_remove_food_item_Callback(hObject, ~, handles) %#ok<DEFNU>
if handles.isSimulationStarted
    answer = questdlg('This will reset simulation, and all progress will be lost. Do you want to continue?', ...
        'Confirm simulation reset','Yes','No','No');
end
if ~handles.isSimulationStarted || strcmpi(answer,'yes')
    i = handles.tab_guilds.listbox_guild_list.Value;
    jj = handles.tab_guilds.listbox_prey.Value;
    removedpreylabel = handles.tab_guilds.listbox_prey.String{jj};
    handles.tab_guilds.listbox_prey.String = ...
        setdiff(handles.tab_guilds.listbox_prey.String,removedpreylabel);
    handles.tab_guilds.listbox_prey.Value = max([1 jj-1]);
    % update edits based on new selection
    listbox_prey_Callback(handles.tab_guilds.listbox_prey, [], handles);
    handles.tab_guilds.listbox_available_prey.String = ...
        union(handles.tab_guilds.listbox_available_prey.String,removedpreylabel);
    handles.tab_guilds.listbox_available_prey.Value = ...
        find(strcmp(handles.tab_guilds.listbox_available_prey.String,removedpreylabel));
    j = find(strcmp({handles.Data.Guilds.label},removedpreylabel));
    handles.Data.communityMatrix(i,j) = 0; %#ok<FNDSB>
    handles.Data = updateGuildInfo(handles.Data);
    handles.Data = addDerivedQuantities(handles.Data);
    handles.isGuildChange = true;
    handles.isNotSaved = true;
    handles.isSimulationStarted = false;
    guidata(hObject, handles)
end
function pushbutton_clearall_Callback(hObject, ~, handles) %#ok<DEFNU>
if handles.isSimulationStarted
    answer = questdlg('This will reset simulation, and all progress will be lost. Do you want to continue?', ...
        'Confirm simulation reset','Yes','No','No');
end
if ~handles.isSimulationStarted || strcmpi(answer,'yes')
    handles.Data.Guilds(1:end) = [];
    handles.Data.communityMatrix = sparse([]);
    handles.Data.B0 = sparse([]);
    handles.d = sparse([]);
    handles.q = sparse([]);
    handles.e = sparse([]);
    handles.y = sparse([]);
    if isfield(handles.Data,'omega')
        handles.Data = rmfield(handles.Data,'omega');
    end
    if isfield(handles.Data,'GuildInfo')
        handles.Data = rmfield(handles.Data,'GuildInfo');
    end
    %%%% TODO: clear all other related variables! %%%%
    handles = updateGuildListbox(handles, 1);
    handles = updatePreyListbox(handles);
    handles = updateAvailablePreyListbox(handles);
    handles = updateEdits(handles);
    handles.isGuildChange = true;
    handles.isNotSaved = true;
    handles.isSimulationStarted = false;
    guidata(hObject, handles)
end
function pushbutton_movedown_Callback(hObject, ~, handles) %#ok<DEFNU>
if handles.isSimulationStarted
    answer = questdlg('This will reset simulation, and all progress will be lost. Do you want to continue?', ...
        'Confirm simulation reset','Yes','No','No');
end
if ~handles.isSimulationStarted || strcmpi(answer,'yes')
    N = length(handles.tab_guilds.listbox_guild_list.String);
    lbtop = handles.tab_guilds.listbox_guild_list.ListboxTop;
    i1 = handles.tab_guilds.listbox_guild_list.Value;
    if i1 < N
        i2 = i1+1;
        handles = swapGuildPositions(handles,i1,i2);
        handles.tab_guilds.listbox_guild_list.ListboxTop = i2-(i1-lbtop);
        handles.isGuildChange = true;
        handles.isNotSaved = true;
        handles.isSimulationStarted = false;
    end
    guidata(hObject, handles)
end
function pushbutton_moveup_Callback(hObject, ~, handles) %#ok<DEFNU>
if handles.isSimulationStarted
    answer = questdlg('This will reset simulation, and all progress will be lost. Do you want to continue?', ...
        'Confirm simulation reset','Yes','No','No');
end
if ~handles.isSimulationStarted || strcmpi(answer,'yes')
    lbtop = handles.tab_guilds.listbox_guild_list.ListboxTop;
    i1 = handles.tab_guilds.listbox_guild_list.Value;
    if i1 > 1
        i2 = i1-1;
        handles = swapGuildPositions(handles,i1,i2);
        handles.tab_guilds.listbox_guild_list.ListboxTop = max(i2-(i1-lbtop),1);
        handles.isGuildChange = true;
        handles.isNotSaved = true;
        handles.isSimulationStarted = false;
    end
    guidata(hObject, handles)
end
function edit_name_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'string',[],[])
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_label_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'labelstring',[],[])
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_intrinsic_growth_rate_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,10)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_metabolic_rate_Callback(hObject, ~, handles) %#ok<DEFNU> % POISTA
if isValidInput(hObject,handles,'real',0,1)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_average_length_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,999)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_lw_a_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,1)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_lw_b_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,10)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_initial_biomass_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,Inf)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_half_saturation_constant_Callback(hObject, ~, handles) %#ok<DEFNU>
if handles.isSimulationStarted
    answer = questdlg('This will reset simulation, and all progress will be lost. Do you want to continue?', ...
        'Confirm simulation reset','Yes','No','No');
end
if ~handles.isSimulationStarted || strcmpi(answer,'yes')
    preystrcell = handles.tab_guilds.listbox_prey.String;
    if ~isempty(preystrcell)
        if isValidInput(hObject,handles,'real',0,Inf)
            handles = saveUIControlValues(handles);
            
            selected_prey_label =  preystrcell{handles.tab_guilds.listbox_prey.Value};
            selected_predator_label = handles.tab_guilds.listbox_guild_list.String{ ...
                handles.tab_guilds.listbox_guild_list.Value};
            allguildlabels = {handles.Data.Guilds.label};
            ipredator = strcmpi(allguildlabels,selected_predator_label);
            iprey = strcmpi(allguildlabels,selected_prey_label);
            handles.Data.B0(ipredator,iprey) = ...
                str2double(handles.tab_guilds.edit_half_saturation_constant.String);
            handles.isGuildChange = true;
            handles.isNotSaved = true;
            handles.isSimulationStarted = false;
        end
    else
        warndlg('No prey selected for this guild. Discarding parameter change.', ...
            'No prey selected')
        hObject.String = '';
    end
    guidata(hObject, handles)
end
function edit_feeding_interference_coefficient_Callback(hObject, ~, handles) %#ok<DEFNU>
if handles.isSimulationStarted
    answer = questdlg('This will reset simulation, and all progress will be lost. Do you want to continue?', ...
        'Confirm simulation reset','Yes','No','No');
end
if ~handles.isSimulationStarted || strcmpi(answer,'yes')
    preystrcell = handles.tab_guilds.listbox_prey.String;
    if ~isempty(preystrcell)
        if isValidInput(hObject,handles,'real',0,1)
            handles = saveUIControlValues(handles);
            
            selected_prey_label = preystrcell{handles.tab_guilds.listbox_prey.Value};
            selected_predator_label = handles.tab_guilds.listbox_guild_list.String{ ...
                handles.tab_guilds.listbox_guild_list.Value};
            allguildlabels = {handles.Data.Guilds.label};
            ipredator = strcmpi(allguildlabels,selected_predator_label);
            iprey = strcmpi(allguildlabels,selected_prey_label);
            handles.Data.d(ipredator,iprey) = ...
                str2double(handles.tab_guilds.edit_feeding_interference_coefficient.String);
            handles.isGuildChange = true;
            handles.isNotSaved = true;
            handles.isSimulationStarted = false;
        end
    else
        warndlg('No prey selected for this guild. Discarding parameter change.', ...
            'No prey selected')
        hObject.String = '';
    end
    guidata(hObject, handles)
end
function edit_basal_respiration_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,1)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_producer_competition_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,2)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_K_mean_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,Inf)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_K_standard_deviation_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,Inf)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_K_autocorrelation_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,1)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_F_q_Callback(hObject, ~, handles) %#ok<DEFNU>
if handles.isSimulationStarted
    answer = questdlg('This will reset simulation, and all progress will be lost. Do you want to continue?', ...
        'Confirm simulation reset','Yes','No','No');
end
if ~handles.isSimulationStarted || strcmpi(answer,'yes')
    preystrcell = handles.tab_guilds.listbox_prey.String;
    if ~isempty(preystrcell)
        if isValidInput(hObject,handles,'real',0,Inf)
            handles = saveUIControlValues(handles);
            
            selected_prey_label =  preystrcell{handles.tab_guilds.listbox_prey.Value};
            selected_predator_label = handles.tab_guilds.listbox_guild_list.String{ ...
                handles.tab_guilds.listbox_guild_list.Value};
            allguildlabels = {handles.Data.Guilds.label};
            ipredator = strcmpi(allguildlabels,selected_predator_label);
            iprey = strcmpi(allguildlabels,selected_prey_label);
            handles.Data.q(ipredator,iprey) = ...
                str2double(handles.tab_guilds.edit_F_q.String);
            handles.isGuildChange = true;
            handles.isNotSaved = true;
            handles.isSimulationStarted = false;
        end
    else
        warndlg('No prey selected for this guild. Discarding parameter change.', ...
            'No prey selected')
        hObject.String = '';
    end
    guidata(hObject, handles)
end
function edit_assimilation_efficiency_Callback(hObject, ~, handles) %#ok<DEFNU>
if handles.isSimulationStarted
    answer = questdlg('This will reset simulation, and all progress will be lost. Do you want to continue?', ...
        'Confirm simulation reset','Yes','No','No');
end
if ~handles.isSimulationStarted || strcmpi(answer,'yes')
    preystrcell = handles.tab_guilds.listbox_prey.String;
    if ~isempty(preystrcell)
        if isValidInput(hObject,handles,'real',0,Inf)
            handles = saveUIControlValues(handles);
            
            selected_prey_label =  preystrcell{handles.tab_guilds.listbox_prey.Value};
            selected_predator_label = handles.tab_guilds.listbox_guild_list.String{ ...
                handles.tab_guilds.listbox_guild_list.Value};
            allguildlabels = {handles.Data.Guilds.label};
            ipredator = strcmpi(allguildlabels,selected_predator_label);
            iprey = strcmpi(allguildlabels,selected_prey_label);
            handles.Data.e(ipredator,iprey) = ...
                str2double(handles.tab_guilds.edit_assimilation_efficiency.String);
            handles.isGuildChange = true;
            handles.isNotSaved = true;
            handles.isSimulationStarted = false;
        end
    else
        warndlg('No prey selected for this guild. Discarding parameter change.', ...
            'No prey selected')
        hObject.String = '';
    end
    guidata(hObject, handles)
end
function edit_maximum_consumption_rate_Callback(hObject, ~, handles) %#ok<DEFNU>
if handles.isSimulationStarted
    answer = questdlg('This will reset simulation, and all progress will be lost. Do you want to continue?', ...
        'Confirm simulation reset','Yes','No','No');
end
if ~handles.isSimulationStarted || strcmpi(answer,'yes')
    preystrcell = handles.tab_guilds.listbox_prey.String;
    if ~isempty(preystrcell)
        if isValidInput(hObject,handles,'real',0,Inf)
            handles = saveUIControlValues(handles);
            
            selected_prey_label =  preystrcell{handles.tab_guilds.listbox_prey.Value};
            selected_predator_label = handles.tab_guilds.listbox_guild_list.String{ ...
                handles.tab_guilds.listbox_guild_list.Value};
            allguildlabels = {handles.Data.Guilds.label};
            ipredator = strcmpi(allguildlabels,selected_predator_label);
            iprey = strcmpi(allguildlabels,selected_prey_label);
            handles.Data.y(ipredator,iprey) = ...
                str2double(handles.tab_guilds.edit_maximum_consumption_rate.String);
            handles.isGuildChange = true;
            handles.isNotSaved = true;
            handles.isSimulationStarted = false;
        end
    else
        warndlg('No prey selected for this guild. Discarding parameter change.', ...
            'No prey selected')
        hObject.String = '';
    end
    guidata(hObject, handles)
end
function edit_invest_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,1)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_hatchery_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,Inf)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_dissolution_rate_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,1)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_fraction_of_exudation_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,1)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)
function edit_activity_respiration_coefficient_Callback(hObject, ~, handles) %#ok<DEFNU>
if isValidInput(hObject,handles,'real',0,1)
    handles = saveUIControlValues(handles);
end
guidata(hObject, handles)

function listbox_guild_list_Callback(hObject, ~, handles) %#ok<DEFNU>
handles = updatePreyListbox(handles);
handles = updateAvailablePreyListbox(handles);
handles = updateEdits(handles);
handles = setGuildEditStatus(handles);
%
preylabels = handles.tab_guilds.listbox_prey.String;
if ~isempty(preylabels)
    selected_prey_label = preylabels{ ...
        handles.tab_guilds.listbox_prey.Value};
    selected_predator_label = hObject.String{hObject.Value};
    allguildlabels = {handles.Data.Guilds.label};
    ipredator = strcmpi(allguildlabels,selected_predator_label);
    iprey = strcmpi(allguildlabels,selected_prey_label);
    dij = handles.Data.d(ipredator,iprey);
    B0ij = handles.Data.B0(ipredator,iprey);
    qij = handles.Data.q(ipredator,iprey);
    eij = handles.Data.e(ipredator,iprey);
    yij = handles.Data.y(ipredator,iprey);
else
    dij = [];
    B0ij = [];
    qij = [];
    eij = [];
    yij = [];
end
handles.tab_guilds.edit_half_saturation_constant.String = num2str(B0ij);
handles.tab_guilds.edit_feeding_interference_coefficient.String = num2str(dij);
handles.tab_guilds.edit_F_q.String = num2str(qij);
handles.tab_guilds.edit_assimilation_efficiency.String = num2str(eij);
handles.tab_guilds.edit_maximum_consumption_rate.String = num2str(yij);

%
guidata(hObject, handles)
function listbox_prey_Callback(hObject, ~, handles)
%
selected_prey_label = hObject.String{hObject.Value};
selected_predator_label = handles.tab_guilds.listbox_guild_list.String{ ...
    handles.tab_guilds.listbox_guild_list.Value};
allguildlabels = {handles.Data.Guilds.label};
ipredator = strcmpi(allguildlabels,selected_predator_label);
iprey = strcmpi(allguildlabels,selected_prey_label);
dij = handles.Data.d(ipredator,iprey);
B0ij = handles.Data.B0(ipredator,iprey);
qij = handles.Data.q(ipredator,iprey);
eij = handles.Data.e(ipredator,iprey);
yij = handles.Data.y(ipredator,iprey);
handles.tab_guilds.edit_half_saturation_constant.String = num2str(B0ij);
handles.tab_guilds.edit_feeding_interference_coefficient.String = num2str(dij);
handles.tab_guilds.edit_F_q.String = num2str(qij);
handles.tab_guilds.edit_assimilation_efficiency.String = num2str(eij);
handles.tab_guilds.edit_maximum_consumption_rate.String = num2str(yij);
%
guidata(hObject, handles)
function listbox_available_prey_Callback(~, ~, ~) %#ok<DEFNU>
function uibuttongroup_type_SelectionChangedFcn(hObject, ~, handles) %#ok<DEFNU>
handles = setGuildEditStatus(handles);
guidata(hObject, handles)
function uibuttongroup_carrying_capacity_SelectionChangedFcn(hObject, ~, handles)
switch hObject.Tag
    case 'radiobutton_K_constant'
        handles.tab_guilds.edit_K_mean.Enable = 'on';
        handles.tab_guilds.edit_K_standard_deviation.Enable = 'off';
        handles.tab_guilds.edit_K_autocorrelation.Enable = 'off';
    case 'radiobutton_K_white_noise'
        handles.tab_guilds.edit_K_mean.Enable = 'on';
        handles.tab_guilds.edit_K_standard_deviation.Enable = 'on';
        handles.tab_guilds.edit_K_autocorrelation.Enable = 'off';
    case 'radiobutton_K_AR1'
        handles.tab_guilds.edit_K_mean.Enable = 'on';
        handles.tab_guilds.edit_K_standard_deviation.Enable = 'on';
        handles.tab_guilds.edit_K_autocorrelation.Enable = 'on';
end
guidata(hObject, handles)


% ----------------------------------------------------------------------- %
% ------------------------------ CreateFcns ----------------------------- %
% ----------------------------------------------------------------------- %

function atn_gui_CreateFcn(~, ~, ~) %#ok<DEFNU>
function edit_fishingmortality_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_numberofyears_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_catchtype_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_axeswidth_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_name_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_label_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_metabolic_rate_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_intrinsic_growth_rate_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function listbox_guild_list_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function listbox_prey_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function listbox_available_prey_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_average_length_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_lw_a_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_lw_b_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_initial_biomass_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_length_season_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_half_saturation_constant_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_feeding_interference_coefficient_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_releaseproportion_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_basal_respiration_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_producer_competition_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_ecosystem_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_K_mean_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_K_standard_deviation_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_K_autocorrelation_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_F_q_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_assimilation_efficiency_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_maximum_consumption_rate_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_invest_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_hatchery_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_dissolution_rate_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_fraction_of_exudation_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_activity_respiration_coefficient_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_biomasstype_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function slider_F_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% ----------------------------------------------------------------------- %
% -------------------------------- Tools -------------------------------- %
% ----------------------------------------------------------------------- %

function status =  isValidInput(hObject,handles,type,low,up)
val = hObject.String;
status = true;
switch type
    case 'int'
        nVal = str2double(val);
        if isempty(nVal) || double(int16(nVal)) ~= nVal || nVal < low || nVal > up ||  isnan(nVal)
            status = false;
            hObject.String = handles.uicontrolvalues.(hObject.Tag);
            errordlg(['Invalid input. Must be integer, >= ' ...
                num2str(low) ' and <= ' num2str(up)])
        end
    case 'real'
        nVal = str2double(val);
        if isempty(nVal) || ~isreal(nVal) || nVal < low || nVal > up || isnan(nVal)
            status = false;
            hObject.String = handles.uicontrolvalues.(hObject.Tag);
            errordlg(['Invalid input. Must be real, >= ' ...
                num2str(low) ' and <= ' num2str(up)])
        end
    case 'labelstring'
        if isempty(val) || ~isvarname(['test_' val])
            status = false;
            hObject.String = handles.uicontrolvalues.(hObject.Tag);
            errordlg('Invalid input. Must be nonempty string containing letters, digits, or underscores')
        end
    case 'string'
        % no restrictions at the moment
    otherwise
        keyboard
end
function handles = saveUIControlValues(handles)
uicv.edit_fishingmortality = handles.tab_simulate.edit_fishingmortality.String;
uicv.edit_numberofyears = handles.tab_simulate.edit_numberofyears.String;
uicv.edit_axeswidth = handles.tab_simulate.edit_axeswidth.String;
uicv.popupmenu_catchtype = handles.tab_simulate.popupmenu_catchtype.Value;
uicv.edit_releaseproportion = handles.tab_simulate.edit_releaseproportion.String;
uicv.popupmenu_ecosystem = handles.tab_simulate.popupmenu_ecosystem.Value;

uicv.edit_name = handles.tab_guilds.edit_name.String;
uicv.edit_label = handles.tab_guilds.edit_label.String;
uicv.edit_intrinsic_growth_rate = handles.tab_guilds.edit_intrinsic_growth_rate.String;
uicv.edit_metabolic_rate = handles.tab_guilds.edit_metabolic_rate.String;
uicv.edit_average_length = handles.tab_guilds.edit_average_length.String;
uicv.edit_lw_a = handles.tab_guilds.edit_lw_a.String;
uicv.edit_lw_b = handles.tab_guilds.edit_lw_b.String;
uicv.edit_initial_biomass = handles.tab_guilds.edit_initial_biomass.String;
uicv.edit_basal_respiration = handles.tab_guilds.edit_basal_respiration.String;
uicv.edit_producer_competition = handles.tab_guilds.edit_producer_competition.String;
uicv.edit_half_saturation_constant = handles.tab_guilds.edit_half_saturation_constant.String;
uicv.edit_feeding_interference_coefficient = handles.tab_guilds.edit_feeding_interference_coefficient.String;
uicv.edit_F_q = handles.tab_guilds.edit_F_q.String;
uicv.edit_length_season = handles.tab_guilds.edit_length_season.String;
uicv.edit_K_mean = handles.tab_guilds.edit_K_mean.String;
uicv.edit_K_standard_deviation = handles.tab_guilds.edit_K_standard_deviation.String;
uicv.edit_K_autocorrelation = handles.tab_guilds.edit_K_autocorrelation.String;
uicv.edit_fraction_of_exudation = handles.tab_guilds.edit_fraction_of_exudation.String;
uicv.edit_dissolution_rate = handles.tab_guilds.edit_dissolution_rate.String;
uicv.edit_hatchery = handles.tab_guilds.edit_hatchery.String;
uicv.edit_activity_respiration_coefficient = handles.tab_guilds.edit_activity_respiration_coefficient.String;
uicv.edit_invest = handles.tab_guilds.edit_invest.String;

% TODO: add more if needed
handles.uicontrolvalues = uicv;
handles.isGuildChange = true;
handles.isNotSaved = true;

function enableUIControls(~, handles, status)
% enableUIControls - used to disable/enable UIControls during computation
handles.tab_simulate.radiobutton_gear_1.Enable = status;
handles.tab_simulate.radiobutton_gear_2.Enable = status;
handles.tab_simulate.radiobutton_gear_3.Enable = status;
handles.tab_simulate.edit_fishingmortality.Enable = status;
handles.tab_simulate.edit_numberofyears.Enable = status;
if handles.Game.gameover
    handles.tab_simulate.pushbutton_gofishing.Enable = 'off';
    handles.tab_simulate.pushbutton_nofishing.Enable = 'off';
else
    handles.tab_simulate.pushbutton_gofishing.Enable = status;
    handles.tab_simulate.pushbutton_nofishing.Enable = status;
end
handles.tab_simulate.pushbutton_reset.Enable = status;
handles.tab_simulate.pushbutton_planktonfigure.Enable = status;
handles.tab_simulate.pushbutton_totalbiomassfigure.Enable = status;
handles.tab_simulate.pushbutton_save_results.Enable = status;
handles.tab_simulate.pushbutton_rewind.Enable = status;
handles.tab_simulate.edit_numberofyears.Enable = status;
if handles.tab_simulate.radiobutton_gear_3.Value
    handles.tab_simulate.edit_releaseproportion.Enable = status;
else
    handles.tab_simulate.edit_releaseproportion.Enable = 'off';
end
handles.tab_simulate.popupmenu_ecosystem.Enable = status;
if strcmp(status,'on')
    handles.tab_simulate.pushbutton_stop.Enable = 'off';
else
    handles.tab_simulate.pushbutton_stop.Enable = 'on';
end
drawnow
function handles = simulateForward(handles)

nYearsFwd = str2double(handles.tab_simulate.edit_numberofyears.String);
handles.Data.nYearsFwd = nYearsFwd;
handles.Data.nGrowthDays = str2double(handles.tab_guilds.edit_length_season.String);
tspan = 0:handles.Data.nGrowthDays;
GuildInfo = handles.Data.GuildInfo;
Guilds = handles.Data.Guilds;
Data = compileOdeData(handles.Data);

Binds = 1:GuildInfo.nGuilds;
Cinds = Binds(end)+(1:GuildInfo.nAdultFishGuilds);
GFinds = Cinds(end)+(1:GuildInfo.nFishGuilds);
Ginds = GFinds(end)+(1:GuildInfo.nFeederGuilds*GuildInfo.nFoodGuilds);
Linds = Ginds(end)+(1:GuildInfo.nFoodGuilds*GuildInfo.nFeederGuilds);

for i = 1:nYearsFwd
    
    if handles.tab_simulate.pushbutton_stop.UserData.simulationStopped
        handles.tab_simulate.text_simulationstatus.String = ...
            ['Status: simulation stopped (' num2str(i) '/' num2str(nYearsFwd) ')'];
        drawnow
        handles.tab_simulate.pushbutton_stop.UserData.simulationStopped = false;
        return
    end
    
    Data.year = i;
    
    handles.tab_simulate.text_simulationstatus.String = ...
        ['Status: running (' num2str(i) '/' num2str(nYearsFwd) ')'];
    drawnow
    
    %% One year at a time with catch
    Binit = handles.Simulation.Binit;
    Cinit = zeros(GuildInfo.nAdultFishGuilds,1);
    GFinit = zeros(GuildInfo.nFishGuilds,1);
    Ginit = zeros(GuildInfo.nFeederGuilds*GuildInfo.nFoodGuilds,1);
    Linit = zeros(GuildInfo.nFoodGuilds*GuildInfo.nFeederGuilds,1);
    
    [~,BC] = ode23(@(t,BC)atn_ode_biomass_catch_gainfish_gain_loss_fast(t,BC,Data), ...
        tspan,[Binit ; Cinit ; GFinit ; Ginit ; Linit]);
    
    if ~isreal(BC) || any(any(isnan(BC)))
        herror = errordlg('Differential eqution solver returned non real biomasses. Check model parameters.');
        uiwait(herror)
        handles = resetSimulation(handles);
        return
    end
    
    %%%%%%% Fix for constraining the biomass of dissolved detritus %%%%%%%%
    bDOC2_max = 1e6;
    BC(end,GuildInfo.iDOC2) = min([BC(end,GuildInfo.iDOC2) bDOC2_max]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.Results.B(:,end+1) = BC(end,Binds)';
    handles.Results.C(:,end+1) = BC(end,Cinds)';
    handles.Results.GF(:,end+1) = BC(end,GFinds)';
    handles.Results.G(GuildInfo.iFeederGuilds,GuildInfo.iFoodGuilds,end+1) = ...
        reshape(BC(end,Ginds)',GuildInfo.nFeederGuilds,GuildInfo.nFoodGuilds);
    handles.Results.L(GuildInfo.iFoodGuilds,GuildInfo.iFeederGuilds,end+1) = ...
        reshape(BC(end,Linds)',GuildInfo.nFoodGuilds,GuildInfo.nFeederGuilds);
    
    ARE = zeros(GuildInfo.nAdultFishGuilds,1);
    for k = 1:GuildInfo.nAdultFishGuilds
        IA = GuildInfo.iAdultFishGuilds(k);
        JA = find(GuildInfo.iFishGuilds == IA);
        age = Data.Guilds(IA).age;
        
        ARE(k) = Data.invest(GuildInfo.iAdultFishGuildsInFish(k))*handles.Results.GF(JA,end)* ...
            BC(end,IA)*Data.P_mat(age+1); %#ok<FNDSB> %biomass of newborns
        if isnan(ARE(k))
            if age < 4
                ARE(k) = 0.001*BC(end,IA);
            else
                ARE(k) = 0.005*BC(end,IA);
            end
        end
    end
    
    handles.Simulation.Binit(GuildInfo.iDetritusGuilds) = ...
        BC(end,GuildInfo.iDetritusGuilds);
    handles.Simulation.Binit(GuildInfo.iProducerGuilds) = ...
        BC(end,GuildInfo.iProducerGuilds);
    handles.Simulation.Binit(GuildInfo.iConsumerGuilds) = ...
        BC(end,GuildInfo.iConsumerGuilds);
    
    for k = 1:GuildInfo.nLarvaeFishGuilds
        IL = GuildInfo.iLarvaeFishGuilds(k);
        gl = Data.Guilds(IL).label;
        gl = gl(1:3);
        
        gal = {Data.Guilds(GuildInfo.iAdultFishGuilds).label};
        newbornBM = 0;
        for j = 1:length(gal)
            if strcmpi(gal{j}(1:3),gl)
                newbornBM = newbornBM + ARE(j);
            end
        end
        iInFish = find(ismember(GuildInfo.iFishGuilds,IL));
        handles.Simulation.Binit(IL) = newbornBM + Data.hatchery(iInFish); %#ok<FNDSB>
    end
    for k = 1:GuildInfo.nJuvenileFishGuilds
        IJ = GuildInfo.iJuvenileFishGuilds(k);
        gl = Data.Guilds(IJ).label;
        gl = gl(1:3);
        
        gal = {Data.Guilds(GuildInfo.iLarvaeFishGuilds).label};
        for j = 1:length(gal)
            if strcmpi(gal{j}(1:3),gl) && gal{j}(4) == '0'
                IL = GuildInfo.iLarvaeFishGuilds(j);
                iInFish = find(ismember(GuildInfo.iFishGuilds,IJ));
                handles.Simulation.Binit(IJ) = BC(end,IL) + Data.hatchery(iInFish); %#ok<FNDSB>
            end
        end
    end
    for k = 1:GuildInfo.nAdultFishGuilds
        IA = GuildInfo.iAdultFishGuilds(k);
        gl = Data.Guilds(IA).label;
        ga = str2double(gl(4));
        gl = gl(1:3);
        
        iJA = [GuildInfo.iJuvenileFishGuilds GuildInfo.iAdultFishGuilds];
        gal = {Data.Guilds(iJA).label};
        for j = 1:length(gal)
            if ga < 4 && strcmpi(gal{j}(1:3),gl) && gal{j}(4) == num2str(ga-1)
                IJA = iJA(j);
                iInFish = find(ismember(GuildInfo.iFishGuilds,IA));
                handles.Simulation.Binit(IA) = BC(end,IJA) + Data.hatchery(iInFish); %#ok<FNDSB>
            elseif ga == 4 && strcmpi(gal{j}(1:3),gl) && gal{j}(4) == num2str(ga-1)
                IJA = iJA(j);
                iInFish = find(ismember(GuildInfo.iFishGuilds,IA));
                handles.Simulation.Binit(IA) = BC(end,IA) + BC(end,IJA) + Data.hatchery(iInFish); %#ok<FNDSB>
            end
            
        end
    end
    
    handles.Results.B(:,end) = handles.Simulation.Binit;
    
    handles = refreshAxes(handles); % Refresh plots at each time step
    
    handles.Results.ARE(:,end+1) = ARE(:);
    handles.Results.allbiomasses = [handles.Results.allbiomasses BC(:,1:GuildInfo.nGuilds)'];
    
    %%% UPDATE THE GAME STATUS %%%
    isGameMode = handles.tab_simulate.checkbox_gamemode.Value;
    
    handles.Game.currentyear = handles.Game.currentyear + 1;
    gameStructDummy = initializeGameStruct();
    handles.Game.history(end+1) = gameStructDummy.history;
    handles.Game.history(end).balance = handles.Game.balance;
    
    if isGameMode
                
        % Calculate profits
        
        % How much money is carried to the next year (no yearly costs yet)
        if any(handles.Data.E)
            handles.Game.history(end).F = handles.Data.F;
            handles.Game.history(end).gear = handles.tab_simulate.uibuttongroup_gear.SelectedObject.String;
            if strcmpi(handles.tab_simulate.uibuttongroup_gear.SelectedObject.String,'Fyke net')
                handles.Game.history(end).prop_released = str2double(handles.tab_simulate.edit_releaseproportion.String);
            else
                handles.Game.history(end).prop_released = NaN;
            end
        else
            handles.Game.history(end).F = 0;
            handles.Game.history(end).gear = '';
            handles.Game.history(end).prop_released = NaN;
        end
        
        fishgls = {Guilds(GuildInfo.iFishGuilds).label};
        S.subs = {1,1:3}; S.type = '()';
        fslabels = unique(cellfun(@(x)subsref(x,S),fishgls,'uniformoutput',false),'stable');
        for ifs = 1:GuildInfo.nFishSpecies
            thisFishSpeciesLabel = fslabels{ifs};
            thisFishGuildInds = GuildInfo.iFishGuilds(strcmpi(fslabels{ifs},cellfun(@(x)subsref(x,S),fishgls,'uniformoutput',false)));
            
            % Calculate profits from catch
            afishgls = {Guilds(GuildInfo.iAdultFishGuilds).label};
            afishgls = cellfun(@(x)subsref(x,S),afishgls,'uniformoutput',false);
            fsinds = find(strcmpi(fslabels{ifs},afishgls));
            thisFishSpeciesCatch = sum(handles.Results.C(fsinds,end)); %#ok<FNDSB>
            thisProfit = thisFishSpeciesCatch * handles.Game.Settings.FISHPRICES{strcmpi(thisFishSpeciesLabel,handles.Game.Settings.FISHPRICES(:,1)),2};
            handles.Game.history(end).profit(ifs) = thisProfit;
            
            handles.Game.history(end).balance = handles.Game.history(end).balance + thisProfit;
            
            % Check fish population status
            thisFishGuildTotalBiomass = sum(handles.Results.B(thisFishGuildInds,end));
            thisFishGuildTotalInitialBiomass = sum([handles.Data.Guilds(thisFishGuildInds).binit]);
            if thisFishGuildTotalBiomass < handles.Game.Settings.FISHPENALTYLIMIT2 * thisFishGuildTotalInitialBiomass
                hwarn = warndlg(['Total biomass of fish species ' fslabels{ifs} ...
                    ' is below ' num2str(handles.Game.Settings.FISHPENALTYLIMIT2*100) ...
                    '% of its equilibrium. That''s a ' num2str(handles.Game.Settings.OVERFISH_PENALTY2) ...
                    ' gold penalty.']);
                uiwait(hwarn)
                % That's a penalty
                handles.Game.history(end).penalty(end+1).type = 'Fish species biomass violation 2';
                handles.Game.history(end).penalty(end).details = thisFishSpeciesLabel;
                handles.Game.history(end).penalty(end).amount = handles.Game.Settings.OVERFISH_PENALTY2;
                handles.Game.history(end).balance = handles.Game.history(end).balance - handles.Game.Settings.OVERFISH_PENALTY2;
            elseif thisFishGuildTotalBiomass < handles.Game.Settings.FISHPENALTYLIMIT * thisFishGuildTotalInitialBiomass
                hwarn = warndlg(['Total biomass of fish species ' fslabels{ifs} ...
                    ' is below ' num2str(handles.Game.Settings.FISHPENALTYLIMIT*100) ...
                    '% of its equilibrium. That''s a ' num2str(handles.Game.Settings.OVERFISH_PENALTY) ...
                    ' gold penalty.']);
                uiwait(hwarn)
                % That's a penalty
                handles.Game.history(end).penalty(end+1).type = 'Fish species biomass violation';
                handles.Game.history(end).penalty(end).details = thisFishSpeciesLabel;
                handles.Game.history(end).penalty(end).amount = handles.Game.Settings.OVERFISH_PENALTY;
                handles.Game.history(end).balance = handles.Game.history(end).balance - handles.Game.Settings.OVERFISH_PENALTY;
            end
            if thisFishGuildTotalBiomass <  handles.Game.Settings.SPECIESCRASHLIMIT * thisFishGuildTotalInitialBiomass;
                herror = errordlg([thisFishSpeciesLabel ' population crashed. Game over!']);
                uiwait(herror)
                % That's a penalty
                handles.Game.history(end).penalty(end+1).type = 'Fish species population crash';
                handles.Game.history(end).penalty(end).details = thisFishSpeciesLabel;
                handles.Game.history(end).penalty(end).amount = [];
                handles.Game.gameover = 1;
                return
            end
        end
        
        
        % Check ecosystem status
        ecoSystemIndicator = (sum(handles.Results.B(:,end)) / sum([handles.Data.Guilds.binit]))^handles.Game.Settings.ECOEXP;
        if ecoSystemIndicator < handles.Game.Settings.ECOSYSINDLIMIT
            hwarn = warndlg(['Ecosystem health dropped below sustainable. That''s a ' ...
                num2str(handles.Game.Settings.ECOSYS_PENALTY) ...
                ' gold penalty.']);
            uiwait(hwarn)
            % That's a penalty
            handles.Game.history(end).penalty(end+1).type = 'Ecosystem abuse';
            handles.Game.history(end).penalty(end).details = '';
            handles.Game.history(end).penalty(end).amount = handles.Game.Settings.ECOSYS_PENALTY;
            handles.Game.history(end).balance = handles.Game.history(end).balance - handles.Game.Settings.ECOSYS_PENALTY;
        end
        
        handles.Game.balance = handles.Game.history(end).balance;
        
        handles.tab_simulate.text_balance.String = ['Balance: ' num2str(round(handles.Game.balance))];
        
        if handles.Game.balance < 0
            herror = errordlg('Balance negative. Game over!');
            uiwait(herror)
            % That's a penalty
            handles.Game.history(end).penalty(end+1).type = 'Bankcruptcy';
            handles.Game.history(end).penalty(end).details = '';
            handles.Game.history(end).penalty(end).amount = [];
            handles.Game.gameover = 1;
            return
        end
        
        iNonFishNonDetritus = setdiff(1:GuildInfo.nGuilds,[GuildInfo.iFishGuilds GuildInfo.iDetritusGuilds]);
        for is = iNonFishNonDetritus
            if handles.Results.B(is,end) < handles.Game.Settings.SPECIESCRASHLIMIT * handles.Data.Guilds(is).binit;
                herror = errordlg([handles.Data.Guilds(is).name ' population crashed. Game over!']);
                uiwait(herror)
                % That's a penalty
                handles.Game.history(end).penalty(end+1).type = 'Plankton species crash';
                handles.Game.history(end).penalty(end).details = Guilds(is).name;
                handles.Game.history(end).penalty(end).amount = [];
                handles.Game.gameover = 1;
                return
            end
        end
        
        if handles.Game.currentyear == handles.Game.Settings.MAXYEAR;
            handles.Game.gameover = 1;
            herror = errordlg('You have reached the end of the game. Game over!');
            uiwait(herror)
            return
        end
        
    end
end
handles.tab_simulate.text_simulationstatus.String = 'Status: ready';
function handles = groupHandles(handles)

% ------ SIMULATE TAB UICONTROLS
tab_simulate.all(1) = handles.uipanel_tab_simulate;
tab_simulate.all(end+1) = handles.axes_biomass;
tab_simulate.all(end+1) = handles.axes_catch;
tab_simulate.all(end+1) = handles.uibuttongroup_gear;
tab_simulate.all(end+1) = handles.radiobutton_gear_1;
tab_simulate.all(end+1) = handles.radiobutton_gear_2;
tab_simulate.all(end+1) = handles.radiobutton_gear_3;
tab_simulate.all(end+1) = handles.uipanel_fishingpressure;
tab_simulate.all(end+1) = handles.edit_fishingmortality;
tab_simulate.all(end+1) = handles.text_fishingmortality;
tab_simulate.all(end+1) = handles.uipanel_simulate;
tab_simulate.all(end+1) = handles.edit_numberofyears;
tab_simulate.all(end+1) = handles.text_numberofyears;
tab_simulate.all(end+1) = handles.pushbutton_gofishing;
tab_simulate.all(end+1) = handles.pushbutton_nofishing;
tab_simulate.all(end+1) = handles.text_axes_biomass;
tab_simulate.all(end+1) = handles.text_axes_catch;
tab_simulate.all(end+1) = handles.text_simulationstatus;
tab_simulate.all(end+1) = handles.pushbutton_reset;
tab_simulate.all(end+1) = handles.popupmenu_catchtype;
tab_simulate.all(end+1) = handles.text_axeswidth;
tab_simulate.all(end+1) = handles.edit_axeswidth;
tab_simulate.all(end+1) = handles.uipanel_biomassplots;
tab_simulate.all(end+1) = handles.pushbutton_save_results;
tab_simulate.all(end+1) = handles.text_releaseproportion;
tab_simulate.all(end+1) = handles.edit_releaseproportion;
tab_simulate.all(end+1) = handles.pushbutton_planktonfigure;
tab_simulate.all(end+1) = handles.axes_bgimage_simulate;
tab_simulate.all(end+1) = handles.popupmenu_ecosystem;
tab_simulate.all(end+1) = handles.pushbutton_rewind;
tab_simulate.all(end+1) = handles.pushbutton_stop;
tab_simulate.all(end+1) = handles.pushbutton_totalbiomassfigure;
tab_simulate.all(end+1) = handles.text_balance;
tab_simulate.all(end+1) = handles.popupmenu_biomasstype;
tab_simulate.all(end+1) = handles.checkbox_gamemode;
tab_simulate.all(end+1) = handles.slider_F;

tab_simulate.uipanel_tab_simulate = handles.uipanel_tab_simulate;
tab_simulate.axes_biomass = handles.axes_biomass;
tab_simulate.axes_catch = handles.axes_catch;
tab_simulate.uibuttongroup_gear = handles.uibuttongroup_gear;
tab_simulate.radiobutton_gear_1 = handles.radiobutton_gear_1;
tab_simulate.radiobutton_gear_2 = handles.radiobutton_gear_2;
tab_simulate.radiobutton_gear_3 = handles.radiobutton_gear_3;
tab_simulate.uipanel_fishingpressure = handles.uipanel_fishingpressure;
tab_simulate.edit_fishingmortality = handles.edit_fishingmortality;
tab_simulate.text_fishingmortality = handles.text_fishingmortality;
tab_simulate.uipanel_simulate = handles.uipanel_simulate;
tab_simulate.edit_numberofyears = handles.edit_numberofyears;
tab_simulate.text_numberofyears = handles.text_numberofyears;
tab_simulate.pushbutton_gofishing = handles.pushbutton_gofishing;
tab_simulate.pushbutton_nofishing = handles.pushbutton_nofishing;
tab_simulate.text_axes_biomass = handles.text_axes_biomass;
tab_simulate.text_axes_catch = handles.text_axes_catch;
tab_simulate.text_simulationstatus = handles.text_simulationstatus;
tab_simulate.pushbutton_reset = handles.pushbutton_reset;
tab_simulate.popupmenu_catchtype = handles.popupmenu_catchtype;
tab_simulate.text_axeswidth = handles.text_axeswidth;
tab_simulate.edit_axeswidth = handles.edit_axeswidth;
tab_simulate.pushbutton_save_results = handles.pushbutton_save_results;
tab_simulate.uipanel_biomassplots = handles.uipanel_biomassplots;
tab_simulate.text_releaseproportion = handles.text_releaseproportion;
tab_simulate.edit_releaseproportion = handles.edit_releaseproportion;
tab_simulate.pushbutton_planktonfigure = handles.pushbutton_planktonfigure;
tab_simulate.axes_bgimage_simulate = handles.axes_bgimage_simulate;
tab_simulate.popupmenu_ecosystem = handles.popupmenu_ecosystem;
tab_simulate.pushbutton_rewind = handles.pushbutton_rewind;
tab_simulate.pushbutton_stop = handles.pushbutton_stop;
tab_simulate.pushbutton_totalbiomassfigure = handles.pushbutton_totalbiomassfigure;
tab_simulate.text_balance = handles.text_balance;
tab_simulate.popupmenu_biomasstype = handles.popupmenu_biomasstype;
tab_simulate.checkbox_gamemode = handles.checkbox_gamemode;
tab_simulate.slider_F = handles.slider_F;

handles = rmfield(handles,{'uipanel_tab_simulate','axes_biomass','axes_catch', ...
    'uibuttongroup_gear','radiobutton_gear_1','radiobutton_gear_2','radiobutton_gear_3', ...
    'uipanel_fishingpressure','edit_fishingmortality', ...
    'text_fishingmortality','uipanel_simulate','edit_numberofyears','pushbutton_gofishing', ...
    'pushbutton_nofishing','text_axes_biomass','text_axes_catch', ...
    'text_simulationstatus', 'pushbutton_reset', 'popupmenu_catchtype', ...
    'text_axeswidth', 'edit_axeswidth', 'text_numberofyears', ...
    'uipanel_biomassplots','pushbutton_save_results', ...
    'text_releaseproportion', 'edit_releaseproportion','pushbutton_planktonfigure', ...
    'axes_bgimage_simulate','popupmenu_ecosystem','pushbutton_rewind', ...
    'pushbutton_stop','pushbutton_totalbiomassfigure','text_balance', ...
    'popupmenu_biomasstype','checkbox_gamemode','slider_F'});
handles.tab_simulate = tab_simulate;


% ------ GUILDS TAB UICONTROLS
tab_guilds.all(1) = handles.uipanel_tab_guilds;
tab_guilds.all(end+1) = handles.uipanel_add_guild;
tab_guilds.all(end+1) = handles.text_name;
tab_guilds.all(end+1) = handles.text_label;
tab_guilds.all(end+1) = handles.text_intrinsic_growth_rate;
tab_guilds.all(end+1) = handles.text_metabolic_rate;
tab_guilds.all(end+1) = handles.edit_name;
tab_guilds.all(end+1) = handles.edit_label;
tab_guilds.all(end+1) = handles.edit_intrinsic_growth_rate;
tab_guilds.all(end+1) = handles.edit_metabolic_rate;
tab_guilds.all(end+1) = handles.pushbutton_add_guild;
tab_guilds.all(end+1) = handles.pushbutton_save_changes;
tab_guilds.all(end+1) = handles.pushbutton_clear;
tab_guilds.all(end+1) = handles.uibuttongroup_type;
tab_guilds.all(end+1) = handles.uipanel_guild_list;
tab_guilds.all(end+1) = handles.listbox_guild_list;
tab_guilds.all(end+1) = handles.pushbutton_delete;
tab_guilds.all(end+1) = handles.pushbutton_import;
tab_guilds.all(end+1) = handles.pushbutton_export;
tab_guilds.all(end+1) = handles.radiobutton_consumer;
tab_guilds.all(end+1) = handles.radiobutton_producer;
tab_guilds.all(end+1) = handles.radiobutton_detritus;
tab_guilds.all(end+1) = handles.listbox_prey;
tab_guilds.all(end+1) = handles.listbox_available_prey;
tab_guilds.all(end+1) = handles.pushbutton_add_food_item;
tab_guilds.all(end+1) = handles.pushbutton_remove_food_item;
tab_guilds.all(end+1) = handles.pushbutton_clearall;
tab_guilds.all(end+1) = handles.radiobutton_fish;
tab_guilds.all(end+1) = handles.text_average_length;
tab_guilds.all(end+1) = handles.edit_average_length;
tab_guilds.all(end+1) = handles.text_length_weight_parameters;
tab_guilds.all(end+1) = handles.edit_lw_a;
tab_guilds.all(end+1) = handles.edit_lw_b;
tab_guilds.all(end+1) = handles.text_lw_a;
tab_guilds.all(end+1) = handles.text_lw_b;
tab_guilds.all(end+1) = handles.edit_initial_biomass;
tab_guilds.all(end+1) = handles.text_initial_biomass;
tab_guilds.all(end+1) = handles.pushbutton_moveup;
tab_guilds.all(end+1) = handles.pushbutton_movedown;
tab_guilds.all(end+1) = handles.text_half_saturation_constant;
tab_guilds.all(end+1) = handles.edit_half_saturation_constant;
tab_guilds.all(end+1) = handles.text_feeding_interference_coefficient;
tab_guilds.all(end+1) = handles.edit_feeding_interference_coefficient;
tab_guilds.all(end+1) = handles.axes_bgimage_guilds;
tab_guilds.all(end+1) = handles.edit_length_season;
tab_guilds.all(end+1) = handles.text_length_season;
tab_guilds.all(end+1) = handles.edit_basal_respiration;
tab_guilds.all(end+1) = handles.text_basal_respiration;
tab_guilds.all(end+1) = handles.edit_producer_competition;
tab_guilds.all(end+1) = handles.text_producer_competition;
tab_guilds.all(end+1) = handles.uibuttongroup_carrying_capacity;
tab_guilds.all(end+1) = handles.radiobutton_K_constant;
tab_guilds.all(end+1) = handles.radiobutton_K_white_noise;
tab_guilds.all(end+1) = handles.radiobutton_K_AR1;
tab_guilds.all(end+1) = handles.text_K_mean;
tab_guilds.all(end+1) = handles.edit_K_mean;
tab_guilds.all(end+1) = handles.text_K_standard_deviation;
tab_guilds.all(end+1) = handles.edit_K_standard_deviation;
tab_guilds.all(end+1) = handles.text_K_autocorrelation;
tab_guilds.all(end+1) = handles.edit_K_autocorrelation;
tab_guilds.all(end+1) = handles.text_F_q;
tab_guilds.all(end+1) = handles.edit_F_q;
tab_guilds.all(end+1) = handles.text_assimilation_efficiency;
tab_guilds.all(end+1) = handles.edit_assimilation_efficiency;
tab_guilds.all(end+1) = handles.text_maximum_consumption_rate;
tab_guilds.all(end+1) = handles.edit_maximum_consumption_rate;
tab_guilds.all(end+1) = handles.text_fraction_of_exudation;
tab_guilds.all(end+1) = handles.edit_fraction_of_exudation;
tab_guilds.all(end+1) = handles.text_dissolution_rate;
tab_guilds.all(end+1) = handles.edit_dissolution_rate;
tab_guilds.all(end+1) = handles.text_hatchery;
tab_guilds.all(end+1) = handles.edit_hatchery;
tab_guilds.all(end+1) = handles.text_invest;
tab_guilds.all(end+1) = handles.edit_invest;
tab_guilds.all(end+1) = handles.text_activity_respiration_coefficient;
tab_guilds.all(end+1) = handles.edit_activity_respiration_coefficient;
tab_guilds.all(end+1) = handles.checkbox_catchable;

tab_guilds.uipanel_tab_guilds = handles.uipanel_tab_guilds;
tab_guilds.uipanel_add_guild = handles.uipanel_add_guild;
tab_guilds.text_name = handles.text_name;
tab_guilds.text_label = handles.text_label;
tab_guilds.text_intrinsic_growth_rate = handles.text_intrinsic_growth_rate;
tab_guilds.text_metabolic_rate = handles.text_metabolic_rate;
tab_guilds.edit_name = handles.edit_name;
tab_guilds.edit_label = handles.edit_label;
tab_guilds.edit_intrinsic_growth_rate = handles.edit_intrinsic_growth_rate;
tab_guilds.edit_metabolic_rate = handles.edit_metabolic_rate;
tab_guilds.pushbutton_add_guild = handles.pushbutton_add_guild;
tab_guilds.pushbutton_save_changes = handles.pushbutton_save_changes;
tab_guilds.pushbutton_clear = handles.pushbutton_clear;
tab_guilds.uibuttongroup_type = handles.uibuttongroup_type;
tab_guilds.uipanel_guild_list = handles.uipanel_guild_list;
tab_guilds.listbox_guild_list = handles.listbox_guild_list;
tab_guilds.pushbutton_delete = handles.pushbutton_delete;
tab_guilds.pushbutton_import = handles.pushbutton_import;
tab_guilds.pushbutton_export = handles.pushbutton_export;
tab_guilds.radiobutton_consumer = handles.radiobutton_consumer;
tab_guilds.radiobutton_producer = handles.radiobutton_producer;
tab_guilds.radiobutton_detritus = handles.radiobutton_detritus;
tab_guilds.listbox_prey = handles.listbox_prey;
tab_guilds.listbox_available_prey = handles.listbox_available_prey;
tab_guilds.pushbutton_add_food_item = handles.pushbutton_add_food_item;
tab_guilds.pushbutton_remove_food_item = handles.pushbutton_remove_food_item;
tab_guilds.pushbutton_clearall = handles.pushbutton_clearall;
tab_guilds.radiobutton_fish = handles.radiobutton_fish;
tab_guilds.text_average_length = handles.text_average_length;
tab_guilds.edit_average_length = handles.edit_average_length;
tab_guilds.text_length_weight_parameters = handles.text_length_weight_parameters;
tab_guilds.edit_lw_a = handles.edit_lw_a;
tab_guilds.edit_lw_b = handles.edit_lw_b;
tab_guilds.text_lw_a = handles.text_lw_a;
tab_guilds.text_lw_b = handles.text_lw_b;
tab_guilds.edit_initial_biomass = handles.edit_initial_biomass;
tab_guilds.text_initial_biomass = handles.text_initial_biomass;
tab_guilds.pushbutton_moveup = handles.pushbutton_moveup;
tab_guilds.pushbutton_movedown = handles.pushbutton_movedown;
tab_guilds.text_half_saturation_constant = handles.text_half_saturation_constant;
tab_guilds.edit_half_saturation_constant = handles.edit_half_saturation_constant;
tab_guilds.text_feeding_interference_coefficient = handles.text_feeding_interference_coefficient;
tab_guilds.edit_feeding_interference_coefficient = handles.edit_feeding_interference_coefficient;
tab_guilds.axes_bgimage_guilds = handles.axes_bgimage_guilds;
tab_guilds.edit_length_season = handles.edit_length_season;
tab_guilds.text_length_season = handles.text_length_season;
tab_guilds.edit_basal_respiration = handles.edit_basal_respiration;
tab_guilds.text_basal_respiration = handles.text_basal_respiration;
tab_guilds.edit_producer_competition = handles.edit_producer_competition;
tab_guilds.text_producer_competition = handles.text_producer_competition;
tab_guilds.uibuttongroup_carrying_capacity = handles.uibuttongroup_carrying_capacity;
tab_guilds.radiobutton_K_constant = handles.radiobutton_K_constant;
tab_guilds.radiobutton_K_white_noise = handles.radiobutton_K_white_noise;
tab_guilds.radiobutton_K_AR1 = handles.radiobutton_K_AR1;
tab_guilds.text_K_mean = handles.text_K_mean;
tab_guilds.edit_K_mean = handles.edit_K_mean;
tab_guilds.text_K_standard_deviation = handles.text_K_standard_deviation;
tab_guilds.edit_K_standard_deviation = handles.edit_K_standard_deviation;
tab_guilds.text_K_autocorrelation = handles.text_K_autocorrelation;
tab_guilds.edit_K_autocorrelation = handles.edit_K_autocorrelation;
tab_guilds.text_F_q = handles.text_F_q;
tab_guilds.edit_F_q = handles.edit_F_q;
tab_guilds.text_assimilation_efficiency = handles.text_assimilation_efficiency;
tab_guilds.edit_assimilation_efficiency = handles.edit_assimilation_efficiency;
tab_guilds.text_maximum_consumption_rate = handles.text_maximum_consumption_rate;
tab_guilds.edit_maximum_consumption_rate = handles.edit_maximum_consumption_rate;
tab_guilds.text_fraction_of_exudation = handles.text_fraction_of_exudation;
tab_guilds.edit_fraction_of_exudation = handles.edit_fraction_of_exudation;
tab_guilds.text_dissolution_rate = handles.text_dissolution_rate;
tab_guilds.edit_dissolution_rate = handles.edit_dissolution_rate;
tab_guilds.text_hatchery = handles.text_hatchery;
tab_guilds.edit_hatchery = handles.edit_hatchery;
tab_guilds.text_invest = handles.text_invest;
tab_guilds.edit_invest = handles.edit_invest;
tab_guilds.text_activity_respiration_coefficient = handles.text_activity_respiration_coefficient;
tab_guilds.edit_activity_respiration_coefficient = handles.edit_activity_respiration_coefficient;
tab_guilds.checkbox_catchable = handles.checkbox_catchable;

handles = rmfield(handles,{'uipanel_tab_guilds','uipanel_add_guild','text_name', ...
    'text_label','text_intrinsic_growth_rate','text_metabolic_rate','edit_name', ...
    'edit_label','edit_intrinsic_growth_rate','edit_metabolic_rate','pushbutton_add_guild', ...
    'pushbutton_save_changes','pushbutton_clear','uibuttongroup_type', ...
    'uipanel_guild_list','listbox_guild_list','pushbutton_delete', ...
    'pushbutton_import','pushbutton_export','radiobutton_consumer', ...
    'radiobutton_producer','radiobutton_detritus','listbox_prey','listbox_available_prey', ...
    'pushbutton_add_food_item','pushbutton_remove_food_item','pushbutton_clearall', ...
    'radiobutton_fish','text_average_length','edit_average_length', ...
    'text_length_weight_parameters','edit_lw_a','edit_lw_b','text_lw_a','text_lw_b', ...
    'edit_initial_biomass','text_initial_biomass','pushbutton_moveup','pushbutton_movedown', ...
    'edit_half_saturation_constant','text_half_saturation_constant', ...
    'edit_feeding_interference_coefficient','text_feeding_interference_coefficient', ...
    'axes_bgimage_guilds','edit_length_season','text_length_season', ...
    'edit_basal_respiration','text_basal_respiration','edit_producer_competition', ...
    'text_producer_competition','uibuttongroup_carrying_capacity', ...
    'radiobutton_K_constant','radiobutton_K_white_noise','radiobutton_K_AR1', ...
    'text_K_mean','edit_K_mean','text_K_standard_deviation','edit_K_standard_deviation', ...
    'text_K_autocorrelation','edit_K_autocorrelation','text_F_q','edit_F_q', ...
    'text_assimilation_efficiency','edit_assimilation_efficiency', ...
    'text_maximum_consumption_rate','edit_maximum_consumption_rate', ...
    'text_fraction_of_exudation','edit_fraction_of_exudation', ...
    'text_dissolution_rate','edit_dissolution_rate','text_hatchery', ...
    'edit_hatchery','text_invest','edit_invest', ...
    'text_activity_respiration_coefficient','edit_activity_respiration_coefficient', ...
    'checkbox_catchable'});


handles.tab_guilds = tab_guilds;

handles = setUIControlPositions(handles);
function handles = resetSimulation(handles)

handles.Game = initializeGameStruct();
handles.Game.Settings = setGameSettings();
handles.tab_simulate.text_balance.String = 'Balance: 0';

GuildInfo = handles.Data.GuildInfo;
handles.Simulation.Binit = vertcat(handles.Data.Guilds.binit);
handles.Results = initializeResultStruct(handles.Data.GuildInfo);

% Delete existing children of axis and remove all handles
handles.tab_simulate.all(ismember(handles.tab_simulate.all, ...
    [handles.tab_simulate.axes_biomass.Children ; handles.tab_simulate.axes_catch.Children])) = [];
handles.plots.BiomassFishGuild = plot([]);
handles.plots.BiomassFishGuildPerSpecies = plot([]);
handles.plots.CatchFishGuildYearly = plot([]);
handles.plots.CatchFishGuildCumulative = plot([]);
handles.plots.CatchFishGuildYearlyPerSpecies = plot([]);
handles.plots.CatchFishGuildCumulativePerSpecies = plot([]);
delete([handles.tab_simulate.axes_biomass.Children ; handles.tab_simulate.axes_catch.Children])

% Delete existing checkboxes from CHECKBOX PANEL
delete(handles.tab_simulate.uipanel_biomassplots.Children)
handles.tab_simulate.plotCheckboxes = [];

% Initialize plots
set(handles.tab_simulate.axes_biomass,'nextplot','add','xlim',[0 handles.minaxeswidth])
set(handles.tab_simulate.axes_catch,'nextplot','add','xlim',[0 handles.minaxeswidth])

% Initialize guild LISTBOX and CHECKBOX PANEL
guildLabels = {handles.Data.Guilds.label};
positions = repmat([10 10 50 15],GuildInfo.nFishGuilds,1);
positions(:,2) = positions(:,2) + (0:20:20*(GuildInfo.nFishGuilds-1))';

% Define line properties
linecolors = [0 1 0
    0.0345 0.7931 1
    0.3103 0.9655 1
    1 0.8621 0.4138];
linestyles = {'--','--','-','-','-'};
markers = {'.','.','.','.','.'};
linewidths = [1 2 1 2 3];

fishgls = {handles.Data.Guilds(GuildInfo.iFishGuilds).label};
S.subs = {1,1:3}; S.type = '()';
fslabels = unique(cellfun(@(x)subsref(x,S),fishgls,'uniformoutput',false),'stable');

% Biomass plots for all fish guilds
for i = 1:GuildInfo.nFishGuilds
    I = GuildInfo.iFishGuilds(i);
    g = handles.Data.Guilds(I);
    li = g.age+1;
    IS = find(strcmpi(g.label(1:3),fslabels));
    
    handles.plots.BiomassFishGuild(i) = plot(handles.tab_simulate.axes_biomass, ...
        0,handles.Data.Guilds(I).binit,'marker',markers{li}, ...
        'color',linecolors(IS,:),'linewidth',linewidths(li), ...
        'linestyle',linestyles{li}); %#ok<FNDSB>
    handles.tab_simulate.all(end+1) = handles.plots.BiomassFishGuild(i);
    
    handles.tab_simulate.plotCheckboxes(i) = uicontrol( ...
        'Parent', handles.tab_simulate.uipanel_biomassplots, ...
        'Style', 'checkbox', ...
        'String', guildLabels{I}, ...
        'Callback', @(hObj,ed)ATN_gui_tabs_10('checkbox_biomassplot_Callback',hObj,ed,guidata(hObj)), ...
        'Value', 1, ...
        'Position', positions(i,:), ...
        'BackgroundColor', handles.plots.BiomassFishGuild(i).Color);
    
end
handles.tab_simulate.uipanel_biomassplots.Units = 'pixels';
handles.tab_simulate.uipanel_biomassplots.Position(4) = ...
    2*positions(1,2) + positions(end,2) + positions(end,4);
handles.tab_simulate.uipanel_biomassplots.Position(3) = ...
    2*positions(1,1) + positions(1,3);
handles.tab_simulate.uipanel_biomassplots.Parent = handles.tab_simulate.uipanel_tab_simulate;

% Catch plots for all adult fish guilds
for i = 1:GuildInfo.nAdultFishGuilds
    I = GuildInfo.iAdultFishGuilds(i);
    g = handles.Data.Guilds(I);
    li = g.age+1;
    IS = find(strcmpi(g.label(1:3),fslabels));
    
    handles.plots.CatchFishGuildYearly(i) = plot(handles.tab_simulate.axes_catch, ...
        0, 0, ...
        'marker', markers{li}, ...
        'color', linecolors(IS,:), ...
        'linewidth', linewidths(li), ...
        'linestyle', linestyles{li});
    handles.tab_simulate.all(end+1) = handles.plots.CatchFishGuildYearly(i);
    
    handles.plots.CatchFishGuildCumulative(i) = plot(handles.tab_simulate.axes_catch, ...
        0, 0, ...
        'marker', markers{li}, ...
        'color', linecolors(IS,:), ...
        'linewidth', linewidths(li), ...
        'linestyle', linestyles{li});
    handles.tab_simulate.all(end+1) = handles.plots.CatchFishGuildCumulative(i);
end

% Catch plots for all fish species
for i = 1:GuildInfo.nFishSpecies
    handles.plots.BiomassFishGuildPerSpecies(i) = plot(handles.tab_simulate.axes_biomass, ...
        0, 0, ...
        'marker', '.', ...
        'linestyle', '-', ...
        'color', linecolors(i,:));
    handles.tab_simulate.all(end+1) = handles.plots.BiomassFishGuildPerSpecies(i);
    
    handles.plots.CatchFishGuildYearlyPerSpecies(i) = plot(handles.tab_simulate.axes_catch, ...
        0, 0, ...
        'marker', '.', ...
        'linestyle', '-', ...
        'color', linecolors(i,:));
    handles.tab_simulate.all(end+1) = handles.plots.CatchFishGuildYearlyPerSpecies(i);
    
    handles.plots.CatchFishGuildCumulativePerSpecies(i) = plot(handles.tab_simulate.axes_catch, ...
        0, 0, ...
        'marker', '.', ...
        'linestyle', '-', ...
        'color', linecolors(i,:));
    handles.tab_simulate.all(end+1) = handles.plots.CatchFishGuildCumulativePerSpecies(i);
end

% Show/hide plots based on the biomasstype popumenu value
switch handles.tab_simulate.popupmenu_biomasstype.Value
    case 1
        set(handles.plots.BiomassFishGuild,'Visible','on')
        set(handles.plots.BiomassFishGuildPerSpecies,'Visible','off')
    case 2
        set(handles.plots.BiomassFishGuild,'Visible','off')
        set(handles.plots.BiomassFishGuildPerSpecies,'Visible','on')
end


% Show/hide plots based on the catchtype popumenu value
switch handles.tab_simulate.popupmenu_catchtype.Value
    case 1
        set(handles.plots.CatchFishGuildYearly,'Visible','on')
        set(handles.plots.CatchFishGuildCumulative,'Visible','off')
        set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','off')
        set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','off')
    case 2
        set(handles.plots.CatchFishGuildYearly,'Visible','off')
        set(handles.plots.CatchFishGuildCumulative,'Visible','on')
        set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','off')
        set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','off')
    case 3
        set(handles.plots.CatchFishGuildYearly,'Visible','off')
        set(handles.plots.CatchFishGuildCumulative,'Visible','off')
        set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','on')
        set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','off')
    case 4
        set(handles.plots.CatchFishGuildYearly,'Visible','off')
        set(handles.plots.CatchFishGuildCumulative,'Visible','off')
        set(handles.plots.CatchFishGuildYearlyPerSpecies,'Visible','off')
        set(handles.plots.CatchFishGuildCumulativePerSpecies,'Visible','on')
end



function handles = refreshAxes(handles)

GuildInfo = handles.Data.GuildInfo;
Guilds = handles.Data.Guilds;
for i = 1:GuildInfo.nFishGuilds
    h = handles.plots.BiomassFishGuild(i);
    IF = GuildInfo.iFishGuilds(i);
    xdata = 0:size(handles.Results.B,2);
    ydata = [Guilds(IF).binit handles.Results.B(IF,:)];
    set(h,'xdata',xdata,'ydata',ydata);
end
fishgls = {Guilds(GuildInfo.iFishGuilds).label};
S.subs = {1,1:3};
S.type = '()';
fslabels = unique(cellfun(@(x)subsref(x,S),fishgls,'uniformoutput',false),'stable');
fishgls = cellfun(@(x)subsref(x,S),fishgls,'uniformoutput',false);
afishgls = {Guilds(GuildInfo.iAdultFishGuilds).label};
afishgls = cellfun(@(x)subsref(x,S),afishgls,'uniformoutput',false);
for i = 1:GuildInfo.nFishSpecies
    afsinds = find(strcmpi(fslabels{i},afishgls));
    fsinds = find(strcmpi(fslabels{i},fishgls));
    h1 = handles.plots.CatchFishGuildYearlyPerSpecies(i);
    h2 = handles.plots.CatchFishGuildCumulativePerSpecies(i);
    h3 = handles.plots.BiomassFishGuildPerSpecies(i);
    xdata = 0:size(handles.Results.C,2);
    ydata1 = [0 sum(handles.Results.C(afsinds,:))]; %#ok<FNDSB>
    ydata2 = cumsum(ydata1);
    ydata3 = [sum([Guilds(GuildInfo.iFishGuilds(fsinds)).binit]) sum(handles.Results.B(GuildInfo.iFishGuilds(fsinds),:))];
    set(h1,'xdata',xdata,'ydata',ydata1);
    set(h2,'xdata',xdata,'ydata',ydata2);
    set(h3,'xdata',xdata,'ydata',ydata3);
end
for i = 1:GuildInfo.nAdultFishGuilds
    h1 = handles.plots.CatchFishGuildYearly(i);
    h2 = handles.plots.CatchFishGuildCumulative(i);
    xdata = 0:size(handles.Results.C,2);
    ydata1 = [0 handles.Results.C(i,:)];
    ydata2 = cumsum(ydata1);
    set(h1,'xdata',xdata,'ydata',ydata1);
    set(h2,'xdata',xdata,'ydata',ydata2);
end
if xdata(end) < handles.minaxeswidth
    handles.tab_simulate.axes_biomass.XLim = [0 handles.minaxeswidth];
    handles.tab_simulate.axes_catch.XLim = [0 handles.minaxeswidth];
else
    axeswidth = str2double(handles.tab_simulate.edit_axeswidth.String); % jne..
    handles.tab_simulate.axes_biomass.XLim = [max([0 xdata(end)-axeswidth]) xdata(end)];
    handles.tab_simulate.axes_catch.XLim = [max([0 xdata(end)-axeswidth]) xdata(end)];
end
drawnow
function handles = setUIControlPositions(handles)

tab_simulate_positions = [ ...
    0 0 204.2 47.7692;          % (uipanel_tab_simulate)
    27.4 17.5 80.4 26.3077;     % (axes_biomass)
    117.6 17.5 80.4 26.3077;    % (axes_catch)
    45.8 3.0769 24.2 6.2308;    % (uibuttongroup_gear)
    0.8 3.3077 21.2 1.7692;     % (radiobutton_gear_1)
    0.8 1.9231 20.6 1.7692;     % (radiobutton_gear_2)
    0.8 0.4615 16.2 1.7692;     % (radiobutton_gear_3)
    71.6 3.0769 66.4 6.2308;    % (uipanel_fishingpressure)
    18.4 3.0769 6.4 1.7;        % (edit_fishingmortality)
    1.1 3.3846 16 1.08;         % (text_fishingmortality)
    139.8 2.7 44 8.3;           % (uipanel_simulate)
    31.4 5.07 10.4 1.7;         % (edit_numberofyears)
    2.2 5.4 28.2 1.08;          % (text_numberofyears)
    2 2.6 20 1.7;               % (pushbutton_gofishing)
    21.8 2.6 20 1.7;            % (pushbutton_nofishing)
    58 45.3 0 NaN;              % (text_axes_biomass)
    151 45.3 0 NaN;             % (text_axes_catch)
    141.8 1 40.2 1.08;          % (text_simulationstatus)
    186 5.9 13.8 2.462;         % (pushbutton_reset)
    170 13 30 1.6154;           % (popupmenu_catchtype)
    27.8 13 32.2 1.08;          % (text_axeswidth)
    61.6 13 6.2 1.7;            % (edit_axeswidth)
    1.8 10.4615 18.2 14.6154;   % (uipanel_biomassplots)
    186 3.354 13.8 2.462;       % (pushbutton_save_results)
    1.1 1.231 39.2 1.08;        % (text_releaseproportion)
    40.8 0.846 6.4 1.7;         % (edit_releaseproportion)
    54 1 17.4 1.7;              % (pushbutton_planktonfigure)
    0 0 1 1;                    % (axes_bgimage_simulate)
    1.8 7.5 40 1.6              % (popupmenu_ecosystem)
    2 0.8 20 1.7;               % (pushbutton_rewind)
    21.8 0.8 20 1.7;            % (pushbutton_stop)
    72 1 22.8 1.7;              % (pushbutton_totalbiomassfigure)
    1 1 36 1.615;               % (text_balance)
    78 13 30 1.6154;            % (popupmenu_biomasstype)
    1 3 17.4 1.8;               % (checkbox_gamemode)
    26 3.2 36.4 1.3;            % (slider_F)
    ];

tab_guilds_positions = [ ...
    0 0 204.2 47.7692;          % (uipanel_tab_guilds)
    7.8 11.3 58.2 33;           % (uipanel_add_guild)
    1.8 31 6.2 1.08;            % (text_name)
    1.8 27.7 6.2 1.08;          % (text_label)
    1.8 17.4 20.2 1.08;         % (text_intrinsic_growth_rate)
    1.8 14.1 20.2 1.08;         % (text_metabolic_rate)
    1.6 29 34.4 1.7;            % (edit_name)
    1.8 25.7 34.2 1.7;          % (edit_label)
    1.8 15.4 10.2 1.7;          % (edit_intrinsic_growth_rate)
    1.8 12.1 10.2 1.7;          % (edit_metabolic_rate)
    1.8 0.3 13.8 1.7;           % (pushbutton_add_guild)
    16 0.3 17.2 1.7;            % (pushbutton_save_changes)
    34 0.3 13.8 1.7;            % (pushbutton_clear)
    1.8 18.385 20.2 7.154;      % (uibuttongroup_type)
    67.6 15.3077 134.4 29.3;    % (uipanel_guild_list)
    2 16.8462 33.8 11.23;       % (listbox_guild_list)
    2 12.1538 34 1.7;           % (pushbutton_delete)
    136 13.6154 13.8 1.7;       % (pushbutton_import)
    154 13.6154 13.8 1.7;       % (pushbutton_export)
    1.8 2.923 14.4 1.7692;      % (radiobutton_consumer)
    1.8 1.615 13.2 1.7692;      % (radiobutton_producer)
    1.8 0.308 13.2 1.7692;      % (radiobutton_detritus)
    38 16.8462 33.8 11.23;      % (listbox_prey)
    98 16.8462 33.8 11.23;      % (listbox_available_prey)
    72.6 23 24.2 1.7;           % (pushbutton_add_food_item)
    72.8 20.5385 24.2 1.7;      % (pushbutton_remove_food_item)
    118 13.6154 13.8 1.7;       % (pushbutton_clearall)
    1.8 4.23 15.6 1.7692;       % (radiobutton_fish)
    23.4 14.1 30.2 1.08;        % (text_average_length)
    23.4 12.1 10.2 1.7;         % (edit_average_length)
    23.4 10.8 26.8 1.08;        % (text_length_weight_parameters)
    27 8.8 10.2 1.7;          	% (edit_lw_a)
    43 8.8 10.2 1.7;            % (edit_lw_b)
    24 9 2 1.08;                % (text_lw_a)
    40 9 2 1.08;                % (text_lw_b)
    23.4 22 10.2 1.7;           % (edit_initial_biomass)
    23.4 24 14 1.08;            % (text_initial_biomass)
    2 14.5 17 1.7;              % (pushbutton_moveup)
    19 14.5 17 1.7;             % (pushbutton_movedown)
    40 15 24 1.08;              % (text_half_saturation_constant)
    40 12.85 10.2 1.7;          % (edit_half_saturation_constant)
    40 11.15 32.2 1.08;         % (text_feeding_interference_coefficient)
    40 9 10.2 1.7;              % (edit_feeding_interference_coefficient)
    0 0 1 1;                    % (axes_bgimage_guilds)
    7.6 7.5 10.2 1.7;           % (edit_length_season)
    7.6 9.5 30.6 1.08;          % (text_length_season)
    23.4 18.7 10.2 1.7;         % (edit_basal_respiration)
    23.4 20.7 28.2 1.08;        % (text_basal_respiration)
    23.4 15.4 10.2 1.7          % (edit_producer_competition)
    23.4 17.4 32.4 1.08;        % (text_producer_competition)
    7.8 1 77.4 5.85;            % (uibuttongroup_carrying_capacity)
    1.5 2.7 14.4 1.7692;        % (radiobutton_K_constant)
    19.4 2.7 17 1.7692;         % (radiobutton_K_white_noise)
    49.6 2.7 11.4 1.7692;       % (radiobutton_K_AR1)
    1.5 1 6 1.08;               % (text_K_mean)
    8.4 0.7 10.2 1.7;           % (edit_K_mean)
    19.4 1 18.8 1.08;           % (text_K_standard_deviation)
    38.4 0.7 10.2 1.7;        	% (edit_K_standard_deviation)
    49.6 1 15.4 1.08;           % (text_K_autocorrelation)
    65.4 0.7 10.2 1.7;          % (edit_K_autocorrelation)
    40 7.3 33.2 1.08;           % (text_F_q)
    40 5.15 10.2 1.7;           % (edit_F_q)
    75 15 22.6 1.08;            % (text_assimilation_efficiency)
    75 12.85 10.2 1.7;          % (edit_assimilation_efficiency)
    75 11.15 27.2 1.08;         % (text_maximum_consumption_rate)
    75 9 10.2 1.7;              % (edit_maximum_consumption_rate)
    1.8 10.8 21.4 1.08;         % (text_fraction_of_exudation)
    1.8 8.8 10.2 1.7;          	% (edit_fraction_of_exudation)
    1.8 7.5 15.4 1.08;          % (text_dissolution_rate)
    1.8 5.5 10.2 1.7;          	% (edit_dissolution_rate)
    1.8 4.2 9.2 1.08;           % (text_hatchery)
    1.8 2.2 10.2 1.7;          	% (edit_hatchery)
    23.4 4.2 32.8 1.08;         % (text_invest)
    23.4 2.2 10.2 1.7;          % (edit_invest)
    23.4 7.5 30.2 1.08;         % (text_activity_respiration_coefficient)
    23.4 5.5 10.2 1.7;          % (edit_activity_respiration_coefficient)
    40 22 13.5 1.4;             % (checkbox_catchable)
    ];


for i = 1:length(handles.tab_simulate.all)
    if strcmpi(handles.tab_simulate.all(i).Type,'text')
        handles.tab_simulate.all(i).Position = tab_simulate_positions(i,1:3);
    else
        handles.tab_simulate.all(i).Position = tab_simulate_positions(i,:);
    end
end

for i = 1:length(handles.tab_guilds.all)
    handles.tab_guilds.all(i).Position = tab_guilds_positions(i,:);
end
function handles = initializeGuildTab(handles)
handles = updateGuildListbox(handles, 1);
handles = updatePreyListbox(handles);
handles = updateAvailablePreyListbox(handles);
handles = setGuildEditStatus(handles);
handles = setParameterEditValues(handles);
function handles = updateGuildListbox(handles, value)
guild_labels = {handles.Data.Guilds.label}';
handles.tab_guilds.listbox_guild_list.String = guild_labels;
handles.tab_guilds.listbox_guild_list.Value = value;
function handles = updatePreyListbox(handles)
glabels = handles.tab_guilds.listbox_guild_list.String;
gvalue = handles.tab_guilds.listbox_guild_list.Value;
if ~isempty(glabels)
    selectedGuildLabel = glabels{gvalue}; % TODO: FIX ERROR IF LIST GETS EMPTY
    gind = find(strcmp(selectedGuildLabel,{handles.Data.Guilds.label}));
    
    % get feeding links for selected guild
    preyinds = find(handles.Data.communityMatrix(gind,:)); %#ok<FNDSB>
    handles.tab_guilds.listbox_prey.String = glabels(preyinds); %#ok<FNDSB>
    handles.tab_guilds.listbox_prey.Value = 1;
else
    handles.tab_guilds.listbox_prey.String = [];
    handles.tab_guilds.listbox_prey.Value = 1;
end
function handles = updateAvailablePreyListbox(handles)
preyitems = handles.tab_guilds.listbox_prey.String;
availablepreyitems = setdiff({handles.Data.Guilds.label}',preyitems);
handles.tab_guilds.listbox_available_prey.String = availablepreyitems;
handles.tab_guilds.listbox_available_prey.Value = 1;
function handles = updateEdits(handles)
ii = handles.tab_guilds.listbox_guild_list.Value;
if ~isempty(handles.tab_guilds.listbox_guild_list.String)
    glabel = handles.tab_guilds.listbox_guild_list.String{ii};
    i = find(strcmp({handles.Data.Guilds.label},glabel));
    g = handles.Data.Guilds(i); %#ok<FNDSB>
    handles.tab_guilds.edit_name.String = g.name;
    handles.tab_guilds.edit_label.String = g.label;
    handles.tab_guilds.edit_intrinsic_growth_rate.String = num2str(g.igr);
    handles.tab_guilds.edit_metabolic_rate.String = num2str(g.mbr);
    handles.tab_guilds.edit_average_length.String = num2str(g.avgl);
    handles.tab_guilds.edit_lw_a.String = num2str(g.lw_a);
    handles.tab_guilds.edit_lw_b.String = num2str(g.lw_b);
    handles.tab_guilds.edit_initial_biomass.String = num2str(g.binit);
    handles.tab_guilds.edit_basal_respiration.String = num2str(g.bx);
    handles.tab_guilds.edit_producer_competition.String = num2str(g.c);
    handles.tab_guilds.edit_fraction_of_exudation.String = num2str(g.s);
    handles.tab_guilds.edit_dissolution_rate.String = num2str(g.diss_rate);
    handles.tab_guilds.edit_hatchery.String = num2str(g.hatchery);
    handles.tab_guilds.edit_activity_respiration_coefficient.String = num2str(g.ax);
    handles.tab_guilds.edit_invest.String = num2str(g.invest);
    gcatchable = 0;
    if isfield(g,'catchable') && ~isempty(g.catchable)
    	gcatchable = g.catchable;
    end
    handles.tab_guilds.checkbox_catchable.Value = gcatchable;
    set(findobj(handles.tab_guilds.uibuttongroup_type.Children, ...
        'String',g.type),'Value',1)
else
    handles.tab_guilds.edit_name.String = '';
    handles.tab_guilds.edit_label.String = '';
    handles.tab_guilds.edit_intrinsic_growth_rate.String = '';
    handles.tab_guilds.edit_metabolic_rate.String = '';
    handles.tab_guilds.edit_average_length.String = '';
    handles.tab_guilds.edit_lw_a.String = '';
    handles.tab_guilds.edit_lw_b.String = '';
    handles.tab_guilds.edit_initial_biomass.String = '';
    handles.tab_guilds.edit_basal_respiration.String = '';
    handles.tab_guilds.edit_producer_competition.String = '';
    handles.tab_guilds.edit_fraction_of_exudation.String = '';
    handles.tab_guilds.edit_dissolution_rate.String = '';
    handles.tab_guilds.edit_hatchery.String = '';
    handles.tab_guilds.edit_activity_respiration_coefficient.String = '';
    handles.tab_guilds.edit_invest.String = '';
    set(handles.tab_guilds.uibuttongroup_type.Children(1),'Value',1)
    
    handles.tab_guilds.edit_half_saturation_constant.String = '';
    handles.tab_guilds.edit_feeding_interference_coefficient.String = '';
    handles.tab_guilds.edit_F_q.String = '';
    handles.tab_guilds.edit_assimilation_efficiency.String = '';
    handles.tab_guilds.edit_maximum_consumption_rate.String = '';
end
function handles = setGuildEditStatus(handles)

switch handles.tab_guilds.uibuttongroup_type.SelectedObject.Tag
    case 'radiobutton_fish'
        set(handles.tab_guilds.edit_intrinsic_growth_rate,'Enable','off')
        set(handles.tab_guilds.edit_metabolic_rate,'Enable','off')
        set(handles.tab_guilds.edit_basal_respiration,'Enable','on')
        set(handles.tab_guilds.edit_producer_competition,'Enable','off')
        set(handles.tab_guilds.edit_average_length,'Enable','on')
        set(handles.tab_guilds.edit_lw_a,'Enable','on')
        set(handles.tab_guilds.edit_lw_b,'Enable','on')
        set(handles.tab_guilds.edit_fraction_of_exudation,'Enable','off')
        set(handles.tab_guilds.edit_dissolution_rate,'Enable','off')
        set(handles.tab_guilds.edit_hatchery,'Enable','on')
        set(handles.tab_guilds.edit_activity_respiration_coefficient,'Enable','on')
        set(handles.tab_guilds.edit_invest,'Enable','on')
    case 'radiobutton_consumer'
        set(handles.tab_guilds.edit_intrinsic_growth_rate,'Enable','off')
        set(handles.tab_guilds.edit_metabolic_rate,'Enable','on')
        set(handles.tab_guilds.edit_basal_respiration,'Enable','on')
        set(handles.tab_guilds.edit_producer_competition,'Enable','off')
        set(handles.tab_guilds.edit_average_length,'Enable','off')
        set(handles.tab_guilds.edit_lw_a,'Enable','off')
        set(handles.tab_guilds.edit_lw_b,'Enable','off')
        set(handles.tab_guilds.edit_fraction_of_exudation,'Enable','off')
        set(handles.tab_guilds.edit_dissolution_rate,'Enable','off')
        set(handles.tab_guilds.edit_hatchery,'Enable','off')
        set(handles.tab_guilds.edit_activity_respiration_coefficient,'Enable','on')
        set(handles.tab_guilds.edit_invest,'Enable','off')
    case 'radiobutton_producer'
        set(handles.tab_guilds.edit_intrinsic_growth_rate,'Enable','on')
        set(handles.tab_guilds.edit_metabolic_rate,'Enable','off')
        set(handles.tab_guilds.edit_basal_respiration,'Enable','off')
        set(handles.tab_guilds.edit_producer_competition,'Enable','on')
        set(handles.tab_guilds.edit_average_length,'Enable','off')
        set(handles.tab_guilds.edit_lw_a,'Enable','off')
        set(handles.tab_guilds.edit_lw_b,'Enable','off')
        set(handles.tab_guilds.edit_fraction_of_exudation,'Enable','on')
        set(handles.tab_guilds.edit_dissolution_rate,'Enable','off')
        set(handles.tab_guilds.edit_hatchery,'Enable','off')
        set(handles.tab_guilds.edit_activity_respiration_coefficient,'Enable','off')
        set(handles.tab_guilds.edit_invest,'Enable','off')
    case 'radiobutton_detritus'
        set(handles.tab_guilds.edit_intrinsic_growth_rate,'Enable','off')
        set(handles.tab_guilds.edit_metabolic_rate,'Enable','off')
        set(handles.tab_guilds.edit_basal_respiration,'Enable','off')
        set(handles.tab_guilds.edit_producer_competition,'Enable','off')
        set(handles.tab_guilds.edit_average_length,'Enable','off')
        set(handles.tab_guilds.edit_lw_a,'Enable','off')
        set(handles.tab_guilds.edit_lw_b,'Enable','off')
        set(handles.tab_guilds.edit_fraction_of_exudation,'Enable','off')
        set(handles.tab_guilds.edit_dissolution_rate,'Enable','on')
        set(handles.tab_guilds.edit_hatchery,'Enable','off')
        set(handles.tab_guilds.edit_activity_respiration_coefficient,'Enable','off')
        set(handles.tab_guilds.edit_invest,'Enable','off')
end
function handles = swapGuildPositions(handles,i1,i2)
g1 = handles.Data.Guilds(i1);
g2 = handles.Data.Guilds(i2);

Cmtx = handles.Data.communityMatrix;
dmtx = handles.Data.d;
B0mtx = handles.Data.B0;
qmtx = handles.Data.q;
emtx = handles.Data.e;
ymtx = handles.Data.y;

Cmtx_new = Cmtx;
dmtx_new = dmtx;
B0mtx_new = B0mtx;
qmtx_new = qmtx;
emtx_new = emtx;
ymtx_new = ymtx;

Cmtx_new(i2,:) = Cmtx(i1,:);
Cmtx_new(i1,:) = Cmtx(i2,:);
Cmtx_new(:,i2) = Cmtx(:,i1);
Cmtx_new(:,i1) = Cmtx(:,i2);
dmtx_new(i2,:) = dmtx(i1,:);
dmtx_new(i1,:) = dmtx(i2,:);
dmtx_new(:,i2) = dmtx(:,i1);
dmtx_new(:,i1) = dmtx(:,i2);
B0mtx_new(i2,:) = B0mtx(i1,:);
B0mtx_new(i1,:) = B0mtx(i2,:);
B0mtx_new(:,i2) = B0mtx(:,i1);
B0mtx_new(:,i1) = B0mtx(:,i2);
qmtx_new(i2,:) = qmtx(i1,:);
qmtx_new(i1,:) = qmtx(i2,:);
qmtx_new(:,i2) = qmtx(:,i1);
qmtx_new(:,i1) = qmtx(:,i2);
emtx_new(i2,:) = emtx(i1,:);
emtx_new(i1,:) = emtx(i2,:);
emtx_new(:,i2) = emtx(:,i1);
emtx_new(:,i1) = emtx(:,i2);
ymtx_new(i2,:) = ymtx(i1,:);
ymtx_new(i1,:) = ymtx(i2,:);
ymtx_new(:,i2) = ymtx(:,i1);
ymtx_new(:,i1) = ymtx(:,i2);

if Cmtx(i1,i2) && ~Cmtx(i2,i1)
    Cmtx_new(i1,i1) = 0;
    Cmtx_new(i2,i1) = 1;
    dmtx_new(i1,i1) = 0;
    dmtx_new(i2,i1) = dmtx(i1,i2);
    B0mtx_new(i1,i1) = 0;
    B0mtx_new(i2,i1) = B0mtx(i1,i2);
    qmtx_new(i1,i1) = 0;
    qmtx_new(i2,i1) = qmtx(i1,i2);
    emtx_new(i1,i1) = 0;
    emtx_new(i2,i1) = emtx(i1,i2);
    ymtx_new(i1,i1) = 0;
    ymtx_new(i2,i1) = ymtx(i1,i2);
elseif ~Cmtx(i1,i2) && Cmtx(i2,i1)
    Cmtx_new(i1,i2) = 1;
    Cmtx_new(i2,i2) = 0;
    dmtx_new(i1,i2) = dmtx(i2,i1);
    dmtx_new(i2,i2) = 0;
    B0mtx_new(i1,i2) = B0mtx(i2,i1);
    B0mtx_new(i2,i2) = 0;
    qmtx_new(i1,i2) = qmtx(i2,i1);
    qmtx_new(i2,i2) = 0;
    emtx_new(i1,i2) = emtx(i2,i1);
    emtx_new(i2,i2) = 0;
    ymtx_new(i1,i2) = ymtx(i2,i1);
    ymtx_new(i2,i2) = 0;
elseif Cmtx(i1,i2) && Cmtx(i2,i1)
    Cmtx_new(i1,i1) = 0;
    Cmtx_new(i2,i1) = 1;
    Cmtx_new(i2,i2) = 0;
    Cmtx_new(i1,i2) = 1;
    dmtx_new(i1,i1) = 0;
    dmtx_new(i2,i1) = dmtx(i1,i2);
    dmtx_new(i2,i2) = 0;
    dmtx_new(i1,i2) = dmtx(i2,i1);
    B0mtx_new(i1,i1) = 0;
    B0mtx_new(i2,i1) = B0mtx(i1,i2);
    B0mtx_new(i2,i2) = 0;
    B0mtx_new(i1,i2) = B0mtx(i2,i1);
    qmtx_new(i1,i1) = 0;
    qmtx_new(i2,i1) = qmtx(i1,i2);
    qmtx_new(i2,i2) = 0;
    qmtx_new(i1,i2) = qmtx(i2,i1);
    emtx_new(i1,i1) = 0;
    emtx_new(i2,i1) = emtx(i1,i2);
    emtx_new(i2,i2) = 0;
    emtx_new(i1,i2) = emtx(i2,i1);
    ymtx_new(i1,i1) = 0;
    ymtx_new(i2,i1) = ymtx(i1,i2);
    ymtx_new(i2,i2) = 0;
    ymtx_new(i1,i2) = ymtx(i2,i1);
end

handles.Data.communityMatrix = Cmtx_new;
handles.Data.d = dmtx_new;
handles.Data.B0 = B0mtx_new;
handles.Data.q = qmtx_new;
handles.Data.e = emtx_new;
handles.Data.y = ymtx_new;
handles.Data.Guilds(i1) = g2;
handles.Data.Guilds(i2) = g1;

handles.Data = updateGuildInfo(handles.Data);
handles.Data = addDerivedQuantities(handles.Data);
handles = updateGuildListbox(handles,i2);

function handles = saveParameterValues(handles)
handles.Data.K.type = handles.tab_guilds.uibuttongroup_carrying_capacity.SelectedObject.String;
switch handles.Data.K.type
    case 'Constant'
        handles.Data.K.mean = str2double(handles.tab_guilds.edit_K_mean.String);
        handles.Data.K.standard_deviation = [];
        handles.Data.K.autocorrelation = [];
    case 'White noise'
        handles.Data.K.mean = str2double(handles.tab_guilds.edit_K_mean.String);
        handles.Data.K.standard_deviation = str2double(handles.tab_guilds.edit_K_standard_deviation.String);
        handles.Data.K.autocorrelation = [];
    case 'AR(1)'
        handles.Data.K.mean = str2double(handles.tab_guilds.edit_K_mean.String);
        handles.Data.K.standard_deviation = str2double(handles.tab_guilds.edit_K_standard_deviation.String);
        handles.Data.K.autocorrelation = str2double(handles.tab_guilds.edit_K_autocorrelation.String);
end
function handles = setParameterEditValues(handles)
handles.tab_guilds.edit_length_season.String = ...
    num2str(handles.Data.nGrowthDays);
switch handles.Data.K.type
    case 'Constant'
        handles.tab_guilds.radiobutton_K_constant.Value = 1;
        handles.tab_guilds.edit_K_mean.String = ...
            num2str(handles.Data.K.mean);
        uibuttongroup_carrying_capacity_SelectionChangedFcn( ...
            handles.tab_guilds.radiobutton_K_constant, [], handles)
    case 'White noise'
        handles.tab_guilds.radiobutton_K_white_noise.Value = 1;
        handles.tab_guilds.edit_K_mean.String = ...
            num2str(handles.Data.K.mean);
        handles.tab_guilds.edit_K_standard_deviation.String = ...
            num2str(handles.Data.K.standard_deviation);
        uibuttongroup_carrying_capacity_SelectionChangedFcn( ...
            handles.tab_guilds.radiobutton_K_white_noise, [], handles)
    case 'AR(1)'
        handles.tab_guilds.radiobutton_K_AR1.Value = 1;
        handles.tab_guilds.edit_K_mean.String = ...
            num2str(handles.Data.K.mean);
        handles.tab_guilds.edit_K_standard_deviation.String = ...
            num2str(handles.Data.K.standard_deviation);
        handles.tab_guilds.edit_K_autocorrelation.String = ...
            num2str(handles.Data.K.autocorrelation);
        uibuttongroup_carrying_capacity_SelectionChangedFcn( ...
            handles.tab_guilds.radiobutton_K_AR1, [], handles)
end

function handles = updateEcosystemPopupmenu(handles,str_sel)
d = dir(fullfile('Data','*.mat'));
ecoStr = cellfun(@(x)x(1:strfind(x,'.mat')-1),{d.name},'uniformoutput',false);
set(handles.tab_simulate.popupmenu_ecosystem,'string',ecoStr)
value = find(strcmpi(ecoStr,str_sel));
if isempty(value) && ~isempty(ecoStr)
    value = 1;
elseif isempty(value) && isempty(ecoStr)
    value = 0;
end
set(handles.tab_simulate.popupmenu_ecosystem,'value',value);

function handles = addTextObjects(handles)

handles.text_axes_biomass = ...
    text(200,10,'Biomasses','parent',handles.axes_bgimage_simulate, ...
    'color','black','fontweight','bold','fontsize',14);
handles.text_axes_catch = ...
    text(530,10,'Catches','parent',handles.axes_bgimage_simulate, ...
    'color','black','fontweight','bold','fontsize',14);

function game = initializeGameStruct()

game.gameover = 0;
game.currentyear = 0;
game.balance = 0;
game.history.balance = 0;
game.history.gear = NaN;
game.history.F = NaN;
game.history.penalty = struct('type',NaN,'details',NaN,'amount',NaN);
game.history.penalty = game.history.penalty([]);
game.history.profit = NaN;
game.history.prop_released = NaN;
