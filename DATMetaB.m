%% DATMetaB v1.0
% DATMetaB Diet Assessment Tool using MetaBarcoding
%   DATMetaB is a pipeline for adjudicating classification results from the
%   Kraken metagenomics software
%   
% USAGE: DATMetaB(varargin)
%
%   varargin            Paired inputs listed below (default)
%     PavianData        use pavian csv files        (true)
%     KrakenData        use kraken report files    (false)
%     IncludeChordata   include chordata           (false)
%     verbose           Verbose output              (true)
%     debug             run code in debug mode     (false)
%     Experiment        name of the analysis  "Stomach_Content"
%     PavianPath        directory of Pavian data
%     filterData        consider Metazoa and Viridiplantae only (false)
%     OutputDir         Output directory
%     OutputOverwrite   overwrite the data in the output directory (false)
%     SaveData          save output data             (true)
%     LoadData          load data from an old analysis (false)
%     LoadDataDir       directory of old data

% 
% EXAMPLE USAGE
%   DATMetaB('PavianPath','/home/PavianPath/',...
%                 'OutputDir','/home/OutputDir/',...
%                 'LoadData',false,'PavianData',true,'filterData',true)
% 
% NOTES
%   * run with Kraken Data is not yet implemented
% 
% See DATMetaB.m file for function list, changelog, and license.
% 
% Copyright (C) 2022 Rahma Amen, Potsdam University, Potsdam, Germany

%   
% CHANGELOG
% v1.00     11NOV2021   beta version
% 
% DISCLAIMERS
%   Opinions, interpretations, conclusions, and recommendations are those
%   of the author and are not necessarily endorsed by Potsdam University.
% 
% 
% DISTRIBUTION
%   Approved for public release; distribution is unlimited.
% 

%% Entry function
% Called by user

% DATMetaB
function DATMetaB(varargin)
% DATMetaB Diet Assessment Tool using MetaBarcoding
% 
% See documentation at beginning of DATMetaB.m file or by typing "help
% DATMetaB" in the MATLAB command prompt
% 
% Copyright (C) 2021 Rahma Amen, Potsdam University, Potsdam, Germany

% Parse varargin to options structure
% Input:  (1) arguments into this function (passed without manipulation)
% Output: (1) structure with assigned and default options
options = processArguments(varargin);

% Get list of files and samples
% Input:  (1) structure with assigned and default options
% Output: (1) options with PAVIAN files and sample names added
options = discoverSamples(options);

% Display info
dispMsg(options,sprintf('Input: %s\n',options.PavianPath));
dispMsg(options,sprintf('Output: %s\n',options.OutputDir));
if ~options.LoadData
    dispMsg(options,sprintf('Samples: %d\n',options.inputNumSamples));
end
fprintf('%s\n',repmat('*',1,100));

% Load/Initiate data structure
if options.LoadData
    f = char(fullfile(options.LoadDataDir, [options.Experiment '.mat']));
    load(f);
    options.inputNumSamples = length(fieldnames(DATA.SAMPLES));
    options.inputSamples = fieldnames(DATA.SAMPLES);
else
    DATA = struct;
end

% Get taxon table
if options.PavianData
    for i = 1:options.inputNumSamples
        % Get Pavian taxon table
        T = readPavianReport(options,i);
        % Add food taxa (class, order, genus, species)
        DATA = addFoodTaxa(DATA,T,i,options);
    end
    
elseif options.KrakenData
    handleError(1,'Input data via Kraken report is''t implemented!',options);
end 

% Calculate RRA for each species indivedually
for i = 1:options.inputNumSamples
    % calculate RRA
    DATA = calculateRRA(DATA,i,options,0.002);
end

% Calculate RRA for all species and save it in one array
DATA = createTableRRA(DATA,options);

% Save data file
if options.SaveData
    f = char(fullfile(options.OutputDir, [options.Experiment '.mat']));
    save(f, 'DATA')
end

end

%% Process functions
% Called by entry function

% processArguments
function options = processArguments(args)
% PROCESSARGUMENTS Parse varargin from DATMetaB main function into options
%
% USAGE: options = processArguments(args)
%   options     Structure of assigned and default options
%   args        Paired inputs listed below (default)
%   brackenLoc  Path containing Bracken scripts
%
% Copyright (C) 2022 Rahma Amen, Potsdam University, Potsdam, Germany

% Set internal options
options = struct;
options.programName = 'DATMetaB';
options.tic = tic; % start timer
options.timestamp = strrep(strrep(char(datetime),':','-'),' ','_');

% Parse valid arguments into options structure
if rem(length(args),2) == 1
    handleError(1,'Unmatched argument detected!',options);
end

valid_args = {'debug'...
              'Experiment'...
              'PavianData'... 
              'PavianPath'... 
              'KrakenData'... 
              'IncludeChordata'... 
              'filterData'... 
              'OutputDir'...
              'OutputOverwrite'...
              'SaveData'...
              'LoadData'...
              'LoadDataDir'...
              'verbose'}; 

for i = 1:2:length(args)
    if ~ischar(args{i})
        handleError(1,'At least one argument is not a string!',options);
    end
    if ~ismember(args{i},valid_args)
        handleError(1,sprintf('"%s" is not a valid argument!',args{i}),options);
    end
    options.(args{i}) = args{i+1};
end
    
% Validate provided options and default remaining

% verbose: Verbose output (true)
opt = 'verbose';
if ~isfield(options,opt)
    options.(opt) = true;
elseif ~islogical(options.(opt))
    error('[FATLERR] %s must be logical!',opt);
end

% PavianData: Use pavian csv files
opt = 'PavianData';
if ~isfield(options,opt)
    options.(opt) = true;
elseif ~islogical(options.(opt))
    handleError(1,sprintf('%s must be logical',opt),options);
end

% KrakenData: Use kraken report files (false)
opt = 'KrakenData';
if ~isfield(options,opt)
    options.(opt) = false;
elseif ~islogical(options.(opt))
    handleError(1,sprintf('%s must be logical',opt),options);
end

% PavianPath: Path containing Pavian data
opt = 'PavianPath';
if ~isfield(options,opt) && options.PavianData
    handleError(1,sprintf('%s is required!',opt),options);
elseif ~isfield(options,opt) && ~options.PavianData
    options.(opt) = '';
elseif isfield(options,opt) && options.PavianData
    if ~ischar(options.(opt))
        handleError(1,sprintf('%s must be a string!',opt),options);
    elseif ~exist(options.(opt),'dir')
        handleError(1,sprintf('%s not found!',opt),options);
    end
end

% outputDir: Output directory (inputDir)
opt = 'OutputDir';
if ~isfield(options,opt)
    handleError(1,sprintf('%s is required!',opt),options);
elseif ~exist(options.(opt),'dir')
        handleError(1,sprintf('%s not found!',opt),options);
end

% IncludeChordata: IncludeChordata (false)
opt = 'IncludeChordata';
if ~isfield(options,opt)
    options.(opt) = true;
elseif ~islogical(options.(opt))
    handleError(1,sprintf('%s must be logical',opt),options);
end

% filterData: filter all the data except Metazoa and Viridiplantae (false)
opt = 'filterData';
if ~isfield(options,opt)
    options.(opt) = false;
elseif ~islogical(options.(opt))
    handleError(1,sprintf('%s must be logical',opt),options);
end

% debug: Allow soft errors to manifest (false)
opt = 'debug';
if ~isfield(options,opt)
    options.(opt) = false;
elseif ~islogical(options.(opt))
    handleError(1,sprintf('%s must be logical',opt),options);
end

% outputOverwrite: Overwrite existing files (false)
opt = 'OutputOverwrite';
if ~isfield(options,opt)
    options.(opt) = false;
elseif ~islogical(options.(opt))
    handleError(1,sprintf('%s must be logical',opt),options);
end

% SaveData: save matlab structure data
opt = 'SaveData';
if ~isfield(options,opt)
    options.(opt) = true;
elseif ~islogical(options.(opt))
    handleError(1,sprintf('%s must be logical',opt),options);
end

% LoadData: Load matlab structure data
opt = 'LoadData';
if ~isfield(options,opt)
    options.(opt) = false;
elseif ~islogical(options.(opt))
    handleError(1,sprintf('%s must be logical',opt),options);
end
if options.PavianData && options.LoadData
    handleError(1,sprintf('%s should be false when PavianData is true!',opt),options);
end

% LoadDataDir: Load data directory
opt = 'LoadDataDir';
if ~isfield(options,opt)
    options.LoadDataDir = options.OutputDir;
elseif ~exist(options.(opt),'dir')
    handleError(1,sprintf('%s not found!',opt),options);
end

% Experiment: Experiment name
opt = 'Experiment';
if ~isfield(options,opt)
    options.(opt) = 'Stomach_Content';
elseif ~ischar(options.(opt))
    handleError(1,sprintf('%s should be character!',opt),options);
end

end % processArguments

% dispMsg
function dispMsg(options,message,source)
% DISPMSG Display provided message if verbose is on
% 
% USAGE: dispMsg(options,message)
%   options     Structure of assigned and default options
%   message     Text to display as status message (ex: 'Running...')
%   source      Source function name (KANTOO)
% 
% Copyright (C) 2022 Rahma Amen, Potsdam University, Potsdam, Germany

    % Check inputs
    if nargin < 3
        source = options.programName;
    elseif ~ischar(source)
        error('source provided to dispMsg must be string!');
    else
        source = upper(source);
    end
    if ~ischar(message)
        error('message provided to dispMsg must be string!');
    end
    if nargin == 0
        error('dispMsg requires the options argument at a minimum!');
    end
    
    % Display provided message
    if options.verbose && ~isempty(message) && ischar(message)
        fprintf('[%s] %s',source,sprintf(message));
    end
end % dispMsg

% handleError
function handleError(stat,ME,options)
% HANDLEERROR Handle errors neatly
% 
% USAGE: handleError(stat,ME,options)
%   stat        Execution status (0=no error,1=fatal error,2=soft error)
%   ME          Error message or MATLAB MExcpetion object
%   options     Structure of assigned and default options
% 
% Copyright (C) 2022 Rahma Amen, Potsdam University, Potsdam, Germany

    % Check inputs
    if nargin < 3
        options.verbose = false;
        tocString = '';
    else
        tocString = sprintf(' (%.2f min)',toc(options.tic)/60);
    end
    if nargin < 2
        ME = '[FATLERR] Undefined error!';
    end
    if nargin == 0
    	error('errorHandle requires the stat argument at a minimum!');
    end
    if ~isfield(options,'verbose')
        options.verbose = true;
    end
    
    % Handle error based on severity and verbosity
    switch stat
        case 0 % no error
            if options.verbose
                fprintf('success!%s\n',tocString);
            end
        case 1 % fatal error
            if options.verbose
                fprintf('failed!%s\n',tocString);
            end
            if isa(ME,'char')
                error('[FATLERR] %s encountered a fatal error during execution:\n%s',options.programName,ME);
            elseif isa(ME,'MException')
                error('[FATLERR] %s encountered a fatal error during execution:\n%s',options.programName,getReport(ME));
            else
                error('Undefined error!');
            end
        otherwise % 2 for warning or any other number for non-fatal error
            if options.verbose
                fprintf('warning!%s\n',tocString);
                if isa(ME,'char')
                    fprintf('%s\n',ME);
                elseif isa(ME,'MException')
                    dispMsg(options,sprintf('%s encountered the following non-fatal error:\n',options.programName),'SOFTERR');
                    fprintf('%s\n',getReport(ME));
                end
            end
    end
end % handleError

% discoverSamples
function options = discoverSamples(options)
% DISCOVERSAMPLES Gather list of FASTQ files and sample names from inputDir
% 
% USAGE: options = makeDirectories(options)
%   options     Structure of assigned and default options with additions
% 
% Copyright (C) 2016-2018 Turner Conrad
    
% Check for Pavian files
dispMsg(options,'Discovering samples...');
if options.PavianData
    fileList = dir(options.PavianPath);
    fileList = {fileList.name}';
    fileList = fileList(3:end);
    options.inputPavian = fileList;
    if isempty(options.inputPavian)
        handleError(1,'No files with inputExt found in inputDir!',options);
    end
    
    % Extract sample names from Pavian names
    options.inputSamples = regexprep(fileList,'_PAVIAN','');
    
    % Annotate number of samples
    options.inputNumSamples = length(options.inputSamples);
    
    handleError(0,'',options);
end
end % discoverSamples

% readPavianReport
function Table = readPavianReport(options,sampleIndex)
% readPavianReport read in the Pavian report and clean it
% 
% USAGE: Table = readPavianReport(options)
%   Table         table containing data for sample i
%   options       Structure of assigned and default options
%   sampleIndex   sample index
% 
% Copyright (C) 2022 Rahma Amen, Potsdam University, Potsdam, Germany
    
% file name of the sample
f = char(fullfile(options.PavianPath, options.inputPavian(sampleIndex)));
Table = readtable(f);

% Clean table
Table.Percentage = [];                                 % Percentage
if options.IncludeChordata
    Table(contains(Table.TaxLineage,'Chordata'),:) = []; % Chordata
end

% if required keep only Viridiplantae and Metazoa
if options.filterData
    Table(~contains(Table.TaxLineage,'Viridiplantae') & ...
          ~contains(Table.TaxLineage,'Metazoa'),:) = []; 
end

end % readPavianReport


% DEV STATUS: in development
function d = addFoodTaxa(d,T,sampleIndex,options)
% addFoodTaxa adds pavian table to the data structure
% 
% USAGE: d = addFoodTaxa(d,T,options)
%   d             data structure containg all the experiment data
%   T             table containing data for sample i
%   sampleIndex   sample index
%   options       Structure of assigned and default options
% 
% Copyright (C) 2022 Rahma Amen, Potsdam University, Potsdam, Germany

% Classes
ids = find(ismember(T.TaxRank,'C'));
ntc = length(ids);
C = cell(ntc,2);
for i=1:ntc
    tName = T.Name{ids(i)};
    C{i,1} = tName(3:end);
    C{i,2} = T.CladeReads(ids(i));
end

% Orders
ids = find(ismember(T.TaxRank,'O'));
nto = length(ids);
O = cell(nto,3);
for i=1:nto
    tName = T.Name{ids(i)};
    x = T.TaxLineage{ids(i)};
    x = split(x,"|");
    if any(contains(x,'c_'))
        x = x{contains(x,'c_')};
    else
        x = 'c_NA';
    end
    O{i,1} = tName(3:end);
    O{i,2} = T.CladeReads(ids(i));
    O{i,3} = x(3:end);
end
% unassigned orders
xn = zeros(1,ntc);
n = nto;
for i=1:ntc
    xn(i) = cell2mat(C(i,2))-sum(cell2mat((O(ismember(O(:,3),C(i,1)),2))));
    if xn(i) > 0
        n = n +1;
        O{n,1} = ['unassigned_' C{i,1}];
        O{n,2} = xn(i);
        O{n,3} = C{i,1};
    end
end

% Family
ids = find(ismember(T.TaxRank,'F'));
ntf = length(ids);
F = cell(ntf,4);
for i=1:ntf
    tName = T.Name{ids(i)};
    x = T.TaxLineage{ids(i)};
    x = split(x,"|");
    if any(contains(x,'c_'))
        xc = x{contains(x,'c_')};
    else
        xc = 'c_NA';
    end
    if any(contains(x,'o_'))
        xo = x{contains(x,'o_')};
    else
        xo = 'o_NA';
    end
    F{i,1} = tName(3:end);
    F{i,2} = T.CladeReads(ids(i));
    F{i,3} = xc(3:end);
    F{i,4} = xo(3:end);
end
% unassigned families
n = ntf;
for i=nto+1:length(O(:,1))
    if contains(O{i,1}, 'unassigned_')
        n = n +1;
        F{n,1} = O{i,1};
        F{n,2} = O{i,2};
        F{n,3} = O{i,3};
        F{n,4} = O{i,1};
    end
end
xn = zeros(1,nto);
for i=1:nto
    xn(i) = cell2mat(O(i,2))-sum(cell2mat((F(ismember(F(:,4),O(i,1)),2))));
    if xn(i) > 0
        n = n + 1;
        F{n,1} = ['unassigned_' O{i,1}];
        F{n,2} = xn(i);
        F{n,3} = O{i,3};
        F{n,4} = O{i,1};
    end
end

% Genus
ids = find(ismember(T.TaxRank,'G'));
ntg = length(ids);
G = cell(ntg,5);
for i=1:ntg
    tName = T.Name{ids(i)};
    x = T.TaxLineage{ids(i)};
    x = split(x,"|");
    if any(contains(x,'c_'))
        xc = x{contains(x,'c_')};
    else
        xc = 'c_NA';
    end
    if any(contains(x,'o_'))
        xo = x{contains(x,'o_')};
    else
        xo = 'o_NA';
    end
    if any(contains(x,'f_'))
        xf = x{contains(x,'f_')};
    else
        xf = 'f_NA';
    end
    G{i,1} = tName(3:end);
    G{i,2} = T.CladeReads(ids(i));
    G{i,3} = xc(3:end);
    G{i,4} = xo(3:end);
    G{i,5} = xf(3:end);
end
% unassigned genus
n = ntg;
for i=ntf+1:length(F(:,1))
    if contains(F{i,1}, 'unassigned_')
        n = n +1;
        G{n,1} = F{i,1};
        G{n,2} = F{i,2};
        G{n,3} = F{i,3};
        G{n,4} = F{i,1};
        G{n,5} = F{i,1};
    end
end

xn = zeros(1,ntf);
for i=1:ntf
    xn(i) = cell2mat(F(i,2))-sum(cell2mat((G(ismember(G(:,5),F(i,1)),2))));
    if xn(i) > 0
        n = n + 1;
        G{n,1} = ['unassigned_' F{i,1}];
        G{n,2} = xn(i);
        G{n,3} = F{i,3};
        G{n,4} = F{i,4};
        G{n,5} = F{i,1};
    end
end

% Speceis
ids = find(ismember(T.TaxRank,'S'));
nts = length(ids);
S = cell(nts,6);
for i=1:nts
    tName = T.Name{ids(i)};
    x = T.TaxLineage{ids(i)};
    x = split(x,"|");
    if any(contains(x,'c_'))
        xc = x{contains(x,'c_')};
    else
        xc = 'c_NA';
    end
    if any(contains(x,'o_'))
        xo = x{contains(x,'o_')};
    else
        xo = 'o_NA';
    end
    if any(contains(x,'f_'))
        xf = x{contains(x,'f_')};
    else
        xf = 'f_NA';
    end
    if any(contains(x,'g_'))
        xg = x{contains(x,'g_')};
    else
        xg = 'g_NA';
    end
    S{i,1} = tName(3:end);
    S{i,2} = T.CladeReads(ids(i));
    S{i,3} = xc(3:end);
    S{i,4} = xo(3:end);
    S{i,5} = xf(3:end);
    S{i,6} = xg(3:end);
end
% unassigned species
n = nts;
for i=ntg+1:length(G(:,1))
    if contains(G{i,1}, 'unassigned_')
        n = n +1;
        S{n,1} = G{i,1};
        S{n,2} = G{i,2};
        S{n,3} = G{i,3};
        S{n,4} = G{i,1};
        S{n,5} = G{i,1};
        S{n,6} = G{i,1};
    end
end
xn = zeros(1,ntg);
for i=1:ntg
    xn(i) = cell2mat(G(i,2))-sum(cell2mat((S(ismember(S(:,6),G(i,1)),2))));
    if xn(i) > 0
        n = n + 1;
        S{n,1} = ['unassigned_' G{i,1}];
        S{n,2} = xn(i);
        S{n,3} = G{i,3};
        S{n,4} = G{i,4};
        S{n,5} = G{i,5};
        S{n,6} = G{i,1};
    end
end

SN = char(options.inputSamples(sampleIndex));

d.SAMPLES.(SN).CLASS   = sortrows(C,2,'descend');
d.SAMPLES.(SN).ORDER   = sortrows(O,2,'descend');
d.SAMPLES.(SN).FAMILY  = sortrows(F,2,'descend');
d.SAMPLES.(SN).GENUS   = sortrows(G,2,'descend');
d.SAMPLES.(SN).SPECIES = sortrows(S,2,'descend');

end

% calculateRRA
function d = calculateRRA(d,sampleIndex,options,threshold)
% calculateRRA calculates the Relative Read Abundance (RRA) for a sample
% sampleIndex. 
% 
% USAGE: d = calculateRRA(d,sampleIndex,options)
%   d             data structure containg all the experiment data
%   sampleIndex   sample index
%   options       Structure of assigned and default options
%   threshold     Threshold to exclude the low-level background noise 
% 
% Copyright (C) 2022 Rahma Amen, Potsdam University, Potsdam, Germany

SN = char(options.inputSamples(sampleIndex));
AB = d.SAMPLES.(SN);

% Classes
x = cell2mat(AB.CLASS(:,2));
n = length(x(x/sum(x)>threshold));
RRA_C = cell(n,2);
for i=1:n
    RRA_C{i,1} = AB.CLASS{i,1};
    RRA_C{i,2} = x(i)/sum(x(1:n));
end
% Orders
x = cell2mat(AB.ORDER(:,2));
n = length(x(x/sum(x)>threshold));
RRA_O = cell(n,2);
for i=1:n
    RRA_O{i,1} = AB.ORDER{i,1};
    RRA_O{i,2} = x(i)/sum(x(1:n));
end

% Families
x = cell2mat(AB.FAMILY(:,2));
n = length(x(x/sum(x)>threshold));
RRA_F = cell(n,2);
for i=1:n
    RRA_F{i,1} = AB.FAMILY{i,1};
    RRA_F{i,2} = x(i)/sum(x(1:n));
end

% Genus
x = cell2mat(AB.GENUS(:,2));
n = length(x(x/sum(x)>threshold));
RRA_G = cell(n,2);
for i=1:n
    RRA_G{i,1} = AB.GENUS{i,1};
    RRA_G{i,2} = x(i)/sum(x(1:n));
end

% Species
x = cell2mat(AB.SPECIES(:,2));
n = length(x(x/sum(x)>threshold));
RRA_S = cell(n,2);
for i=1:n
    RRA_S{i,1} = AB.SPECIES{i,1};
    RRA_S{i,2} = x(i)/sum(x(1:n));
end

d.SAMPLES.(SN).RRA.RRA_CLASS = RRA_C;
d.SAMPLES.(SN).RRA.RRA_ORDER = RRA_O;
d.SAMPLES.(SN).RRA.RRA_FAMILY = RRA_F;
d.SAMPLES.(SN).RRA.RRA_GENUS = RRA_G;
d.SAMPLES.(SN).RRA.RRA_SPECIES = RRA_S;

end

% createTableRRA
function d = createTableRRA(d,options)
% createTableRRA creates RRA table for for all sample
% 
% USAGE: d = createTableRRA(d,options)
%   d             data structure containg all the experiment data
%   options       Structure of assigned and default options
% 
% Copyright (C) 2022 Rahma Amen, Potsdam University, Potsdam, Germany

species = fieldnames(d.SAMPLES);
SN = char(options.inputSamples(1));

levels = {'CLASS','ORDER','FAMILY','GENUS','SPECIES'};

for ilev=1:length(levels)
    lev = ['RRA_' levels{ilev}];
    taxaName = d.SAMPLES.(SN).RRA.(lev)(:,1);
    for i=2:length(species)
        SN = char(options.inputSamples(i));
        taxaNamei = d.SAMPLES.(SN).RRA.(lev)(:,1);
        C = union(taxaName,taxaNamei);
        taxaName = C;
    end
    
    rra = zeros(length(species),length(C));
    for ic=1:length(species)
        SN = char(options.inputSamples(ic));
        taxaName = d.SAMPLES.(SN).RRA.(lev)(:,1);
        taxaRRA = cell2mat(d.SAMPLES.(SN).RRA.(lev)(:,2));
        for ir=1:length(taxaName)
            idr = find(ismember(C,taxaName(ir)));
            rra(ic,idr) = taxaRRA(ir);
        end
    end
    
    % Save data to data structure
    d.RRA.(levels{ilev}).RRA = rra;
    d.RRA.(levels{ilev}).TAXA = C;
    
    % save data to csv file
    char_rem = {'/';'''';'-';':';'.';' ';'(';')';',';'<';'>'};
    for i=1:length(char_rem);C = strrep(C,char_rem{i},''); end
    T = array2table(rra);
    T.Properties.VariableNames = C;
    f = char(fullfile(options.OutputDir, [lev '.csv']));
    writetable(T,f)
end

end