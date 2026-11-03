function systems = buildSensitivitySystems(opts)
% buildSensitivitySystems  Collect the benchmark and SARS-CoV-2 records analyzed by runSensitivity.
%
% CALLS     generateCase1Series, defineCOVIDRegionConfigs, preprocessCOVIDPeriod
% CALLED BY runSensitivity, runSampleSize
%
% INPUT
%   opts.Benchmarks     cell array of ids from {'1a','1bi','1bii','1biii','1c'}
%   opts.COVID.Regions  cell array of region ids from defineCOVIDRegionConfigs, {} for none
%   opts.COVID.Periods  cell array from {'full','preVacc','postVacc'}
% OUTPUT
%   systems   struct array with id, label, group ('benchmark' or 'covid'), Mu, Beta
systems = struct('id',{},'label',{},'group',{},'Mu',{},'Beta',{});
benchLabel = containers.Map({'1a','1bi','1bii','1biii','1c'}, ...
    {'Case 1(a)','Case 1(b)(i)','Case 1(b)(ii)','Case 1(b)(iii)','Case 1(c)'});
for k = 1:numel(opts.Benchmarks)
    id = opts.Benchmarks{k};
    [Mu, Beta] = generateCase1Series(id);
    systems(end+1) = struct('id',id,'label',benchLabel(id),'group','benchmark','Mu',Mu,'Beta',Beta);
end
regions = defineCOVIDRegionConfigs();
perDefs = {'full','Full period'; 'preVacc','Pre-mass-vaccination'; 'postVacc','Post-mass-vaccination'};
for r = 1:numel(regions)
    if ~any(strcmp(regions(r).id, opts.COVID.Regions)), continue; end
    S = load(regions(r).dataFile);
    Nt = numel(S.Mu); brk = regions(r).periodBreakIndex;
    ranges = struct('full',[1 Nt],'preVacc',[1 brk-1],'postVacc',[brk Nt]);
    for q = 1:size(perDefs,1)
        pid = perDefs{q,1};
        if ~any(strcmp(pid, opts.COVID.Periods)), continue; end
        per = preprocessCOVIDPeriod(S.Mu(:), S.Beta(:), ranges.(pid)(1), ranges.(pid)(2), pid, perDefs{q,2});
        systems(end+1) = struct('id',[regions(r).id '_' pid], ...
            'label',[regions(r).label ', ' lower(perDefs{q,2})],'group','covid', ...
            'Mu',per.Mu,'Beta',per.Beta);
    end
end
end
