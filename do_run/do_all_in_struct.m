function do_all_in_struct(whichPts)

%% Parameters
overwrite = 1;
also_validate = 0;

%% Get various path locations
locations = cceps_files; % Need to make a file pointing to you own path
script_folder = locations.script_folder;
results_folder = locations.results_folder;
mkdir(results_folder);

data_folder = locations.data_folder;


% add paths
addpath(genpath(script_folder));

%% Load pt file
pt = load([script_folder,'support/clinical_info/pt.mat']);
pt = pt.pt;

if isempty(whichPts)
    whichPts = 1:length(pt);
end

for p = whichPts
    
    fname = pt(p).ccep.file.name;
    
    if exist([results_folder,'results_',fname,'.mat'],'file') ~= 0
        if overwrite == 0
            fprintf('\nAlready did %s, skipping\n',fname);
            continue
        end
    end
    
    if strcmp(pt(p).ccep.file.ann,'empty')
        fprintf('\nNo annotations for %s, skipping\n',fname);
        continue
    end
    
    fprintf('\nDoing %s\n',fname);
    out = cceps_struct(pt,p);
    
    if also_validate
        fprintf('\nValidating...\n');
        random_rejections_keeps(out)
    end
    
end



end