clear all; clc;
addpath(genpath('.')); 
locations = cceps_files;
addpath(genpath(locations.freesurfer_matlab));
%% parcellate data

atlas_table = atlas_def;
coords = load('../data/elecs.mat');
waveforms = {'N1','N2'};

for subj = locations.subjects
    subj = char(subj); 
    load([locations.results_folder,sprintf('/results_%s_CCEP.mat',subj)]);
    if isfield(out,'parcellation')
        out = rmfield(out,'parcellation');    
    end
    out = ADD_COORDINATES_BIPOLAR(out,coords);
    for j_wave = 1:length(waveforms)
        wave = char(waveforms(j_wave));
        for atlas = 1:size(atlas_table,1) % loop through atlases
            parc_fname = char(atlas_table.atlas_file_name(atlas));
            parc_sname = char(atlas_table.atlas_short_name(atlas));
            fprintf([repmat('%',1,20),'\n','Parcellating CCEP networks for %s in %s atlas \n',repmat('%',1,20),'\n'],subj,parc_sname);           
            nifti = load_nifti(fullfile('../data','nifti',[parc_fname,'.nii.gz']));
            hemi_idx = 1:max(nifti.vol,[],'all'); % don't limit search to one hemisphere for now - in theory you could if you felt that electrode names were more accurate than MNI registration
            [A_parc,allElectrodesROIAssignments] = CCEPS_TO_PARCELS_MNI(out.locs_bipolar,out.network(j_wave).A,nifti,hemi_idx);
            out.parcellation.(wave).(parc_sname) = A_parc;
            out.parcellation.elec_parcels.(parc_sname) = allElectrodesROIAssignments;
            out.parcellation.parcels_n_elecs.(parc_sname) = cell2mat(cellfun(@(x) nansum(allElectrodesROIAssignments==x),num2cell(1:length(hemi_idx))','Uni',false));
            out.parcellation.n_elecs_parcellated.(parc_sname) = sum(~isnan(allElectrodesROIAssignments));
            
            % do parcellation with very tight threshold for purpose of being conservative with WM/GM locations
            [~,allElectrodesROIAssignments_tight] = CCEPS_TO_PARCELS_MNI(out.locs_bipolar,out.network(j_wave).A,nifti,hemi_idx,1);
            out.parcellation_tight.elec_parcels.(parc_sname) = allElectrodesROIAssignments_tight;
        end
    end
    save([locations.results_folder,sprintf('/results_%s_CCEP.mat',subj)],'out');    
end