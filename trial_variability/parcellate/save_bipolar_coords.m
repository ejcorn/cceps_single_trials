%% this script computes bipolar coordinates for each patient and stores in one struct
clear all; close all; clc
addpath(genpath('.'));
%% define patient and save directories

locations = cceps_files;
mkdir(fullfile(locations.results_folder,'coords'));
coords = load('../data/elecs.mat');

locs_bipolar = struct();
pt = 'HUP213';
for pt = locations.subjects
    pt = char(pt);
    disp(pt);
    load([locations.results_folder,'/results_',pt,'_CCEP.mat']);    
    
    out = ADD_COORDINATES_BIPOLAR(out,coords);
    locs_bipolar.(pt) = out.locs_bipolar;    
end

save(fullfile(locations.results_folder,'coords','BipolarCoordinates.mat'),'locs_bipolar');