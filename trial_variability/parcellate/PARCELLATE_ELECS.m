function A_Parcellated = PARCELLATE_ELECS(A,allROIElectrodeAssignments,count,dim)

% INPUT:
% A: NxN electrode x electrode matrix
% allROIElectrodeAssignments: Mx1 cell for M anatomic parcels, each element
% contains electrode indices
% count: whether to sum or take mean
% dim: which dimension to use for mean/sum
%
% OUTPUTS:
% CCEPs: MxN if dim = 1, NxM if dim = 2, matrix of partially parcellated
% (either by stim or recording electrodes)

if count
    fun = @nansum;
elseif ~count
    fun = @nanmean;
end

nparc = length(allROIElectrodeAssignments);
nElecs = length(A);

if dim == 1
    A_Parcellated = nan(nparc,nElecs);
else
    A_Parcellated = nan(nElecs,nparc);
end

allROIElectrodeAssignments = cell(nparc,1);
for ROI = 1:nparc
    if ~isempty(allROIElectrodeAssignments{ROI}) % if there are electrodes in an ROI
        % assign mean expression across those electrodes to the corresponding ROI
        % average over rows
        if dim == 1
            A_Parcellated(ROI,:) = ...
                fun(A(allROIElectrodeAssignments{ROI},:),dim);
        else
            A_Parcellated(:,ROI) = ...
                fun(A(:,allROIElectrodeAssignments{ROI}),dim);
        end
    end
end

