function idx_poly = polyphase_indices(idx_mono,n_phases)
% idx_poly = polyphase_indices(idx_mono,n_phases)
%
% INPUT
% - idx_mono    monophase indices
% - n_phases    number of phases
%
% OUTPUT
% - idx_poly    polyphase indices

idx_poly = cell(length(idx_mono),1);

for i=1:size(idx_mono,1)
    idx_poly{i} = (idx_mono(i,:)-1)*n_phases + (1:n_phases)';
end

idx_poly = cell2mat(idx_poly);

end