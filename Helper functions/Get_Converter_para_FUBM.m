function Grid_para = Get_Converter_para_FUBM(idx1,linedata,Grid_para,mpc)

    n_ph = Grid_para.n_ph;
    n_ac = Grid_para.n_ac;
    AFE_type = strings;

    pos_ac = zeros(length(idx1.vscdc_pq)+length(idx1.vscdc_vq),2);
    pos_dc = zeros(length(idx1.vscdc_pq)+length(idx1.vscdc_vq),2);

    idx_dc = [(idx1.vscdc_vq);(idx1.vscdc_pq)];
    idx_ac = [(idx1.vscac_vq);(idx1.vscac_pq)];

    for i = 1: length(idx1.vscdc_pq)+length(idx1.vscdc_vq)

        ind_dc = [find(linedata(:,2) == idx_dc(i)),1];
        ind_ac = [find(linedata(:,2) == idx_ac(i)),1];

        pos_dc(i,:) = [idx_dc(i),linedata(ind_dc(1),ind_dc(2))];
        pos_ac(i,:) = [idx_ac(i),linedata(ind_ac(1),ind_ac(2))];

        pos_ac3 = polyphase_indices(pos_ac,n_ph);
        pos_dc3 = (n_ph-1)*n_ac + pos_dc;

        if find(idx_dc(i) == idx1.vscdc_pq)
            AFE_type(i) = "PQ";
        elseif find(idx_dc(i) == idx1.vscdc_vq)
            AFE_type(i) = "VdcQ";
        end

    end

    Grid_para.AFE_type = AFE_type;
    Grid_para.n_AFE = length(idx1.vscdc_pq) + length(idx1.vscdc_vq);
    Grid_para.pos_ac3 = pos_ac3;
    Grid_para.pos_dc3 = pos_dc3;
end