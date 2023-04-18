function idx = Get_multiphase_Node_indices(idx1,Grid_para)

    n_ph = Grid_para.n_ph;
    n_ac = Grid_para.n_ac;

    idx.slack = polyphase_indices(idx1.slack,n_ph);
    idx.pqac = polyphase_indices(idx1.pqac,n_ph);
    idx.pvac = polyphase_indices(idx1.pvac,n_ph);

    idx.pdc = (n_ph-1)*n_ac + idx1.pdc;
    idx.vdc = (n_ph-1)*n_ac + idx1.vdc;

    idx.vscac_pq = polyphase_indices(idx1.vscac_pq,n_ph);
    idx.vscac_vq = polyphase_indices(idx1.vscac_vq,n_ph);

    idx.vscdc_pq = (n_ph-1)*n_ac + idx1.vscdc_pq;
    idx.vscdc_vq = (n_ph-1)*n_ac + idx1.vscdc_vq;
end