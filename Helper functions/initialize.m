function [Grid_para,idx,Filter_para,simulation_para] = initialize(nb_phases)

    %% Base values
    Grid_para.A_b = 1e5;
    Grid_para.V_b= 400;
    Grid_para.Y_b = Grid_para.A_b/Grid_para.V_b^2; 

    Grid_para.Vdc_b = 800;
    Grid_para.Adc_b = Grid_para.A_b;
    Grid_para.Ydc_b = Grid_para.Adc_b/Grid_para.Vdc_b^2; 

    %% Set the Grid parameters
    Grid_para.n_dc = 8;
    Grid_para.n_ac = 18;
    Grid_para.n_ph = nb_phases;
    Grid_para.n_nodes = Grid_para.n_ac*Grid_para.n_ph + Grid_para.n_dc;
    Grid_para.V_b =  Grid_para.V_b;
    Grid_para.Y_b =  Grid_para.Y_b;

    [Yac, YYL, YL, YT, YYT, I_b, Ampacities, y_lx, y_tx, A, linedata_ac]  = Ymatrix('linedata_AC.txt',Grid_para.A_b,Grid_para.V_b,[]);
    [Ydc, YYLdc, YLdc, YT_dc, YYTdc, I_bdc, Ampacitiesdc, y_ihdc, y_idc, Adc, linedata_dc]  = Ymatrix('linedata_DC.txt',Grid_para.Adc_b,Grid_para.Vdc_b,[]);

    dc_idx_min = min(linedata_dc(:,1:2),[],"all");
    dc_idx_max = max(linedata_dc(:,1:2),[],"all");

    Ydc = Ydc(dc_idx_min:dc_idx_max,dc_idx_min:dc_idx_max)/2; %remove empty rows and collumns
    Yac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),Yac,'UniformOutput',false));
    Grid_para.YY = blkdiag(Yac,Ydc);
    Grid_para.G = real(Grid_para.YY);
    Grid_para.B = imag(Grid_para.YY);

    %% Set the nodes types
    idx1.slack = 1;
    idx1.pqac = [2:14]';
    idx1.pvac = []';

    idx1.pdc = [23:26]';
    idx1.vdc = []';

    idx1.vscac_pq = [17,15]';
    idx1.vscac_vq = [18,16]';

    idx1.vscdc_pq = [21,19]';
    idx1.vscdc_vq = [22,20]';

    idx = Get_multiphase_Node_indices(idx1,Grid_para);
    linedata = [linedata_ac;linedata_dc];
    Grid_para = Get_Converter_para(idx1,linedata,Grid_para);

    %% Set the filter parameters
    Filter_para.R = 0.008*Grid_para.Y_b; %checked
    Filter_para.X = 0.04*Grid_para.Y_b;  %checked
    Filter_para.IGBT_piecewise = [  0                   0
                                    0.04926559627563	0.7
                                    2.30625399327864	0.8
                                    15.7793399043317	0.85
                                    107.547461516782	0.8999
                                    735.837403888342	0.9499
                                    1588.01477341768	0.9699];
    Filter_para.Exclude_losses = 1;

    %% Initialize load flow
    if Grid_para.n_ph == 3
        simulation_para.E_0 = [repmat([1; exp(-1i*pi*2/3);  exp(1i*pi*2/3) ], Grid_para.n_ac,1 ) ; ones(Grid_para.n_dc,1)];
    elseif Grid_para.n_ph == 1
        simulation_para.E_0 = [repmat(1, Grid_para.n_ac,1 ) ; ones(Grid_para.n_dc,1)];
    end
    
    simulation_para.tol = 1e-7;
    simulation_para.n_max = 100;

end