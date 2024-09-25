function [num_Hf17, num_Hf27, num_Hf37, num_Hf47, num_Hf57, num_Hf67, num_Hf71, num_Hf72, num_Hf73, num_Hf74, num_Hf75, num_Hf76, num_Hf77] = d2Pfdp_dxvtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart)
%D2PFDP_DXVTMA2PERT  Computes 2nd derivatives of Droop Control "from" w.r.t. vtmaVa, vtmaVm, vtmaBeqz, vtmaBeqv, vtmaSh, vtmaqtma, Vavtma, Vmvtma, Beqzvtma, Beqvvtma, Shvtma, qtmavtma, vtmavtma (Finite Differences Method).
%
%   The derivatives will be taken with respect to polar or cartesian coordinates
%   depending on the 7th argument. So far only Polar has been coded
%   
%   [HfvtmaVa, HfvtmaVm, HfvtmaBz, HfvtmaBv, HfvtmaSh, Hfvtmaqtma, HfVavtma, HfVmvtma, HfBzvtma, HfBvvtma, HfShvtma, Hfqtmavtma, Hfvtmavtma] = D2PFDP_DXVTMA2(BASEMVA, BUS, BRANCH, V, LAM, PERT, 0)
%
%   Returns 13 matrices containing the partial derivatives w.r.t. Va, Vm,
%   Beq, Sh and ma of the product of a vector LAM with the 1st partial 
%   derivatives of the complex branch power flows.
%
%   [HfvtmaVa, HfvtmaVm, HfvtmaBz, HfvtmaBv, HfvtmaSh, Hfvtmaqtma, HfVavtma, HfVmvtma, HfBzvtma, HfBvvtma, HfShvtma, Hfqtmavtma, Hfvtmavtma] = D2PFDP_DXVTMA2(BASEMVA, BUS, BRANCH, V, LAM, PERT, 1)
%
%   Not Coded Yet
%
%   Examples:
%   [HfvtmaVa, HfvtmaVm, HfvtmaBz, HfvtmaBv, HfvtmaSh, Hfvtmaqtma, HfVavtma, HfVmvtma, HfBzvtma, HfBvvtma, HfShvtma, Hfqtmavtma, Hfvtmavtma] = d2Pfdp_dxvtma2Pert(baseMVA, bus, branch, V, lam, pert, 0)
%                                                               
%       Here the output matrices correspond to:
%           HvtmaVa   = d/dvtma (dPfdp_dVa.'   * mu)
%           HvtmaVm   = d/dvtma (dPfdp_dVm.'   * mu)
%           HvtmaBz   = d/dvtma (dPfdp_dBeqz.' * mu)
%           HvtmaBv   = d/dvtma (dPfdp_dBeqv.' * mu)
%           HvtmaSh   = d/dvtma (dPfdp_dSh.'   * mu)
%           Hvtmaqtma = d/dvtma (dPfdp_dqtma.' * mu)
%           HVavtma   = d/dVa   (dPfdp_dvtma.' * mu)
%           HVmvtma   = d/dVm   (dPfdp_dvtma.' * mu)
%           HBzvtma   = d/dBeqz (dPfdp_dvtma.' * mu)
%           HBvvtma   = d/dBeqv (dPfdp_dvtma.' * mu)
%           HShvtma   = d/dSh   (dPfdp_dvtma.' * mu)
%           Hqtmavtma = d/dqtma (dPfdp_dvtma.' * mu)
%           Hvtmavtma = d/dvtma (dPfdp_dvtma.' * mu)
%
%   For more details on the derivations behind the derivative code used
%   in MATPOWER information, see:
%
%   [TN2]  R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
%          their Derivatives using Complex Matrix Notation", MATPOWER
%          Technical Note 2, February 2010.
%             http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf
%   [TN4]  B. Sereeter and R. D. Zimmerman, "AC Power Flows and their
%          Derivatives using Complex Matrix Notation and Cartesian
%          Coordinate Voltages," MATPOWER Technical Note 4, April 2018.
%             http://www.pserc.cornell.edu/matpower/
%                                           TN4-OPF-Derivatives-Cartesian.pdf
%   [TNX]  A. Alvarez-Bustos, "AC/DC Power Flows and their Derivatives using
%          FUBM Complex Matrix Notation" MATPOWER Technical Note x, Month 20XX.
%             http://www.pserc.cornell.edu/matpower/
%                                           TNX-OPF-Derivatives-FUBM.pdf   %%AAB- Technical note to be written
                                          
%   ABRAHAM ALVAREZ BUSTOS
%   This code is a modification of MATPOWER code to include
%   The Flexible Universal Branch Model (FUBM) for MATPOWER
%   For more info about the model, email:
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk

%   MATPOWER
%   Copyright (c) 2008-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3, KDP] = idx_brch;%<<AAB-extra fields for FUBM
%% default input args
if nargin < 7
    vcart = 0;      %% default to polar coordinates
end

%% constants
nl = length(branch(:,1)); %% number of lines
nb = length(V);  %% number of buses
Vm = bus(:, VM);
Va = bus(:, VA) * pi/180;
%Changing Perturbation to Degrees for Theta_sh since inside makeYbus it gets changed to radians.
pertDeg = (pert*180)/pi;
[stat, Cf, Ct, k2, tap, Ys, Bc, Beq] = getbranchdata(branch, nb); %AAB- Gets the requested data from branch
YttB = Ys + 1j*Bc/2 + 1j*Beq;

%% identifier of AC/DC grids
%%AAB---------------------------------------------------------------------- 
iBeqz = find ((branch(:,CONV)==1 | branch(:,CONV)==3 ) & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq

%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360)& (branch(:, CONV)~=3) & (branch(:, CONV)~=4)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift
iQtma = find (branch(:,QT)~=0 & branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap

%% Find elements with Voltage Droop Control and slope
iPfdp = find( (branch(:,VF_SET)~=0) & (branch(:,KDP)~=0) & (branch(:, BR_STATUS)~=0) & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360) & (branch(:, CONV)==3 | branch(:, CONV)==4) ); %AAB- Find branch locations of the branch elements with Pf-Vdc Droop Control [nPfdp,1] (VSCIII)
nPfdp = length(iPfdp); %AAB- Number of VSC with Voltage Droop Control by theta_shift

Kdp = sparse(zeros(nl,1));
Kdp(iPfdp) = branch(iPfdp,KDP); %Droop Control Slope 

%% Calculation of derivatives
if vcart
    error('d2Pfdp_dxvtma2Pert: Derivatives of Flow Limit equations w.r.t ma using Finite Differences Method in cartasian have not been coded yet')    

else %AAB- Polar Version
    %Auxiliary 1st Derivatives
    diagV = sparse(1:nb, 1:nb, V, nb, nb); %AAB- diagV = sparse(diag(V))
    dVm=diag(V./abs(V)); %dV_Vm
    dVa=1j*diagV; %dV_Va  
       
    %Make Ybus, Yf, Yt
    [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

    %Derivatives of Voltage Magnitude w.r.t. Voltage magnitude
    dVmf_dVm = sparse(zeros(nl,nb));                          % Initialize for speed [nl,nb]
    fdp = branch(iPfdp, F_BUS);                               % List of "from" buses with Voltage Droop Control [nPfdp, 1]
    Cfdp = sparse(1:nPfdp, fdp, ones(nPfdp, 1), nPfdp, nb);   % connection matrix for line & from buses with Voltage Droop Control [nPfdp, nb]
    dVmf_dVm(iPfdp,:)=Cfdp;        
    
    %Sbr 1st Derivatives 
    [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
    [dSf_dBeqz, dSt_dBeqz] = dSbr_dBeq(branch, V, 1, vcart);
    [dSf_dBeqv, dSt_dBeqv] = dSbr_dBeq(branch, V, 2, vcart);
    [dSf_dPfsh, dSt_dPfsh] = dSbr_dsh(branch, V, 1, vcart);
    [dSf_dQtma, dSt_dQtma] = dSbr_dma(branch, V, 2, vcart);
    [dSf_dVtma, dSt_dVtma] = dSbr_dma(branch, V, 4, vcart);    
    
    %Pfdp 1st Derivatives
    dPfdp_dV1 = -real(dSf_dV1);
    dPfdp_dV2 = -real(dSf_dV2) + Kdp.*( dVmf_dVm );
    dPfdp_dBeqz = -real(dSf_dBeqz);
    dPfdp_dBeqv = -real(dSf_dBeqv);
    dPfdp_dPfsh = -real(dSf_dPfsh);
    dPfdp_dQtma = -real(dSf_dQtma); 
    dPfdp_dVtma = -real(dSf_dVtma);     
       
    %Selector of active Beqz 
    BeqzAux1 = sparse( zeros(nl,1) );       %AAB- Vector of zeros for the seclector
    BeqzAux1(iBeqz) = 1;                    %AAB- Fill the selector with 1 where Beq is active
    diagBeqzsel = sparse( diag(BeqzAux1) ); %AAB- Beq Selector [nl,nBeqx]
    BeqzAux2 = sparse( zeros(nl,nl));       %AAB- Beq second derivative Selector [nl, nl], the second derivative is zero 
      
    %Selector of active Beq 
    BeqvAux1 = sparse( zeros(nl,1) );       %AAB- Vector of zeros for the seclector
    BeqvAux1(iBeqv) = 1;                    %AAB- Fill the selector with 1 where Beq is active
    diagBeqvsel = sparse( diag(BeqvAux1) ); %AAB- Beq Selector [nl,nBeqx]
    BeqvAux2 = sparse( zeros(nl,nl));       %AAB- Beq second derivative Selector [nl, nl], the second derivative is zero 
        
    %Selector of active Theta_sh 
    ShSel = sparse( zeros(nl,1) );        %AAB- Vector of zeros for the selector
    ShSel(iPfsh) = 1;                     %AAB- Fill the selector with 1 where Theta_sh is active controlling Pf
    diagShSel = sparse( diag(ShSel));     %AAB- Diagonal of the selector for derivative w.r.t. ShSh, size [nl,nl]
    diagYssh = sparse( diag(ShSel.*Ys) ); %AAB- Theta_shift selector multilied by the series addmitance Ys, size [nl,nl]
 
    %Selector of active Qt ma/tap 
    QtmaSel = sparse( zeros(nl,1) );          %AAB- Vector of zeros for the selector
    QtmaSel(iQtma) = 1;                       %AAB- Fill the selector with 1 where ma is active controlling Qt
    diagQtmaSel = sparse( diag(QtmaSel));     %AAB- Diagonal of the selector for derivative w.r.t. Qtma, size [nl,nl]
    diagYsQtma = sparse( diag(QtmaSel.*Ys) ); %AAB- Qtma selector multilied by the series addmitance Ys, size [nl,nl]
    
    %Selector of active Vt ma/tap 
    VtmaSel = sparse( zeros(nl,1) );          %AAB- Vector of zeros for the selector
    VtmaSel(iVtma) = 1;                       %AAB- Fill the selector with 1 where ma is active controlling Vt
    diagVtmaSel = sparse( diag(VtmaSel));     %AAB- Diagonal of the selector for derivative w.r.t. Vtma, size [nl,nl]
    diagYsVtma = sparse( diag(VtmaSel.*Ys) ); %AAB- Vtma selector multilied by the series addmitance Ys, size [nl,nl]
     
    %Dimensionalize (Allocate for computational speed)  
    d2Pfdp_dVtmaVa   = sparse( zeros(nb,   nVtma) ); 
    d2Pfdp_dVtmaVm   = sparse( zeros(nb,   nVtma) );
    d2Pfdp_dVtmaBeqz = sparse( zeros(nBeqz,nVtma) );
    d2Pfdp_dVtmaBeqv = sparse( zeros(nBeqv,nVtma) );
    d2Pfdp_dVtmaPfsh = sparse( zeros(nPfsh,nVtma) );   
    d2Pfdp_dVtmaQtma = sparse( zeros(nQtma,nVtma) );  
    
    d2Pfdp_dVaVtma   = sparse( zeros(nVtma,nb   ) ); 
    d2Pfdp_dVmVtma   = sparse( zeros(nVtma,nb   ) );
    d2Pfdp_dBeqzVtma = sparse( zeros(nVtma,nBeqz) );
    d2Pfdp_dBeqvVtma = sparse( zeros(nVtma,nBeqv) );
    d2Pfdp_dPfshVtma = sparse( zeros(nVtma,nPfsh) );   
    d2Pfdp_dQtmaVtma = sparse( zeros(nVtma,nQtma) );     
    
    d2Pfdp_dVtma2    = sparse( zeros(nVtma,nVtma) );

    %VtmaVa num_G71
    for k = 1:nb
        V1p = V;
        V1p(k) = Vm(k) * exp(1j * (Va(k) + pert));  %% perturb Va
        [dSf_dVtma_PertVa, dSt_dVtma_PertVa] = dSbr_dma(branch, V1p, 4, vcart); %dSbr_dVtmaPertVa
        dPfdp_dVtma_PertVa = -real(dSf_dVtma_PertVa);  
        d2Pfdp_dVaVtma(:, k) = (dPfdp_dVtma_PertVa - dPfdp_dVtma).' * lam / pert; %VtmaVa from, size of [nVtma, nb] 
    end
    %VtmaVm num_G72
    for k = 1:nb
        V2p = V;
        V2p(k) = (Vm(k) + pert) * exp(1j * Va(k));  %% perturb Vm
        [dSf_dVtma_PertVm, dSt_dVtma_PertVm] = dSbr_dma(branch, V2p, 4, vcart); %dSbr_dVtmaPertVm
        dPfdp_dVtma_PertVm = -real(dSf_dVtma_PertVm);  
        d2Pfdp_dVmVtma(:, k) = (dPfdp_dVtma_PertVm - dPfdp_dVtma).' * lam / pert; %VtmaVm from, size of [nVtma, nb] 
    end
    %VtmaBeqz num_G73
    for k=1:nBeqz 
        PertSel=diagBeqzsel(:,iBeqz(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqz
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dVtmaPertBeqz evaluated in x+pert
        [dSf_dVtma_PertBeqz, dSt_dVtma_PertBeqz] = dSbr_dma(branch_Pert, V, 4, vcart); %dSbr_dVtmaPertBeqz
        dPfdp_dVtma_PertBeqz = -real(dSf_dVtma_PertBeqz);         
        %2nd Derivatives of Pfdp w.r.t. VtmaBeqz
        d2Pfdp_dBeqzVtma(:, k) = (dPfdp_dVtma_PertBeqz - dPfdp_dVtma).' * lam / pert;  %VtmaBeqz from, size of [nVtma, nBeqz]
    end
    %VtmaBeqv num_G74
    for k=1:nBeqv 
        PertSel=diagBeqvsel(:,iBeqv(k)); %AAB- Selects the column of diagBeqsel representing only the active Beqv
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Beq in the Perturbed branch (One vsc at a time)
        branch_Pert(:,BEQ) = branch(:,BEQ) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dVtmaPertBeqv evaluated in x+pert
        [dSf_dVtma_PertBeqv, dSt_dVtma_PertBeqv] = dSbr_dma(branch_Pert, V, 4, vcart); %dSbr_dVtmaPertBeqv
        dPfdp_dVtma_PertBeqv = -real(dSf_dVtma_PertBeqv);
        %2nd Derivatives of Pfdp w.r.t. VtmaBeqv
        d2Pfdp_dBeqvVtma(:, k) = (dPfdp_dVtma_PertBeqv - dPfdp_dVtma).' * lam / pert;  %VtmaBeqv from, size of [nVtma, nBeqv]
    end
    %VtmaPfsh num_G75
    for k=1:nPfsh 
        PertSel=diagShSel(:,iPfsh(k)); %AAB- Selects the column of diagShsel representing only the active Pfsh
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing Theta_sh in the Perturbed branch (One Pfsh at a time)
        branch_Pert(:,SHIFT) = branch(:,SHIFT) + (pertDeg.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dVtmaPertPfsh evaluated in x+pert
        [dSf_dVtma_PertPfsh, dSt_dVtma_PertPfsh] = dSbr_dma(branch_Pert, V, 4, vcart); %dSbr_dVtmaPertPfsh
        dPfdp_dVtma_PertPfsh = -real(dSf_dVtma_PertPfsh);
        %2nd Derivatives of Sbr w.r.t. VtmaPfsh
        d2Pfdp_dPfshVtma(:, k) = (dPfdp_dVtma_PertPfsh - dPfdp_dVtma).' * lam / pert;  %VtmaPfsh from, size of [nVtma, nPfsh] 
    end    
    %VtmaQtma num_G76
    for k=1:nQtma 
        PertSel=diagQtmaSel(:,iQtma(k)); %AAB- Selects the column of diagQtmasel representing only the active Qtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Qtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dVtmaPertQtma evaluated in x+pert
        [dSf_dVtma_PertQtma, dSt_dVtma_PertQtma] = dSbr_dma(branch_Pert, V, 4, vcart); %dSbr_dVtmaPertQtma
        dPfdp_dVtma_PertQtma = -real(dSf_dVtma_PertQtma);        
        %2nd Derivatives of Sbr w.r.t. VtmaQtma
        d2Pfdp_dQtmaVtma(:, k) = (dPfdp_dVtma_PertQtma - dPfdp_dVtma).' * lam / pert;  %VtmaQtma from, size of [nVtma, nQtma] 
    end 
    
    %VxVtma num_G17 num_G27
    for k=1:nVtma 
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmasel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        [Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dVxPertQtma evaluated in x+pert
        [dSf_dV1_PertVtma, dSf_dV2_PertVtma, dSt_dV1_PertVtma, dSt_dV2_PertVtma, Sf, St] = dSbr_dV(branch_Pert, Yf_Pert, Yt_Pert, V, vcart);
        
        %dPfdp_dVxPertVtma evaluated in x+pert
        dPfdp_dV1_PertVtma = -real(dSf_dV1_PertVtma);
        dPfdp_dV2_PertVtma = -real(dSf_dV2_PertVtma) + Kdp.*( dVmf_dVm );
        
        %2nd Derivatives of Pfdp w.r.t. VtmaVx
        d2Pfdp_dVtmaVa(:, k) = (dPfdp_dV1_PertVtma - dPfdp_dV1).' * lam / pert;  %VaVtma from, size of [nb, nVtma] 
        d2Pfdp_dVtmaVm(:, k) = (dPfdp_dV2_PertVtma - dPfdp_dV2).' * lam / pert;  %VmVtma  to , size of [nb, nVtma] 
    end
    %BeqzVtma num_G37
    for k=1:nVtma 
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagQtmasel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dBeqzPertVtma evaluated in x+pert
        [dSf_dBeqz_PertVtma, dSt_dBeqz_PertVtma] = dSbr_dBeq(branch_Pert, V, 1, vcart); %dSbr_dBeqzPertVtma
        dPfdp_dBeqz_PertVtma = -real(dSf_dBeqz_PertVtma);         
        %2nd Derivatives of Pfdp w.r.t. BeqzVtma
        d2Pfdp_dVtmaBeqz(:, k) = (dPfdp_dBeqz_PertVtma - dPfdp_dBeqz).' * lam / pert;  %BeqzVtma from, size of [nBeqz, nVtma] 
    end
    %BeqvVtma num_G47
    for k=1:nVtma 
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmasel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dBeqvPertVtma evaluated in x+pert
        [dSf_dBeqv_PertVtma, dSt_dBeqv_PertVtma] = dSbr_dBeq(branch_Pert, V, 2, vcart); %dSbr_dBeqvPertVtma
        dPfdp_dBeqv_PertVtma = -real(dSf_dBeqv_PertVtma);        
        %2nd Derivatives of Sbr w.r.t. BeqvVtma
        d2Pfdp_dVtmaBeqv(:, k) = (dPfdp_dBeqv_PertVtma - dPfdp_dBeqv).' * lam / pert;  %BeqvVtma from, size of [nBeqv, nVtma] 
    end
    %PfshVtma num_G57
    for k=1:nVtma
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmaSel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);
        %dSbr_dPfshPertVtma evaluated in x+pert
        [dSf_dPfsh_PertVtma, dSt_dPfsh_PertVtma] = dSbr_dsh(branch_Pert, V, 1, vcart); %dSbr_dPfshPertVtma
        dPfdp_dPfsh_PertVtma = -real(dSf_dPfsh_PertVtma);
        %2nd Derivatives of Sbr w.r.t. PfshVtma   
        d2Pfdp_dVtmaPfsh(:, k) = (dPfdp_dPfsh_PertVtma - dPfdp_dPfsh).' * lam / pert;  %PfshVtma from, size of [nPfsh , nVtma] 
    end
    %QtmaVtma num_G67
    for k=1:nVtma
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmaSel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %dSbr_dQtmaPertVtma evaluated in x+pert
        [dSf_dQtma_PertVtma, dSt_dQtma_PertVtma] = dSbr_dma(branch_Pert, V, 2, vcart); %dSbr_dQtmaPertVtma
        dPfdp_dQtma_PertVtma = -real(dSf_dQtma_PertVtma);
        %2nd Derivatives of Pfdp w.r.t. QtmaVtma   
        d2Pfdp_dVtmaQtma(:, k) = (dPfdp_dQtma_PertVtma - dPfdp_dQtma).' * lam / pert;  %QtmaVtma from, size of [nQtma , nVtma] 
    end
    %Vtma2 num_G77
    for k=1:nVtma
        PertSel=diagVtmaSel(:,iVtma(k)); %AAB- Selects the column of diagVtmaSel representing only the active Vtma
        %Restoring perturbated branch to the original one
        branch_Pert = branch;
        %Perturbing ma in the Perturbed branch (One Vtma at a time)
        branch_Pert(:,TAP) = branch(:,TAP) + (pert.* PertSel); 
        %Make Ybus, Yf, Yt Perturbated
        %[Ybus_Pert, Yf_Pert, Yt_Pert] = makeYbus(baseMVA, bus, branch_Pert);        
        %dSbr_dVtmaPertVtma evaluated in x+pert
        [dSf_dVtma_PertVtma, dSt_dVtma_PertVtma] = dSbr_dma(branch_Pert, V, 4, vcart); %dSbr_dVtmaPertVtma
        dPfdp_dVtma_PertVtma = -real(dSf_dVtma_PertVtma);        
        %2nd Derivatives of Pfdp w.r.t. Vtma2   
        d2Pfdp_dVtma2(:, k) = (dPfdp_dVtma_PertVtma - dPfdp_dVtma).' * lam / pert;  %Vtma2 from, size of [nVtma , nVtma] 
    end
end

%Assigning the partial derivatives with their respective outputs. 
num_Hf17 = sparse(d2Pfdp_dVtmaVa);
num_Hf27 = sparse(d2Pfdp_dVtmaVm);
num_Hf37 = sparse(d2Pfdp_dVtmaBeqz);
num_Hf47 = sparse(d2Pfdp_dVtmaBeqv);
num_Hf57 = sparse(d2Pfdp_dVtmaPfsh);
num_Hf67 = sparse(d2Pfdp_dVtmaQtma);

num_Hf71 = sparse(d2Pfdp_dVaVtma);
num_Hf72 = sparse(d2Pfdp_dVmVtma);
num_Hf73 = sparse(d2Pfdp_dBeqzVtma);
num_Hf74 = sparse(d2Pfdp_dBeqvVtma);
num_Hf75 = sparse(d2Pfdp_dPfshVtma);
num_Hf76 = sparse(d2Pfdp_dQtmaVtma);

num_Hf77 = sparse(d2Pfdp_dVtma2);
