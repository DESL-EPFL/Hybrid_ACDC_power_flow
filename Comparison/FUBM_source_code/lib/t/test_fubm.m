function test_fubm
%% TEST_FUBM Tests FUBM 
clear all;
clc; 
addpath '/Users/willem/Documents/phd/Git clones/matpower-fubm-master/lib'
t_jacobian_fubm;
t_hessian_fubm;
t_pf_acdc_fubm;
% t_opf_fubm_knitro;
% t_opf_fubm_mips;
% t_opf_fubm_ipopt;
% t_opf_fubm_fmincon;

