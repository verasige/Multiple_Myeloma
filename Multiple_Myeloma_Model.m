%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Co-culture Myeloma Model
% Elias Siguenza
% Institute of Metabolism and Systems Research (IMSR)
% School of Medicine and Dentistry / School of Mathematics
% University of Birmingham, England, United Kingdom
% 16 December 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following MATLAB script describes the interactions of the metabolic
% reactions in the simplified system. Please read the README.txt file
% for mire details. In short, the model requires an optimisation package
% called Gurobi, for which an academic licence is being used. 
% The model is compiled and ran using the COBRA toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initialise the Cobra Toolbox and do not update: (false)
clc; clear; close all; clc; initCobraToolbox(false); clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reactions of the model:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reactions are tagged to identify to which cell they belong:

% b_ = bone marrow cell
% [c] = bone marrow cell cytoplasm 
% [m] = bone marrow cell mitochondrion 
% [i] = bone marrow intermitochondrial membrane

% p_ = plasma cell
% [x] = plasma cell cytoplasm
% [l] = plasma cell mitochondrion 
% [n] = plasma intermitochondrial membrane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Extracellular Fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EX_glc = '1 glc_D[e] <=>';  % D-Glucose
EX_h   = '1 h[e] <=>';      % Protons (H+)
EX_h2o = '1 h2o[e] <=>';    % Water
EX_lac = '1 lac_L[e] <=>';  % L-Lactate
EX_co2 = '1 co2[e] <=>';    % Carbon Dioxide
EX_o2  = '1 o2[e] <=>';     % Oxygen
EX_glu = '1 glu_L[e] <=>';  % Glutamate
EX_gln = '1 gln_L[e] <=>';  % Glutamine
EX_hco3= '1 hco3[e] <=>';   % Bicarbonate
EX_asp = '1 asp_L[e] <=>';  % L-Aspartate
EX_urea= '1 urea[e] <=>';   % Urea
EX_hdca= '1 hdca[e] <=>';   % Palmitic Acid
EX_Lcys= '1 Lcystin[e] <=>';% L-Cystine
EX_cys = '1 cys_L[e] <=>';  % L-Cysteine
EX_oaa = '1 oaa[e] <=>';    % Oxaloacetate
EX_pyr = '1 pyr[e] <=>';    % Pyruvate
EX_pro = '1 pro_L[e] <=>';  % L-Proline
EX_val = '1 val_L[e] <=>';  % L-Valine
EX_crn = '1 crn[e] <=>';    % L-Carnitine
EX_drib= '1 drib[e] <=>';   % Deoxyribose
EX_met = '1 met_L[e]';      % Methionine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Extracellular Reactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eCA = '1 h2o[e] + 1 co2[e] <=> 1 hco3[e] + 1 h[e]'; % Carbonic Anhydrases
eSOD= '2 h[e] + 2 o2s[e] -> 1 h2o2[e] + 1 o2[e]';   % Superoxide Dimutases

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Monocarboxylate transporter 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_MCT1 = '1 pyr[e] <=> 1 pyr[x]';
b_MCT1 = '1 pyr[c] <=> 1 pyr[e]';   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Monocarboxylate transporter 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_MCT4 = '1 lac_L[e] <=> 1 lac_L[x]'; 
b_MCT4 = '1 lac_L[c] <=> 1 lac_L[e]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Other PM Transport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_xCT  = '1 glu_L[c] + 1 Lcystin[e] <=> 1 glu_L[e] + 1 Lcystin[c]'; % Cystine/ Glutamate
b_aspt = '1 asp_L[e] + 1 h[e] <=> 1 asp_L[c] + 1 h[c]';
b_Nhe  = '1 h[c] <=> 1 h[e]';
b_AQP2 = '1 h2o[c] <=> 1 h2o[e]';
b_AQP3 = '1 h2o[e] + 1 h2o2[e] <=> 1 h2o[c] + 1 h2o2[c]';
b_o2t  = '1 o2[e] <=> 1 o2[c]';
b_co2t = '1 co2[e] <=> 1 co2[c]';
b_ureat= '1 urea[c] <=> 1 urea[e]';
b_palt = '1 hdca[e] <=> 1 hdca[c]';
b_oaat = '1 oaa[e] <=> 1 oaa[c]';
b_cyst = '1 cys_L[e] <=> 1 cys_L[c]';
b_glnt = '1 gln_L[e] -> 1 gln_L[c]';
b_crnt = '1 crn[e] <=> 1 crn[c]';
b_dribt= '1 drib[e] <=> 1 drib[c]';
b_VALt = '1 val_L[e] <=> val_L[c]';
b_FUMt = '1 fum[c] <=>';

p_xCT  = '1 glu_L[x] + 1 Lcystin[e] <=> 1 glu_L[e] + 1 Lcystin[x]'; % Cystine/ Glutamate
p_aspt = '1 asp_L[e] + 1 h[e] <=> 1 asp_L[x] + 1 h[x]';
p_Nhe  = '1 h[x] <=> 1 h[e]';
p_AQP2 = '1 h2o[e] <=> 1 h2o[x]';
p_AQP3 = '1 h2o[e] + 1 h2o2[e] <=> 1 h2o[x] + 1 h2o2[x]';
p_o2t  = '1 o2[e] <=> 1 o2[x]';
p_co2t = '1 co2[e] <=> 1 co2[x]';
p_ureat= '1 urea[x] <=> 1 urea[e]';
p_palt = '1 hdca[e] <=> 1 hdca[x]';
p_oaat = '1 oaa[e] <=> 1 oaa[x]';
p_cyst = '1 cys_L[e] <=> 1 cys_L[x]';
p_glnt = '1 gln_L[e] -> 1 gln_L[x]';
p_crnt = '1 crn[e] <=> 1 crn[x]';
p_dribt= '1 drib[e] <=> 1 drib[x]';
p_VALt = '1 val_L[e] <=> val_L[x]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Demand Reactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are necessary reactions that assume other pathways, not involved in
% the model are demanding ATP/GTP/AMP etc...
b_DM_ATP  = '1 atp[c] + 1 h2o[c] <=> 1 adp[c] + 1 pi[c] + 1 h[c]';
b_DM_NADH = '1 h[c] + 1 nadh[c] <=> 1 nad[c]';
b_DM_AMP  = '1 amp[c] + 1 ppi[c] <=> 1 atp[c]';
b_DM_GMP  = '1 gmp[c] + 1 ppi[c] <=> 1 gtp[c]';
b_sink_coa= '1 coa[c] <=>';
b_NOX     = '1 nadph[c] + 1 o2[e] <=> 1 o2s[e] + 1 h[c] + 1 nadp[c]';
b_CA      = '1 h2o[c] + 1 co2[c] <=> 1 hco3[c] + h[c]';

p_DM_ATP = '1 atp[x] + 1 h2o[x] <=> 1 adp[x] + 1 pi[x] + 1 h[x]';
p_DM_NADH= '1 h[x] + 1 nadh[x] <=> 1 nad[x]';
p_DM_AMP = '1 amp[x] + 1 ppi[x] <=> 1 atp[x]';
p_DM_GMP = '1 gmp[x] + 1 ppi[x] <=> 1 gtp[x]';
p_sink_coa='1 coa[x] <=>';

p_NOX = '1 nadph[x] + 2 o2[e] <=> 2 o2s[e] + 1 h[x] + 1 nadp[x]';
p_CA  = '1 h2o[x] + 1 co2[x] <=> 1 hco3[x] + h[x]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Mito Membrane Transport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_AQPtm = '1 h2o[c] <=> 1 h2o[m]';
b_Nhetm = '1 h[c] <=> 1 h[m]';
b_CO2tm = '1 co2[c] <=> 1 co2[m]';
b_O2tm  = '1 o2[c] <=> 1 o2[m]';
b_mCA   = '1 h2o[c] + 1 co2[c] <=> 1 h[m] + 1 hco3[m]';
b_AQP8m = '1 h2o[c] + 1 nh4[c] <=> 1 h2o[m] + 1 nh4[m]';% Soria et al. (2011) AQP8 water & amonia
b_OOAtm = '1 oaa[c] <=> 1 oaa[m]';
b_VALtm = '1 val_L[c] <=> 1 val_L[m]';
b_GLUtm = '1 glu_L[c] + 1 h[c] <=> 1 glu_L[m] + 1 h[m]';
b_ALAtm = '1 ala_L[c] + 1 h[c] <=> 1 ala_L[m] + 1 h[m]';
b_FUMtm = '1 fum[m] <=>';
b_SERtm = '1 ser_L[c] <=> 1 ser_L[m]';
b_GLYtm = '1 gly[c] <=> 1 gly[m]';
b_PROtm = '1 pro_L[c] <=> 1 pro_L[m]';
b_sinkm_coa= '1 coa[m] <=>';

p_AQPtm = '1 h2o[x] <=> 1 h2o[l]';
p_Nhetm = '1 h[x] <=> 1 h[l]';
p_CO2tm = '1 co2[x] <=> 1 co2[l]';
p_O2tm  = '1 o2[x] <=> 1 o2[l]';
p_mCA   = '1 h2o[x] + 1 co2[x] <=> 1 h[l] + 1 hco3[l]';
p_AQP8m = '1 h2o[x] + 1 nh4[x] <=> 1 h2o[l] + 1 nh4[l]';% Soria et al. (2011) AQP8 water & amonia
p_OOAtm = '1 oaa[x] <=> 1 oaa[l]';
p_VALtm = '1 val_L[x] <=> 1 val_L[l]';
p_GLUtm = '1 glu_L[x] + 1 h[x] <=> 1 glu_L[l] + 1 h[l]';
p_ALAtm = '1 ala_L[x] + 1 h[x] <=> 1 ala_L[l] + 1 h[l]';
p_FUMtm = '1 fum[l] + 1 pi[l] <=>';
p_SERtm = '1 ser_L[x] <=> 1 ser_L[l]';
p_GLYtm = '1 gly[x] <=> 1 gly[l]';
p_PROtm = '1 pro_L[x] <=> 1 pro_L[l]';
p_fumtm = '1 fum[l] <=> 1 fum[x]';
p_fumt  = '1 fum[x] ->';
p_sinkm_coa = '1 coa[l] <=>';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Glycolysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_GLUT= '1 glc_D[e] <=> 1 glc_D[c]'; 
b_HEX = '1 atp[c] + 1 glc_D[c] -> 1 adp[c] + 1 g6p[c] + 1 h[c]';
b_PGI = '1 g6p[c] -> 1 f6p[c]';
b_PFK = '1 atp[c] + 1 f6p[c] -> 1 adp[c] + 1 fdp[c] + 1 h[c]';
b_FBA = '1 fdp[c] -> 1 dhap[c] + 1 g3p[c]';
b_TPI = '1 dhap[c] -> 1 g3p[c]';
b_GAPD= '1 g3p[c] + 1 nad[c] + 1 pi[c] -> 1 13dpg[c] + 1 h[c] + 1 nadh[c]';
b_PGK = '1 13dpg[c] + 1 adp[c] -> 1 3pg[c] + 1 atp[c]';
b_PGM = '1 3pg[c] -> 1 2pg[c]';
b_ENO = '1 2pg[c] -> 1 h2o[c] + 1 pep[c]';
b_PYK = '1 adp[c] + 1 h[c] + 1 pep[c] -> 1 atp[c] + 1 pyr[c]';
b_LDH = '1 h[c] + 1 nadh[c] + 1 pyr[c] <=> 1 lac_L[c] + 1 nad[c]';

p_GLUT= '1 glc_D[e] <=> 1 glc_D[x]'; 
p_HEX = '1 atp[x] + 1 glc_D[x] -> 1 adp[x] + 1 g6p[x] + 1 h[x]';
p_PGI = '1 g6p[x] -> 1 f6p[x]';
p_PFK = '1 atp[x] + 1 f6p[x] -> 1 adp[x] + 1 fdp[x] + 1 h[x]';
p_FBA = '1 fdp[x] -> 1 dhap[x] + 1 g3p[x]';
p_TPI = '1 dhap[x] -> 1 g3p[x]';
p_GAPD= '1 g3p[x] + 1 nad[x] + 1 pi[x] -> 1 13dpg[x] + 1 h[x] + 1 nadh[x]';
p_PGK = '1 13dpg[x] + 1 adp[x] -> 1 3pg[x] + 1 atp[x]';
p_PGM = '1 3pg[x] -> 1 2pg[x]';
p_ENO = '1 2pg[x] -> 1 h2o[x] + 1 pep[x]';
p_PYK = '1 adp[x] + 1 h[x] + 1 pep[x] -> 1 atp[x] + 1 pyr[x]';
p_LDH = '1 h[x] + 1 nadh[x] + 1 pyr[x] <=> 1 lac_L[x] + 1 nad[x]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Phospho-Pentose pathway
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_G6PDH = '1 g6p[x] + 1 nadp[x] -> 1 6pgl[x] + 1 h[x] + 1 nadph[x]';  
p_PGL = '1 6pgl[x] + 1 h2o[x] -> 1 6pgc[x] + 1 h[x]';
p_GND = '1 6pgc[x] + 1 nadp[x] -> 1 co2[x] + 1 nadph[x] + 1 ru5p_D[x]';
p_RPI = '1 ru5p_D[x] <=> 1 r5p[x]';
p_PPM = '1 r1p[x] <=> 1 r5p[x]';
p_PRPP= '1 atp[x] + 1 r5p[x] <=> 1 amp[x] + 1 h[x] + 1 prpp[x] + 1 ppi[x]';
p_TKT1= '1 r5p[x] + 1 xu5p_D[x] <=> 1 g3p[x] + 1 s7p[x]';
p_TALA= '1 g3p[x] + 1 s7p[x] <=> 1 e4p[x] + 1 f6p[x]';
p_TKT2= '1 e4p[x] + 1 f6p[x] <=> 1 xu5p_D[x] + 1 g3p[x]';

b_G6PDH = '1 g6p[c] + 1 nadp[c] -> 1 6pgl[c] + 1 h[c] + 1 nadph[c]';  
b_PGL = '1 6pgl[c] + 1 h2o[c] -> 1 6pgc[c] + 1 h[c]';
b_GND = '1 6pgc[c] + 1 nadp[c] -> 1 co2[c] + 1 nadph[c] + 1 ru5p_D[c]';
b_RPI = '1 ru5p_D[c] <=> 1 r5p[c]';
b_PPM = '1 r1p[c] <=> 1 r5p[c]';
b_PRPP= '1 atp[c] + 1 r5p[c] <=> 1 amp[c] + 1 h[c] + 1 prpp[c] + 1 ppi[c]';
b_TKT1= '1 r5p[c] + 1 xu5p_D[c] <=> 1 g3p[c] + 1 s7p[c]';
b_TALA= '1 g3p[c] + 1 s7p[c] <=> 1 e4p[c] + 1 f6p[c]';
b_TKT2= '1 e4p[c] + 1 f6p[c] <=> 1 xu5p_D[c] + 1 g3p[c]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% L-Aalanine Synthesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_ALATA = '1 glu_L[x] + 1 pyr[x] <=> 1 akg[x] + 1 ala_L[x]';
p_ALA_L = '1 ala_L[x] -> 1 bio[x]';

b_ALATA = '1 glu_L[c] + 1 pyr[c] <=> 1 akg[c] + 1 ala_L[c]';
b_ALA_L = '1 ala_L[c] -> 1 bio[c]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% L-Serine Synthesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_PGCD = '1 3pg[x] + 1 nad[x] -> 1 3php[x] + 1 h[x] + 1 nadh[x]';
p_PSERT= '1 3php[x] + 1 glu_L[x] -> 1 akg[x] + 1 pser_L[x]';
p_PSP  = '1 h2o[x] + 1 pser_L[x] -> 1 ser_L[x]';
p_GHMT = '1 ser_L[x] + 1 thf[x] <=> 1 gly[x] + 1 h2o[x] + 1 10fthf[x]';

b_PGCD = '1 3pg[c] + 1 nad[c] -> 1 3php[c] + 1 h[c] + 1 nadh[c]';
b_PSERT= '1 3php[c] + 1 glu_L[c] -> 1 akg[c] + 1 pser_L[c]';
b_PSP  = '1 h2o[c] + 1 pser_L[c] -> 1 ser_L[c]';
b_GHMT = '1 ser_L[c] + 1 thf[c] <=> 1 gly[c] + 1 h2o[c] + 1 10fthf[c]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Serine derived lipids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_SERPT= '1 h[x] + 1 pmtcoa[x] + 1 ser_L[x] -> 1 3dsphgn[x] + 1 co2[x] + 1 coa[x]';
p_DSPHR= '1 3dsphgn[x] + 1 h[x] + 1 nadph[x] -> 1 nadp[x] + 1 sphgn[x]';
p_SLIP = '1 sphgn[x] -> 1 bio[x]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Glutamine Metabolism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_GLUN = '1 gln_L[x] + 1 h2o[x] -> 1 glu_L[x] + 1 nh4[x]';
p_GLNtm= '1 gln_L[x] <=> 1 gln_L[l]';
p_GLUNm= '1 gln_L[l] + 1 h2o[l] -> 1 glu_L[l] + 1 nh4[l]';
p_GLUDm= '1 glu_L[l] + 1 h2o[l] + 1 nad[l] <=> 1 akg[l] + 1 h[l] + 1 nadh[l] + 1 nh4[l]';

b_GLUN = '1 gln_L[c] + 1 h2o[c] -> 1 glu_L[c] + 1 nh4[c]';
b_GLNtm= '1 gln_L[c] <=> 1 gln_L[m]';
b_GLUNm= '1 gln_L[m] + 1 h2o[m] -> 1 glu_L[m] + 1 nh4[m]';
b_GLUDm= '1 glu_L[m] + 1 h2o[m] + 1 nad[m] <=> 1 akg[m] + 1 h[m] + 1 nadh[m] + 1 nh4[m]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Valine Metabolism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_VALTAm = '1 akg[l] + 1 val_L[l] <=> 1 3mob[l] + 1 glu_L[l]';
p_OVID2m = '1 3mob[l] + 1 coa[l] + 1 nad[l] -> 1 co2[l] + 1 ibcoa[l] + 1 nadh[l]';
p_RE3190M= '1 accoa[l] + 1 ibcoa[l] + 1 nad[l]  + 1 fad[l] <=>  1 h[l] + 1 nadh[l] + 1 fadh2[l] + 1 ppcoa[l]';
p_MMCDm  = '1 co2[l] + 1 ppcoa[l] -> 1 h[l] + 1 mmcoa_S[l]';
p_MMEm   = '1 mmcoa_S[l] <=> 1 mmcoa_R[l]';
p_MMMm   = '1 mmcoa_R[l] <=> 1 succoa[l]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Metathione Metabolism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_ADOMET = '1 met_L[c] -> 1 hcys_L[c]';
b_METS   = '1 5mthf[c] + 1 hcys_L[c] -> 1 h[c] + 1 met_L[c] + 1 thf[c]';
b_MS     = '1 amet[c] + 1 thf[c] -> 1 5mthf[c] + 1 ahcys[c]';
b_HMR4710= '1 ahcys[c] + 1 h[c] -> 1 amet[c]';
b_CYSTS  = '1 hcys_L[c] + 1 ser_L[c] -> 1 cyst_L[c] + 1 h2o[c]';
b_r0210  = '1 cyst_L[c] -> 1 cys_L[c]';
b_CYSt   = '1 cys_L[c] <=> 1 cys_L[e]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Purine Synthesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_GLNS  = '1 glu_L[x] + 1 atp[x] <=> 1 adp[x] + 1 gln_L[x] + 1 h[x] + 1 pi[x]';
p_GLUPRT= '1 gln_L[x] + 1 h2o[x] + 1 prpp[x] -> 1 glu_L[x] + 1 pram[x]';
p_PRAGS = '1 pram[x] + 1 atp[x] + 1 gly[x] <=> 1 h[x] + 1 pi[x] + 1 gar[x] + 1 adp[x]';
p_GARFT = '1 10fthf[x] + 1 gar[x] <=> 1 h[x] + 1 thf[x] + 1 fgam[x]';
p_PRFGS = '1 atp[x] + 1 fgam[x] + 1 gln_L[x] + 1 h2o[x] -> 1 adp[x] + 1 fpram[x] + 1 glu_L[x] + 1 h[x] + 1 pi[x]';
p_r0666 = '1 atp[x] + 1 fpram[x] <=> 1 adp[x] + 1 air[x] + 1 h[x] + 1 pi[x]';
p_AIRC  = '1 air[x] + 1 co2[x] <=> 1 5aizc[x] + 1 h[x]';
p_PRASCS= '1 5aizc[x] + 1 asp_L[x] + 1 atp[x] <=> 1 25aics[x] + 1 adp[x] + 1 h[x] + 1 pi[x]';
p_ADSL2 = '1 25aics[x] -> 1 aicar[x] + 1 fum[x]';
p_AICART= '1 10fthf[x] + 1 aicar[x] <=> 1 fprica[x] + 1 thf[x]';
p_IMPC  = '1 fprica[x] <=> 1 h2o[x] + 1 imp[x]'; 
p_FTHFL = '1 atp[x] + 1 thf[x] <=> 1 10fthf[x] + 1 adp[x] + 1 pi[x]';

b_GLNS  = '1 glu_L[c] + 1 atp[c] <=> 1 adp[c] + 1 gln_L[c] + 1 h[c] + 1 pi[c]';
b_GLUPRT= '1 gln_L[c] + 1 h2o[c] + 1 prpp[c] -> 1 glu_L[c] + 1 pram[c]';
b_PRAGS = '1 pram[c] + 1 atp[c] + 1 gly[c] <=> 1 h[c] + 1 pi[c] + 1 gar[c] + 1 adp[c]';
b_GARFT = '1 10fthf[c] + 1 gar[c] <=> 1 h[c] + 1 thf[c] + 1 fgam[c]';
b_PRFGS = '1 atp[c] + 1 fgam[c] + 1 gln_L[c] + 1 h2o[c] -> 1 adp[c] + 1 fpram[c] + 1 glu_L[c] + 1 h[c] + 1 pi[c]';
b_r0666 = '1 atp[c] + 1 fpram[c] <=> 1 adp[c] + 1 air[c] + 1 h[c] + 1 pi[c]';
b_AIRC  = '1 air[c] + 1 co2[c] <=> 1 5aizc[c] + 1 h[c]';
b_PRASCS= '1 5aizc[c] + 1 atp[c] <=> 1 25aics[c] + 1 adp[c] + 1 h[c] + 1 pi[c]';
b_ADSL2 = '1 25aics[c] -> 1 aicar[c] + 1 fum[c]';
b_AICART= '1 10fthf[c] + 1 aicar[c] <=> 1 fprica[c] + 1 thf[c]';
b_IMPC  = '1 fprica[c] <=> 1 h2o[c] + 2 imp[c]'; 
b_FTHFL = '1 atp[c] + 1 thf[c] <=> 1 10fthf[c] + 1 adp[c] + 1 pi[c]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Nucleotide interconversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_DRBK  = '1 atp[c] + 1 drib[c] -> 1 2dr5p[c] + 1 adp[c] + 1 h[c] ';
b_DRPP  = '1 2dr5p[c] <=> 1 2dr1p[c]';
%%%%%%%%%%%%%%%%%%%%%%%%%%% Adenine:
b_ADSS  = '1 asp_L[c] + 1 imp[c] -> 1 dcamp[c] + 1 h[c] + 1 pi[c]';
b_ADSL  = '1 dcamp[c] -> 1 amp[c] + 1 fum[c]';
b_PUNP0 = '1 amp[c] + 1 h2o[c] -> 1 r1p[c] + 1 ade[c]';
b_PUNP2 = '1 2dr1p[c] + 1 ade[c] <=> 1 dad_2[c] + 1 pi[c]';
b_dA    = '1 dad_2[c] -> 1 bio[c]';
%%%%%%%%%%%%%%%%%%%%%%%%%%% Guanine:
b_IMPD  = '1 h2o[c] + 1 imp[c] + 1 nad[c] -> 1 h[c] + 1 nadh[c] + 1 xmp[c]';
b_GMPS  = '1 atp[c] + 1 gln_L[c] + 1 h2o[c] + 1 xmp[c] -> 1 amp[c] + 1 glu_L[c] + 1 gmp[c] + 1 h[c] + 1 ppi[c]';
b_PUNP3 = '1 gmp[c] + 1 h2o[c] -> 1 r1p[c] + 1 gua[c]';
b_PUNP4 = '1 2dr1p[c] + 1 gua[c] <=> 1 dgsn[c] + 1 pi[c]';
b_dG    = '1 dgsn[c] -> 1 bio[c]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_DRBK  = '1 atp[x] + 1 drib[x] -> 1 2dr5p[x] + 1 adp[x] + 1 h[x]';
p_DRPP  = '1 2dr5p[x] <=> 1 2dr1p[x]';
%%%%%%%%%%%%%%%%%%%%%%%%%%% Adenine:
p_ADSS  = '1 asp_L[x] + 1 imp[x] -> 1 dcamp[x] + 1 h[x]';
p_ADSL  = '1 dcamp[x] -> 1 amp[x] + 1 fum[x]';
p_PUNP0 = '1 amp[x] + 1 h2o[x] -> 1 r1p[x] + 1 ade[x]';
p_PUNP2 = '1 2dr1p[x] + 1 ade[x] <=> 1 dad_2[x] + 1 pi[x]';
p_dA    = '1 dad_2[x] -> 1 bio[x]';
%%%%%%%%%%%%%%%%%%%%%%%%%%% Guanine:
p_IMPD = '1 h2o[x] + 1 imp[x] -> 1 h[x] + 1 xmp[x]';
p_GMPS = '1 atp[x] + 1 gln_L[x] + 1 h2o[x] + 1 xmp[x] -> 1 amp[x] + 1 glu_L[x] + 1 gmp[x] + 1 h[x] + 1 ppi[x]';
p_PUNP3= '1 gmp[x] + 1 h2o[x] -> 1 r1p[x] + 1 gua[x]';
p_PUNP4= '1 2dr1p[x] + 1 gua[x] <=> 1 dgsn[x] + 1 pi[x]';
p_dG   = '1 dgsn[x] -> 1 bio[x]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Pyrimidine Synthesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_CBPS  = '1 atp[x] + 1 gln_L[x] + 1 h2o[x] + 1 hco3[x] -> 1 adp[x] + 1 cbp[x] + 1 glu_L[x] + 1 h[x] + 1 pi[x]';
p_ASPCT = '1 asp_L[x] + 1 cbp[x] <=> 1 cbasp[x] + 1 h[x] + 1 pi[x]';
p_DHORD = '1 cbasp[x] + 1 h[x] <=> 1 dhor_S[x] + 1 h2o[x]';
p_DHORD9= '1 dhor_S[x] + 1 q10[l] -> 1 orot[x] + 1 q10h2[l]';
p_ORPT  = '1 orot[x] + 1 prpp[x] <=> 1 orot5p[x]';
p_OMPDC = '1 orot5p[x] + 1 h[x] -> 1 co2[x] + 1 ump[x]';
p_URIK  = '1 adp[x] + 1 pi[x] + 1 h[x] + 1 ump[x] <=> 1 atp[x] + 1 uri[x]';
p_PYNP  = '1 pi[x] + 1 uri[x] <=> 1 r1p[x] + 1 ura[x]';
p_CSND  = '1 nh4[x] + 1 ura[x] -> 1 csn[x] + 1 h[x] + 1 h2o[x]';
p_dC    = '1 csn[x] -> 1 bio[x]';

b_CBPS  = '1 atp[c] + 1 gln_L[c] + 1 h2o[c] + 1 hco3[c] -> 1 adp[c] + 1 cbp[c] + 1 glu_L[c] + 1 h[c] + 1 pi[c]';
b_ASPCT = '1 asp_L[c] + 1 cbp[c] <=> 1 cbasp[c] + 1 h[c] + 1 pi[c]';
b_DHORD = '1 cbasp[c] + 1 h[c] <=> 1 dhor_S[c] + 1 h2o[c]';
b_DHORD9= '1 dhor_S[c] + 1 q10[m] -> 1 orot[c] + 1 q10h2[m]';
b_ORPT  = '1 orot[c] + 1 prpp[c] <=> 1 orot5p[c]';
b_OMPDC = '1 orot5p[c] + 1 h[c] -> 1 co2[c] + 1 ump[c]';
b_URIK  = '1 adp[c] + 1 pi[c] + 1 h[c] + 1 ump[c] <=> 1 atp[c] + 1 uri[c]';
b_PYNP  = '1 pi[c] + 1 uri[c] <=> 1 r1p[c] + 1 ura[c]';
b_CSND  = '1 nh4[c] + 1 ura[c] -> 1 csn[c] + 1 h[c] + 1 h2o[c]';
b_dC    = '1 csn[c] -> 1 bio[c]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Myeloma Biomass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_BioMass= '1 bio[x] ->';
p_r1p = '1 r1p[x] ->';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% TCA Cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Pyruvate Metabolism
b_PYRtm = '1 h[c] + 1 pyr[c] -> 1 h[m] + 1 pyr[m]';
b_PDHm  = '1 coa[m] + 1 nad[m]  + 1 pyr[m] -> 1 accoa[m] + 2 co2[m] + 1 nadh[m]';
b_PCm   = '1 atp[m] + 1 hco3[m] + 1 pyr[m] -> 1 adp[m] + 1 h[m] + 1 oaa[m] + 1 pi[m]';
b_CSm   = '1 accoa[m] + 1 h2o[m] + 1 oaa[m] -> 1 cit[m] + 1 coa[m] + 1 h[m]';
b_ACONTm= '1 cit[m] -> 1 icit[m]';
b_ICDHm = '1 icit[m] + 1 nad[m] -> 1 akg[m] + 1 co2[m] + 1 h[m] + 1 nadh[m]';
b_AKGDm = '1 coa[m] + 1 akg[m] + 1 nad[m] -> 1 nadh[m] + 1 succoa[m]';
b_SUCOASm='1 adp[m] + 1 pi[m] + 1 succoa[m] -> 1 atp[m] + 1 succ[m]';                  
b_SUCDm = '1 fad[m] + 1 succ[m] -> 1 fadh2[m] + 1 fum[m]'; 
b_FUMm  = '1 fum[m] + 1 h2o[m] -> 1 mal_L[m]';                           
b_MDHm  = '1 mal_L[m] + 1 nad[m] -> 1 h[m] + 1 nadh[m] + 1 oaa[m]';

p_PYRtm = '1 h[x] + 1 pyr[x] -> 1 h[l] + 1 pyr[l]';
p_PDHm  = '1 coa[l] + 1 nad[l]  + 1 pyr[l] -> 1 accoa[l] + 2 co2[l] + 1 nadh[l]';
p_PCm   = '1 atp[l] + 1 hco3[l] + 1 pyr[l] -> 1 adp[l] + 1 h[l] + 1 oaa[l] + 1 pi[l]';
p_CSm   = '1 accoa[l] + 1 h2o[l] + 1 oaa[l] -> 1 cit[l] + 1 coa[l] + 1 h[l]';
p_ACONTm= '1 cit[l] -> 1 icit[l]';
p_ICDHm = '1 icit[l] + 1 nad[l] -> 1 akg[l] + 1 co2[l] + 1 h[l] + 1 nadh[l]';
p_AKGDm = '1 coa[l] + 1 akg[l] + 1 nad[l] -> 1 nadh[l] + 1 succoa[l]';
p_SUCOASm='1 adp[l] + 1 pi[l] + 1 succoa[l] -> 1 atp[l] + 1 succ[l]';                  
p_SUCDm = '1 fad[l] + 1 succ[l] -> 1 fadh2[l] + 1 fum[l]'; 
p_FUMm  = '1 fum[l] + 1 h2o[l] -> 1 mal_L[l]';                           
p_MDHm  = '1 mal_L[l] + 1 nad[l] -> 1 h[l] + 1 nadh[l] + 1 oaa[l]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Fatty Acid Beta-Oxidation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_FACOAL= '1 atp[c] + 1 coa[c] + 1 hdca[c] -> 1 amp[c] + 1 pmtcoa[c] + 1 ppi[c]';
b_CPT1  = '1 crn[c] + 1 pmtcoa[c] <=> 1 coa[c] + 1 pmtcrn[c]';
b_CTRt  = '1 pmtcrn[c] <=> 1 pmtcrn[m]';
b_CPT2  = '1 coa[m] + 1 pmtcrn[m] <=> 1 pmtcoa[m]';
b_FAOX  = '4 coa[m] + 4 fad[m] + 4 h2o[m] + 4 nad[m] + 1 pmtcoa[m] -> 4 accoa[m] + 4 fadh2[m] + 4 h[m] + 4 nadh[m] + 1 occoa[m]';
b_FAOXC = '3 coa[m] + 3 fad[m] + 3 h2o[m] + 3 nad[m] + 1 occoa[m] -> 3 accoa[m] + 3 fadh2[m] + 3 h[m] + 3 nadh[m]';

p_FACOAL= '1 atp[x] + 1 coa[x] + 1 hdca[x] -> 1 amp[x] + 2 pmtcoa[x] + 1 ppi[x]';
p_CPT1  = '1 crn[x] + 1 pmtcoa[x] <=> 1 coa[x] + 1 pmtcrn[x]';
p_CTRt  = '1 pmtcrn[x] <=> 1 pmtcrn[l]';
p_CPT2  = '1 coa[l] + 1 pmtcrn[l] <=> 1 pmtcoa[l]';
p_FAOX  = '4 coa[l] + 4 fad[l] + 4 h2o[l] + 4 nad[l] + 1 pmtcoa[l] -> 4 accoa[l] + 4 fadh2[l] + 4 h[l] + 4 nadh[l] + 1 occoa[l]';
p_FAOXC = '3 coa[l] + 3 fad[l] + 3 h2o[l] + 3 nad[l] + 1 occoa[l] -> 3 accoa[l] + 3 fadh2[l] + 3 h[l] + 3 nadh[l]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ETC/ ROS-Detox:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_Cplx1 = '1 nadh[m] + 1 q10[m] + 5 h[m] + 0.3 o2[m] -> 1 nad[m] + 1 q10h2[m] + 4 h[i] + 0.05 o2s[m]'; % .068
b_Cplx2 = '1 fadh2[m] + 1 q10[m] -> 1 fad[m] + 1 q10h2[m]';
b_Cplx3 = '1 q10h2[m] + 2 h[m] + 2 ficytC[m] + 0.3 o2[m] -> 1 q10[m] + 6 h[i] + 2 focytC[m] + 0.05 o2s[m]'; % .068
b_Cplx4 = '4 focytC[m] + 1 h[m] + 0.3 o2[m] -> 4 ficytC[m] + 1 h2o[m] + 0.3 o2s[m] + 4 h[i]';
b_ATPSm = '4 h[i] + 1 adp[m] + 1 pi[m] -> 1 atp[m] + 1 h[m]';
b_ATPtm = '1 atp[m] + 1 adp[c] + 1 pi[c] + 1 h[c] <=> 1 atp[c] + 1 adp[m] + 1 pi[m] + 1 h[m]';

p_Cplx1 = '1 nadh[l] + 1 q10[l] + 5 h[l] + 0.3 o2[l] -> 1 nad[l] + 1 q10h2[l] + 4 h[n] + 0.05 o2s[l]'; % .068
p_Cplx2 = '1 fadh2[l] + 1 q10[l] -> 1 fad[l] + 1 q10h2[l]';
p_Cplx3 = '1 q10h2[l] + 2 h[l] + 2 ficytC[l] + 0.3 o2[l] -> 1 q10[l] + 6 h[n] + 2 focytC[l] + 0.05 o2s[l]'; % .068
p_Cplx4 = '4 focytC[l] + 1 h[l] + 0.3 o2[l] -> 4 ficytC[l] + 2 h2o[l] + 0.3 o2s[l] + 4 h[n]';
p_ATPSm = '4 h[n] + 1 adp[l] + 1 pi[l] -> 1 atp[l] + 1 h[l]';
p_ATPtm = '1 atp[l] + 1 adp[x] + 1 pi[x] + 1 h[x] <=> 1 atp[x] + 1 adp[l] + 1 pi[l] + 1 h[l]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ROS Detox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_SOD =   '1 h[m] + 1 o2s[m] -> 1 h2o2[m] + 1 o2[m]';
b_CAT =   '2 h2o2[m] -> 2 h2o[m] + 1 o2[m]';

p_SOD =   '1 h[l] + 1 o2s[l] -> 1 h2o2[l] + 1 o2[l]';
p_CAT =   '2 h2o2[l] -> 2 h2o[l] + 1 o2[l]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Urea Cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_ARGSS = '1 asp_L[x] + 1 atp[x] + 1 citr_L[x] -> 1 amp[x] + 1 argsuc[x] + 1 h[x] + 1 ppi[x]';
p_ARGSL = '1 argsuc[x] <=> 1 arg_L[x] + 1 fum[x]';
p_ARGN  = '1 arg_L[x] + 1 h2o[x] -> 1 orn[x] + 1 urea[x]';
p_ORNtm = '1 h[l] + 1 orn[x] <=> 1 h[x] + 1 orn[l]';
p_OCBTm = '1 cbp[l] + 1 orn[l] -> 1 citr_L[l] + 1 h[l] + 1 pi[l]';
p_CITRtm= '1 citr_L[l] <=> 1 citr_L[x]';
p_CBPSm = '1 atp[l] + 1 hco3[l] + 1 nh4[l] -> 1 adp[l] + 1 cbp[l] + 1 h[l] + 1 pi[l]';

b_ARGSS = '1 asp_L[c] + 1 atp[c] + 1 citr_L[c] -> 1 amp[c] + 1 argsuc[c] + 1 h[c] + 1 ppi[c]';
b_ARGSL = '1 argsuc[c] <=> 1 arg_L[c] + 1 fum[c]';
b_ARGN  = '1 arg_L[c] + 1 h2o[c] -> 1 orn[c] + 1 urea[c]';
b_ORNtm = '1 h[m] + 1 orn[c] <=> 1 h[c] + 1 orn[m]';
b_OCBTm = '1 cbp[m] + 1 orn[m] -> 1 citr_L[m] + 1 h[m]';
b_CITRtm= '1 citr_L[m] <=> 1 citr_L[c]';
b_CBPSm = '1 atp[m] + 1 hco3[m] + 1 nh4[m] -> 1 adp[m] + 1 cbp[m] + 1 h[m] + 1 pi[m]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Malate-Aspartate Shuttle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_MDH    = '1 h[x] + 1 nadh[x] + 1 oaa[x] <=> 1 mal_L[x] + 1 nad[x]';
p_AKGMALm= '1 akg[l] + 1 mal_L[x] <=> 1 akg[x] + 1 mal_L[l]';
p_ASPTA  = '1 akg[x] + 1 asp_L[x] <=> 1 glu_L[x] + 1 oaa[x]';
p_ASPTAm = '1 glu_L[l] + 1 oaa[l] <=> 1 akg[l] + 1 asp_L[l]';
p_ASPGLU = '1 asp_L[l] + 1 glu_L[x] + 1 h[x] <=> 1 asp_L[x] + 1 glu_L[l] + 1 h[l]';

b_MDH    = '1 h[c] + 1 nadh[c] + 1 oaa[c] <=> 1 mal_L[c] + 1 nad[c]';
b_AKGMALm= '1 akg[m] + 1 mal_L[c] <=> 1 akg[c] + 1 mal_L[m]';
b_ASPTA  = '1 akg[c] + 1 asp_L[c] <=> 1 glu_L[c] + 1 oaa[c]';
b_ASPTAm = '1 glu_L[m] + 1 oaa[m] <=> 1 akg[m] + 1 asp_L[m]';
b_ASPGLU = '1 asp_L[m] + 1 glu_L[c] + 1 h[c] <=> 1 asp_L[c] + 0 glu_L[m] + 1 h[m]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Glutathione Metabolism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_Cys    = '1 h[x] + 1 Lcystin[x] + 1 nadh[x] -> 1 cys_L[x] + 1 nad[x]';
p_GLUCYS = '1 atp[x] + 1 cys_L[x] + 1 glu_L[x] -> 1 adp[x] + 1 glucys[x] + 1 h[x] + 1 pi[x]';
p_GTHS   = '1 atp[x] + 1 glucys[x] + 1 gly[x] -> 1 adp[x] + 1 gthrd[x] + 1 h[x] + 1 pi[x]';
p_GTHP   = '1 gthrd[x] + 1 h2o2[x] -> 1 gthox[x] + 1 h2o[x]';
p_GTHO   = '2 gthox[x] + 1 h[x] + 1 nadph[x] -> 2 gthrd[x] + 1 nadp[x]';
p_GTHRDt = '1 atp[x] + 1 gthrd[x] + 1 h2o[x] -> 1 adp[x] + 1 gthrd[l] + 1 h[x] + 1 pi[x]';
p_GTHPm  = '2 gthrd[l] + 1 h2o2[l] -> 2 gthox[l] + 1 h2o[l]'; % Glutathionne Peroxidase, 
p_GTHOm  = '1 gthox[l] + 1 h[l]  -> 1 gthrd[l]';
p_r2520  = '1 akg[x] + 2 gthrd[l] <=> 1 akg[l] + 1 gthrd[x]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% 5-ALA/HEME  Metabolism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_ALASm   = '1 gly[m] + 1 h[m] + 1 succoa[m] -> 1 5aop[m] + 1 co2[m] + 0 coa[m]';
b_AOPtm   = '1 5aop[m] <=> 1 5aop[c]';
b_PPBNGS  = '1 5aop[c] -> 1 h[c] + 1 h2o[c] + 1 ppbng[c]';
b_HMBS    = '1 h2o[c] + 4 ppbng[c] -> 1 hmbil[c] + 1 nh4[c]';
b_UPP3S   = '1 hmbil[c] -> 1 h2o[c] + 1 uppg3[c]';
b_UPPDC1  = '1 h[c] + 1 uppg3[c] -> 4 co2[c] + 1 cpppg3[c]';
b_CPPPGO  = '1 cpppg3[c] + 1 h[c] + 1 o2[c] -> 1 co2[c] + 1 h2o[c] + 1 pppg9[c]';
b_PPPG9tm = '1 pppg9[c] <=> 1 pppg9[m]';
b_PPPGOm  = '1.5 o2[m] + 1 pppg9[m] -> 3 h2o[m] + 1 ppp9[m]';
b_FCLTm   = '1 ppp9[m] -> 1 h[m] + 1 pheme[m]';
b_PHEMEtm = '1 pheme[m] <=> 1 pheme[c]';
b_FCLTc   = '1 ppp9[c] -> 1 h[c] + 1 pheme[c]';
b_PHEMEe  = '1 pheme[c] -> 1 bio[c]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% BMMSC Biomass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_bio = '1 bio[c] ->';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Construct the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Reaction Names:
%%%%%%%%%%%%%%%%%%%%%%%%%%%
reactionIDs = {'EX_glc','EX_h','EX_h2o','EX_lac','EX_co2',... % Extracellular
'EX_o2','EX_gln','EX_hco3','EX_asp','EX_urea','EX_hdca',... % Extracellular
'EX_Lcys','EX_cys','EX_oaa','EX_pyr','EX_glu','EX_pro','EX_val','EX_crn','EX_drib','EX_met',...
'eCA','eSOD',... % Extracellular
'b_MCT1','p_MCT1','b_MCT4','p_MCT4',...% Monocarboxylate transporters
'b_xCT','p_xCT',... % Cystine/Glutamate transporters
'b_aspt','b_Nhe','b_AQP2','b_AQP3','b_o2t','b_co2t',...
'b_ureat','b_palt','b_oaat','b_cyst','b_glnt','b_crnt',...
'b_dribt','b_VALt','b_DM_ATP','b_DM_NADH','b_DM_AMP',...
'b_DM_GMP','p_aspt','p_Nhe','p_AQP2','p_AQP3','p_o2t',...
'p_co2t','p_ureat','p_palt','p_oaat','p_cyst','p_glnt',...
'p_crnt','p_dribt','p_VALt','p_sink_coa','p_DM_ATP','p_DM_NADH',...
'p_DM_AMP','p_DM_GMP',...
'b_AQPtm','b_Nhetm','b_CO2tm','b_O2tm','b_mCA','b_AQP8m',...
'b_OOAtm','b_sinkm_coa','b_VALtm','b_GLUtm','b_ALAtm','b_FUMtm',...
'b_SERtm','b_GLYtm','b_PROtm','p_AQPtm','p_Nhetm','p_CO2tm','p_O2tm',...
'p_mCA','p_OOAtm','p_sinkm_coa','p_VALtm','p_GLUtm','p_ALAtm','p_FUMtm',...
'p_SERtm','p_GLYtm','p_PROtm',...
'b_GLUT','b_HEX','b_PGI','b_PFK','b_FBA','b_TPI','b_GAPD','b_PGK',... %BM Glycolysis
'b_PGM','b_ENO','b_PYK','b_LDH',... % BM Glycolysis
'p_PYRtm','p_PDHm','p_PCm','p_CSm','p_ACONTm','p_ICDHm','p_AKGDm',... % PC TCA Cycle
'p_SUCOASm','p_SUCDm','p_FUMm','p_MDHm',... % PC TCA Cycle
'p_FACOAL','p_CPT1','p_CTRt','p_CPT2','p_FAOX','p_FAOXC',... % PC FA beta-oxidation
'p_Cplx1','p_Cplx2','p_Cplx3','p_Cplx4','p_ATPSm','p_ATPtm','p_SOD','p_CAT',... % PC ETC
'p_GLUT','p_HEX','p_PGI','p_PFK','p_FBA','p_TPI','p_GAPD','p_PGK',... % PC Glycolysis
'p_PGM','p_ENO','p_PYK','p_LDH',... % PC Glycolysis
'p_G6PDH','p_PGL','p_GND','p_RPI','p_PPM','p_PRPP',... % PC PP Pathway
'p_TKT1','p_TALA','p_TKT2',...% PC PP Pathway
'p_ALATA','p_ALA_L',... % PC Alanine Synthesis
'p_PGCD','p_PSERT','p_PSP','p_GHMT',... % PC Serine Synthesis
'p_BioMass','p_r1p',... % PC OBJECTIVE
'p_SERPT','p_DSPHR','p_SLIP',... % PC Serine Derived Lipids
'p_GLUN','p_GLNtm','p_GLUNm','p_GLUDm',... % PC Glutamine Metabolism
'p_NOX',... % PC NADPH oxidase
'p_GLNS','p_GLUPRT','p_PRAGS','p_GARFT','p_PRFGS','p_r0666',... % PC Purine Synth
'p_AIRC','p_PRASCS','p_ADSL2','p_AICART','p_IMPC','p_FTHFL',... % PC Purine Synth
'p_DRBK', 'p_DRPP','p_ADSS', 'p_ADSL', 'p_PUNP0', 'p_PUNP2', 'p_dA',... % PC Nucleotide Interconv.
'p_IMPD', 'p_GMPS', 'p_PUNP3', 'p_PUNP4', 'p_dG',... % PC Nucleotide Interconv.
'p_CA',...
'p_CBPS', 'p_ASPCT', 'p_DHORD', 'p_DHORD9', 'p_ORPT', 'p_OMPDC',... % PC Pyrimidine Synth
'p_URIK', 'p_PYNP', 'p_CSND', 'p_dC',... % PC Pyrimidine Synth
'p_ARGSS', 'p_ARGSL', 'p_ARGN', 'p_ORNtm',... % PC Urea Cycle 
'p_OCBTm', 'p_CITRtm', 'p_CBPSm','p_fumtm','p_fumt',... % PC Urea Cycle
'p_MDH', 'p_AKGMALm', 'p_ASPTA', 'p_ASPTAm', 'p_ASPGLU',... % PC Asp-Mal Shuttle
'p_Cys', 'p_GLUCYS', 'p_GTHS', 'p_GTHP', 'p_GTHO', 'p_GTHRDt', 'p_GTHPm', 'p_GTHOm', 'p_r2520',... % PC Glutathione
'p_VALTAm', 'p_OVID2m', 'p_RE3190M', 'p_MMCDm', 'p_MMEm', 'p_MMMm',... % PC Valine
'b_PYRtm', 'b_PDHm', 'b_PCm', 'b_CSm', 'b_ACONTm', 'b_ICDHm', 'b_AKGDm', 'b_SUCOASm', 'b_SUCDm', 'b_FUMm', 'b_MDHm',... % BM TCA Cycle
'b_Cplx1', 'b_Cplx2', 'b_Cplx3', 'b_Cplx4', 'b_ATPSm', 'b_ATPtm','b_SOD','b_CAT',... % BM ETC 
'b_FACOAL', 'b_CPT1', 'b_CTRt', 'b_CPT2', 'b_FAOX', 'b_FAOXC','b_sink_coa',... % BM FA Oxidation
'b_G6PDH','b_PGL','b_GND','b_RPI','b_PPM','b_PRPP','b_TKT1','b_TALA','b_TKT2',... % BM PP Pathway
'b_PGCD','b_PSERT','b_PSP','b_GHMT',... % BM Serine
'b_GLNS', 'b_GLUPRT', 'b_PRAGS', 'b_GARFT', 'b_PRFGS', 'b_r0666',... % BM Purine
'b_AIRC', 'b_PRASCS', 'b_ADSL2', 'b_AICART', 'b_IMPC', 'b_FTHFL',... % BM Purine
'b_DRBK', 'b_DRPP', 'b_ADSS', 'b_ADSL', 'b_PUNP0', 'b_PUNP2', 'b_dA',... % BM Nucleotide interconv
'b_IMPD', 'b_GMPS', 'b_PUNP3', 'b_PUNP4', 'b_dG',... % BM Nucleotide interconv
'b_NOX','b_CA','b_ALATA','b_ALA_L',... % BM Alanine
'b_ARGSS', 'b_ARGSL', 'b_ARGN', 'b_ORNtm', 'b_OCBTm', 'b_CITRtm', 'b_CBPSm',... % Urea Cycle
'b_GLUN','b_GLNtm','b_GLUNm','b_GLUDm',... % BM Glutamine Metabolism
'b_MDH','b_AKGMALm','b_ASPTA','b_ASPTAm','b_ASPGLU',... % BM Asp-Mal Shuttle
'b_CBPS','b_ASPCT','b_DHORD','b_DHORD9','b_ORPT','b_OMPDC',... % BM Pyrimidine
'b_URIK','b_PYNP','b_CSND','b_dC',...% BM Pyrimidine
'b_ADOMET', 'b_METS', 'b_MS', 'b_HMR4710', 'b_CYSTS','b_r0210', 'b_CYSt',... % BM Metathione
'b_ALASm', 'b_AOPtm', 'b_PPBNGS', 'b_HMBS', 'b_UPP3S', 'b_UPPDC1',... % BM Heme 5-ALA
'b_CPPPGO', 'b_PPPG9tm', 'b_PPPGOm', 'b_FCLTm', 'b_PHEMEtm', 'b_FCLTc',... % BM Heme 5-ALA
'b_PHEMEe','b_FUMt','b_bio'... % BM Heme 5-ALA
}; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These ought to be modified in the future, as more data is available
%%%%%%%%%%%%%%%%%%%%%%%%%%%
reactionNames = reactionIDs;
grRule = reactionIDs;
geneNames = reactionIDs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reaction Order
%%%%%%%%%%%%%%%%%%%%%%%%%%%
reactionForm = {EX_glc,EX_h,EX_h2o,EX_lac,EX_co2,... % Extracellular
EX_o2,EX_gln,EX_hco3,EX_asp,EX_urea,EX_hdca,... % Extracellular
EX_Lcys,EX_cys,EX_oaa,EX_pyr,EX_glu,EX_pro,EX_val,EX_crn,EX_drib,EX_met,...
eCA,eSOD,... % Extracellular
b_MCT1,p_MCT1,b_MCT4,p_MCT4,...
b_xCT,p_xCT,b_aspt,b_Nhe,b_AQP2,b_AQP3,b_o2t,b_co2t,...
b_ureat,b_palt,b_oaat,b_cyst,b_glnt,b_crnt,b_dribt,b_VALt,...
b_DM_ATP,b_DM_NADH,b_DM_AMP,b_DM_GMP,p_aspt,p_Nhe,p_AQP2,...
p_AQP3,p_o2t,p_co2t,p_ureat,p_palt,p_oaat,p_cyst,p_glnt,...
p_crnt,p_dribt,p_VALt,p_sink_coa,p_DM_ATP,p_DM_NADH,p_DM_AMP,p_DM_GMP,...
b_AQPtm,b_Nhetm,b_CO2tm,b_O2tm,b_mCA,b_AQP8m,b_OOAtm,b_sinkm_coa,...
b_VALtm,b_GLUtm,b_ALAtm,b_FUMtm,b_SERtm,b_GLYtm,b_PROtm,p_AQPtm,...
p_Nhetm,p_CO2tm,p_O2tm ,p_mCA,p_OOAtm,p_sinkm_coa,p_VALtm,p_GLUtm,...
p_ALAtm,p_FUMtm,p_SERtm,p_GLYtm,p_PROtm,...
b_GLUT,b_HEX,b_PGI,b_PFK,b_FBA,b_TPI,b_GAPD,b_PGK,b_PGM,b_ENO,b_PYK,b_LDH,...
p_PYRtm,p_PDHm,p_PCm,p_CSm,p_ACONTm,p_ICDHm,p_AKGDm,p_SUCOASm,p_SUCDm,...
p_FUMm,p_MDHm,...
p_FACOAL,p_CPT1,p_CTRt,p_CPT2,p_FAOX,p_FAOXC,...
p_Cplx1,p_Cplx2,p_Cplx3,p_Cplx4,p_ATPSm,p_ATPtm,p_SOD,p_CAT,...
p_GLUT,p_HEX,p_PGI,p_PFK,p_FBA,p_TPI,p_GAPD,p_PGK,p_PGM,p_ENO,p_PYK,p_LDH,...
p_G6PDH,p_PGL,p_GND,p_RPI,p_PPM,p_PRPP,p_TKT1,p_TALA,p_TKT2,...
p_ALATA,p_ALA_L,...
p_PGCD,p_PSERT,p_PSP,p_GHMT,...
p_BioMass,p_r1p,... % PC OBJECTIVE
p_SERPT,p_DSPHR,p_SLIP,... 
p_GLUN,p_GLNtm,p_GLUNm,p_GLUDm,...
p_NOX,...
p_GLNS,p_GLUPRT,p_PRAGS,p_GARFT,p_PRFGS,p_r0666,p_AIRC,p_PRASCS,...
p_ADSL2,p_AICART,p_IMPC,p_FTHFL,...
p_DRBK, p_DRPP, p_ADSS, p_ADSL, p_PUNP0, p_PUNP2, p_dA,...
p_IMPD, p_GMPS, p_PUNP3, p_PUNP4, p_dG,...
p_CA,...
p_CBPS,p_ASPCT,p_DHORD,p_DHORD9,p_ORPT,p_OMPDC,p_URIK,p_PYNP,p_CSND,p_dC,...
p_ARGSS, p_ARGSL, p_ARGN, p_ORNtm, p_OCBTm, p_CITRtm, p_CBPSm,p_fumtm,p_fumt...
p_MDH, p_AKGMALm, p_ASPTA, p_ASPTAm, p_ASPGLU,...
p_Cys, p_GLUCYS, p_GTHS, p_GTHP, p_GTHO, p_GTHRDt, p_GTHPm, p_GTHOm, p_r2520,...
p_VALTAm, p_OVID2m, p_RE3190M, p_MMCDm, p_MMEm, p_MMMm,...
b_PYRtm, b_PDHm, b_PCm, b_CSm, b_ACONTm, b_ICDHm, b_AKGDm, b_SUCOASm, b_SUCDm, b_FUMm, b_MDHm,...
b_Cplx1, b_Cplx2, b_Cplx3, b_Cplx4, b_ATPSm, b_ATPtm,b_SOD,b_CAT,...
b_FACOAL,b_CPT1,b_CTRt,b_CPT2,b_FAOX,b_FAOXC,b_sink_coa,...
b_G6PDH,b_PGL,b_GND,b_RPI,b_PPM,b_PRPP,b_TKT1,b_TALA,b_TKT2,...
b_PGCD,b_PSERT,b_PSP,b_GHMT,...
b_GLNS,b_GLUPRT,b_PRAGS,b_GARFT,b_PRFGS,b_r0666,b_AIRC,b_PRASCS,...
b_ADSL2,b_AICART,b_IMPC,b_FTHFL,...
b_DRBK,b_DRPP,b_ADSS,b_ADSL,b_PUNP0,b_PUNP2,b_dA,...
b_IMPD,b_GMPS,b_PUNP3,b_PUNP4,b_dG,...
b_NOX,b_CA,b_ALATA,b_ALA_L,...
b_ARGSS,b_ARGSL,b_ARGN,b_ORNtm,b_OCBTm,b_CITRtm,b_CBPSm,...
b_GLUN,b_GLNtm,b_GLUNm,b_GLUDm,...
b_MDH,b_AKGMALm,b_ASPTA,b_ASPTAm,b_ASPGLU,...
b_CBPS,b_ASPCT,b_DHORD,b_DHORD9,b_ORPT,b_OMPDC,b_URIK,b_PYNP,b_CSND,b_dC,...
b_ADOMET, b_METS, b_MS, b_HMR4710, b_CYSTS, b_r0210, b_CYSt,...
b_ALASm, b_AOPtm, b_PPBNGS, b_HMBS, b_UPP3S, b_UPPDC1, b_CPPPGO, b_PPPG9tm,...
b_PPPGOm, b_FCLTm, b_PHEMEtm, b_FCLTc, b_PHEMEe,b_FUMt,b_bio,...
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is the reaction reversible? (1 = yes, 0 = no)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
rev = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,...
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,...
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,...
1,0,0,0,0,0,0,0,0,0,0,1,... % BM Glycolysis
0,0,0,0,0,0,0,0,0,0,0,... % PC TCA Cycle
0,1,1,1,0,0,... % PC FA beta oxidation
0,0,0,0,0,0,0,0,... % PC ETC
1,0,0,0,0,0,0,0,0,0,0,0,... % PC Glycolysis
0,0,0,0,0,0,0,0,0,... % PC PP Pathway
0,0,... % PC Alanine Synth
0,0,0,0,... % PC Serine Synth
0,0,... % PC Objective
0,0,0,... % PC Sphingolipids
0,0,0,0,... % PC Glutamine Metabolism
1,... % PC NADH Oxidoreductase
1,0,1,1,0,1,1,1,0,1,1,1,... % PC Purine Synth
0,0,0,0,0,0,... % PC Nucleotide Interconv
0,0,0,0,0,0,... % PC Nucleotide Interconv
1,...
0,0,0,0,0,0,0,0,0,0,... % PC Pyrimidine Synth
0,0,0,0,0,0,0,... % PC Urea Cycle
0,0,0,0,0,... % PC Asp-Mal-Shuttle
0,0,0,0,0,0,0,0,0,... % PC Glutathionne
0,0,0,0,0,0,... % PC Valinne
0,0,0,0,0,0,0,0,0,0,0,... % BM TCA Cycle
0,0,0,0,0,0,0,0,0,0,0,... % BM ETC ROS
0,0,0,0,0,0,1,... % BM FA Oxidation
0,0,0,0,0,0,0,0,0,... % BM PP Pathway
0,0,0,0,... % BM Serine
0,0,0,0,0,0,0,0,0,0,0,0,... % BM Purine & Nucleotide Int.
0,0,0,0,0,0,0,0,0,0,0,0,... % BM Purine & Nucleotide Int.
1,1,0,0,...
0,0,0,0,0,0,0,0,... % BM Urea Cycle
0,0,0,0,... % BM Glutamine Metabolism
0,0,0,0,0,... % BM Mal-Asp Shuttle
0,0,0,0,0,0,0,0,0,0,... % BM Pyrimidine
0,0,0,0,0,0,0,... % BM Metathione
0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,... % BM Heme 5-ALA 
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimisation's Lower Bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb = [-1000 % -> EX_glc
0 % -> EX_h
-1000 % -> EX_h2o
0 % -> EX_lac
-1000 % -> EX_co2
-1000 % -> EX_o2
-1000 % -> EX_gln
0 % -> EX_hco3
0 % -> EX_asp
0 % -> EX_urea
-1000 % -> EX_hdca
-1000 % -> EX_Lcys
0 % -> EX_cys
0 % -> EX_oaa
0 % -> EX_pyr
0 % -> EX_glu
0 % -> EX_pro
-1000 % -> EX_val
-1000 % -> EX_crn
-1000 % -> EX_Drib
0     % -> EX_met
0 % -> eCA
0 % -> eSOD
0 % -> b_MCT1 %%%%%%%%%%%%%%
0 % -> p_MCT1 % MCT's
0 % -> b_MCT4  MCT's
-0.1 % -> p_MCT4 %%%%%%%%%%%%%%
0 % -> b_xCT
0 % -> p_xCT
-1000 % ->  b_aspt
0 % ->  b_Nhe
0 % ->  b_AQP2
0 % ->  b_AQP3
0 % ->  b_o2t
-1000 % ->  b_co2t
0 % ->  b_ureat
0 % ->  b_palt
0 % ->  b_oaat
0 % ->  b_cyst
0 % ->  b_glnt
0 % ->  b_crnt
-1000 % ->  b_dribt
0 % ->  b_VALt
-1000 % -> p_sink_coa
0 % ->  b_DM_ATP
0 % ->  b_DM_NADH
0 % ->  b_DM_AMP
0 % ->  b_DM_GMP
0 % ->  p_aspt
-1000 % ->  p_Nhe
-1000 % ->  p_AQP2
-1000 % ->  p_AQP3
0 % ->  p_o2t
-1000 % ->  p_co2t
0 % ->  p_ureat
-1000 % ->  p_palt
-1000 % ->  p_oaat
0 % ->  p_cyst
0 % ->  p_glnt
0 % ->  p_crnt
0 % ->  p_dribt
0 % ->  p_VALt
0 % ->  p_DM_ATP
0 % ->  p_DM_NADH
0 % ->  p_DM_AMP
0 % ->  p_DM_GMP
0 % -> b_AQPtm
-1000 % -> b_Nhetm
-1000 % -> b_CO2tm
0 % -> b_O2tm
0 % -> b_mCA
-1000 % -> b_AQP8m
-1000 % -> b_OOAtm
-1000 % -> b_sinkm_coa
0 % -> b_VALtm
0 % -> b_GLUtm
0 % -> b_ALAtm
0 % -> b_FUMtm
0 % -> b_SERtm
0 % -> b_GLYtm
0 % -> b_PROtm
-1000 % -> p_AQPtm
-1000 % -> p_Nhetm
-1000 % -> p_CO2tm
15.5 % -> p_O2tm 
0 % -> p_mCA
-1000 % -> p_OOAtm
-1000 % -> p_sinkm_coa
0 % -> p_VALtm
0 % -> p_GLUtm
0 % -> p_ALAtm
0 % -> p_FUMtm
0 % -> p_SERtm
0 % -> p_GLYtm
0 % -> p_PROtm
1 % -> b_GLUT %%%%%%%%%%%%%% 
0 % -> b_HEX  BM GLycolysis
0 % -> b_PGI  BM GLycolysis
0 % -> b_PFK  BM GLycolysis
0 % -> b_FBA  BM GLycolysis
0 % -> b_TPI  BM GLycolysis
0 % -> b_GAPD BM GLycolysis
0 % -> b_PGK  BM GLycolysis
0 % -> b_PGM  BM GLycolysis
0 % -> b_ENO  BM GLycolysis
0 % -> b_PYK  BM GLycolysis
0 % -> b_LDH  %%%%%%%%%%%%%%
0 % -> p_PYRtm PC TCA Cycle
0 % -> p_PDHm  PC TCA Cycle
0.1 % -> p_PCm   PC TCA Cycle
0 % -> p_CSm   PC TCA Cycle
0 % -> p_ACONTm PC TCA Cycle
0 % -> p_ICDHm PC TCA Cycle
0 % -> p_AKGDm PC TCA Cycle
0 % -> p_SUCOASm PC TCA Cycle
0 % -> p_SUCDm PC TCA Cycle
0.1 % -> p_FUMm  PC TCA Cycle
0 % -> p_MDHm %%%%%%%%%%%%%%
0.1 % -> p_FACOAL PC FA Oxidation
0.1 % -> p_CPT1   PC FA Oxidation
0 % -> p_CTRt   PC FA Oxidation
0 % -> p_CPT2   PC FA Oxidation
0 % -> p_FAOX   PC FA Oxidation
0 % -> p_FAOXC %%%%%%%%%%%%%%
0 % -> p_Cplx1  PC ETC
0 % -> p_Cplx2  PC ETC
0 % -> p_Cplx3  PC ETC
0 % -> p_Cplx4  PC ETC
0 % -> p_ATPSm  PC ETC
0 % -> p_ATPtm  PC ETC
0 % -> p_SOD    ROS DETOX
0 % -> p_CAT  %%%%%%%%%%%%%%
0 % -> p_GLUT PC Glycolysis
0 % -> p_HEX  PC Glycolysis
0 % -> p_PGI  PC Glycolysis
0.1 % -> p_PFK  PC Glycolysis
0 % -> p_FBA  PC Glycolysis
0 % -> p_TPI  PC Glycolysis
0 % -> p_GAPD PC Glycolysis
0 % -> p_PGK  PC Glycolysis
0 % -> p_PGM  PC Glycolysis
0 % -> p_ENO  PC Glycolysis
0 % -> p_PYK  PC Glycolysis
0 % -> p_LDH  %%%%%%%%%%%%%%
1 % -> p_G6PDH PC PP Pathway
0 % -> p_PGL   PC PP Pathway
0 % -> p_GND   PC PP Pathway
0 % -> p_RPI   PC PP Pathway
0 % -> p_PPM   PC PP Pathway
0 % -> p_PRPP  PC PP Pathway
0 % -> p_TKT1  PC PP Pathway
0 % -> p_TALA  PC PP Pathway
0 % -> p_TKT2 %%%%%%%%%%%%%%
1 % -> p_ALATA PC Alanine Pathway
0 % -> p_ALA_L %%%%%%%%%%%%%%
0 % -> p_PGCD  PC Serine Pathway
0 % -> p_PSERT PC Serine Pathway
0 % -> p_PSP   PC Serine Pathway
0 % -> p_GHMT  %%%%%%%%%%%%%%
0 % -> BIOMASS OBJECTIVE PLASMA CELL
-1000 % -> p_r1p
0 % -> p_SERPT PC Serine Lipds
0 % -> p_DSPHR PC Serine Lipds
0 % -> p_SLIP %%%%%%%%%%%%%%
0 % -> p_GLUN  PC Glutmine Met
0 % -> p_GLNtm PC Glutmine Met
0 % -> p_GLUNm PC Glutmine Met
0 % -> p_GLUDm %%%%%%%%%%%%%%
0 % -> p_NOX  PC NADPH
0 % -> p_GLNS  %%%%%%%%%%%%%%
1 % -> p_GLUPRT PC Purine Synth
0 % -> p_PRAGS  PC Purine Synth
0 % -> p_GARFT  PC Purine Synth
0 % -> p_PRFGS  PC Purine Synth
0 % -> p_r0666  PC Purine Synth
0 % -> p_AIRC   PC Purine Synth
0 % -> p_PRASCS PC Purine Synth
0 % -> p_ADSL2  PC Purine Synth
0 % -> p_AICART PC Purine Synth
0 % -> p_IMPC  PC Purine Synth
0 % -> p_FTHFL %%%%%%%%%%%%%%
0 % -> p_DRBK  PC Nucleotide Interconv
0 % -> p_DRPP  PC Nucleotide Interconv
0.5 % -> p_ADSS  PC Nucleotide Interconv
0 % -> p_ADSL  PC Nucleotide Interconv
0 % -> p_PUNP0 PC Nucleotide Interconv
0 % -> p_PUNP2 PC Nucleotide Interconv
0 % -> p_dA    PC Nucleotide Interconv
0.5 % -> p_IMPD  PC Nucleotide Interconv
0 % -> p_GMPS  PC Nucleotide Interconv
0 % -> p_PUNP3 PC Nucleotide Interconv
0 % -> p_PUNP4 PC Nucleotide Interconv
0 % -> p_dG    %%%%%%%%%%%%%%
-1000 % -> p_CA
0.1 % -> p_CBPS %%%%%%%%%%%%%%
0 % -> p_ASPCT PC Pyrimidine Synth
0 % -> p_DHORD PC Pyrimidine Synth
0 % -> p_DHORD9 PC Pyrimidine Synth
0 % -> p_ORPT  PC Pyrimidine Synth
0 % -> p_OMPDC PC Pyrimidine Synth
0 % -> p_URIK  PC Pyrimidine Synth
0 % -> p_PYNP  PC Pyrimidine Synth
0 % -> p_CSND  PC Pyrimidine Synth
0 % -> p_dC %%%%%%%%%%%%%%
0.1 % -> p_ARGSS PC Urea
0 % -> p_ARGSL PC Urea
0 % -> p_ARGN  PC Urea
0 % -> p_ORNtm PC Urea
0 % -> p_OCBTm PC Urea
0 % -> p_CITRtm PC Urea
0 % -> p_CBPSm %%%%%%%%%%%%%%
0 % -> p_fumtm
0 % -> p_fumt %%%%%%%%%%%%%%
0 % -> p_MDH      PC Asp-Mal Shuttle
0 % -> p_AKGMALm  PC Asp-Mal Shuttle
0 % -> p_ASPTA    PC Asp-Mal Shuttle
0 % -> p_ASPTAm   PC Asp-Mal Shuttle
0 % -> p_ASPGLU %%%%%%%%%%%%%%
0.1 % - > p_Cys    PC Glutathione
0 % - > p_GLUCYS PC Glutathione
0 % - > p_GTHS   PC Glutathione
0 % - > p_GTHP   PC Glutathione
0 % - > p_GTHO   PC Glutathione
0 % - > p_GTHRDt PC Glutathione
0 % - > p_GTHPm  PC Glutathione
0.1 % - > p_GTHOm  PC Glutathione
0 % - > p_r2520  %%%%%%%%%%%%%%
0.1 % -> p_VALTAm  PC VALINE
0 % -> p_OVID2m  PC VALINE
0 % -> p_RE3190M PC VALINE
0 % -> p_MMCDm   PC VALINE
0 % -> p_MMEm    PC VALINE
0 % -> p_MMMm   %%%%%%%%%%%%%%
0 % -> b_PYRtm  BM TCA Cycle
0 % -> b_PDHm   BM TCA Cycle
0 % -> b_PCm    BM TCA Cycle
0 % -> b_CSm    BM TCA Cycle
0 % -> b_ACONTm BM TCA Cycle
0 % -> b_ICDHm  BM TCA Cycle
0 % -> b_AKGDm  BM TCA Cycle
0 % -> b_SUCOASm  BM TCA Cycle
0 % -> b_SUCDm  BM TCA Cycle
0 % -> b_FUMm   BM TCA Cycle
0 % -> b_MDHm   %%%%%%%%%%%%%%
0 % -> b_Cplx1  BM ETC 
0 % -> b_Cplx2  BM ETC
0 % -> b_Cplx3  BM ETC
0 % -> b_Cplx4  BM ETC
0 % -> b_ATPSm  BM ETC
0 % -> b_ATPtm  BM ETC
0 % -> b_SOD    BM ROS
0 % -> b_CAT    %%%%%%%%%%%%%%
0.1 % -> b_FACOAL   BM FA Oxidation
0 % -> b_CPT1   BM FA Oxidation
0 % -> b_CTRt   BM FA Oxidation
0 % -> b_CPT2   BM FA Oxidation
0 % -> b_FAOX   BM FA Oxidation
0 % -> b_FAOXC   BM FA Oxidation
-1000 % -> b_sink_coa
0.1 % - > b_G6PDH %%%%%%%%%%%%%%
0 % - > b_PGL  BM PP Pathway
0 % - > b_GND  BM PP Pathway
0 % - > b_RPI  BM PP Pathway
0 % - > b_PPM  BM PP Pathway
0 % - > b_PRPP BM PP Pathway
0 % - > b_TKT1 BM PP Pathway
0 % - > b_TALA BM PP Pathway
0 % - > b_TKT20 %%%%%%%%%%%%%%
0 % -> b_PGCD   BM Serine
0 % -> b_PSERT  BM Serine
0 % -> b_PSP    BM Serine
0 % -> b_GHMT  %%%%%%%%%%%%%%
0 % -> b_GLNS   BM Purine & Nucleotide Int.
0.1 % -> b_GLUPRT BM Purine & Nucleotide Int.
0 % -> b_PRAGS  BM Purine & Nucleotide Int.
0 % -> b_GARFT  BM Purine & Nucleotide Int.
0 % -> b_PRFGS  BM Purine & Nucleotide Int.
0 % -> b_r0666  BM Purine & Nucleotide Int.
0 % -> b_AIRC   BM Purine & Nucleotide Int.
0 % -> b_PRASCS BM Purine & Nucleotide Int.
0 % -> b_ADSL2  BM Purine & Nucleotide Int.
0 % -> b_AICART BM Purine & Nucleotide Int.
0 % -> b_IMPC   BM Purine & Nucleotide Int.
0 % -> b_FTHFL  BM Purine & Nucleotide Int.
0 % -> b_DRBK   BM Purine & Nucleotide Int.
0 % -> b_DRPP   BM Purine & Nucleotide Int.
0.1 % -> b_ADSS   BM Purine & Nucleotide Int.
0 % -> b_ADSL   BM Purine & Nucleotide Int.
0 % -> b_PUNP0  BM Purine & Nucleotide Int.
0 % -> b_PUNP2  BM Purine & Nucleotide Int.
0 % -> b_dA     BM Purine & Nucleotide Int.
0.1 % -> b_IMPD BM Purine & Nucleotide Int.
0 % -> b_GMPS   BM Purine & Nucleotide Int.
0 % -> b_PUNP3  BM Purine & Nucleotide Int.
0 % -> b_PUNP4  BM Purine & Nucleotide Int.
0 % -> b_dG    %%%%%%%%%%%%%%
0 % -> b_NOX    BM NADPH Oxidase
0 % -> b_CA     BM Carbonic Anhydrases
0 % -> b_ALATA  BM Alanine
0 % -> b_ALA    BM Alanine
0.1 % -> b_ARGSS %%%%%%%%%%%%%%
0 % -> b_ARGSL BM Urea Cycle
0 % -> b_ARGN BM Urea Cycle
0 % -> b_ORNtm BM Urea Cycle
0 % -> b_OCBTm BM Urea Cycle
0 % -> b_CITRtm BM Urea Cycle
0.1 % -> b_CBPSm %%%%%%%%%%%%%%
0 % -> b_GLUN BM Glutamine Metabolism
0 % -> b_GLNtm BM Glutamine Metabolism
0 % -> b_GLUNm BM Glutamine Metabolism
0 % -> b_GLUDm BM Glutamine Metabolism
0 % -> b_MDH  %%%%%%%%%%%%%%
0 % -> b_AKGMALm  BM Asp-Mal Shuttle
0 % -> b_ASPTA  BM Asp-Mal Shuttle
0 % -> b_ASPTAm  BM Asp-Mal Shuttle
0 % -> b_ASPGLU  BM Asp-Mal Shuttle
0.1 % -> b_CBPS  %%%%%%%%%%%%%%
0 % -> b_ASPCT  BM Pyrimidines
0 % -> b_DHORD  BM Pyrimidines
0 % -> b_DHORD9 BM Pyrimidines
0 % -> b_ORPT  BM Pyrimidines
0 % -> b_OMPDC BM Pyrimidines
0 % -> b_URIK  BM Pyrimidines
0 % -> b_PYNP  BM Pyrimidines
0 % -> b_CSND  BM Pyrimidines
0 % -> b_dC    %%%%%%%%%%%%%%
0.1 % -> b_ADOMET  BM Metathione
0 % -> b_METS  	 BM Metathione
0 % -> b_MS   	 BM Metathione
0 % -> b_HMR4710 BM Metathione
0 % -> b_CYSTS   BM Metathione
0 % -> b_r0210   BM Metathione
0 % -> b_CYSt    BM Metathione
0.1 % -> b_ALASm %%%%%%%%%%%%%%
0 % -> b_AOPtm   BM Heme 5-ALA 
0 % -> b_PPBNGS  BM Heme 5-ALA 
0 % -> b_HMBS    BM Heme 5-ALA 
0 % -> b_UPP3S   BM Heme 5-ALA 
0 % -> b_UPPDC1  BM Heme 5-ALA 
0 % -> b_CPPPGO  BM Heme 5-ALA 
0 % -> b_PPPG9tm BM Heme 5-ALA 
0 % -> b_PPPGOm  BM Heme 5-ALA 
0 % -> b_FCLTm   BM Heme 5-ALA 
0 % -> b_PHEMEtm BM Heme 5-ALA 
0 % -> b_FCLTc   BM Heme 5-ALA 
0 % -> b_PHEMEe  BM Heme 5-ALA
0 % -> Fumt
0 % -> b_bio
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimisation's Upper Bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%
ub = ones(length(lb),1)*1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the model object
%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = createModel(reactionIDs,reactionNames,reactionForm);
model.lb = lb;
model.ub = ub;
model.grRules = grRule';
model.genes = geneNames';
model.rev = rev;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change specific bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not accumulate pyruvate in the extracellular compartment:
model = changeRxnBounds(model, 'EX_pyr', 0, 'l');
model = changeRxnBounds(model, 'EX_pyr', 0, 'u');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a test function, here I am disabling it:
model = changeRxnBounds(model, 'p_DM_NADH', 0, 'u');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not accumulate cysteine in the extracellular compartment:
model = changeRxnBounds(model, 'EX_cys', 0, 'u');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hypoxia (Px O2 ~ 2%)
model = changeRxnBounds(model, 'p_o2t', 15.5, 'u');
model = changeRxnBounds(model, 'b_o2t', 15.5, 'u');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Bounds
model = changeRxnBounds(model, 'p_PRAGS', 1, 'u'); 
model = changeRxnBounds(model, 'p_PGM', 2.7, 'u');
model = changeRxnBounds(model, 'p_PGK', 3.90, 'u');
model = changeRxnBounds(model, 'p_MCT1', 20.526, 'u');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HIF1 This is an approximation and will change when the model is guided
% by genomic data.
model = changeRxnBounds(model, 'p_PDHm', 1, 'l');
model = changeRxnBounds(model, 'p_PDHm', 1, 'u');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the model's optimisation objectives
model.c=model.c;
model = changeObjective(model,{'p_BioMass','b_PYK'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimise
%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = {optimizeCbModel(model,'max'),optimizeCbModel(model,'min')};

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fluxes = [F{1}.v,F{2}.v];
Reactions = model.rxns;
clc
T = table(Reactions,Fluxes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Script
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Create an SBML file for visualisation and knockout simulation purposes
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The resulting file can be then uploaded to fluxer, where the simulations
% are performed best.
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This line of code will be added, not removed, once we fit genomic data:
model=rmfield(model,'rxnGeneMat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This line of code is to avoid execution of example in non 
% gui-environments:
if usejava('desktop')
writeCbModel(model, 'fileName', 'MM_Test.sbml','format','sbml')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%



