#ifndef __PC3MATCHING_H__
#define __PC3MATCHING_H__

static const int pc3_zed_bins  = 10;
static const int pc3_pt_bins   = 15;
static const int pc3_cent_bins = 10;
static const int pc3_num_det   = 8;
static const int pc3_func_deg_mean_I   = 7;
static const int pc3_func_deg_sigma_I  = 8;
static const int pc3_func_deg_mean_II  = 7;
static const int pc3_func_deg_sigma_II = 7;

static double pc3_zed_mean_I[pc3_cent_bins][pc3_num_det][pc3_zed_bins][pc3_func_deg_mean_I];
static double pc3_zed_mean_II[pc3_cent_bins][pc3_num_det][pc3_zed_bins][pc3_func_deg_mean_II];
static double pc3_zed_sigma_I[pc3_cent_bins][pc3_num_det][pc3_zed_bins][pc3_func_deg_sigma_I];
static double pc3_zed_sigma_II[pc3_cent_bins][pc3_num_det][pc3_zed_bins][pc3_func_deg_sigma_II];
static double point_lines_I[pc3_cent_bins][pc3_num_det][pc3_zed_bins][pc3_pt_bins-1][4];
static double point_lines_II[pc3_cent_bins][pc3_num_det][pc3_zed_bins][pc3_pt_bins-1][4];

void  pc3_init_fit_pars_I();
void  pc3_init_fit_pars_II();

float pc3_sdphi_func_I(float charge, short dcarm, int iCent, int iZed, float pt, float pc3dphi);
float pc3_sdphi_func_II(float charge, short dcarm, int iCent, int iZed, float pt, float pc3dphi);
float pc3_sdz_func_I(float charge, short dcarm, int iCent, int iZed, float pt, float pc3dz);
float pc3_sdz_func_II(float charge, short dcarm, int iCent, int iZed, float pt, float pc3dz);

float pc3_eval_mean_I(float pt, int iDet, int iCent, int iZed);
float pc3_eval_mean_II(float pt, int iDet, int iCent, int iZed);
float pc3_eval_sigma_I(float pt, int iDet, int iCent, int iZed);
float pc3_eval_sigma_II(float pt, int iDet, int iCent, int iZed);


#endif /* __PC3_MATCHING_H__ */



