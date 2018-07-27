#ifndef __EMCMATCHING_H__
#define __EMCMATCHING_H__


enum type1_enum{DPHI, DZ};
enum type2_enum{INITIAL, INTERMEDIATE, FINAL};

const static int NCHARGE = 2;
const static int NARMS = 2;
const static int NZED = 10;
const static int NCENT = 1;
const static int NARMSECTS = 8;
const static int N_HIGH_PT_TEST_BINS = 0;
const static int NPT = 13;


int GetPtBin(float pT);
float GetPtBinCenter(int pt_bin);
void LoadMean();
void LoadSigma();
void LoadAllEMCFiles();
void LoadPoints();

float CalcMean(float pt, int icharge, int iarmsect, int ized, int icent, int STAGE, int DPHI_DZ);
float CalcSigma(float pt, int icharge, int iarmsect, int ized, int icent, int STAGE, int DPHI_DZ);
float CalcEMC(float pt, int icharge, int iarmsect, int ized, int icent, float dphi_dz, int STAGE, int DPHI_DZ);

	 
static float mean_points_emc[NCHARGE][NARMSECTS][NZED][NCENT][NPT][3][2];
static float sigma_points_emc[NCHARGE][NARMSECTS][NZED][NCENT][NPT][3][2];
static float mean0[NARMSECTS][NPT];
static float mean1[NARMSECTS][NPT];
static float mean2[NARMSECTS][NPT];
static float sigma0[NARMSECTS][NPT];
static float sigma1[NARMSECTS][NPT];
static float sigma2[NARMSECTS][NPT];



#endif //__EMCMATCHING_H__
