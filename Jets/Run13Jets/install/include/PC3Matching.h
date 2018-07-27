#ifndef __PC3MATCHING_H__
#define __PC3MATCHING_H__

//#include "CommonVar.h"
//#include "EmcMatching.h"


//enum type1_enum{DPHI, DZ};
//enum type2_enum{INITIAL, INTERMEDIATE, FINAL};

const static int NCHARGE_PC3 = 2;
const static int NARMS_PC3 = 2;
const static int NZED_PC3 = 10;
const static int NCENT_PC3 = 1;
const static int NARMSECTS_PC3 = 8;
const static int N_HIGH_PT_TEST_BINS_PC3 = 0;
const static int NPT_PC3 = 13;


int GetPtBinPC3(float pT);
float GetPtBinCenterPC3(int pt_bin);
void LoadAllPC3Files();
void LoadPointsPC3();

float CalcMeanPC3(float pt, int icharge, int iarmsect, int ized, int icent, int STAGE, int DPHI_DZ);
float CalcSigmaPC3(float pt, int icharge, int iarmsect, int ized, int icent, int STAGE, int DPHI_DZ);
float CalcPC3(float pt, int icharge, int iarmsect, int ized, int icent, float dphi_dz, int STAGE, int DPHI_DZ);
  
static float mean_points_pc3[NCHARGE_PC3][NARMS_PC3][NZED_PC3][NCENT_PC3][NPT_PC3][3][2];
static float sigma_points_pc3[NCHARGE_PC3][NARMS_PC3][NZED_PC3][NCENT_PC3][NPT_PC3][3][2];

#endif // __PC3MATCHING_H__
