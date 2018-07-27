#ifndef __EMCMATCHING_H__
#define __EMCMATCHING_H__

#include <iostream>
#include <fstream>

static const int NCHARGE = 2;
static const int NSECT = 8;
static const int NZED = 10;
static const int NCENT = 6;
static const int NPT = 14;

static float initialdPhiMeanParam_0[NCHARGE][NSECT][NZED][NCENT];
static float initialdPhiMeanParam_1[NCHARGE][NSECT][NZED][NCENT];
static float initialdPhiMeanParam_2[NCHARGE][NSECT][NZED][NCENT];

static float initialdPhiSigmaParam_0[NCHARGE][NSECT][NZED][NCENT];
static float initialdPhiSigmaParam_1[NCHARGE][NSECT][NZED][NCENT];
static float initialdPhiSigmaParam_2[NCHARGE][NSECT][NZED][NCENT];

static float finaldPhiMeanParam_0[NCHARGE][NSECT][NZED][NCENT];
static float finaldPhiMeanParam_1[NCHARGE][NSECT][NZED][NCENT];
static float finaldPhiMeanParam_2[NCHARGE][NSECT][NZED][NCENT];

static float finaldPhiSigmaParam_0[NCHARGE][NSECT][NZED][NCENT];

static float mean0[NSECT][NPT];
static float mean1[NSECT][NPT];
static float mean2[NSECT][NPT];

static float sigma0[NSECT][NPT];
static float sigma1[NSECT][NPT];
static float sigma2[NSECT][NPT];

static float initialdZMeanParam_0[NCHARGE][NSECT][NZED][NCENT];

static float initialdZSigmaParam_0[NCHARGE][NSECT][NZED][NCENT];
static float initialdZSigmaParam_1[NCHARGE][NSECT][NZED][NCENT];
static float initialdZSigmaParam_2[NCHARGE][NSECT][NZED][NCENT];

static float finaldZMeanParam_0[NCHARGE][NSECT][NZED][NCENT];

static float finaldZSigmaParam_0[NCHARGE][NSECT][NZED][NCENT];
static float finaldZSigmaParam_1[NCHARGE][NSECT][NZED][NCENT];
static float finaldZSigmaParam_2[NCHARGE][NSECT][NZED][NCENT];

int GetPtBin(float);

void LoadEmcMatchingParameters();

void LoadInitialdPhiMeanParameters();
void LoadInitialdPhiSigmaParameters();

void LoadFinaldPhiMeanParameters();
void LoadFinaldPhiSigmaParameters();

void LoadMean();
void LoadSigma();

void LoadInitialdZMeanParameters();
void LoadInitialdZSigmaParameters();

void LoadFinaldZMeanParameters();
void LoadFinaldZSigmaParameters();

float CalculateInitialEmcsdPhi(int, int, int, int, float, float);
float CalculateFinalEmcsdPhi(int, int, int, int, float, float);

float CalculateCorrectedEmcdZ(float, float, float, int);

float CalculateInitialEmcsdZ(int, int, int, int, float, float);
float CalculateFinalEmcsdZ(int, int, int, int, float, float);

#endif
