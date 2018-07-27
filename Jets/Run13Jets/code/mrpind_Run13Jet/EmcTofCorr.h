#ifndef __EMCTOFCORR_H__
#define __EMCTOFCORR_H__

const static int NARMSECTS_TOF = 8;
const static int NIYVAL = 48;
const static int NIZVAL = 96;

void ReadTOFMap(int fillnumber);
float TofCorrection(int armsect, int iy, int iz);

static float tofmap[NARMSECTS_TOF][NIYVAL][NIZVAL];

#endif //__EMCTOFCORR_H__
