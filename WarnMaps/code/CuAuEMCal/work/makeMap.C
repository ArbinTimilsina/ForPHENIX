#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH3D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>

using namespace std;

#define NARMSECT  8
#define YPOS_PBGL 48
#define YPOS_PBSC 36
#define ZPOS_PBGL 96
#define ZPOS_PBSC 72
#define NECOREBIN 8

#define MASK_EDGE            0x01
#define MASK_AROUND_DEAD     0x02
#define MASK_AROUND_HOT      0x04
#define MASK_DEAD            0x08
#define MASK_HOT             0x10
#define MASK_AROUND_UNCALIB  0x20
#define MASK_UNCALIB         0x40

#define SIGMACUT_HOT 3.5
#define SIGMACUT_DEAD 3.5

int IsPbGl(int armsect)
{
  return ((armsect == 4 || armsect == 5) ? 1 : 0);
}


int IsValidYZ(int as, int y, int z)
{
  int ret = 0;
  if (IsPbGl(as))
    {
      if (y >= 0 && y < YPOS_PBGL && z >= 0 && z < ZPOS_PBGL)
        {
	  ret = 1;
        }
    }
  else
    {
      if (y >= 0 && y < YPOS_PBSC && z >= 0 && z < ZPOS_PBSC)
        {
	  ret = 1;
        }
    }
  return ret;
}

void makeMap()
{
  gStyle->SetOptStat(0);

  const char *armsect[8] = {"Sector: 0(PbSc)", "Sector: 1(PbSc)", "Sector: 2(PbSc)", "Sector: 3(PbSc)",
			    "Sector: 4(PbGl)", "Sector: 5(PbGl)", "Sector: 6(PbSc)", "Sector: 7(PbSc)"
  };

  TFile *file = new TFile("/direct/phenix+hhj/arbint/RootFiles/CuAuEMCal.root", "r");

  int warnmap[NARMSECT][YPOS_PBGL][ZPOS_PBGL];

  //Get the 3D histogram
  TH3D *hEcore[NARMSECT];
  for (int ias = 0; ias < NARMSECT; ias++)
    {
      char hname[256];
      sprintf(hname, "hEmc3D_Armsect_%i", ias);
      hEcore[ias] = (TH3D*)file->Get(hname);
    }

  //Read in mean and sigma
  float mean[NARMSECT][NECOREBIN];
  float sigma[NARMSECT][NECOREBIN];

  cout << "Reading the mean and sigma " << endl;

  ifstream fmean;
  fmean.open("meanSigma.txt");

  if (!fmean)
    {
      cout << "ERROR!!! Mean and sigma list file could not be located!" << endl;
    }


  for (int ias = 0; ias < NARMSECT; ias++)
    {
      for (int iec = 0; iec < NECOREBIN; iec++)
        {
	  fmean >> mean[ias][iec]
		>> sigma[ias][iec];
        }
    }

  fmean.close();

  //Set around edge and dead, and hot
  memset(warnmap, 0, sizeof(warnmap));
  for (int ias = 0; ias < NARMSECT; ias++)
    {
      for (int iz = 0; iz < ZPOS_PBGL; iz++ )
        {
	  for (int iy = 0; iy < YPOS_PBGL; iy++ )
            {
	      if (IsValidYZ(ias, iy, iz))
                {
		  //Check edege
		  if (IsPbGl(ias))
                    {
		      if (iy == 0 || iy == YPOS_PBGL - 1 ||
			  iz == 0 || iz == ZPOS_PBGL - 1)
                        {
			  warnmap[ias][iy][iz] |= MASK_EDGE;
                        }
                    }
		  else
                    {
		      if (iy == 0 || iy == YPOS_PBSC - 1 ||
			  iz == 0 || iz == ZPOS_PBSC - 1)
                        {
			  warnmap[ias][iy][iz] |= MASK_EDGE;
                        }
                    }

		  //Check if hot- done at all ecore bins- starting at 500 MeV
		  int ebin = 0;
		  while(ebin < NECOREBIN-2)
                    {
		  //Check if dead- done at 500 MeV
		  int nhitDead = 0;
		  for (int iec = 0; iec < NECOREBIN; iec++)
                    {
		      nhitDead += hEcore[ias]->GetBinContent(hEcore[ias]->GetXaxis()->FindBin(iz),
							     hEcore[ias]->GetYaxis()->FindBin(iy),
							     iec + 1);
                    }
		  float limitDead = mean[ias][0] - SIGMACUT_DEAD * sigma[ias][0];

		  if (nhitDead < limitDead)
                    {
		      warnmap[ias][iy][iz] |= MASK_DEAD;
                    }

		      float nhitHot = 0.0;
		      for (int iec = ebin; iec < NECOREBIN; iec++)
                        {
			  nhitHot += hEcore[ias]->GetBinContent(hEcore[ias]->GetXaxis()->FindBin(iz),
								hEcore[ias]->GetYaxis()->FindBin(iy),
								iec + 1);
                        }

		      float limitHot = mean[ias][ebin] + SIGMACUT_HOT * sigma[ias][ebin];

		      if (nhitHot > limitHot)
                        {
			  warnmap[ias][iy][iz] |= MASK_HOT;
                        }
		      ebin++;
                    }

                }
            }
        }
    }

  //Exceptional Hot tower
  warnmap[5][17][45] |= MASK_HOT;

  //Set uncalibrated
  cout << "Reading the uncalibrated list" << endl;

  ifstream fcalib;
  fcalib.open("uncalibrated.txt");

  if (!fcalib)
    {
      cout << "ERROR!!! Uncalibrated list file could not be located!" << endl;
    }

  string line;
  while (getline(fcalib, line))
    {
      istringstream one_line(line);
      int as, y, z;
      if (one_line >> as >> y >> z)
        {
	  warnmap[as][y][z] |= MASK_UNCALIB;
        }
    }

  fcalib.close();

  //Set around hot/dead/uncalib
  for (int ias = 0; ias < NARMSECT; ias++)
    {
      for (int iz = 0; iz < ZPOS_PBGL; iz++ )
        {
	  for (int iy = 0; iy < YPOS_PBGL; iy++ )
            {
	      if (IsValidYZ(ias, iy, iz))
                {
		  for (int dy = -1; dy <= 1; dy++)
                    {
		      for (int dz = -1; dz <= 1; dz++)
                        {
			  if ((! (dy == 0 && dz == 0)) && // except itself
			      IsValidYZ(ias, iy + dy, iz + dz))
                            {
			      if (warnmap[ias][iy + dy][iz + dz] & MASK_HOT)
                                {
				  warnmap[ias][iy][iz] |= MASK_AROUND_HOT;
                                }
			      if (warnmap[ias][iy + dy][iz + dz] & MASK_DEAD)
                                {
				  warnmap[ias][iy][iz] |= MASK_AROUND_DEAD;
                                }
			      if (warnmap[ias][iy + dy][iz + dz] & MASK_UNCALIB)
                                {
				  warnmap[ias][iy][iz] |= MASK_AROUND_UNCALIB;
                                }
                            }
                        }
                    }
                }
            }
        }
    }


  //Write the reasult
  cout << "Writing out the warnmap" << endl;

  ofstream myfileWarn ("/direct/phenix+u/arbint/Jets/Analysis/warnmap/warnmapCuAu.txt");
  if (myfileWarn.is_open())
    {
      for (int ias = 0; ias < NARMSECT; ias++)
        {
	  for (int iz = 0; iz < ZPOS_PBGL; iz++ )
            {
	      for (int iy = 0; iy < YPOS_PBGL; iy++ )
                {
		  if (IsValidYZ(ias, iy, iz))
                    {
		      myfileWarn << ias << " "
				 << iy << " "
				 << iz << " "
				 << warnmap[ias][iy][iz] << "\n";
                    }
                }
            }
        }
    }

  myfileWarn.close();

  // current date/time based on current system
  time_t now = time(0);

  // convert now to string form
  char* dt = ctime(&now);

  //Write out info for record
  ofstream myfileInfo ("/direct/phenix+u/arbint/Jets/Analysis/warnmap/infoCuAuMap");
  if (myfileInfo.is_open())
    {
      myfileInfo << "\n" << "Warn map made on: " << dt << "\n"
		 << "Sigma cut used for hot: " << SIGMACUT_HOT << "\n"
		 << "Sigma cut used for dead: " << SIGMACUT_DEAD;
    }
  myfileInfo.close();

  //Plot the warnmap
  TH2F* hMap[NARMSECT];
  for (int ias = 0; ias < NARMSECT; ias++)
    {
      char hname[256];
      sprintf(hname, "hMap_sector_%i", ias);

      char htitle[256];
      sprintf(htitle, "%s", armsect[ias]);

      int ny, nz;
      if (IsPbGl(ias))
        {
	  ny = YPOS_PBGL;
	  nz = ZPOS_PBGL;
        }
      else
        {
	  ny = YPOS_PBSC;
	  nz = ZPOS_PBSC;
        }

      hMap[ias] = new TH2F(hname, htitle, nz, 0, nz, ny, 0, ny);

      const int ncont = 5; // good, hot, dead, uncalib, around/edge
      hMap[ias]->GetZaxis()->SetRangeUser(1, ncont);
      gStyle->SetNumberContours(ncont);

      for (int iz = 0; iz < ZPOS_PBGL; iz++ )
        {
	  for (int iy = 0; iy < YPOS_PBGL; iy++ )
            {
	      if (IsValidYZ(ias, iy, iz))
                {
		  int status;
		  if      (warnmap[ias][iy][iz] & MASK_HOT)
                    {
		      status = 4;
                    }
		  else if (warnmap[ias][iy][iz] & MASK_DEAD)
                    {
		      status = 3;
                    }
		  else if (warnmap[ias][iy][iz] & MASK_UNCALIB)
                    {
		      status = 2;
                    }
		  else if (warnmap[ias][iy][iz] & MASK_AROUND_HOT     ||
			   warnmap[ias][iy][iz] & MASK_AROUND_DEAD    ||
			   warnmap[ias][iy][iz] & MASK_AROUND_UNCALIB ||
                             warnmap[ias][iy][iz] & MASK_EDGE )
                    {
                        status = 1;
                    }
                    else
                    {
                        status = 0;
                    }
                    hMap[ias]->Fill(iz, iy, status);
                }
            }
        }
    }

    int palette[4];
    palette[0] = 16;
    palette[1] = 7;
    palette[2] = 1;
    palette[3] = 2;

    gStyle->SetPalette(4, palette);
    gStyle->SetOptStat(0);

    TCanvas* canvas[2];
    for (int iarm = 0; iarm < 2; iarm++)
    {
        canvas[iarm] = new TCanvas(Form("Canvas_arm%d", iarm), "", 500, 700);
        canvas[iarm]->Divide(1, 4);
        for (int isector = 0; isector < 4; isector++)
        {
            int ias = iarm * 4 + isector;
            canvas[iarm]->cd(4 - isector);
            hMap[ias]->SetStats(kFALSE);
            hMap[ias]->Draw("col");
        }
        canvas[iarm]->SaveAs(Form("/direct/phenix+hhj/arbint/plots/JetAnalyzer/WarnMaps/CuAu_Map_Arm%d.png", iarm));
    }
}
