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
#include "/direct/phenix+u/arbint/histoMakeup/histoMakeup.C"

using namespace std;

#define NARMSECT  8
#define YPOS_PBGL 48
#define YPOS_PBSC 36
#define ZPOS_PBGL 96
#define ZPOS_PBSC 72

#define MASK_EDGE            0x01
#define MASK_AROUND_DEAD     0x02
#define MASK_AROUND_HOT      0x04
#define MASK_DEAD            0x08
#define MASK_HOT             0x10
#define MASK_AROUND_UNCALIB  0x20
#define MASK_UNCALIB         0x40

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

void drawMap(int what = 0)
{
  gStyle->SetOptStat(0);

  const char *armsect[8] = {"Sector: 0 (PbSc)", "Sector: 1 (PbSc)", "Sector: 2 (PbSc)", "Sector: 3 (PbSc)",
			    "Sector: 4 (PbGl)", "Sector: 5 (PbGl)", "Sector: 6 (PbSc)", "Sector: 7 (PbSc)"
  };


  int warnmap[NARMSECT][YPOS_PBGL][ZPOS_PBGL];
  //Read the warnmap
  cout << "Reading the mean and sigma " << endl;

  ifstream fmap;
  if(what == 0)
    {
      fmap.open("warnmapCuAu.txt");
    }
  else if(what == 1)
    {
      fmap.open("warnmapPP.txt");
    }
  else
    {
      cout << "Chose 0 for CuAu, 1 for PP" << endl;
      exit(1);
    }

  if (!fmap)
    {
      cout << "ERROR!!! Warnmap for p+p not found!" << endl;
    }

  string line;
  memset(warnmap, 0, sizeof(warnmap));
  while (getline(fmap, line))
    {
      istringstream one_line(line);
      int as, y, z, status;
      if (one_line >> as >> y >> z >> status)
        {
	  if (IsValidYZ(as, y, z))
            {
	      warnmap[as][y][z] = status;
            }
        }
    }

  fmap.close();

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
      canvas[iarm] = new TCanvas(Form("Canvas_arm%d", iarm), "", 360, 720);
      canvas[iarm]->Divide(1, 4);
      for (int isector = 0; isector < 4; isector++)
        {
	  int ias = iarm * 4 + isector;
	  canvas[iarm]->cd(4 - isector);
	  hMap[ias]->SetStats(kFALSE);

	  hMap[ias]->GetXaxis()->SetTitle("Z Position");
	  hMap[ias]->GetXaxis()->SetTitleSize(0.04);
	  hMap[ias]->GetXaxis()->SetTitleOffset(1.2);
	  hMap[ias]->GetXaxis()->SetTitleFont(42);

	  hMap[ias]->GetYaxis()->SetTitle("Y Position");
	  hMap[ias]->GetYaxis()->SetTitleSize(0.04);
	  hMap[ias]->GetYaxis()->SetTitleOffset(0.5);
	  hMap[ias]->GetYaxis()->SetTitleFont(42);

	  hMap[ias]->SetTitle("");

	  TPaveText *title = getHistoTitle();
	  if(ias == 3 || ias == 7)
            {
	      if(what == 0)
                {
		  title->AddText(Form("%s                                        For Run 12 Cu+Au @ 200 GeV", armsect[ias]));
                }
                else
                {
                    title->AddText(Form("%s                                                For Run 12 p+p @ 200 GeV", armsect[ias]));
                }
            }
            else
            {
                title->AddText(Form("%s", armsect[ias]));
            }
            hMap[ias]->Draw("col");
            title->Draw("Same");
        }

        if(what == 0)
        {
            canvas[iarm]->SaveAs(Form("/direct/phenix+hhj/arbint/plots/JetAnalyzer/WarnMaps/CuAu_Map_Arm%d.png", iarm));
        }
        else
        {
            canvas[iarm]->SaveAs(Form("/direct/phenix+hhj/arbint/plots/JetAnalyzer/WarnMaps/PP_Map_Arm%d.png", iarm));
        }
    }
}

