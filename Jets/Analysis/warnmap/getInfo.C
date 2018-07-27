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
#include <iomanip>

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

void getInfo(int what = 0)
{
  gStyle->SetOptStat(0);

  const char *armsect[8] = {"Sector: 0(PbSc)", "Sector: 1(PbSc)", "Sector: 2(PbSc)", "Sector: 3(PbSc)",
			    "Sector: 4(PbGl)", "Sector: 5(PbGl)", "Sector: 6(PbSc)", "Sector: 7(PbSc)"
  };

  //Read the warnmap
  cout << "Reading the warnmap" << endl;

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
      cout << "ERROR!!! WarnMap file could not be located!" << endl;
    }

  int dead[NARMSECT];
  int hot[NARMSECT];
  int edge[NARMSECT];
  int uncalibrated[NARMSECT];
  int dead_uncalibrated[NARMSECT];
  int hot_uncalibrated[NARMSECT];
  int all[NARMSECT];
  int around[NARMSECT];

  for (int ias = 0; ias < NARMSECT; ias++)
    {
      dead[ias] = 0;
      hot[ias] = 0;
      edge[ias] = 0;
      uncalibrated[ias] = 0;
      dead_uncalibrated[ias] = 0;
      hot_uncalibrated[ias] = 0;
      all[ias] = 0;
      around[ias] = 0;
    }

  if(what == 0)
    {
      ofstream myfileInfo("/direct/phenix+u/arbint/Jets/Analysis/warnmap/CuAuMapInfo.txt");
    }
  else
    {
      ofstream myfileInfo("/direct/phenix+u/arbint/Jets/Analysis/warnmap/PPMapInfo.txt");
    }

  if (myfileInfo.is_open())
    {
      string line;
      while (getline(fmap, line))
        {
	  istringstream one_line(line);
	  int as, y, z, map;
	  if (one_line >> as >> y >> z >> map)
            {
	      if (map & MASK_DEAD)
                {
		  dead[as]++;
                }
	      if (map & MASK_HOT)
                {
		  hot[as]++;
                }
	      if (map & MASK_UNCALIB)
                {
		  uncalibrated[as]++;
                }
	      if (map & MASK_EDGE)
                {
		  edge[as]++;
                }
	      if ((map & MASK_UNCALIB) && (map & MASK_DEAD))
                {
		  dead_uncalibrated[as]++;
                }
	      if ((map & MASK_UNCALIB) && (map & MASK_HOT))
                {
		  hot_uncalibrated[as]++;
                }
	      if ((map & MASK_DEAD) || (map & MASK_HOT) || (map & MASK_UNCALIB))
                {
		  all[as]++;
                }
	      if ((map & MASK_AROUND_DEAD) || (map & MASK_AROUND_HOT) || (map & MASK_AROUND_UNCALIB))
                {
		  around[as]++;
                }
            }
        }

      fmap.close();

      myfileInfo << "--------------------------------------------------------------------------------------------------------------------------------" << "\n";
      myfileInfo << setw(15) << "Sector" << setw(16) << "Dead" << setw(17) << "Hot" << setw(18) << "Uncalib"
		 << setw(19) << "Dead and Uncalib" << setw(20) << "Hot and Uncalib" << setw(21) << "Dead||Hot||Uncalib" << "\n";
      myfileInfo << "--------------------------------------------------------------------------------------------------------------------------------" << "\n";

      int deadd = 0;
      int hott = 0;
      int uncalibratedd = 0;
      int edgee = 0;
      int dead_uncalibratedd = 0;
      int hot_uncalibratedd = 0;
      int alll = 0;
      int aroundd = 0;
      for (int ias = 0; ias < NARMSECT; ias++)
        {
	  myfileInfo << setw(15) << armsect[ias] << setw(16) << dead[ias] << setw(17) << hot[ias] << setw(18) << uncalibrated[ias]
		     << setw(19) << dead_uncalibrated[ias] << setw(20) << hot_uncalibrated[ias] << setw(21) << all[ias] << "\n";
	  deadd = deadd + dead[ias];
	  hott = hott + hot[ias];
	  uncalibratedd = uncalibratedd + uncalibrated[ias];
	  edgee = edgee + edge[ias];
	  dead_uncalibratedd = dead_uncalibratedd + dead_uncalibrated[ias];
	  hot_uncalibratedd = hot_uncalibratedd + hot_uncalibrated[ias];
	  alll = alll + all[ias];
	  aroundd = aroundd + around[ias];
        }
      myfileInfo << "--------------------------------------------------------------------------------------------------------------------------------" << "\n";
      myfileInfo << setw(15) << "Total: " << setw(16) << deadd << setw(17) << hott << setw(18) << uncalibratedd
		 << setw(19) << dead_uncalibratedd << setw(20) << hot_uncalibratedd << setw(21) << alll << "\n";
      myfileInfo << "--------------------------------------------------------------------------------------------------------------------------------" << "\n";




      int nTotal[NARMSECT + 1] = {2592, 2592, 2592, 2592, 4608, 4608, 2592, 2592, 24768};

      myfileInfo << "--------------------------------------------------------------------" << "\n";
      myfileInfo << "--------------------------------------------------------------------" << "\n";
      myfileInfo << setw(15) << "Sector" << setw(16) << "Total towers" << setw(17) << "Inactive" << setw(18) << "% Inactive" << "\n";
      myfileInfo << "--------------------------------------------------------------------" << "\n";

      for (int ias = 0; ias < NARMSECT; ias++)
        {
	  myfileInfo << setw(15) << armsect[ias] << setw(16) << nTotal[ias] << setw(17) << all[ias] << setw(18) << setprecision(3) << (float)all[ias] * 100 / nTotal[ias] << "\n";
        }
      myfileInfo << "--------------------------------------------------------------------" << "\n";
      myfileInfo << setw(15) << "Total" << setw(16) << nTotal[8] << setw(17) << alll << setw(18) << (float)alll * 100 / nTotal[8] << "\n";
      myfileInfo << "--------------------------------------------------------------------" << "\n";


      myfileInfo << "-----------------------------------------------------------------------------" << "\n";
      myfileInfo << "---------------------------------------------------------------------------------------" << "\n";
      myfileInfo << setw(15) << "Sector" << setw(16) << "Total towers" << setw(17) << "Edge" << setw(18) << "Around" << setw(19) << "% Semi-inactive" << "\n";
      myfileInfo << "---------------------------------------------------------------------------------------" << "\n";

        for (int ias = 0; ias < NARMSECT; ias++)
        {
            myfileInfo << setw(15) << armsect[ias] << setw(16) << nTotal[ias] << setw(17) << edge[ias] << setw(18) << around[ias] << setw(19) << (float)(edge[ias] + around[ias]) * 100 / nTotal[ias] << "\n";
        }
        myfileInfo << "---------------------------------------------------------------------------------------" << "\n";
        myfileInfo << setw(15) << "Total" << setw(16) << nTotal[8] << setw(17) << edgee << setw(18) << aroundd << setw(19) << (float)(edgee + aroundd) * 100 / nTotal[8] << "\n";
        myfileInfo << "---------------------------------------------------------------------------------------" << "\n";
    }

    cout << "Writing the warn map information" << endl;
    myfileInfo.close();
}


