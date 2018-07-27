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
#include <ctime>
#include <string>
#include <TROOT.h>
#include <cmath>

using namespace std;

void convertArmSect()
{

  ifstream originalMap;
  originalMap.open("/direct/phenix+u/arbint/Jets/Analysis/warnmap/warnmapPPOriginal.txt");

  if (!originalMap)
    {
      cout << "ERROR!!! File " << originalMap << " could not be located!" << endl;
    }

  ofstream newMap ("/direct/phenix+u/arbint/Jets/Analysis/warnmap/warnmapPP.txt");
  if (newMap.is_open())
    {
      string originalLine;
      while (getline(originalMap, originalLine))
        {
	  istringstream original_one_line(originalLine);
	  int as, y, z, status;
	  if (original_one_line >> as >> y >> z >> status)
            {
	      int ass = -1;
	      if(as == 4)
                {
		  ass = 7;
                }
	      else if(as == 5)
                {
		  ass = 6;
                }
	      else if(as == 6)
                {
		  ass = 5;
                }
	      else if(as == 7)
                {
		  ass = 4;
                }
	      else
                {
		  ass = as;
                }

	      newMap   << ass << " "
                       << y << " "
                       << z << " "
                       << status << "\n";
            }
        }
        originalMap.close();
    }
    newMap.close();
}
