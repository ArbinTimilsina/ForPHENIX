#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <PHCentralTrack.h>
#include <PHGlobal.h>
#include <getClass.h>

#include <TMath.h>
#include <TH2.h>

#include "AlphaVsPhi.h"

using namespace std;

AlphaVsPhi::AlphaVsPhi(std::string outfilename)
    : SubsysReco("AlphaVsPhi"),
      outfname(outfilename)
{
    outfile = new TFile(outfname.c_str(), "RECREATE");

    nEvents = 0;

    hAlphaVsPhi_East = new TH2F("hAlphaVsPhi_East", "Alpha vs. phi, east", 360, 2.2, 3.7, 100, -0.01, 0.01);
    hAlphaVsPhi_West = new TH2F("hAlphaVsPhi_West", "Alpha vs. phi, west", 360, -0.55, 0.95, 100, -0.01, 0.01);

}

void AlphaVsPhi::SetCuAu(bool what = true)
{
    cout << "***********************************************************************" << endl;
    if(what)
        {
            cout << endl;
            cout << "Running over Cu+Au data setup" << endl << endl;
        }
    else
        {
            cout << endl;
            cout << "Running over p+p data setup" << endl << endl;
        }

    isCuAu = what;
    cout << "***********************************************************************" << endl << endl;
}

int AlphaVsPhi::process_event(PHCompositeNode *topNode)
{
    nEvents++;
    if (nEvents % 1000 == 0)
        {
                    cout << "Events analyzed = " << nEvents << endl;
        }

    PHCentralTrack *track = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");

    for (unsigned int i = 0; i < track->get_npart(); i++)
        {

            int quality = track->get_quality(i);
            float phi = track->get_phi(i);
            float alpha = track->get_alpha(i);
            int arm = track->get_dcarm(i);
            int n0 = track->get_n0(i);
            float zed = track->get_zed(i);

            bool qualityTrack = quality == 31 || quality == 63;
            bool passN0Cut = n0 < 1;
            bool passZedCut = abs(zed) < 75;

            bool goodTrack = qualityTrack && passN0Cut && passZedCut;
 
	    if(isCuAu)
                {
		    alpha = getNewAlpha(alpha, phi, arm);
		}

           if (goodTrack)
                {
                    if(arm == 0)
                        {
                            hAlphaVsPhi_East->Fill(phi, alpha);
                        }
                    else
                        {
                            hAlphaVsPhi_West->Fill(phi, alpha);
                        }
                }
        }

    return 0;
}

int AlphaVsPhi::End(PHCompositeNode *topNode)
{
    outfile->Write();
    outfile->Close();

    return 0;
}
