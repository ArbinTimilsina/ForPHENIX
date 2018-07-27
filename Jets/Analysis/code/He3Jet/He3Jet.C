//General PHENIX tools
#include <getClass.h>
#include <PHCompositeNode.h>
#include <phool.h>
#include <RunHeader.h>

//Fun4All tools
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <Fun4AllReturnCodes.h>

//FastJet tools
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

//My source file
#include "He3Jet.h"

//For ERT trigger
static const int ERT_BIT_MASK_A = 0x00001000;
static const int ERT_BIT_MASK_B = 0x00002000;
static const int ERT_BIT_MASK_C = 0x00004000;

using namespace std;
using namespace findNode;
using namespace fastjet;

//================================ Constructor ================================
//Here we can initiate some variables
He3Jet::He3Jet(std::string outfilename)
    : SubsysReco("He3Jet"),
      verbo(1),
      outfname(outfilename),
      writeTree(true)
{
    outfile = new TFile(outfname.c_str(), "RECREATE");

    nTotalEvents = 0;
    nGoodEvents = 0;
    nJets = 0;

    InitTrees(writeTree);
    InitHistograms();

    return;
}

int He3Jet::Init(PHCompositeNode *topNode)
{
    if (verbo)
        {
            cout << ">>>>>>>>>>>>>>>>>>>  Init called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
        }

    //Load and read the warnmap
    EmcMap::ReadWarnMap("warnmap.txt");

    return 0;
}

int He3Jet::InitRun(PHCompositeNode *topNode)
{
    if (verbo)
        {
            cout << ">>>>>>>>>>>>>>>>>>>> InitRun called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
        }

    RunHeader *run_header = getClass<RunHeader> (topNode, "RunHeader");
    if (!run_header)
        {
            cout << "No RunHeader! No sense continuing" << endl;
            exit(1);
        }

    runNumber = run_header->get_RunNumber();

    triggerHelper = new TriggerHelper(topNode);

    return EVENT_OK;
}

int He3Jet::process_event(PHCompositeNode *topNode)
{
    //He3+Au runs
    if(runNumber < 415370 || runNumber > 416893)
        {
            return DISCARDEVENT;
        }


    //Make trigger selection
    TrigLvl1 *lvl1trigger = getClass<TrigLvl1>(topNode, "TrigLvl1");

    int ertFired = 0;
    if((lvl1trigger->get_lvl1_trigscaled()&ERT_BIT_MASK_A) ||
       (lvl1trigger->get_lvl1_trigscaled()&ERT_BIT_MASK_B) ||
       (lvl1trigger->get_lvl1_trigscaled()&ERT_BIT_MASK_C))
        {
            ertFired = 1;
        }

    if(!ertFired)
        {
            return DISCARDEVENT;
        }

    //Get the data I need...
    phglobal                               = getClass<PHGlobal>                (topNode, "PHGlobal");
    if (nTotalEvents == 1 && !phglobal)
        {
            cout << "No PHGlobal!  No sense continuing" << endl;
            exit(1);
        }

    evtheader                              = getClass<EventHeader>             (topNode, "EventHeader");
    if (nTotalEvents == 1 && !evtheader)
        {
            cout << "No EventHeader! No sense continuing" << endl;
            exit(1);
        }

    //Vertex
    zVertex = phglobal->getBbcZVertex();
    hVertex->Fill(zVertex);

    //These are the events for the trigger selection and good runs
    nTotalEvents++;

    // Informational message...
    if (nTotalEvents % 1000 == 0 && verbosity)
        {
            if (verbo)
                {
                    cout << "Events analyzed for run " << runNumber << " = " << nTotalEvents << endl;
                }
        }

    //Vertex cut
    if (fabs(zVertex) > VERTEX_CUT)
        {
            return DISCARDEVENT;
        }
    nGoodEvents++;

    //Get tracks from this event
    GetTracks(topNode, tracks_list);

    //Get clusters from this event
    GetClusters(topNode, clusters_list, zVertex);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Make anti-kt jet here
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    float R = 0.3;
    fastjet::JetDefinition *antikt = new fastjet::JetDefinition(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best);

    std::vector<fastjet::PseudoJet> jetParticles_all;
    jetParticles_all.clear();

    unsigned int indexTotal = tracks_list.size() + clusters_list.size();
    unsigned int indexCharged = tracks_list.size();

    float particlePt[indexTotal];
    fill(particlePt, particlePt + indexTotal / sizeof(float), -999.9);
    bool particleErtTrigger[indexTotal];
    fill(particleErtTrigger, particleErtTrigger + indexTotal / sizeof(bool), false);

    int index = 0;
    for (unsigned int h = 0; h < tracks_list.size(); h++)
        {
            fastjet::PseudoJet pseudoCharged(tracks_list[h].px,
                                             tracks_list[h].py,
                                             tracks_list[h].pz,
                                             tracks_list[h].mom);
            pseudoCharged.set_user_index(index);
            particlePt[index] = tracks_list[h].pT;
            particleErtTrigger[index] = false;
            jetParticles_all.push_back(pseudoCharged);
            index++;
        }

    for (unsigned int n = 0; n < clusters_list.size(); n++)
        {
            fastjet::PseudoJet pseudoNeutral(clusters_list[n].px,
                                             clusters_list[n].py,
                                             clusters_list[n].pz,
                                             clusters_list[n].energy);
            pseudoNeutral.set_user_index(index);
            particlePt[index] = clusters_list[n].pT;
            particleErtTrigger[index] = clusters_list[n].ertFired;
            jetParticles_all.push_back(pseudoNeutral);
            index++;
        }

    fastjet::ClusterSequence jetAll(jetParticles_all, *antikt);
    std::vector<fastjet::PseudoJet> fastAll = jetAll.inclusive_jets();
    for (unsigned int n = 0; n < fastAll.size(); n++)
        {
            fastjet::PseudoJet aFastJet = fastAll[n];

            float chargedPt = 0.0;

            float jetPt  = aFastJet.perp();
            float jetEta = aFastJet.pseudorapidity();
            float jetPhi = phiReduce(aFastJet.phi());

            vector<fastjet::PseudoJet> constituents = jetAll.constituents(aFastJet);
            unsigned int nconst = constituents.size();

            bool ertTriggerFired = false;
            for (unsigned int iconst = 0; iconst < nconst; iconst++)
                {
                    unsigned int indx = constituents[iconst].user_index();

                    if (indx < indexCharged)//Charged particles
                        {
                            chargedPt += particlePt[indx];
                        }
                    else  //Neutral particles
                        {
                            if(particleErtTrigger[indx] == true)
                                {
                                    ertTriggerFired = true;
                                }
                        }
                }

            int jetArm = 1;
            if (jetPhi > 1.57)
                {
                    jetArm = 0;
                }

            if(!ertTriggerFired)
                {
                    continue;
                }

            float jetCf = chargedPt / jetPt;
            if((jetPt > 6.0) && ((float)nconst >= 3.0))
                {
                    hCf->Fill(jetCf);
                }

            //Jet selection
            bool passJetLevelCuts = (jetPt > 6.0) && ((float)nconst >= 3.0);
            if(!passJetLevelCuts)
                {
                    continue;
                }
            hJets->Fill(jetPt);

            jets temp;
            temp.arm           = jetArm;
            temp.pT            = jetPt;
            temp.eta           = jetEta;
            temp.phi           = jetPhi;
            temp.nc            = (float)nconst;
            temp.cf            = jetCf;
            jets_list.push_back(temp);

            nJets++;
        }
    delete antikt;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Fill Trees
    FillTrees(writeTree);

    // any other return code might lead to aborting the event or analysis
    return 0;
}

int He3Jet::ResetEvent(PHCompositeNode *topNode)
{
    tracks_list.clear();
    clusters_list.clear();
    jets_list.clear();

    return EVENT_OK;
}

int He3Jet::EndRun(const int runNumber)
{
    if (verbo)
        {
            cout << ">>>>>>>>>>>>>>>>>>>  EndRun called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
        }

    if (verbo)
        {
            std::ios_base::fmtflags originalFlags = std::cout.flags();

            cout.precision(3);
            cout << "+++++++++++++  Statistics:     +++++++++++++++++++" << endl;
            cout << "Run number:                                            " << runNumber << endl;
            cout << "Number of total events for this run:                   " << nTotalEvents << endl;
            cout << "Number of good events for this run:                   " << nGoodEvents << endl;

            cout << "++++++++++  Statistics: Jets   +++++++++++++++++" << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "Number of total Jets (Anti-kt, R = 0.3) for this run:  " << setw(5) << nJets << endl;
            std::cout.flags(originalFlags);
        }

    return EVENT_OK;
}

int He3Jet::End(PHCompositeNode *topNode)
{
    if (verbo)
        {
            cout << ">>>>>>>>>>>>>>>>>>>  End called <<<<<<<<<<<<<<<<<<<<<<<" << endl;
        }

    hEvents->SetBinContent(1, nTotalEvents);
    hEvents->SetBinContent(2, nGoodEvents);

    outfile->Write();
    outfile->Close();

    delete triggerHelper;

    return EVENT_OK;
}

void He3Jet::GetTracks(PHCompositeNode *topNode, std::vector<tracks>& track_list)
{
    PHCentralTrack  *phcentral = getClass<PHCentralTrack>(topNode, "PHCentralTrack");
    if (!phcentral)
        {
            cout << "No PHCentral!  No sense continuing" << endl;
            exit(1);
        }

    for (unsigned int i = 0; i < phcentral->get_npart(); i++)
        {
            float mom           = phcentral->get_mom(i); // Magnitude of the momentum.
            float theta         = phcentral->get_the0(i); //The track's theta direction at the vertex
            float pT            = mom * sin(theta);
            float eta           = (float)(-log(tan(0.5 * theta)));
            float phi           = phcentral->get_phi0(i); // The track's phi direction at the vertex
            float phiDC         = phcentral->get_phi(i); // The phi coordinate at which the track crosses the drift chamber reference radius
            float zedDC         = phcentral->get_zed(i); //The Z coordinate at which the track crosses PC1
            float alpha         = phcentral->get_alpha(i); // This is the inclination of the track w.r.t. an infinite momentum track
            float energy        = phcentral->get_ecore(i); // EMC "shower core" energy

            int charge          = phcentral->get_charge(i);
            int quality         = phcentral->get_quality(i); // Quality of the Drift Chamber Tracks
            int n0              = phcentral->get_n0(i); //The number of phototubes that fired in the normally sized ring area
            int arm             = phcentral->get_dcarm(i); //Arm containing the track (East=0, West=1)
            int sector          = phcentral->get_sect(i); //EMC sector of the associated cluster
            int emcid           = phcentral->get_emcid(i); // Index of the emc cluster used

            float pc3dphi       = phcentral->get_pc3dphi(i);
            float pc3sdphi      = phcentral->get_pc3sdphi(i);
            float pc3dz         = phcentral->get_pc3dz(i);
            float pc3sdz        = phcentral->get_pc3sdz(i);
            float emcdphi       = phcentral->get_emcdphi(i);
            float emcsdphi      = phcentral->get_emcsdphi(i);
            float emcdz         = phcentral->get_emcdz(i);
            float emcsdz        = phcentral->get_emcsdz(i);

            float beta          = phcentral->get_beta(i);

            float px = pT * cos(phi);
            float py = pT * sin(phi);
            float pz = mom * cos(theta);
            float eT = energy * sin(theta);

            int armsect = -99;
            if (arm == 1)
                {
                    armsect = sector;
                }
            else
                {
                    armsect = 4 + sector;
                }

            float board = -9999.99;

            //Track selection
            bool passPt = (pT > TRACK_MIN_PT_CUT);
            bool passZedDC = fabs(zedDC) < 75.0;
            bool passQuality = (quality == 63) || (quality == 31);

            bool minRequirement = passPt && passZedDC && passQuality;
            if(!minRequirement)
                {
                    continue;
                }

            hPC3dphi->Fill(pc3dphi);
            hPC3dz->Fill(pc3dz);

            hEMCdphi->Fill(emcdphi);
            hEMCdz->Fill(emcdz);

            bool pc3Matching = (fabs(pc3dphi) < 0.00804) && (fabs(pc3dz) < 6.015);
            bool emcMatching = (fabs(emcdphi) < 0.0216) && (fabs(emcdz) < 10.905);
            bool passMatching = pc3Matching || emcMatching;

            if(!passMatching)
                {
                    continue;
                }

            tracks temp;
            temp.mom            = mom;
            temp.theta          = theta;
            temp.pT             = pT;
            temp.eta            = eta;
            temp.phi            = phi;
            temp.phiDC          = phiDC;
            temp.zedDC          = zedDC;
            temp.alpha          = alpha;
            temp.energy         = energy;
            temp.board          = board;

            temp.charge         = charge;
            temp.quality        = quality;
            temp.n0             = n0;
            temp.arm            = arm;
            temp.armsect        = armsect;
            temp.emcid          = emcid;

            temp.px             = px;
            temp.py             = py;
            temp.pz             = pz;
            temp.eT             = eT;

            temp.pc3dphi        = pc3dphi;
            temp.pc3dz          = pc3dz;
            temp.pc3sdphi       = pc3sdphi;
            temp.pc3sdz         = pc3sdz;
            temp.emcdphi        = emcdphi;
            temp.emcdz          = emcdz;
            temp.emcsdphi       = emcsdphi;
            temp.emcsdz         = emcsdz;

            track_list.push_back(temp);
        }
}

void He3Jet::GetClusters(PHCompositeNode *topNode, std::vector<clusters>& cluster_list, float vertex)
{
    //The number of towers is 48(y) x 96(z) x 2(PbGl) + 36(y) x 72(z) x 6(PbSc) = 9216+15552 = 24768.
    //Tower ID is a unique ID of all towers in eight sectors (0~24767).
    emcClusterContainer *emcclustercontainer = getClass<emcClusterContainer>(topNode, "emcClusterContainer");
    if (!emcclustercontainer)
        {
            cout << "No emcClusterContainer!  No sense continuing" << endl;
            exit(1);
        }

    //For ERT trigger
    ertOut = getClass<ErtOut>(topNode, "ErtOut");
    if (!ertOut)
        {
            cout << "No ertOut!  No sense continuing" << endl;
            exit(1);
        }

    int Nclus = emcclustercontainer->size();
    for (int iclus = 0; iclus < Nclus; iclus++)
        {
            emcClusterContent* clus = emcclustercontainer->getCluster(iclus);
            int arm         = clus->arm(); //In EMCal convention, West Arm is 0 and East Arm is 1, and thus armsector 0...7 are W0...W3 E0...E3
            int sector      = clus->sector();
            int emcid       = clus->id();
            float energy    = clus->ecore();
            float prob      = clus->prob_photon();
            float chi2      = clus->chi2();

            //need to construct cluster theta and pT
            float x         = clus->x();
            float y         = clus->y();
            float z         = clus->z() - vertex;
            float theta     = clus->theta();//acos(z / sqrt(x * x + y * y + z * z));

            TVector3 v3;
            TLorentzVector v4;
            v3.SetXYZ(x, y, z);
            v3 = energy * v3.Unit();
            v4.SetT(energy);
            v4.SetVect(v3);

            float pT        = v4.Pt();
            float phi       = phiReduce(v4.Phi());
            float px        = v4.Px();
            float py        = v4.Py();
            float pz        = v4.Pz();
            float eta       = v4.Eta();
            float eT        = v4.Et();

            int yTowerPos   = clus->iypos();
            int zTowerPos   = clus->izpos();

            int armsect     = (arm * 4) + sector; //Arm for emcal is different

            int towerId         = EmcMap::GetTowerID(armsect, yTowerPos, zTowerPos);
            bool passIsValid    = EmcMap::IsValidYZ(armsect, yTowerPos, zTowerPos);
            bool passNotBad     = !EmcMap::IsBad(armsect, yTowerPos, zTowerPos);

            //ERT information; ertBit: 0=4x4a, 1=4x4b, 2=4x4c
            int ertSM = (arm == 1 && sector < 2)  ?  (yTowerPos / 12) * 8 + zTowerPos / 12  :  (yTowerPos / 12) * 6 + zTowerPos / 12;
            int ert4x4A = ertOut->get_ERTbit(0, arm, sector, ertSM);
            int ert4x4B = ertOut->get_ERTbit(1, arm, sector, ertSM);
            int ert4x4C = ertOut->get_ERTbit(2, arm, sector, ertSM);

            bool ertFired = ert4x4A || ert4x4B || ert4x4C;

            //Cluster selection
            bool passEnergy = energy > CLUSTER_MIN_ENERGY_CUT;
            bool passEverything = passIsValid && passNotBad && passEnergy;

            if(!passEverything)
                {
                    continue;
                }

            //Additional hot towers- temp soln
            if(towerId == 756 || towerId == 757 || towerId == 828 || towerId == 829 ||
	       towerId == 2417 || towerId == 2418 || towerId == 2457 || towerId == 2490 ||

	       towerId == 2764 || towerId == 2782 || towerId == 2783 || towerId == 2784 ||
	       towerId == 2805 || towerId == 2839 || towerId == 2853 || towerId == 2854 ||
	       towerId == 2855 ||

	       towerId == 5361 || towerId == 5512 || towerId == 5722 || towerId == 5723 ||
	       towerId == 5794 || towerId == 5795 || towerId == 5800 || towerId == 5935 ||
	       towerId == 5948 || towerId == 6019 || towerId == 6023 || towerId == 6024 ||
	       towerId == 6111 ||

	       towerId == 8749 || towerId == 8790 || towerId == 8823 || towerId == 8826 ||
	       towerId == 9176 || towerId == 9549 || towerId == 9762 ||

	       towerId == 10543 || towerId == 11322 || towerId == 11323 || towerId == 11324 ||
	       towerId == 11395 || towerId == 11396 || towerId == 11397 || towerId == 11468 ||
	       towerId == 11469 || towerId == 12863 ||

	       towerId == 13073 || towerId == 13077 || towerId == 13187 || towerId == 13258 ||
	       towerId == 13259 || towerId == 13260 || towerId == 13330 || towerId == 13331 ||
	       towerId == 13332 || towerId == 13591 || towerId == 13592 || towerId == 13594 ||

	       towerId == 17011 || towerId == 17714 || towerId == 17714 || towerId == 18080 ||
	       towerId == 18609 || towerId == 18784 || towerId == 19071 || towerId == 19171 ||
	       towerId == 19349 ||

	       towerId == 20480 || towerId == 21296 || towerId == 22113 || towerId == 22415 ||
	       towerId == 23035 || towerId == 24269 || towerId == 24380 || towerId == 24427 ||
	       towerId == 24521 || towerId == 24522 || towerId == 24523 || towerId == 24618 ||
	       towerId == 24619 || towerId == 24620 || towerId == 24623 ||

	       towerId == 7696 || towerId == 7697 || towerId == 7701 || towerId == 24524)
                {
                    continue;
                }

            hTowerE->Fill(towerId, energy);
            hClustersTower->Fill(towerId);
            hSectorHits[armsect]->Fill(zTowerPos, yTowerPos);

            //Convert EMC arm to DC convention
            if(arm == 0)
                {
                    arm = 1;
                }
            else
                {
                    arm = 0;
                }

            clusters temp;
            temp.arm             = arm;
            temp.sector          = sector;
            temp.armsect         = armsect;
            temp.emcid           = emcid;
            temp.yTowerPos       = yTowerPos;
            temp.zTowerPos       = zTowerPos;
            temp.towerId         = towerId;

            temp.ertFired        = ertFired;

            temp.energy          = energy;
            temp.theta           = theta;
            temp.pT              = pT;
            temp.eT              = eT;
            temp.phi             = phi;
            temp.px              = px;
            temp.py              = py;
            temp.pz              = pz;
            temp.eta             = eta;
            temp.chi2            = chi2;
            temp.prob            = prob;

            cluster_list.push_back(temp);
        }
}

void He3Jet::InitHistograms()
{
    hVertex = new TH1F("hVertex", "Vertex distribution", 300, -75, 75);
    hEvents = new TH1F("hEvents", "Number of events", 5, 0, 5);

    hPC3dphi = new TH1F("hPC3dphi", "PC3dphi distribution", 1600, -0.1, 0.1);
    hPC3dz = new TH1F("hPC3dz", "PC3dz distribution", 800, -50, 50);

    hEMCdphi = new TH1F("hEMCdphi", "EMCdphi distribution", 1600, -0.1, 0.1);
    hEMCdz = new TH1F("hEMCdz", "EMCdz distribution", 800, -50, 50);

    hTowerE = new TH2F("hTowerE", "Ecore vs Tower-id", NTOWER, 0, NTOWER, 40, 0, 40.0);
    hClustersTower = new TH1F("hClustersTower", "Hits vs Tower-id", NTOWER, 0, NTOWER);
    for (int ias = 0; ias < NARMSECT; ias++)
        {
            int ny, nz;
            if (EmcMap::IsPbGl(ias))
                {
                    ny = YPOS_PBGL;
                    nz = ZPOS_PBGL;
                }
            else
                {
                    ny = YPOS_PBSC;
                    nz = ZPOS_PBSC;
                }
            hSectorHits[ias] = new TH2F(Form("hSecorHits_%i", ias), Form("Sector %i", ias), nz, 0, nz, ny, 0, ny);
        }

    hCf = new TH1F("hCf", "Charged fraction", 100, -0.02, 1.02);
    hJets = new TH1F("hJets", "Reco jet pT", NPTBINS, PTBINS);
}

void He3Jet::InitTrees(bool writeTrees)
{
    if(writeTrees)
        {
            //tTracks = new TTree("tTracks", "Central arm tracks");
            //tTracks->Branch("tTracks", &tracks_list);

            //tClusters = new TTree("tClusters", "EMCal clusters");
            //tClusters->Branch("tClusters", &clusters_list);

            tJets = new TTree("tJets", "Anti-kt reconstructed jets");
            tJets->Branch("tJets", &jets_list);
        }
}

void He3Jet::FillTrees(bool writeTrees)
{
    if(writeTrees)
        {
            //tTracks->Fill();
            //tClusters->Fill();
            tJets->Fill();
        }
}















































