void getDeadCh(int what = 0)
{
    gSystem->Load("libpad.so");

    TFile *infile;
    int runNumber = -1;
    double HITS_THRESHOLD = -99.99;
    if(what == 0)
        {
            infile = new TFile("/direct/phenix+hhj/arbint/ExtraFiles/Run_372524_PadCal.root", "READ");
            runNumber = 372524;
            HITS_THRESHOLD = 20.0;
        }
    else if(what == 1)
        {
            infile = new TFile("/direct/phenix+hhj/arbint/ExtraFiles/Run_360934_PadCal.root", "READ");
            runNumber = 360934;
            HITS_THRESHOLD = 10.0;
        }
    else
        {
            cout << "Choose 0 for Cu+Au, 1 for p+p" << endl;
            exit(1);
        }

    PadAddressObject *addressObj = new PadAddressObject;
    TH2F *padFemPadxPadz[80] = {0};

    TH2F *pad1w = new TH2F("pad1w", "PC1.W padx vs padz", 216, 0, 216, 160, 0, 160);
    TH2F *pad1e = new TH2F("pad1e", "PC1.E padx vs padz", 216, 0, 216, 160, 0, 160);
    TH2F *deadpad1w = new TH2F("deadpad1w", "PC1.W padx vs padz deadmap", 216, 0, 216, 160, 0, 160);
    TH2F *deadpad1e = new TH2F("deadpad1e", "PC1.E padx vs padz deadmap", 216, 0, 216, 160, 0, 160);

    TString hname;
    for (int ihist = 0; ihist < 80; ihist++)
        {
            hname = "padFemPadxPadz";
            hname += ihist;
            padFemPadxPadz[ihist] = (TH2F*)infile->Get( hname );
        }

    const int NPADZ = 108;
    const int NPADX = 20;

    ofstream badchfile("tempdead");

    //  TH1F *hpadhits = new TH1F("hpadhits", "num pad hits", 10000, 0, 10000);
    TH1F *hpadhits = new TH1F("hpadhits", "num pad hits", 1000, 0, 1000);
    hpadhits->SetLineColor(2);
    int nbad_padchannels = 0;

    for (int ihist = 0; ihist < 80; ihist++)
        {
            if ( padFemPadxPadz[ihist] == 0 )
                {
                    continue;
                }

            int packet_id = ihist + 4001;
            if ( ihist > 48 )
                {
                    packet_id += 16;
                }
            cout << "pktid " << ihist << "\t" << packet_id << endl;

            for (int ipadz = 0; ipadz < 108; ipadz++)
                {
                    for (int ipadx = 0; ipadx < 20; ipadx++)
                        {
                            Double_t npadhits = padFemPadxPadz[ihist]->GetBinContent( ipadz + 1, ipadx + 1 );
                            hpadhits->Fill( npadhits );
                            int ch_id = addressObj->getChannelid(0, ipadz, ipadx);

                            int side = addressObj->getSide(packet_id);
                            int gpadz = side * 108 + ipadz;
                            int gpadx = ((packet_id - 4001) % 8) * 20 + ipadx;

                            if ( npadhits < HITS_THRESHOLD && packet_id != 4010 )
                                {
                                    badchfile << packet_id << "\t" << ch_id << "\t" << 1 << endl;
                                    nbad_padchannels++;

                                    if ( packet_id <= 4016 )
                                        {
                                            deadpad1w->Fill( gpadz, gpadx );
                                        }
                                    else if ( packet_id > 4016 && packet_id <= 4032 )
                                        {
                                            deadpad1e->Fill( gpadz, gpadx );
                                        }
                                }

                            // fill 2d hist of whole sector
                            if ( packet_id <= 4016 )
                                {
                                    pad1w->Fill( gpadz, gpadx, npadhits );
                                }
                            else if ( packet_id > 4016 && packet_id <= 4032 )
                                {
                                    pad1e->Fill( gpadz, gpadx, npadhits );
                                }
                        }
                }
        }

    badchfile.close();

    // deadpad file must have the number of bad channels on 1st line
    TString dead = "pad_deadch_";
    dead += runNumber;
    dead += ".dat";
    ofstream final_badchfile(dead);
    final_badchfile << nbad_padchannels << endl;
    final_badchfile.close();

    TString command = "cat tempdead >> " + dead;
    gSystem->Exec(command);

    TCanvas *ac = new TCanvas("ac", "ac", 550, 425);
    hpadhits->Draw();
    gPad->SetLogy(1);

    TCanvas *bc = new TCanvas("c_pad1w", "pad1w", 425, 550);
    gStyle->SetPalette(1);
    pad1w->Draw("colz");
    gPad->SetLogz(1);

    TCanvas *cc = new TCanvas("c_pad1e", "pad1e", 425, 550);
    pad1e->Draw("colz");
    gPad->SetLogz(1);

    TCanvas *dc = new TCanvas("c_deadpad1w", "pad1w dead", 425, 550);
    deadpad1w->Draw("colz");
    gPad->SetLogz(1);

    TCanvas *ec = new TCanvas("c_deadpad1e", "pad1e dead", 425, 550);
    deadpad1e->Draw("colz");
    gPad->SetLogz(1);
}




