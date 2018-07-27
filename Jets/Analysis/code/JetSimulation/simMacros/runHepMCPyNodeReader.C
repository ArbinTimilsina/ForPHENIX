void runHepMCPyNodeReader(char* outfile = "JetSimWithoutDetector.root",
                          const char* sHijingFile = "sHijing.dat",
                          const char *outputName = "phhijing.root",
                          int nevents = 0)
{
    gSystem->Load("libfun4all.so");        // framework + reco modules
    gSystem->Load("libsimreco.so");     // framework + reco modules

    gSystem->Load("libHepMCPyNodeReader.so");


    recoConsts *rc = recoConsts::instance();
    rc->set_IntFlag("RUNNUMBER", 0);

    /////////////////////////////////////////////////////////////////
    //  Server...
    Fun4AllServer *se = Fun4AllServer::instance();

    SubsysReco *sync = new SyncSimreco();
    se->registerSubsystem(sync);

    HepMCPyNodeReader* hep = new HepMCPyNodeReader();
    hep->Verbosity(0);
    se->registerSubsystem(hep);

    Fun4AllInputManager *in = new Fun4AllHepMCInputManager("DSTIN");
    se->registerInputManager( in );
    in->fileopen(sHijingFile);

    TString OUTPUT = outputName;
    Fun4AllDstOutputManager *dst_output_mgr  = new Fun4AllDstOutputManager("PHHIJING", OUTPUT.Data());
    dst_output_mgr->AddNode("Sync");
    dst_output_mgr->AddNode("PHHijingHeader");
    dst_output_mgr->AddNode("PHHijing");

    se->registerOutputManager(dst_output_mgr);

    if(nevents > 0)
        {
            cout << "Running over " << nevents << " Events" << endl;
        }
    else
        {
            cout << "Running over all Events" << endl;
        }

    se->run(nevents);  // run over all events

    cout << "Calling Fun4AllServer::End()" << endl;
    se->End();

    cout << "Done.  Enjoy your histos." << endl;
}



