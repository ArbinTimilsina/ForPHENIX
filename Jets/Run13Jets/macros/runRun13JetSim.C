void runRun13JetSim(const float minPt = 5.0, const int nevents = 0,
                    const char *pythiaFile = "phpythia.root",
                    const char *simDstFile = "SimDST.root")
{
    gSystem->Load("libfun4all.so");	// framework + reco modules
    gSystem->Load("libPHPythiaEventGen.so");
    gSystem->Load("libPHPythia.so");
    gSystem->Load("librecal.so");
    gSystem->Load("libsimreco.so");	// framework + reco modules

    //For fastjet package
    gSystem->Load("libCGAL.so");
    gSystem->Load("libfastjet.so");
    gSystem->Load("libfastjettools.so");
    gSystem->Load("libsiscone.so");
    gSystem->Load("libsiscone_spherical.so");

    gSystem->Load("libRun13JetTriggerPythia.so");
    gSystem->Load("libRun13Jet.so");
    gSystem->Load("libRun13JetSim.so");

    recoConsts *rc = recoConsts::instance();
    rc->set_IntFlag("RUNNUMBER", 0);

    Fun4AllServer *se = Fun4AllServer::instance();

    MasterRecalibratorManager *mr = new MasterRecalibratorManager();
    se->registerSubsystem(mr);

    Run13JetSim *Run13JetSim = new Run13JetSim(minPt);
    se->registerSubsystem(Run13JetSim);

    //Input
    rc->set_CharFlag("PHPYTHIA_TOPNODE", "PHPYTHIA");
    Fun4AllDstInputManager *in1 = new Fun4AllNoSyncDstInputManager("DSTin1", "DST" , "PHPYTHIA");
    se->registerInputManager(in1);
    in1->AddFile(pythiaFile);

    //PISA+Reco TopNode name
    rc->set_CharFlag("PISARECO_TOPNODE", "PISARECO");
    Fun4AllInputManager *in2 = new Fun4AllNoSyncDstInputManager("DSTin2", "DST", "PISARECO");
    se->registerInputManager(in2);
    in2->AddFile(simDstFile);

    if(nevents > 0)
        {
            cout << "Running over " << nevents << " Events" << endl;
        }
    else
        {
            cout << "Running over all Events" << endl;
        }

    se->run(nevents);

    cout << "Calling Fun4AllServer::End()" << endl;
    se->End();

    cout << "Done.  Enjoy your histos." << endl;
}




