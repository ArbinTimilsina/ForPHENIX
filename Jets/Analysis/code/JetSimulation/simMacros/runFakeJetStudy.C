void runFakeJetStudy(const char* sHijingfile = "sHijing.dat",
		     const char *simDstFile = "SimDST.root",
		     int nevents = 0)
{
    gSystem->Load("libfun4all.so");
    gSystem->Load("libfun4allfuncs.so");
    gSystem->Load("librecal.so");
    gSystem->Load("libert.so");

    //For the fastjet package
    gSystem->Load("libCGAL.so");
    gSystem->Load("libfastjet.so");
    gSystem->Load("libfastjettools.so");
    gSystem->Load("libsiscone.so");
    gSystem->Load("libsiscone_spherical.so");

    //For Gaussian Filter
    gSystem->Load("libjetbase.so");
    gSystem->Load("libjetevent.so");
    gSystem->Load("libjetrec.so");

    gSystem->Load("libphhepmc.so");

    gSystem->Load("libJetAnalyzer.so");
    gSystem->Load("libFakeJetStudy.so");

    recoConsts *rc = recoConsts::instance();
    rc->set_IntFlag("RUNNUMBER", 0);

    /////////////////////////////////////////////////////////////////
    //  Server...
    Fun4AllServer *se = Fun4AllServer::instance();

    /////////////////////////////////////////////////////////////////
    //   analyzing compact CNTs need the master recalibrator
    //   to create the CNT in memory
    MasterRecalibratorManager *mr = new MasterRecalibratorManager();
    se->registerSubsystem(mr);

    /////////////////////////////////////////////////////////////////
    //  Reconstruction Modules...
    SubsysReco *FakeJetStudy = new FakeJetStudy();
    FakeJetStudy->Verbosity(1);
    se->registerSubsystem(FakeJetStudy);

    /////////////////////////////////////////////////////////////////
    //  Input Managers...
    /////////////////////////////////////////////////////////////////

    //PISA+Reco TopNode name
    rc->set_CharFlag("PISARECO_TOPNODE", "PISARECO");

    Fun4AllInputManager *in1 = new Fun4AllNoSyncDstInputManager("DSTin1", "DST", "PISARECO");
    se->registerInputManager(in1);
    in1->AddFile(simDstFile);

    Fun4AllHepMCInputManager *in2 = new Fun4AllHepMCInputManager("HEPMC");
    se->registerInputManager(in2);
    in2->fileopen(sHijingfile);

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


