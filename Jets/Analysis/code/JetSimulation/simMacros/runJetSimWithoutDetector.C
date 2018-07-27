void runJetSimWithoutDetector(char* outfile = "JetSimWithoutDetector.root",
                              const char* sHijingfile = "sHijing.dat",
                              const char* pythiafile = "phpythia.root",
                              int nevents = 0)
{
    gSystem->Load("libfun4all.so");
    gSystem->Load("libfun4allfuncs.so");
    gSystem->Load("librecal.so");
    gSystem->Load("libert.so");

    gSystem->Load("libPHPythiaEventGen.so");
    gSystem->Load("libPHPythia.so");

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

    gSystem->Load("libJetSimWithoutDetector.so");

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
    SubsysReco *JetSimWithoutDetector = new JetSimWithoutDetector(outfile);
    JetSimWithoutDetector->Verbosity(1);
    se->registerSubsystem(JetSimWithoutDetector);

    /////////////////////////////////////////////////////////////////
    //  Input Managers...
    Fun4AllNoSyncDstInputManager *in1 = new Fun4AllNoSyncDstInputManager("DSTin1", "DST");
    se->registerInputManager(in1);
    in1->AddFile(pythiafile);

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


