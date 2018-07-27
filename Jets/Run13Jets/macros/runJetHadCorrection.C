void runJetHadCorrection(char* outfile = "JetHadCorrection.root",
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

    gSystem->Load("libJetHadCorrection.so");

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
    SubsysReco *JetHadCorrection = new JetHadCorrection();
    JetHadCorrection->Verbosity(1);
    se->registerSubsystem(JetHadCorrection);

    /////////////////////////////////////////////////////////////////
    //  Input Managers...
    Fun4AllDstInputManager *in = new Fun4AllDstInputManager("DSTin1", "DST");
    se->registerInputManager(in);
    in->AddFile(pythiafile);

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




