void runJetSimAcceptance(
			 const int nevents = 0,
			 const bool isCuAu = true,
			 const char *simDstFile = "SimDST.root",
			 )
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

    //For Gaussian Filter
    gSystem->Load("libjetbase.so");
    gSystem->Load("libjetevent.so");
    gSystem->Load("libjetrec.so");


    gSystem->Load("libJetAnalyzer.so");
    gSystem->Load("libJetSimAcceptance.so");

    recoConsts *rc = recoConsts::instance();
    rc->set_IntFlag("RUNNUMBER", 0);

    Fun4AllServer *se = Fun4AllServer::instance();

    MasterRecalibratorManager *mr = new MasterRecalibratorManager();
    se->registerSubsystem(mr);

    JetSimAcceptance *JetSimAcceptance = new JetSimAcceptance();
    JetSimAcceptance->SetCuAu(isCuAu);
    se->registerSubsystem(JetSimAcceptance);

    Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
    se->registerInputManager(in1);
    in1->AddFile(simDstFile);

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


