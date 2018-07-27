void runJetConstituentPerformance(
				  const float R                = 0.2,
				  const float minPt            = 10.0,
				  const float ncPythia         = 3.0,
				  const float minCfPythia      = 0.0,
				  const float maxCfPythia      = 1.0,
				  const float minChargedDeltaR = 0.03,
				  const float minNeutralDeltaR = 0.03,
				  const int nevents            = 0,
				  const bool perfectEMCal      = false,
				  const bool phenixParticle    = false,
				  const char *pythiaFile       = "phpythia.root",
				  const char *simDstFile       = "SimDST.root")
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

    gSystem->Load("libJetTriggerPythia.so");
    gSystem->Load("libJetAnalyzer.so");
    gSystem->Load("libJetConstituentPerformance.so");

    recoConsts *rc = recoConsts::instance();
    rc->set_IntFlag("RUNNUMBER", 0);

    Fun4AllServer *se = Fun4AllServer::instance();

    MasterRecalibratorManager *mr = new MasterRecalibratorManager();
    se->registerSubsystem(mr);

    JetConstituentPerformance *JetConstituentPerformance = new JetConstituentPerformance(R, minPt,
											 ncPythia, minCfPythia, maxCfPythia,
											 minChargedDeltaR, minNeutralDeltaR,
											 perfectEMCal, phenixParticle);
    se->registerSubsystem(JetConstituentPerformance);

    //PYTHIA TopNode name
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









