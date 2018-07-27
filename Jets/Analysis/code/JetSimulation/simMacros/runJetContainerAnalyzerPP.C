void runJetContainerAnalyzerPP(
			       const float R             = 0.2,
			       const float minPt         = 6.59,
			       const float ncPythia      = 0.0,
			       const float minCfPythia   = 0.0,
			       const float maxCfPythia   = 1.0,
			       const int nevents         = 0,
			       const char *pythiaFile    = "phpythia.root",
			       const char *simDstFile    = "SimDST.root")
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

    gSystem->Load("libJetTriggerPythia.so");
    gSystem->Load("libJetAnalyzer.so");
    gSystem->Load("libJetContainerAnalyzer.so");

    recoConsts *rc = recoConsts::instance();
    rc->set_IntFlag("RUNNUMBER", 0);

    Fun4AllServer *se = Fun4AllServer::instance();

    MasterRecalibratorManager *mr = new MasterRecalibratorManager();
    se->registerSubsystem(mr);

    JetContainerAnalyzer *jetContainerAnalyzer = new JetContainerAnalyzer(R, minPt,
									  ncPythia, minCfPythia, maxCfPythia);
    jetContainerAnalyzer->SetCuAu(false);
    se->registerSubsystem(jetContainerAnalyzer);

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

    Fun4AllDstOutputManager *dst_output_mgr  = new Fun4AllDstOutputManager("PHPYTHIA", "JetContainerAnalyzer.root");
    dst_output_mgr->AddNode("EventContainerNode");
    dst_output_mgr->AddNode("JetContainerNode");
    dst_output_mgr->AddNode("ParticleContainerNode");
    se->registerOutputManager(dst_output_mgr);

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

    cout << "Done!!" << endl;
}





