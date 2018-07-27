void runCentralityStudy(
			const int nevents = 0,
			const char *simDstFile = "SimDST.root",
			)
{
    gSystem->Load("libfun4all.so");	// framework + reco modules
    gSystem->Load("libPHPythiaEventGen.so");
    gSystem->Load("libPHPythia.so");
    gSystem->Load("librecal.so");
    gSystem->Load("libsimreco.so");	// framework + reco modules

    gSystem->Load("libCentralityStudy.so");

    recoConsts *rc = recoConsts::instance();
    rc->set_IntFlag("RUNNUMBER", 0);

    Fun4AllServer *se = Fun4AllServer::instance();

    MasterRecalibratorManager *mr = new MasterRecalibratorManager();
    se->registerSubsystem(mr);

    CentralityStudy *CentralityStudy = new CentralityStudy();
    se->registerSubsystem(CentralityStudy);

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


