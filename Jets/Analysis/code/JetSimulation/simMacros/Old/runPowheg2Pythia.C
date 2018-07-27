void runPowheg2Pythia(const char* powhegFile = "pwgevents.lhe",
                      char* pythiaFile = "phpythia.root",
                      const int nevents = 0)
{
    gSystem->Load("libfun4all.so"); // framework + reco modules
    gSystem->Load("libPowheg2Pythia.so");
    gSystem->Load("libsimreco.so");

    Fun4AllServer *se = Fun4AllServer::instance();

    SubsysReco *sync = new SyncSimreco();
    se->registerSubsystem(sync);

    Powheg2Pythia *ana = new Powheg2Pythia(powhegFile);
    se->registerSubsystem(ana);

    Fun4AllDummyInputManager *in1 = new Fun4AllDummyInputManager("DSTin1", "DST");
    se->registerInputManager(in1);

    Fun4AllDstOutputManager *out1 = new Fun4AllDstOutputManager("DSTin1", pythiaFile);
    out1->AddNode("Sync");
    out1->AddNode("PHPythiaHeader");
    out1->AddNode("PHPythia");
    se->registerOutputManager(out1);

    // run over all events
    se->run(nevents);
    se->End();
}

