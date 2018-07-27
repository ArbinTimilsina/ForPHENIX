void generateHijing(
		    const int nevents = 1000,
		    const char *outputname = "phhijing.root",
		    )
{
    gSystem->Load("libfun4all.so");	// framework + reco modules
    gSystem->Load("libsimreco.so");
    gSystem->Load("libPHHijing.so");

    recoConsts *rc = recoConsts::instance();
    rc->set_IntFlag("RUNNUMBER", 0);

    Fun4AllServer *se = Fun4AllServer::instance();

    SubsysReco *sync = new SyncSimreco();
    se->registerSubsystem(sync);

    PHHijing *phhijing = new PHHijing();
    se->registerSubsystem(phhijing);

    Fun4AllDummyInputManager *in1 = new Fun4AllDummyInputManager("DSTin1", "DST");
    se->registerInputManager(in1);

    // DST output manager
    TString OUTPUT = outputname;
    Fun4AllDstOutputManager *dst_output_mgr  = new Fun4AllDstOutputManager("PHHIJING", OUTPUT.Data());
    dst_output_mgr->AddNode("Sync");
    dst_output_mgr->AddNode("PHHijingHeader");
    dst_output_mgr->AddNode("PHHijing");

    se->registerOutputManager(dst_output_mgr);

    se->run(nevents);
    se->End();
}



