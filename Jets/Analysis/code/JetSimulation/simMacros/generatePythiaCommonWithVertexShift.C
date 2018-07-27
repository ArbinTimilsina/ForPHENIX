void generatePythiaCommonWithVertexShift(
					 const int nevents = 10,
					 const char *outputname = "phpythia.root"
					 )
{
    gSystem->Load("libfun4all.so");	// framework + reco modules
    gSystem->Load("libPHPythiaEventGen.so");
    gSystem->Load("libPHPythia.so");
    gSystem->Load("libsimreco.so");	// framework + reco modules

    recoConsts *rc = recoConsts::instance();
    rc->set_IntFlag("RUNNUMBER", 0);

    Fun4AllServer *se = Fun4AllServer::instance();

    SubsysReco *sync = new SyncSimreco();
    se->registerSubsystem(sync);

    PHPythia *phpythia = new PHPythia();
    se->registerSubsystem(phpythia);

    //Set vertex shift
    PHPyVertexShift *vtxShift = new PHPyVertexShift("PHPyVertexShift");
    vtxShift->SetVtxFile("vertex.txt");
    vtxShift->SetVtxFileVersion(PHPyVertexShift::VERSION1);
    se->registerSubsystem(vtxShift);

    Fun4AllDummyInputManager *in1 = new Fun4AllDummyInputManager("DSTin1", "DST");
    se->registerInputManager(in1);

    // DST output manager
    TString OUTPUT = outputname;
    Fun4AllDstOutputManager *dst_output_mgr  = new Fun4AllDstOutputManager("PHPYTHIA", OUTPUT.Data());
    dst_output_mgr->AddNode("Sync");
    dst_output_mgr->AddNode("PHPythiaHeader");
    dst_output_mgr->AddNode("PHPythia");

    se->registerOutputManager(dst_output_mgr);

    se->run(nevents);
    se->End();
}




