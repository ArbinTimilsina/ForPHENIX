void generateData(
		  const float R          = 0.2,
		  const float nc         = 3.0,
		  const float minPt      = 10.4,
		  const float minCf      = 0.2,
		  const float maxCf      = 0.7,
		  const int nevents      = 0,
		  const int nskipevents  = 1,
		  const char *outputname = "DataDST.root"
		  )
{
    gSystem->Load("libfun4all.so");
    gSystem->Load("librecal.so");
    gSystem->Load("libsimreco.so");

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
    gSystem->Load("libJetTriggerData.so");

    Fun4AllServer *se = Fun4AllServer::instance();

    MasterRecalibratorManager *mr = new MasterRecalibratorManager();
    se->registerSubsystem(mr);

    SubsysReco *JetTriggerData = new JetTriggerData(R, nc, minPt, minCf, maxCf, nevents);
    se->registerSubsystem(JetTriggerData);

    Fun4AllDstInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
    se->registerInputManager(in1);
    in1->AddListFile("filelist.list");
    in1->run(1);
    in1->skip(nskipevents);

    // DST output manager
    TString OUTPUT = outputname;
    Fun4AllDstOutputManager *manager  = new Fun4AllDstOutputManager("RUN12CUAU", OUTPUT.Data());
    se->registerOutputManager(manager);

    se->run(nevents * 1e6);
    se->End();
}





