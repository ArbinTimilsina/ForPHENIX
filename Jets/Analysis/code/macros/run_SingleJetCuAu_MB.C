void run_SingleJetCuAu_MB(char* outfile = "JetAnalyzerCuAu_MB.root", int nevents = 100000)
{
  gSystem->Load("libfun4all.so");
  gSystem->Load("librecal.so");
  gSystem->Load("libsimreco.so");

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
  
  JetAnalyzer *JetAnalyzer = new JetAnalyzer(outfile);
  JetAnalyzer->Verbosity(1);
  JetAnalyzer->SetData(true);
  JetAnalyzer->SetCuAu(true);
  JetAnalyzer->SetMB(true);
  se->registerSubsystem(JetAnalyzer);

  /////////////////////////////////////////////////////////////////
  //  Input Managers...
  Fun4AllDstInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  se->registerInputManager(in1);
  in1->AddFile("CNT_MB_run12CuAu_200GeV_CA_pro99-0000374003-9000.root");
  //in1->AddFile("CNT_MB_run12CuAu_200GeV_CA_pro99-0000375428-9000.root");

  //in1->AddFile("CNT_MB_run12CuAu_200GeV_CA_pro99-0000373780-9025.root");

  //Fun4AllDstInputManager *in2 = new Fun4AllDstInputManager("DSTin2", "DST_EVE");
  //se->registerInputManager(in2);
  //in2->AddFile("DST_EVE_MB_run12CuAu_200GeV_CA_pro99-0000374003-9000.root");

  if(nevents>0){
    cout << "Running over " << nevents << " Events" << endl;
  }
  else{
    cout << "Running over all Events" << endl;
  }

  se->run(nevents);  // run over all events                                                                        

  cout << "Calling Fun4AllServer::End()" << endl;
  se->End();

  cout << "Done.  Enjoy your histos." << endl;
}
