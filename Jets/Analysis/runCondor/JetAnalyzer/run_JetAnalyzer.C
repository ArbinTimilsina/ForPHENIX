void run_JetAnalyzer(char* outfile = "JetAnalyzer.root", int nevents = 0)
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
  in1->AddListFile("filelist.list");



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
