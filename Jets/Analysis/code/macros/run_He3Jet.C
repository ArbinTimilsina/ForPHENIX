void run_He3Jet(char* outfile = "He3Jet.root", int nevents = 500000)
{
  gSystem->Load("libfun4all.so");
  gSystem->Load("librecal.so");
  gSystem->Load("libsimreco.so");

  gSystem->Load("libCGAL.so");
  gSystem->Load("libfastjet.so");
  gSystem->Load("libfastjettools.so");
  gSystem->Load("libsiscone.so");
  gSystem->Load("libsiscone_spherical.so");

  gSystem->Load("libHe3Jet.so");

  //  gSystem->ListLibraries();

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

  He3Jet *He3Jet= new He3Jet(outfile);
  He3Jet->Verbosity(1);
  se->registerSubsystem(He3Jet);

  /////////////////////////////////////////////////////////////////
  //  Input Managers...
  Fun4AllDstInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  se->registerInputManager(in1);
  in1->AddFile("CNT_ERT_run14HeAu_200GeV_CA_pro102-0000415751-9000.root");

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
