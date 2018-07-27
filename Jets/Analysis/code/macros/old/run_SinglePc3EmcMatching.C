void run_SinglePc3EmcMatching(char* outfile = "Pc3EmcMatching.root", int nevents = 0)
{
  gSystem->Load("libfun4all.so");
  gSystem->Load("librecal.so");
  gSystem->Load("libert.so");
  gSystem->Load("libPc3EmcMatching.so");

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
  
  SubsysReco *Pc3EmcMatching = new Pc3EmcMatching(outfile);
  Pc3EmcMatching->Verbosity(1);
  se->registerSubsystem(Pc3EmcMatching);

  /////////////////////////////////////////////////////////////////
  //  Input Managers...
  Fun4AllDstInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  se->registerInputManager(in1);

  in1->AddFile("CNT_run12_online_central-0000373828-0000.root");

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
