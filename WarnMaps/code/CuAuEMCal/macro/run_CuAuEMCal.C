void run_CuAuEMCal(char* outfile = "CuAuEMCal.root", int nevents = 0)
{
  gSystem->Load("libfun4all.so");
  gSystem->Load("librecal.so");
  gSystem->Load("libert.so");
  gSystem->Load("libsimreco.so");

  gSystem->Load("/direct/phenix+u/arbint/WarnMaps/install/lib/libEMCalMap.so");

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
  
  SubsysReco *EMCalMap = new EMCalMap(outfile);
  EMCalMap->Verbosity(1);
  se->registerSubsystem(EMCalMap);

  /////////////////////////////////////////////////////////////////
  //  Input Managers...
  Fun4AllDstInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  se->registerInputManager(in1);
  in1->AddFile("CNT_MB_run12CuAu_200GeV_pro95-0000373139-9005.root");

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
