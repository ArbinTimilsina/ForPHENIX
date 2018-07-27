void runAlphaVsPhi_CuAu(char* outfile = "/phenix/hhj/arbint/RootFiles/AlphaVsPhi/AlphaVsPhi_CuAu.root", int nevents = 0)
{
    gSystem->Load("libfun4all.so");
    gSystem->Load("librecal.so");
    gSystem->Load("libsimreco.so");

    gSystem->Load("libAlphaVsPhi.so");

    Fun4AllServer *se = Fun4AllServer::instance();

    MasterRecalibratorManager *mr = new MasterRecalibratorManager();
    se->registerSubsystem(mr);

    AlphaVsPhi *AlphaVsPhi = new AlphaVsPhi(outfile);
    AlphaVsPhi->Verbosity(1);
    AlphaVsPhi->SetCuAu(true);
    se->registerSubsystem(AlphaVsPhi);

    Fun4AllDstInputManager *in = new Fun4AllDstInputManager("DSTin1", "DST");
    se->registerInputManager(in);
    in->AddFile("CNT_MB_run12CuAu_200GeV_ZF_ana337-0000372527-9000.root");

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
