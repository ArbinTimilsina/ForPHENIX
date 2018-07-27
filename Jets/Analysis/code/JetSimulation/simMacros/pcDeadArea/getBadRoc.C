void getBadRoc(int what = 0)
{
    int runNumber = -1;

    if(what == 0)
        {
            runNumber = 372524;
        }
    else if(what == 1)
        {
            runNumber = 360934;
        }
    else
        {
            cout << "Choose 0 for Cu+Au, 1 for p+p" << endl;
            exit(1);
        }

    gSystem->Load("libfun4allfuncs.so");
    gSystem->Load("libpad.so");

    RunToTime *runtotime = RunToTime::instance();
    PHTimeStamp *TS = runtotime->getBeginTime(runNumber);

    PadCalibrationObject* PCO = new PadCalibrationObject();
    PCO->setTimeStamp(*TS);
    PCO->FetchBadROCObjy();

    PCO->setDebugLevel(2);
    PCO->print();
    TString name = "pad_deadroc_";
    name += runNumber;
    name += ".dat";
    PCO->PutBadROCToFile( name );

    delete PCO;

    runtotime->DisconnectDB();
}



