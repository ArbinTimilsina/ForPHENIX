void pisaToDST(Int_t nEvents = 0, int runNumber = 393483, char *filein = "PISAEvent.root", char *dstout = "SimDST.root")
{
    cout << "Running over " << nEvents << " events" << endl;
    cout << "Inupt file:  " << filein << endl;
    cout << "Output file: " << dstout << endl;

    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Running pisaToDST.C" << endl;
    cout << "Number of events is: " << nEvents << endl;
    cout << "Reference Run number is: " << runNumber << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    // Load libraries
    gSystem->Load("libfun4all");
    gSystem->Load("libmutoo_subsysreco");
    gSystem->Load("libfun4allfuncs");
    gSystem->Load("libsimreco");

    gROOT->ProcessLine(".L pisaToDST_IOManager.C");

    ///////////////////////////////////////////
    // recoConsts setup
    //////////////////////////////////////////
    recoConsts *rc = recoConsts::instance();

    // 2 means PISA-To-DST
    rc->set_IntFlag("SIMULATIONFLAG", 2);

    // Disable embedding
    rc->set_IntFlag("EMBEDFLAG", 0);

    // Reference Run number
    rc->set_IntFlag("RUNNUMBER", runNumber);

    // Requested by EMCal
    rc->set_IntFlag("EVALUATIONFLAG", 1);

    // Needed for TofSimreco
    rc->set_FloatFlag("TOFTIMINGRESOLUTION", 0.100);

    /*
      Flags to abort event if required number of GEANT hits is not present in the subsystem
      Defaults are all 0 except for the Drift Chamber (default setting for DC is 3 wire plane hits)
    */
    rc->set_IntFlag("DCHREQFLAG", 3);
    rc->set_IntFlag("PC1REQFLAG", 0);
    rc->set_IntFlag("PC2REQFLAG", 0);
    rc->set_IntFlag("PC3REQFLAG", 0);
    rc->set_IntFlag("TOFREQFLAG", 0);

    // Not yet operational
    rc->set_IntFlag("EMCREQFLAG", 0);

    //For AFS files, 0 uses AFS files, 1 uses local files
    //This flag is used by DchSimreco and not by PadSimreco
    rc->set_IntFlag("AFSABSENT", 1);

    //For PC dead area
    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Loading PC dead area" << endl;
    rc->set_CharFlag("PADDEADCHFILE", "pad_deadch.dat");
    rc->set_CharFlag("PADDEADROCFILE", "pad_deadroc.dat");
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    // Don't know what the defult values are- so set everything to 0.0
    // simVertexFlag = 0 means that the BBC Z0 value will be used
    // simVertexFlag = 1 means that the same simZ0Vertex centroid value is used for all events
    // simVertexFlag = 2 means that the Z0 centroid is taken from the PISA event header for each event
    // The centroid values are modified by the Width values which are Gaussian sigma values
    rc->set_IntFlag("SIMVERTEXFLAG", 2);
    rc->set_FloatFlag("SIMZ0VERTEX", 0.0);
    rc->set_FloatFlag("SIMT0VERTEX", 0.0);
    rc->set_FloatFlag("SIMZ0VERTEXWIDTH", 0.0);
    rc->set_FloatFlag("SIMT0VERTEXWIDTH", 0.0);

    ///////////////////////////////////////////
    // Make the Server
    //////////////////////////////////////////
    Fun4AllServer *se = Fun4AllServer::instance();

    ///////////////////////////////////////////
    // Activate the subsystems
    //////////////////////////////////////////


    // Run header and trigger setting
    se->registerSubsystem( new HeadSimreco() );
    se->registerSubsystem( new TrigSimreco() );

    // BBC simReco
    se->registerSubsystem(new BbcSimreco("BBC"));

    //This should happen after BBCSimreco
    //The input vertex is smeared for BBC vertex resolution
    //https://www.phenix.bnl.gov/WWW/offline/wikioff/index.php/VERTEX_PISA_TO_DST
    VtxSimreco* vtx_sim = new VtxSimreco();
    vtx_sim->UseVtx( VtxSimreco::PISA );
    vtx_sim->SmearZ( true );
    vtx_sim->ZVertexSigma(0.5);//The vertex is smeared along z, with a Gaussian distribution of width 5mm
    vtx_sim->UseXY(false);
    vtx_sim->OverwriteBBC(true);
    se->registerSubsystem(vtx_sim);

    // t0
    T0Simreco* t0_sim = new T0Simreco();
    t0_sim->T0Sigma(0.04);
    se->registerSubsystem( t0_sim );

    // Pad chambers
    se->registerSubsystem(new PadSimreco("PAD"));

    // Pad Vertexing Code, maybe needed for multiple vertices
    se->registerSubsystem(new PadVtxReco("PADVTX"));

    // The VtxReco works unchanged for both real and simulation events
    se->registerSubsystem(new VtxReco("VTX"));

    // The T0Reco works unchanged for both real and simulation events
    se->registerSubsystem(new T0Reco());

    // As of January 2, 2004 the Dch has uninitialized variable warnings from Valgrind
    // There are also log file output warning messages
    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Loading DC dead area" << endl;
    se->registerSubsystem( new DchSimreco("DCH") );
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

    // Time expansion chamber
    se->registerSubsystem( new TecSimreco("TEC"));

    // Time of flight detector
    se->registerSubsystem(new TofSimreco("TOF"));

    // Tof west
    se->registerSubsystem(new TfwSimreco("TFW"));

    // HBD
    //se->registerSubsystem(new HbdSimreco("HBD"));

    // RICH
    se->registerSubsystem(new CrkSimreco("CRK"));

    // Aerogel
    se->registerSubsystem(new AccSimreco("ACC"));
    se->registerSubsystem(new AccReco());

    //EMCal
    //Feb 21, 2014: According to Gabor "high" should be 8 times "low": maybe the correct numbers are 0.001 and 0.008
    rc->set_FloatFlag("EMCTOWERLOWGAIN", 0.0015625);
    rc->set_FloatFlag("EMCTOWERHIGHGAIN", 0.0125);
    se->registerSubsystem( new EmcReco3() );

    // The CglReco works unchanged for both real and simulation events
    se->registerSubsystem(new CglReco("CGL"));

    //Aerogel cluster  (Needs to be after cglRec)
    se->registerSubsystem(new AccclusterReco());

    //  This is the class which makes the RICH Ring data structure
    se->registerSubsystem( new RingReco() );

    // This is the class which makes the Central Tracks nanoDST output
    // 22 corresponds to the version used in pro.78 for Run7 Au+Au
    se->registerSubsystem(new CentraltrackReco( 22 ));

    //  This is the class which makes the GlobalEvent data on the nanoDST output
    se->registerSubsystem(new GlobalReco());

    // This is the class which checks for charged particles going into EMCal
    se->registerSubsystem(new ChargedvetoReco());

    //added the DC based global evaluation module
    se->registerSubsystem( new McEvalSimreco() );

    // dumper
    if( false )
        {
            gSystem->Load( "libphnodedump" );
            se->registerSubsystem( new Dumper() );
        }

    ///////////////////////////////////////////
    // InputManager
    ///////////////////////////////////////////
    Fun4AllInputManager *input_manager = new Fun4AllPisaInputManager("PisaIn", "TOP");
    se->registerInputManager(input_manager);

    ///////////////////////////////////////////
    // OutputManagers Set up functions
    ///////////////////////////////////////////
    if( dstout )
        {
            DST_IOManager(dstout, se);
        }

    ///////////////////////////////////////////
    // Open input file
    se->fileopen(input_manager->Name(), filein);

    // Process input events
    gBenchmark->Start("eventLoop");
    se->run(nEvents);
    se->End();
    gBenchmark->Show("eventLoop");

    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Reconstruction completed!!!" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;
}







