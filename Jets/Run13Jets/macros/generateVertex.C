#include <TROOT.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TRandom3.h>

using namespace std;

void generateVertex(int nevents = 10)
{
    ofstream vertexFile;
    vertexFile.open("vertex.txt");

    int events = 0;
    while(events < nevents)
        {
            TRandom3 *myRandom3 = new TRandom3();
            myRandom3->SetSeed(0);
            float zvertex = myRandom3->Gaus(0.0, 17);

            if(fabs(zvertex) <= 30.0)
                {
                    vertexFile << events << " " << zvertex << " " << 0.0 << "\n";
                    events++;
                }
        }
    vertexFile.close();
}


