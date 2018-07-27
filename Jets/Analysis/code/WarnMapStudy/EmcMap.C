#include "EmcMap.h"
#include "TOAD.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void ReadWarnMap()
{
  TOAD *toad_loader = new TOAD("WarnMapStudy");
  string emcmap= toad_loader->location("warnmap.txt");
  cout << "Loading EMCal warnmap from file " << emcmap << endl;

  ifstream fmap;
  fmap.open(emcmap.c_str());

  if (!fmap)
    {
      cout << "ERROR!!! File " << emcmap << " could not be located!" << endl;
    }

  string line;
  memset(warnmap,0,sizeof(warnmap));
  while (getline(fmap, line))
    {
      istringstream one_line(line);
      int as, y, z, status;
      if (one_line >> as >> y >> z >> status)
        {
	    if (IsValidYZ(as, y, z)) {
	      warnmap[as][y][z] = status;
	    }
	}
    }

  fmap.close();
  delete toad_loader;
  cout << "Successfully EMCal warnmap from  file " << emcmap << endl << endl;
}

int IsBad(int as, int y, int z)
{
  return ( IsHot(as, y, z) || IsDead(as, y, z) || IsUncalib(as, y, z) ||
           IsAroundHot(as, y, z) || IsAroundDead(as, y, z) ||
           IsAroundUncalib(as, y, z) || IsEdge(as, y, z));
};

int IsHot(int as, int y, int z)
{
  return (warnmap[as][y][z] & MASK_HOT);
};

int IsDead(int as, int y, int z)
{
  return (warnmap[as][y][z] & MASK_DEAD);
};

int IsUncalib(int as, int y, int z)
{
  return (warnmap[as][y][z] & MASK_UNCALIB);
};
int IsAroundHot(int as, int y, int z)
{
  return (warnmap[as][y][z] & MASK_AROUND_HOT);
};

int IsAroundDead(int as, int y, int z)
{
  return (warnmap[as][y][z] & MASK_AROUND_DEAD);
};

int IsAroundUncalib(int as, int y, int z)
{
  return (warnmap[as][y][z] & MASK_AROUND_UNCALIB);
};

int IsEdge(int as, int y, int z)
{
  return (warnmap[as][y][z] & MASK_EDGE);
};

int IsPbGl(int armsect)
{
  return ((armsect == 4 || armsect == 5) ? 1 : 0);
}

int IsValidYZ(int as, int y, int z)
{
  int ret = 0;
  if (IsPbGl(as))
    {
      if (y >= 0 && y < YPOS_PBGL && z >= 0 && z < ZPOS_PBGL) ret = 1;
    }
  else
    {
      if (y >= 0 && y < YPOS_PBSC && z >= 0 && z < ZPOS_PBSC)   ret = 1;
    }
  return ret;
}

int GetTowerID(int as, int y, int z)
{
  int towerid = -1;
  if (0 <= as && as <= 3) { towerid = 2592*as + 72*y + z; }
  else if (as == 4 || as == 5) { towerid = 15552 + 4608*(as-4) + 96*y + z; }
  else if (as == 6 || as == 7) { towerid = 10368 + 2592*(as-6) + 72*y + z; }
  return towerid;
}
