#ifndef __EMCMAP_H__
#define __EMCMAP_H__

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <TROOT.h>
#include "TOAD.h"
#include <cmath>

namespace EmcMap
{
#define NARMSECT 8
#define YPOS_PBGL 48
#define YPOS_PBSC 36
#define ZPOS_PBGL 96
#define ZPOS_PBSC 72
#define NTOWER 24768

#define MASK_EDGE            0x01
#define MASK_AROUND_DEAD     0x02
#define MASK_AROUND_HOT      0x04
#define MASK_DEAD            0x08
#define MASK_HOT             0x10
#define MASK_AROUND_UNCALIB  0x20
#define MASK_UNCALIB         0x40

    static int warnmap[NARMSECT][YPOS_PBGL][ZPOS_PBGL];

    void ReadWarnMap(std::string);

    int IsBad(int, int, int);
    int IsHot(int, int, int);
    int IsDead(int, int, int);
    int IsUncalib(int, int, int);
    int IsAroundHot(int, int, int);
    int IsAroundDead(int, int, int);
    int IsAroundUncalib(int, int, int);
    int IsPbGl(int);
    int IsEdge(int, int, int);
    int IsValidYZ(int, int, int);
    int GetTowerID(int, int, int);

    inline void ReadWarnMap(std::string mapFile)
    {
	TOAD *toad_loader = new TOAD("JetAnalyzer");
	std::string emcmap = toad_loader->location(mapFile);

	std::cout << std::endl;
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "Loading EMCal warnmap from file " << emcmap << std::endl;
	std::cout << "***********************************************************************" << std::endl << std::endl;

	std::ifstream fmap;
	fmap.open(emcmap.c_str());

	if (!fmap)
	    {
		std::cout << "ERROR!!! File " << emcmap << " could not be located!" << std::endl;
	    }

	std::string line;
	memset(warnmap, 0, sizeof(warnmap));
	while (getline(fmap, line))
	    {
		std::istringstream one_line(line);
		int as, y, z, status;
		if (one_line >> as >> y >> z >> status)
		    {
			if (IsValidYZ(as, y, z))
			    {
				warnmap[as][y][z] = status;
			    }
		    }
	    }

	fmap.close();
	delete toad_loader;
    }

    inline int IsBad(int as, int y, int z)
    {
	return ( IsHot(as, y, z) || IsDead(as, y, z) || IsUncalib(as, y, z) ||
		 IsAroundHot(as, y, z) || IsAroundDead(as, y, z) ||
		 IsAroundUncalib(as, y, z) || IsEdge(as, y, z));
    };

    inline int IsHot(int as, int y, int z)
    {
	return (warnmap[as][y][z] & MASK_HOT);
    };

    inline int IsDead(int as, int y, int z)
    {
	return (warnmap[as][y][z] & MASK_DEAD);
    };

    inline int IsUncalib(int as, int y, int z)
    {
	return (warnmap[as][y][z] & MASK_UNCALIB);
    };

    inline int IsAroundHot(int as, int y, int z)
    {
	return (warnmap[as][y][z] & MASK_AROUND_HOT);
    };

    inline int IsAroundDead(int as, int y, int z)
    {
	return (warnmap[as][y][z] & MASK_AROUND_DEAD);
    };

    inline int IsAroundUncalib(int as, int y, int z)
    {
	return (warnmap[as][y][z] & MASK_AROUND_UNCALIB);
    };

    inline int IsEdge(int as, int y, int z)
    {
	return (warnmap[as][y][z] & MASK_EDGE);
    };

    inline int IsPbGl(int armsect)
    {
	return ((armsect == 4 || armsect == 5) ? 1 : 0);
    }

    inline int IsValidYZ(int as, int y, int z)
    {
	int ret = 0;
	if (IsPbGl(as))
	    {
		if (y >= 0 && y < YPOS_PBGL && z >= 0 && z < ZPOS_PBGL)
		    {
			ret = 1;
		    }
	    }
	else
	    {
		if (y >= 0 && y < YPOS_PBSC && z >= 0 && z < ZPOS_PBSC)
		    {
			ret = 1;
		    }
	    }
	return ret;
    }

    inline int GetTowerID(int as, int y, int z)
    {
	int towerid = -1;
	if (0 <= as && as <= 3)
	    {
		towerid = 2592 * as + 72 * y + z;
	    }
	else if (as == 4 || as == 5)
	    {
		towerid = 15552 + 4608 * (as - 4) + 96 * y + z;
        }
    else if (as == 6 || as == 7)
        {
            towerid = 10368 + 2592 * (as - 6) + 72 * y + z;
        }
    return towerid;


}

}

#endif

