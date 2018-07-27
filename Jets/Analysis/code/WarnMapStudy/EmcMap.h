#ifndef __EMCMAP_H__
#define __EMCMAP_H__

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <TROOT.h>

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

void ReadWarnMap();

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

#endif
