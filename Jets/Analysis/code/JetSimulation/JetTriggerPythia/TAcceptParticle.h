#ifndef __TACCEPTPARTICLE_H__
#define __TACCEPTPARTICLE_H__

#include <TLorentzVector.h>

class TAcceptParticle
{

private:
  double line1(double phi);

public:
  TAcceptParticle();
  ~TAcceptParticle(){;}

  int acceptParticle(TLorentzVector *part, int charge, double vertex, 
		     bool &phiAccept, bool &zedAccept, bool &bigzedAccept, bool &deadArea );

  void SetFullTwoPi(bool in){fullTwoPi = in;}

protected:

  double RadiusDC;

  double PhiBotE;
  double PhiTopE;
  double PhiBotW;
  double PhiTopW;

  double Sl1;
  double Sl2;

  bool fullTwoPi; 

};

#endif  /* __TACCEPTPARTICLE_HH__ */


