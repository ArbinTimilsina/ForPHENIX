namespace TrackQualityCuAu
{
  inline float getBoard(float phi, int arm)
    {
      float board = -9999;
      if (arm == 0) board = (3.72402 - phi + 0.008047 * cos(phi + 0.87851)) / 0.01963496;
      if (arm == 1) board = (0.573231 + phi - 0.0046 * cos(phi + 0.05721)) / 0.01963496;

      return board;
    }

  inline bool inBrokenX1(float phi, float alpha, float zed, int arm)
    {
      float board = getBoard(phi,arm);

      //NE
      if (zed > 0 && arm == 0)
	{
	  if ((alpha > (0.30*board) - 10.65) && (alpha < (0.30*board) - 9.48)) return true;
	  if ((alpha > (0.29*board) - 7.44) && (alpha < (0.32*board) - 7.82)) return true;
	}

      //SE
      if (zed < 0 && arm == 0)
	{
	  if ((alpha > (0.30*board) - 10.65) && (alpha < (0.30*board) - 9.48)) return true;
	}

      //NW
      if (zed > 0 && arm == 1)
	{
	  if ((alpha > (-0.31*board) + 6.28) && (alpha < (-0.33*board) + 7.19)) return true;
	  if ((alpha > (-0.34*board) + 13.03) && (alpha < (-0.34*board) + 14.85)) return true;
	}

      //SW
      if (zed < 0 && arm == 1)
	{
	  if ((alpha > (-0.30*board) + 5.15) && (alpha < (-0.30*board) + 5.49)) return true;
	  if ((alpha > (-0.30*board) + 6.22) && (alpha < (-0.34*board) + 7.31)) return true;
	  if ((alpha > (-0.33*board) + 12.56) && (alpha < (-0.33*board) + 14.07)) return true;
	  if ((alpha > (-0.32*board) + 17.69) && (alpha < (-0.32*board) + 18.13)) return true;
	  if ((alpha > (-0.32*board) + 18.34) && (alpha < (-0.35*board) + 20.39)) return true;
	}

      return false;
    }


  inline bool inBrokenX2(float phi, float alpha, float zed, int arm)
    {
      float board = getBoard(phi,arm);

      //NE
      if (zed > 0 && arm == 0)
	{
	  if ((alpha > (-0.46*board) + 17.48) && (alpha < (-0.48*board) + 19.11)) return true;
	}

      //SE
      if (zed < 0 && arm == 0)
	{
	  if ((alpha > (-0.35*board) + 12.63) && (alpha < (-0.90*board) + 32.39)) return true;
	  if ((alpha > (-0.44*board) + 16.84) && (alpha < (-0.49*board) + 19.47)) return true;
	  if ((alpha > (-0.44*board) + 23.95) && (alpha < (-0.52*board) + 28.74)) return true;
	}

      //NW
      if (zed > 0 && arm == 1)
	{
	  if ((alpha > (0.48*board) - 4.89) && (alpha < (0.57*board) - 5.14)) return true;
	  if ((alpha > (0.44*board) - 11.53) && (alpha < (0.54*board) - 13.52)) return true;
	  if ((alpha > (0.46*board) - 15.61) && (alpha < (0.55*board) - 18.11)) return true;
	  if ((alpha > (0.49*board) - 20.42) && (alpha < (0.70*board) - 25.95)) return true;
	  if ((alpha > (0.46*board) - 26.17) && (alpha < (0.55*board) - 30.55)) return true;
	}

      //SW
      if (zed < 0 && arm == 1)
	{
	  if ((alpha > (0.44*board) - 8.02) && (alpha < (0.51*board) - 8.70)) return true;
	  if ((alpha > (0.46*board) - 12.10) && (alpha < (0.54*board) - 13.45)) return true;
	  if ((alpha > (0.46*board) - 15.85) && (alpha < (0.52*board) - 17.03)) return true;
	  if ((alpha > (0.46*board) - 19.03) && (alpha < (0.46*board) - 17.18)) return true;
	  if ((alpha > (0.43*board) - 24.73) && (alpha < (0.51*board) - 28.23)) return true;
	}

      return false;
    }


  inline bool inBrokenUV(float phi, float alpha, float zed, int arm)
    {
      float board = getBoard(phi,arm);

      //NE
      if (zed > 0 && arm == 0)
	{
	  //Nothing
	}

      //SE
      if (zed < 0 && arm == 0)
	{
	  if ((alpha > (0.25*board) - 14.52) && (alpha > (-0.19*board) + 11.27)) return true;
	}

      //NW
      if (zed > 0 && arm == 1)
	{
	  if (board > 36.0 && board < 43.0) return true;
	}

      //SW
      if (zed < 0 && arm == 1)
	{
	  if (board > 37.0 && board < 43.0) return true;
	}

      return false;
    }


  inline bool passQualityMask(int quality, float phi, float alpha, float zed, int arm)
    {
      //(quality & 1) == 0 ==> no X1 bit
      //(quality & 1) == 1 ==> X1 bit

      if((quality & 1) == 0 && (quality & 2) == 0){
	return false;
      }

      if ((quality & 1) == 0 && !inBrokenX1(phi, alpha, zed, arm))
	{
	  return false;
	}

      if ((quality & 2) == 0 && !inBrokenX2(phi, alpha, zed, arm))
	{
	  return false;
	}

      if ((quality & 16) == 0)
        {
	  //Require PC to be found
          return false;
        }

      if ((quality & 32) == 0 && (quality & 12) == 0)
	{
	  // PC not unique- so UV needs to be unique
	  return false;
	}

      if ((quality & 12) == 0 && !inBrokenUV(phi, alpha, zed, arm))
	{
	  return false;
	}

      //To match with simulation
      float board = getBoard(phi,arm);
      if ((arm ==1) && (alpha > (0.47*board) - 19.29) && (alpha > (-0.41*board) + 15.39)){
	return false;
      }
      return true;
    }

};
