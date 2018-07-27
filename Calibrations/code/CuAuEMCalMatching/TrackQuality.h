namespace TrackQuality
{

  inline bool inBrokenX1(float board, float alpha, float zed, int arm)
    {
      //NE
      if (zed > 0 && arm == 0)
	{
	  //Nothing
	}

      //SE
      if (zed < 0 && arm == 0)
	{
	  //Nothing
	}

      //NW
      if (zed > 0 && arm == 1)
	{
	  if ((alpha > (-0.32*board) + 6.80) && (alpha < (-0.32*board) + 7.12)) return true;
	  if ((alpha > (-0.32*board) + 12.31) && (alpha < (-0.32*board) + 14.05)) return true;
	}

      //SW
      if (zed < 0 && arm == 1)
	{
	  if ((alpha > (-0.30*board) + 5.27) && (alpha < (-0.30*board) + 5.63)) return true;
	  if ((alpha > (-0.30*board) + 6.36) && (alpha < (-0.30*board) + 6.64)) return true;
	  if ((alpha > (-0.30*board) + 11.57) && (alpha < (-0.30*board) + 13.16)) return true;
	  if ((alpha > (-0.30*board) + 16.71) && (alpha < (-0.30*board) + 17.03)) return true;
	  if ((alpha > (-0.30*board) + 17.31) && (alpha < (-0.30*board) + 17.64)) return true;
	}

      return false;
    }


  inline bool inBrokenX2(float board, float alpha, float zed, int arm)
    {
      //NE
      if (zed > 0 && arm == 0)
	{
	  if ((alpha > (-0.44*board) + 17.06) && (alpha < (-0.44*board) + 17.62)) return true;
	}

      //SE
      if (zed < 0 && arm == 0)
	{
	  if ((alpha > (-0.45*board) + 17.13) && (alpha < (-0.45*board) + 17.77)) return true;
	  if ((alpha > (-0.45*board) + 24.28) && (alpha < (-0.45*board) + 24.70)) return true;
	}

      //NW
      if (zed > 0 && arm == 1)
	{
	  if ((alpha > (0.49*board) - 5.16) && (alpha < (0.49*board) - 4.71)) return true;
	  if ((alpha > (0.49*board) - 12.88) && (alpha < (0.49*board) - 12.47)) return true;
	  if ((alpha > (0.49*board) - 16.81) && (alpha < (0.49*board) - 16.41)) return true;
	  if ((alpha > (0.49*board) - 20.30) && (alpha < (0.49*board) - 18.65)) return true;
	  if ((alpha > (0.49*board) - 27.98) && (alpha < (0.49*board) - 27.47)) return true;
	}

      //SW
      if (zed < 0 && arm == 1)
	{
	  if ((alpha > (0.45*board) - 8.43) && (alpha < (0.45*board) - 7.85)) return true;
	  if ((alpha > (0.45*board) - 11.94) && (alpha < (0.45*board) - 11.56)) return true;
	  if ((alpha > (0.45*board) - 15.62) && (alpha < (0.45*board) - 15.11)) return true;
	  if ((alpha > (0.45*board) - 18.81) && (alpha < (0.45*board) - 17.13)) return true;
	  if ((alpha > (0.45*board) - 25.94) && (alpha < (0.45*board) - 25.42)) return true;
	}

      return false;
    }


  inline bool inBrokenUV(float board, float alpha, float zed, int arm)
    {
      //NE
      if (zed > 0 && arm == 0)
	{
	  //Nothing
	}

      //SE
      if (zed < 0 && arm == 0)
	{
	  //Nothing
	}

      //NW
      if (zed > 0 && arm == 1)
	{
	  if (board > 37.0 && board < 43.0) return true;
	}

      //SW
      if (zed < 0 && arm == 1)
	{
	  if ((alpha > (0.57*board) - 23.79 ) && (alpha > (-0.91*board) + 34.69)) return true;
	}

      return false;
    }


  inline bool passQualityMask(int quality, float board, float alpha, float zed, int arm)
    {
      //(quality & 1) == 0 ==> no X1 bit
      //(quality & 1) == 1 ==> X1 bit

      if ((quality & 1) == 0 && !inBrokenX1(board, alpha, zed, arm))
	{
	  return false;
	}

      if ((quality & 2) == 0 && !inBrokenX2(board, alpha, zed, arm))
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

      if ((quality & 12) == 0 && !inBrokenUV(board, alpha, zed, arm))
	{
	  return false;
	}

      return true;
    }

};
