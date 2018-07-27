#include <iostream>

namespace TrackQualityPP
{
  enum arms {EAST, WEST}; //this is DC (PHENIX) convention, not EMC convention
  static const unsigned int X1_USED = 1;
  static const unsigned int X2_USED = 2;
  static const unsigned int UV_FOUND_AND_UNIQUE = 12;
  static const unsigned int PC_FOUND = 16;
  static const unsigned int PC_UNIQUE = 32;
  static const unsigned int PC_FOUND_AND_UNIQUE = 48;

  static const float TOP_SLICE_ALPHA = 0.21;
  static const float BOTTOM_SLICE_ALPHA = -0.21;
  static const float delta_y = TOP_SLICE_ALPHA - BOTTOM_SLICE_ALPHA;

  static std::vector<float> x1_NE_slopes;
  static std::vector<float> x1_NE_intercepts;
  static std::vector<float> x1_NW_slopes;
  static std::vector<float> x1_NW_intercepts;
  static std::vector<float> x1_SE_slopes;
  static std::vector<float> x1_SE_intercepts;
  static std::vector<float> x1_SW_slopes;
  static std::vector<float> x1_SW_intercepts;

  static std::vector<float> x2_NE_slopes;
  static std::vector<float> x2_NE_intercepts;
  static std::vector<float> x2_NW_slopes;
  static std::vector<float> x2_NW_intercepts;
  static std::vector<float> x2_SE_slopes;
  static std::vector<float> x2_SE_intercepts;
  static std::vector<float> x2_SW_slopes;
  static std::vector<float> x2_SW_intercepts;

  inline std::vector<float> CalcSlopes(std::vector<float> top, std::vector<float> bot)
  {
    std::vector<float> slopes_temp;
    if(top.size() != bot.size())
      {
        std::cout << "TrackQuality.h : Top and Bottom vector sizes not equal... crashing..." << std::endl;
      }

    for(UInt_t i = 0; i < top.size(); i++)
      {
	float slope = delta_y / (top.at(i) - bot.at(i));
        slopes_temp.push_back(slope);
        //std::cout << "Top-side board: " << top.at(i) << ", slope: " << slope << std::endl;
      }

    return slopes_temp;

  }

  inline std::vector<float> CalcIntercept(std::vector<float> top, std::vector<float> slope)
  {
    std::vector<float> intercepts_temp;
    if(top.size() != slope.size())
      {
        std::cout << "TrackQuality.h : Top and Bottom vector sizes not equal... crashing..." << std::endl;
      }

    for(UInt_t i = 0; i < top.size(); i++)
      {
        float intercept = TOP_SLICE_ALPHA - slope.at(i) * top.at(i);
        intercepts_temp.push_back(intercept);
        //std::cout << "Top-side board: " << top.at(i) << ", intercept: " << intercept << std::endl;
      }

    return intercepts_temp;

  }


  inline void InitTrackQuality()
  {
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "Executing InitTrackQuality() in TrackQualityPP.h" << std::endl;
    std::cout << "***********************************************************************" << std::endl<<std::endl;

    //RunGroup1
    //x1_NE top    25.9-26.3
    //x1_NE bottom 24.3-24.8
    std::vector<float> x1_NE_top;
    std::vector<float> x1_NE_bot;
    x1_NE_top.push_back(32.3);
    x1_NE_top.push_back(36.2);
    x1_NE_top.push_back(25.9);
    x1_NE_top.push_back(26.3);
    x1_NE_bot.push_back(30.9);
    x1_NE_bot.push_back(34.8);
    x1_NE_bot.push_back(24.3);
    x1_NE_bot.push_back(24.8);
    // std::cout << std::endl << "x1_NE: " << std::endl;
    x1_NE_slopes = CalcSlopes(x1_NE_top, x1_NE_bot);
    x1_NE_intercepts = CalcIntercept(x1_NE_top, x1_NE_slopes);

    //x1 NW top    10.3-23.0, 30.5-34.6, 37.3-42.6, 56.4-57.4, 68.9-69.4
    //x1 NW bottom 11.9-24.4, 31.9-35.7, 38.8-43.8, 55.8-56.4, 69.7-70.4
    std::vector<float> x1_NW_top;
    std::vector<float> x1_NW_bot;
    x1_NW_top.push_back(10.3);
    x1_NW_top.push_back(23.0);
    x1_NW_top.push_back(30.5);
    x1_NW_top.push_back(34.6);
    x1_NW_top.push_back(37.3);
    x1_NW_top.push_back(42.6);
    x1_NW_bot.push_back(11.9);
    x1_NW_bot.push_back(24.4);
    x1_NW_bot.push_back(31.9);
    x1_NW_bot.push_back(35.7);
    x1_NW_bot.push_back(38.8);
    x1_NW_bot.push_back(43.8);
    // std::cout << std::endl << "x1_NW: " << std::endl;
    x1_NW_slopes = CalcSlopes(x1_NW_top, x1_NW_bot);
    x1_NW_intercepts = CalcIntercept(x1_NW_top, x1_NW_slopes);

    //x1 SE top    --
    //x1 SE bottom --
    std::vector<float> x1_SE_top;
    std::vector<float> x1_SE_bot;
    x1_SE_top.push_back(32.3);
    x1_SE_top.push_back(36.2);
    x1_SE_bot.push_back(30.9);
    x1_SE_bot.push_back(34.8);

    // std::cout << std::endl << "x1_SE: " << std::endl;
    x1_SE_slopes = CalcSlopes(x1_SE_top, x1_SE_bot);
    x1_SE_intercepts = CalcIntercept(x1_SE_top, x1_SE_slopes);

    //x1 SW top    10.3-23.0, 30.5-34.6, 37.3-42.6, 56.4-57.4, 69.9-69.4
    //x1 SW bottom 11.9-24.4, 31.9-35.7, 38.8-43.8, 55.8-56.4, 69.7-70.4
    std::vector<float> x1_SW_top;
    std::vector<float> x1_SW_bot;
    x1_SW_top.push_back(10.3);
    x1_SW_top.push_back(23.0);
    x1_SW_top.push_back(30.5);
    x1_SW_top.push_back(34.6);
    x1_SW_top.push_back(37.3);
    x1_SW_top.push_back(42.6);
    x1_SW_bot.push_back(11.9);
    x1_SW_bot.push_back(24.4);
    x1_SW_bot.push_back(31.9);
    x1_SW_bot.push_back(35.7);
    x1_SW_bot.push_back(38.8);
    x1_SW_bot.push_back(43.8);
    // std::cout << std::endl << "x1_SW: " << std::endl;
    x1_SW_slopes = CalcSlopes(x1_SW_top, x1_SW_bot);
    x1_SW_intercepts = CalcIntercept(x1_SW_top, x1_SW_slopes);

    //x2 NE top    38.0-39.6
    //x2 NE bottom 39.0-40.3
    std::vector<float> x2_NE_top;
    std::vector<float> x2_NE_bot;
    x2_NE_top.push_back(38.0);
    x2_NE_top.push_back(39.6);
    x2_NE_bot.push_back(39.0);
    x2_NE_bot.push_back(40.3);
    // std::cout << std::endl << "x2_NE: " << std::endl;
    x2_NE_slopes = CalcSlopes(x2_NE_top, x2_NE_bot);
    x2_NE_intercepts = CalcIntercept(x2_NE_top, x2_NE_slopes);

    //x2 NW top    7.5-12.7, 25.6-26.6, 33.4-34.6, 37.5-41.5, 56.4-57.4, 59.4-63.5, 67.3-71.3
    //x2 NW bottom 6.8-11.2, 24.8-25.6, 32.6-33.7, 37.2-40.7, 55.3-56.5, 58.6-62.5, 66.6-70.6
    std::vector<float> x2_NW_top;
    std::vector<float> x2_NW_bot;
    x2_NW_top.push_back(7.5);
    x2_NW_top.push_back(12.7);
    x2_NW_top.push_back(25.6);
    x2_NW_top.push_back(26.6);
    x2_NW_top.push_back(33.4);
    x2_NW_top.push_back(34.6);
    x2_NW_top.push_back(37.5);
    x2_NW_top.push_back(41.5);
    x2_NW_top.push_back(56.4);
    x2_NW_top.push_back(57.4);
    x2_NW_top.push_back(59.4);
    x2_NW_top.push_back(63.5);
    x2_NW_top.push_back(67.3);
    x2_NW_top.push_back(71.3);
    x2_NW_bot.push_back(6.8);
    x2_NW_bot.push_back(11.2);
    x2_NW_bot.push_back(24.8);
    x2_NW_bot.push_back(25.6);
    x2_NW_bot.push_back(32.6);
    x2_NW_bot.push_back(33.7);
    x2_NW_bot.push_back(37.2);
    x2_NW_bot.push_back(40.7);
    x2_NW_bot.push_back(55.3);
    x2_NW_bot.push_back(56.5);
    x2_NW_bot.push_back(58.6);
    x2_NW_bot.push_back(62.5);
    x2_NW_bot.push_back(66.6);
    x2_NW_bot.push_back(70.6);
    // std::cout << std::endl << "x2_NW: " << std::endl;
    x2_NW_slopes = CalcSlopes(x2_NW_top, x2_NW_bot);
    x2_NW_intercepts = CalcIntercept(x2_NW_top, x2_NW_slopes);

    //x2 SE top    38.1-39.6, 54.1-55.1
    //x2 SE bottom 39.1-40.4, 55.1-55.9
    std::vector<float> x2_SE_top;
    std::vector<float> x2_SE_bot;
    x2_SE_top.push_back(38.1);
    x2_SE_top.push_back(39.6);
    x2_SE_top.push_back(54.1);
    x2_SE_top.push_back(55.1);
    x2_SE_bot.push_back(39.1);
    x2_SE_bot.push_back(40.4);
    x2_SE_bot.push_back(55.1);
    x2_SE_bot.push_back(55.9);
    // std::cout << std::endl << "x2_SE: " << std::endl;
    x2_SE_slopes = CalcSlopes(x2_SE_top, x2_SE_bot);
    x2_SE_intercepts = CalcIntercept(x2_SE_top, x2_SE_slopes);

    //x2 SW top    7.6-12.7, 25.6-26.6, 33.4-34.6, 37.3-41.5, 56.4-57.4, 59.4-63.5, 67.3-71.3
    //x2 SW bottom 6.8-11.2, 24.8-25.6, 32.6-33.7, 37.2-40.7, 55.3-56.5, 58.6-62.5, 66.6-70.6
    std::vector<float> x2_SW_top;
    std::vector<float> x2_SW_bot;
    x2_SW_top.push_back(7.6);
    x2_SW_top.push_back(12.7);
    x2_SW_top.push_back(25.6);
    x2_SW_top.push_back(26.6);
    x2_SW_top.push_back(33.4);
    x2_SW_top.push_back(34.6);
    x2_SW_top.push_back(37.3);
    x2_SW_top.push_back(41.5);
    x2_SW_top.push_back(56.4);
    x2_SW_top.push_back(57.4);
    x2_SW_top.push_back(59.4);
    x2_SW_top.push_back(63.5);
    x2_SW_top.push_back(67.3);
    x2_SW_top.push_back(71.3);
    x2_SW_bot.push_back(6.8);
    x2_SW_bot.push_back(11.2);
    x2_SW_bot.push_back(24.8);
    x2_SW_bot.push_back(25.6);
    x2_SW_bot.push_back(32.6);
    x2_SW_bot.push_back(33.7);
    x2_SW_bot.push_back(37.2);
    x2_SW_bot.push_back(40.7);
    x2_SW_bot.push_back(55.3);
    x2_SW_bot.push_back(56.5);
    x2_SW_bot.push_back(58.6);
    x2_SW_bot.push_back(62.5);
    x2_SW_bot.push_back(66.6);
    x2_SW_bot.push_back(70.6);
    // std::cout << std::endl << "x2_SW: " << std::endl;
    x2_SW_slopes = CalcSlopes(x2_SW_top, x2_SW_bot);
    x2_SW_intercepts = CalcIntercept(x2_SW_top, x2_SW_slopes);


    //uv top
    //uv bottom
  }

inline float getBoard(float phi, int arm)
{
    float board = -9999;
    if (arm == 0) board = (3.72402 - phi + 0.008047 * cos(phi + 0.87851)) / 0.01963496;
    if (arm == 1) board = (0.573231 + phi - 0.0046 * cos(phi + 0.05721)) / 0.01963496;

    return board;
}

  inline bool inBrokenX1(float phi, float alpha, float zed, int arm)
  {
    float board = getBoard(phi, arm);
    std::vector<float> temp_vector_slopes;
    std::vector<float> temp_vector_intercepts;

    if (zed > 0 && arm == WEST)
      {
        // std::cout << "TrackQualiy.h inBrokenX1 North West" << std::endl;
        temp_vector_slopes = x1_NW_slopes;
        temp_vector_intercepts = x1_NW_intercepts;
      }
    if (zed > 0 && arm == EAST)
      {
        // std::cout << "TrackQualiy.h inBrokenX1 North East" << std::endl;
        temp_vector_slopes = x1_NE_slopes;
        temp_vector_intercepts = x1_NE_intercepts;
      }
    if (zed < 0 && arm == WEST)
      {
        // std::cout << "TrackQualiy.h inBrokenX1 South West" << std::endl;
        temp_vector_slopes = x1_SW_slopes;
        temp_vector_intercepts = x1_SW_intercepts;
      }
    if (zed < 0 && arm == EAST)
      {
        // std::cout << "TrackQualiy.h inBrokenX1 South East" << std::endl;
        temp_vector_slopes = x1_SE_slopes;
        temp_vector_intercepts = x1_SE_intercepts;
      }

    if(temp_vector_slopes.size() < 1) return false;

    for(UInt_t i = 0; i < temp_vector_slopes.size() - 1; i++)
      {
        if(temp_vector_slopes.at(i) > 0)
	  {
            if(alpha < temp_vector_slopes.at(i)*board + temp_vector_intercepts.at(i)
	       && alpha > temp_vector_slopes.at(i + 1)*board + temp_vector_intercepts.at(i + 1))
	      {
                return true;
	      }
	  }//positive slope
        else if(temp_vector_slopes.at(i) < 0)
	  {
            if(alpha > temp_vector_slopes.at(i)*board + temp_vector_intercepts.at(i)
	       && alpha < temp_vector_slopes.at(i + 1)*board + temp_vector_intercepts.at(i + 1))
	      {
                return true;
	      }
	  }//negative slope
        else
	  {
            std::cout << "Slope == 0, returning false... but this is bad (TrackQuality.h)" << std::endl;
            return false;
	  }

        i++;

      }//loop over vector of slopes

    return false;
  }//inBrokenX1

  inline bool inBrokenX2(float phi, float alpha, float zed, int arm)
  {
    float board = getBoard(phi, arm);

    std::vector<float> temp_vector_slopes;
    std::vector<float> temp_vector_intercepts;

    if (zed > 0 && arm == WEST)
      {
        temp_vector_slopes = x2_NW_slopes;
        temp_vector_intercepts = x2_NW_intercepts;
      }
    if (zed > 0 && arm == EAST)
      {
        temp_vector_slopes = x2_NE_slopes;
        temp_vector_intercepts = x2_NE_intercepts;
      }
    if (zed < 0 && arm == WEST)
      {
        temp_vector_slopes = x2_SW_slopes;
        temp_vector_intercepts = x2_SW_intercepts;
      }
    if (zed < 0 && arm == EAST)
      {
        temp_vector_slopes = x2_SE_slopes;
        temp_vector_intercepts = x2_SE_intercepts;
      }

    if(temp_vector_slopes.size() < 1) return false;

    for(UInt_t i = 0; i < temp_vector_slopes.size() - 1; i++)
      {
        if(temp_vector_slopes.at(i) > 0)
	  {
            if(alpha < temp_vector_slopes.at(i)*board + temp_vector_intercepts.at(i)
	       && alpha > temp_vector_slopes.at(i + 1)*board + temp_vector_intercepts.at(i + 1))
	      {
                return true;
	      }
	  }//positive slope
        else if(temp_vector_slopes.at(i) < 0)
	  {
            if(alpha > temp_vector_slopes.at(i)*board + temp_vector_intercepts.at(i)
	       && alpha < temp_vector_slopes.at(i + 1)*board + temp_vector_intercepts.at(i + 1))
	      {
                return true;
	      }
	  }//negative slope
        else
	  {
            std::cout << "Slope == 0, returning false... but this is bad (TrackQuality.h)" << std::endl;
            return false;
	  }

        i++;

      }//loop over vector of slopes

    return false;
  }//inBrokenX2


  inline bool inBrokenUV(float phi, float alpha, float zed, int arm)
  {
    float board = getBoard(phi, arm);

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
    if (zed > 0 && arm == WEST)
      {
        if ((alpha > (-0.4 * board) + 4.34 ) && (alpha > (0.4 * board) - 4.9)) return true;
        if(board > 37. && board < 42.) return true;
        if ((alpha < (0.2857 * board) - 70293) && (alpha < (-0.8 * board) + 19.4)) return true;
      }

    //SW
    if (zed < 0 && arm == WEST)
      {
        if ((alpha > (-0.294 * board) + 3.265 ) && (alpha > (0.3875 * board) - 4.7528)) return true;
        if ((alpha < (0.4667 * board) - 8.0 ) && (alpha > (0.4667 * board) - 8.517)) return true;
        if(board > 38.6 && board < 40.5) return true;
        if(board > 40.5 && (alpha > (0.5 * board) - 20.6 ) && board < 42.0) return true;
      }

    return false;
  }//inBrokenUV

  inline bool passQualityMask(int quality, float phi, float alpha, float zed, int arm)
{
    //(quality & 1) == 0 ==> no X1 bit
    //(quality & 1) == 1 ==> X1 bit

    if((quality & 1) == 0 && (quality & 2) == 0)
    {
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
