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

  static std::vector<float> uv_NE_slopes;
  static std::vector<float> uv_NE_intercepts;
  static std::vector<float> uv_NW_slopes;
  static std::vector<float> uv_NW_intercepts;
  static std::vector<float> uv_SE_slopes;
  static std::vector<float> uv_SE_intercepts;
  static std::vector<float> uv_SW_slopes;
  static std::vector<float> uv_SW_intercepts;

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
    //x1_NE top    25.77-26.18
    //x1_NE bottom 24.36-24.75
    std::vector<float> x1_NE_top;
    std::vector<float> x1_NE_bot;
    x1_NE_top.push_back(25.77);
    x1_NE_top.push_back(26.18);
    
    x1_NE_bot.push_back(24.36);
    x1_NE_bot.push_back(24.75);

    x1_NE_slopes = CalcSlopes(x1_NE_top, x1_NE_bot);
    x1_NE_intercepts = CalcIntercept(x1_NE_top, x1_NE_slopes);

    //x1 NW top    20.26-20.88 | 22.58-26.77 | 37.36-38.43
    //x1 NW bottom 21.53-22.11 | 23.95-27.75 | 38.70-39.71
    std::vector<float> x1_NW_top;
    std::vector<float> x1_NW_bot;
    x1_NW_top.push_back(20.26);
    x1_NW_top.push_back(20.88);
    x1_NW_top.push_back(22.58);
    x1_NW_top.push_back(26.77);
    x1_NW_top.push_back(37.36);
    x1_NW_top.push_back(38.43);

    x1_NW_bot.push_back(21.53);
    x1_NW_bot.push_back(22.11);
    x1_NW_bot.push_back(23.95);
    x1_NW_bot.push_back(27.75);
    x1_NW_bot.push_back(38.70);
    x1_NW_bot.push_back(39.71);

    // std::cout << std::endl << "x1_NW: " << std::endl;
    x1_NW_slopes = CalcSlopes(x1_NW_top, x1_NW_bot);
    x1_NW_intercepts = CalcIntercept(x1_NW_top, x1_NW_slopes);

    //x1 SE top    -- NOTHING TO DO YET1
    //x1 SE bottom --
    //std::vector<float> x1_SE_top;
    //std::vector<float> x1_SE_bot;
    //x1_SE_top.push_back(32.3);
    //x1_SE_top.push_back(36.2);
    //x1_SE_bot.push_back(30.9);
    //x1_SE_bot.push_back(34.8);

    // std::cout << std::endl << "x1_SE: " << std::endl;
    //x1_SE_slopes = CalcSlopes(x1_SE_top, x1_SE_bot);
    //x1_SE_intercepts = CalcIntercept(x1_SE_top, x1_SE_slopes);

    //x1 SW top    16.56-17.63 | 20.22-20.91 | 22.52-26.92 | 37.28-38.45 | 54.42-55.27 | 56.43-57.23
    //x1 SW bottom 18.02-18.73 | 21.55-22.18 | 23.90-27.97 | 38.71-39.79 | 55.73-56.58 | 57.86-58.58
    std::vector<float> x1_SW_top;
    std::vector<float> x1_SW_bot;
    x1_SW_top.push_back(16.56);
    x1_SW_top.push_back(17.63);
    x1_SW_top.push_back(20.22);
    x1_SW_top.push_back(20.91);
    x1_SW_top.push_back(22.52);
    x1_SW_top.push_back(26.92);
    x1_SW_top.push_back(37.28);
    x1_SW_top.push_back(38.45);
    x1_SW_top.push_back(54.42);
    x1_SW_top.push_back(55.27);
    x1_SW_top.push_back(56.43);
    x1_SW_top.push_back(57.23);

    x1_SW_bot.push_back(18.02);
    x1_SW_bot.push_back(18.73);
    x1_SW_bot.push_back(21.55);
    x1_SW_bot.push_back(22.18);
    x1_SW_bot.push_back(23.90);
    x1_SW_bot.push_back(27.97);
    x1_SW_bot.push_back(38.71);
    x1_SW_bot.push_back(39.79);
    x1_SW_bot.push_back(55.73);
    x1_SW_bot.push_back(56.58);
    x1_SW_bot.push_back(57.86);
    x1_SW_bot.push_back(58.58);

    // std::cout << std::endl << "x1_SW: " << std::endl;
    x1_SW_slopes = CalcSlopes(x1_SW_top, x1_SW_bot);
    x1_SW_intercepts = CalcIntercept(x1_SW_top, x1_SW_slopes);

    //x2 NE top    38.2-39.5
    //x2 NE bottom 39.1-40.3
    std::vector<float> x2_NE_top;
    std::vector<float> x2_NE_bot;
    x2_NE_top.push_back(38.2);
    x2_NE_top.push_back(39.5);
    x2_NE_bot.push_back(39.1);
    x2_NE_bot.push_back(40.3);
    // std::cout << std::endl << "x2_NE: " << std::endl;
    x2_NE_slopes = CalcSlopes(x2_NE_top, x2_NE_bot);
    x2_NE_intercepts = CalcIntercept(x2_NE_top, x2_NE_slopes);

    //x2 NW top    9.65-10.52 | 25.56-26.5  | 33.62-34.44 | 38.11-41.45 | 56.50-57.27 | 59.61-67.81
    //x2 NW bottom 8.88-9.52  | 24.79-25.57 | 32.81-33.51 | 37.36-40.55 | 55.68-56.39 | 58.81-66.83
    std::vector<float> x2_NW_top;
    std::vector<float> x2_NW_bot;
    x2_NW_top.push_back(9.65);
    x2_NW_top.push_back(10.52);
    x2_NW_top.push_back(25.56);
    x2_NW_top.push_back(26.5);
    x2_NW_top.push_back(33.62);
    x2_NW_top.push_back(34.44);
    x2_NW_top.push_back(38.11);
    x2_NW_top.push_back(41.45);
    x2_NW_top.push_back(56.50);
    x2_NW_top.push_back(57.27);
    x2_NW_top.push_back(59.61);
    x2_NW_top.push_back(67.81);

    x2_NW_bot.push_back(8.88);
    x2_NW_bot.push_back(9.52);
    x2_NW_bot.push_back(24.79);
    x2_NW_bot.push_back(25.57);
    x2_NW_bot.push_back(32.81);
    x2_NW_bot.push_back(33.51);
    x2_NW_bot.push_back(37.36);
    x2_NW_bot.push_back(40.55);
    x2_NW_bot.push_back(55.68);
    x2_NW_bot.push_back(56.39);
    x2_NW_bot.push_back(58.81);
    x2_NW_bot.push_back(66.83);
    
    // std::cout << std::endl << "x2_NW: " << std::endl;
    x2_NW_slopes = CalcSlopes(x2_NW_top, x2_NW_bot);
    x2_NW_intercepts = CalcIntercept(x2_NW_top, x2_NW_slopes);

    //x2 SE top    38.20-39.44 | 54.20-54.90
    //x2 SE bottom 39.11-40.19 | 55.10-55.77
    std::vector<float> x2_SE_top;
    std::vector<float> x2_SE_bot;
    x2_SE_top.push_back(38.20);
    x2_SE_top.push_back(39.44);
    x2_SE_top.push_back(54.20);
    x2_SE_top.push_back(54.90);

    x2_SE_bot.push_back(39.11);
    x2_SE_bot.push_back(40.19);
    x2_SE_bot.push_back(55.10);
    x2_SE_bot.push_back(55.77);
    // std::cout << std::endl << "x2_SE: " << std::endl;
    x2_SE_slopes = CalcSlopes(x2_SE_top, x2_SE_bot);
    x2_SE_intercepts = CalcIntercept(x2_SE_top, x2_SE_slopes);

    //x2 SW top    17.58-18.52 | 25.53-26.55 | 33.58-34.45 | 38.10-41.46 | 56.58-57.40 | 59.52-67.68
    //x2 SW bottom 16.85-17.64 | 24.71-25.60 | 32.80-33.57 | 37.33-40.53 | 55.75-56.54 | 58.65-66.79
    std::vector<float> x2_SW_top;
    std::vector<float> x2_SW_bot;
    x2_SW_top.push_back(17.58);
    x2_SW_top.push_back(18.52);
    x2_SW_top.push_back(25.53);
    x2_SW_top.push_back(26.55);
    x2_SW_top.push_back(33.58);
    x2_SW_top.push_back(34.45);
    x2_SW_top.push_back(38.10);
    x2_SW_top.push_back(41.46);
    x2_SW_top.push_back(56.58);
    x2_SW_top.push_back(57.40);
    x2_SW_top.push_back(59.52);
    x2_SW_top.push_back(67.68);
   
    x2_SW_bot.push_back(16.85);
    x2_SW_bot.push_back(17.64);
    x2_SW_bot.push_back(24.71);
    x2_SW_bot.push_back(25.60);
    x2_SW_bot.push_back(32.80);
    x2_SW_bot.push_back(33.57);
    x2_SW_bot.push_back(37.33);
    x2_SW_bot.push_back(40.53);
    x2_SW_bot.push_back(55.75);
    x2_SW_bot.push_back(56.54);
    x2_SW_bot.push_back(58.65);
    x2_SW_bot.push_back(66.79);
   
    // std::cout << std::endl << "x2_SW: " << std::endl;
    x2_SW_slopes = CalcSlopes(x2_SW_top, x2_SW_bot);
    x2_SW_intercepts = CalcIntercept(x2_SW_top, x2_SW_slopes);


    //uv NE top    NOTHING
    //uv NE bottom 
    // std::vector<float> uv_NE_top;
    //std::vector<float> uv_NE_bot;

    //uv NW top    25.58-26.63 | 36.79-41.88 | 63.98-68.89
    //uv NW bottom 24.77-25.61 | 36.98-42.24 | 63.89-67.89
    std::vector<float> uv_NW_top;
    std::vector<float> uv_NW_bot;
    uv_NW_top.push_back(25.58);
    uv_NW_top.push_back(26.63);
    uv_NW_top.push_back(36.79);
    uv_NW_top.push_back(41.88);
    uv_NW_top.push_back(63.98);
    uv_NW_top.push_back(68.89);

    uv_NW_bot.push_back(24.77);
    uv_NW_bot.push_back(25.61);
    uv_NW_bot.push_back(36.98);
    uv_NW_bot.push_back(42.24);
    uv_NW_bot.push_back(63.89);
    uv_NW_bot.push_back(67.89);

    uv_NW_slopes = CalcSlopes(uv_NW_top, uv_NW_bot);
    uv_NW_intercepts = CalcIntercept(uv_NW_top, uv_NW_slopes);

    //uv SE top    57.05-58.10 | 72.86-74.13
    //uv SE bottom 56.90-57.84 | 72.75-73.85
    std::vector<float> uv_SE_top;
    std::vector<float> uv_SE_bot;
    uv_SE_top.push_back(57.05);
    uv_SE_top.push_back(58.10);
    uv_SE_top.push_back(72.86);
    uv_SE_top.push_back(74.13);

    uv_SE_bot.push_back(56.90);
    uv_SE_bot.push_back(57.84);
    uv_SE_bot.push_back(72.75);
    uv_SE_bot.push_back(73.85);

    uv_SE_slopes = CalcSlopes(uv_SE_top, uv_SE_bot);
    uv_SE_intercepts = CalcIntercept(uv_SE_top, uv_SE_slopes);

    //uv SW top    25.44-26.62 | 37.38-40.49 | 66.97-67.89
    //uv SW bottom 24.65-25.70 | 37.87-40.68 | 66.35-67.21
    std::vector<float> uv_SW_top;
    std::vector<float> uv_SW_bot;
    uv_SW_top.push_back(25.44);
    uv_SW_top.push_back(26.62);
    uv_SW_top.push_back(37.38);
    uv_SW_top.push_back(40.49);
    uv_SW_top.push_back(66.97);
    uv_SW_top.push_back(67.89);

    uv_SW_bot.push_back(24.65);
    uv_SW_bot.push_back(25.70);
    uv_SW_bot.push_back(37.87);
    uv_SW_bot.push_back(40.68);
    uv_SW_bot.push_back(66.35);
    uv_SW_bot.push_back(67.21);

    uv_SW_slopes = CalcSlopes(uv_SW_top, uv_SW_bot);
    uv_SW_intercepts = CalcIntercept(uv_SW_top, uv_SW_slopes);

  }//end inittrack...()

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
    std::vector<float> temp_vector_slopes;
    std::vector<float> temp_vector_intercepts;
    
    if (zed > 0 && arm == WEST)
      {
        // std::cout << "TrackQualiy.h inBrokenUV North West" << std::endl;
        temp_vector_slopes = uv_NW_slopes;
        temp_vector_intercepts = uv_NW_intercepts;
      }
    if (zed > 0 && arm == EAST)
      {
        // std::cout << "TrackQualiy.h inBrokenUV North East" << std::endl;
        temp_vector_slopes = uv_NE_slopes;
        temp_vector_intercepts = uv_NE_intercepts;
      }
    if (zed < 0 && arm == WEST)
      {
        // std::cout << "TrackQualiy.h inBrokenUV South West" << std::endl;
        temp_vector_slopes = uv_SW_slopes;
        temp_vector_intercepts = uv_SW_intercepts;
      }
    if (zed < 0 && arm == EAST)
      {
        // std::cout << "TrackQualiy.h inBrokenUV South East" << std::endl;
        temp_vector_slopes = uv_SE_slopes;
        temp_vector_intercepts = uv_SE_intercepts;
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
  }//inBrokenUV

  inline bool passQualityMask(int quality, float phi, float alpha, float zed, int arm)
  {
    //(quality & 1) == 0 ==> no X1 bit
    //(quality & 1) == 1 ==> X1 bit

    if((quality & 1) == 0 && (quality & 2) == 0)
      {//no x1 or x2 bit used
        return false;
      }

    if ((quality & 1) == 0 && !inBrokenX1(phi, alpha, zed, arm))
      {//no x1 bit and not in brokenx1 region
        return false;
      }

    if ((quality & 2) == 0 && !inBrokenX2(phi, alpha, zed, arm))
      {//no x2 bit and not in broken x2 region
        return false;
      }

    if ((quality & 16) == 0)
      {//no PC1 bit
        return false;
      }

    if ((quality & 32) == 0 && (quality & 12) == 0)
      {//no PC unique bit && no UV unique bit
        // PC not unique- so UV needs to be unique
        return false;
      }

    if ((quality & 12) == 0 && !inBrokenUV(phi, alpha, zed, arm))
      {//no UV unique bit && not in region of weak UV acceptance
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
