#ifndef CENTBIN
#define CENTBIN

#include <string>
#include "TMath.h"

class CentralityBin{
public:
  
  CentralityBin();
  ~CentralityBin();

  int getHiBinFromhiHF( const Double_t hiHF, int variation = 0);
  int getHiBinFromhiHFSides( const Double_t hiHFP, const Double_t hiHFM, int variation = 0);


private:

  std::string mode;
  static const Int_t nBins = 200;

  const Double_t binTableData[nBins+1] = {0, 10.5072, 11.2099, 11.8364, 12.478, 13.1194, 13.7623, 14.4081, 15.0709, 15.7532, 16.4673, 17.1881, 17.923, 18.673, 19.4865, 20.3033, 21.1536, 22.0086, 22.9046, 23.8196, 24.7924, 25.8082, 26.8714, 27.9481, 29.0828, 30.2757, 31.5043, 32.8044, 34.1572, 35.6142, 37.1211, 38.6798, 40.3116, 42.0398, 43.8572, 45.6977, 47.6312, 49.6899, 51.815, 54.028, 56.3037, 58.7091, 61.2024, 63.8353, 66.5926, 69.3617, 72.2068, 75.2459, 78.3873, 81.5916, 84.9419, 88.498, 92.1789, 95.9582, 99.8431, 103.739, 107.78, 111.97, 116.312, 120.806, 125.46, 130.269, 135.247, 140.389, 145.713, 151.212, 156.871, 162.729, 168.762, 174.998, 181.424, 188.063, 194.907, 201.942, 209.19, 216.683, 224.37, 232.291, 240.43, 248.807, 257.416, 266.256, 275.348, 284.668, 294.216, 304.053, 314.142, 324.488, 335.101, 345.974, 357.116, 368.547, 380.283, 392.29, 404.564, 417.122, 429.968, 443.116, 456.577, 470.357, 484.422, 498.78, 513.473, 528.479, 543.813, 559.445, 575.411, 591.724, 608.352, 625.344, 642.686, 660.361, 678.371, 696.749, 715.485, 734.608, 754.068, 773.846, 794.046, 814.649, 835.608, 856.972, 878.719, 900.887, 923.409, 946.374, 969.674, 993.435, 1017.62, 1042.21, 1067.28, 1092.72, 1118.64, 1144.96, 1171.71, 1198.98, 1226.67, 1254.82, 1283.46, 1312.65, 1342.21, 1372.27, 1402.85, 1433.93, 1465.49, 1497.62, 1530.29, 1563.49, 1597.22, 1631.49, 1666.37, 1701.8, 1737.75, 1774.35, 1811.51, 1849.29, 1887.75, 1926.79, 1966.6, 2006.97, 2047.99, 2089.71, 2132.1, 2175.23, 2219.17, 2263.72, 2309.2, 2355.43, 2402.47, 2450.33, 2499.05, 2548.66, 2599.16, 2650.59, 2703.03, 2756.32, 2810.75, 2866.27, 2922.91, 2980.54, 3039.47, 3099.53, 3160.98, 3223.66, 3287.71, 3353.18, 3420.34, 3489.13, 3559.72, 3632.06, 3706.18, 3782.42, 3860.78, 3941.42, 4024.52, 4110.27, 4199.4, 4292.8, 4394.49, 4519.52, 5199.95};
  const Double_t binTableDataDown[nBins+1] = {0, 10.5071, 11.2094, 11.8357, 12.4763, 13.117, 13.7597, 14.4049, 15.0671, 15.7491, 16.4622, 17.1812, 17.9144, 18.6674, 19.4797, 20.2963, 21.1435, 21.9974, 22.8928, 23.8068, 24.7805, 25.7931, 26.8556, 27.9308, 29.0638, 30.2582, 31.4795, 32.7816, 34.1349, 35.5834, 37.0941, 38.6474, 40.2782, 42.0035, 43.8112, 45.6576, 47.5758, 49.6381, 51.6667, 53.7353, 55.8903, 58.1259, 60.4528, 62.8712, 65.3859, 67.9968, 70.7065, 73.5231, 76.4519, 79.4922, 82.6461, 85.9264, 89.3269, 92.8562, 96.5212, 100.322, 104.262, 108.344, 112.585, 116.971, 121.521, 126.225, 131.09, 136.127, 141.328, 146.721, 152.284, 158.014, 163.935, 170.054, 176.372, 182.878, 189.602, 196.532, 203.653, 211.017, 218.599, 226.387, 234.418, 242.667, 251.16, 259.886, 268.852, 278.071, 287.498, 297.2, 307.184, 317.409, 327.894, 338.66, 349.686, 360.996, 372.607, 384.508, 396.669, 409.133, 421.86, 434.906, 448.258, 461.916, 475.906, 490.16, 504.74, 519.663, 534.911, 550.453, 566.322, 582.525, 599.08, 615.968, 633.211, 650.805, 668.76, 687.048, 705.707, 724.774, 744.163, 763.9, 783.999, 804.528, 825.432, 846.746, 868.429, 890.523, 913.007, 935.952, 959.211, 982.919, 1007.08, 1031.63, 1056.62, 1082.08, 1107.96, 1134.24, 1160.99, 1188.22, 1215.91, 1244.06, 1272.69, 1301.85, 1331.45, 1361.51, 1392.07, 1423.18, 1454.77, 1486.93, 1519.57, 1552.81, 1586.55, 1620.87, 1655.79, 1691.26, 1727.27, 1763.93, 1801.12, 1838.97, 1877.47, 1916.61, 1956.45, 1996.89, 2038.04, 2079.84, 2122.35, 2165.52, 2209.53, 2254.24, 2299.83, 2346.19, 2393.31, 2441.28, 2490.16, 2539.86, 2590.57, 2642.16, 2694.74, 2748.23, 2802.81, 2858.47, 2915.33, 2973.2, 3032.28, 3092.56, 3154.24, 3217.19, 3281.45, 3347.18, 3414.6, 3483.65, 3554.56, 3627.2, 3701.66, 3778.25, 3856.97, 3937.98, 4021.48, 4107.62, 4197.21, 4291.05, 4393.19, 4518.6, 5199.95};
  const Double_t binTableDataUp[nBins+1] = {0, 10.5075, 11.2107, 11.838, 12.4797, 13.1213, 13.7641, 14.4124, 15.0745, 15.7577, 16.473, 17.1939, 17.9297, 18.6812, 19.4958, 20.3143, 21.1648, 22.0218, 22.9159, 23.8328, 24.8059, 25.8204, 26.89, 27.9702, 29.1042, 30.3022, 31.528, 32.8347, 34.1896, 35.6439, 37.1542, 38.7172, 40.3518, 42.091, 43.9053, 45.7415, 47.6853, 49.7457, 51.8755, 54.0983, 56.3594, 58.7848, 61.2861, 63.9228, 66.6825, 69.4421, 72.297, 75.3547, 78.4967, 81.6977, 85.0755, 88.6211, 92.3058, 96.1071, 99.9975, 104.065, 108.272, 112.512, 116.906, 121.601, 126.465, 131.482, 136.866, 142.229, 147.786, 153.546, 159.571, 165.586, 171.902, 178.419, 185.063, 191.856, 199.055, 206.261, 213.999, 221.719, 229.671, 237.84, 246.088, 254.828, 263.883, 272.907, 282.236, 291.925, 301.519, 311.477, 321.691, 332.153, 342.892, 353.878, 365.161, 376.742, 388.577, 400.684, 413.075, 425.746, 438.711, 451.989, 465.556, 479.45, 493.608, 508.077, 522.891, 538.003, 553.415, 569.151, 585.216, 601.601, 618.354, 635.422, 652.84, 670.599, 688.699, 707.161, 726.014, 745.185, 764.687, 784.557, 804.838, 825.489, 846.537, 867.951, 889.752, 911.955, 934.588, 957.52, 980.912, 1004.73, 1028.94, 1053.57, 1078.67, 1104.17, 1130.07, 1156.39, 1183.2, 1210.47, 1238.17, 1266.38, 1295.02, 1324.16, 1353.71, 1383.77, 1414.35, 1445.41, 1477, 1509.09, 1541.74, 1574.88, 1608.59, 1642.83, 1677.66, 1713.07, 1748.98, 1785.47, 1822.63, 1860.33, 1898.72, 1937.73, 1977.42, 2017.71, 2058.62, 2100.25, 2142.57, 2185.56, 2229.38, 2273.91, 2319.2, 2365.33, 2412.22, 2459.94, 2508.52, 2557.98, 2608.35, 2659.61, 2711.86, 2765, 2819.23, 2874.58, 2930.97, 2988.46, 3047.12, 3106.95, 3168.15, 3230.6, 3294.37, 3359.58, 3426.47, 3494.95, 3565.21, 3637.21, 3711.03, 3786.91, 3864.85, 3945.11, 4027.8, 4113.06, 4201.73, 4294.72, 4395.9, 4520.5, 5199.95}; 
  const Double_t binTableMC_Drum5F[nBins+1] = {0,12.2187,13.0371,13.7674,14.5129,15.2603,16.0086,16.7623,17.5335,18.3283,19.1596,19.9989,20.8532,21.7297,22.6773,23.6313,24.6208,25.6155,26.6585,27.7223,28.8632,30.041,31.2865,32.5431,33.8655,35.2539,36.6912,38.2064,39.7876,41.4818,43.2416,45.0605,46.9652,48.9918,51.1,53.2417,55.5094,57.9209,60.3817,62.9778,65.6099,68.4352,71.3543,74.4154,77.6252,80.8425,84.1611,87.7395,91.3973,95.1286,99.0571,103.185,107.482,111.929,116.45,121.178,126.081,130.995,136.171,141.612,147.298,153.139,159.419,165.633,172.114,178.881,185.844,192.845,200.244,207.83,215.529,223.489,231.878,240.254,249.319,258.303,267.508,277.037,286.729,296.845,307.458,317.882,328.787,340.074,351.295,362.979,375.125,387.197,399.604,412.516,425.683,439.001,452.667,466.816,481.007,495.679,510.588,526.138,541.782,557.641,574.141,591.071,608.379,626.068,643.616,661.885,680.288,699.449,718.925,738.968,758.983,779.459,800.376,821.638,843.555,865.771,888.339,911.031,934.979,958.56,982.582,1007.02,1031.9,1057.81,1084.01,1111.71,1138.21,1165.72,1193.73,1221.65,1251.51,1281.23,1311.01,1341.1,1372.4,1404.29,1436.52,1468.65,1501.91,1535.56,1569.69,1604.69,1640.65,1676.05,1712.62,1749.28,1787.43,1825.89,1866.07,1906.58,1947.84,1989.66,2031.4,2072.8,2115.32,2159.5,2205.23,2252.68,2298.58,2345.65,2393.36,2442.87,2491.45,2541.04,2592.81,2645.52,2699.1,2753.29,2807.93,2864.37,2922.6,2979.42,3038.68,3098.72,3159.29,3221.66,3285.9,3350.95,3415.81,3482.69,3552.62,3623.61,3694.63,3767.25,3840.28,3917.04,3993.66,4073.36,4154.33,4238.13,4322.21,4409.83,4498.89,4589.72,4681.56,4777.09,4877.95,4987.05,5113.04,5279.58, 6242.82};

};

CentralityBin::CentralityBin(){
}

CentralityBin::~CentralityBin(){

}

int CentralityBin::getHiBinFromhiHFSides( const Double_t hiHFP, const Double_t hiHFM, int variation){
  return getHiBinFromhiHF( hiHFP + hiHFM, variation);
}

int CentralityBin::getHiBinFromhiHF(const Double_t hiHF, int variation)
{
  if(variation<=2 && hiHF >= 5199.95) return 0;
  if(variation == 3 && hiHF >= 6242.82) return 0;

  Int_t binPos = -1;

  if( variation == 0){
    for(int i = 0; i < nBins; ++i){
      if(hiHF >= binTableData[i] && hiHF < binTableData[i+1]){
        binPos = i;
        break;
      }
    }
  }
  if( variation == 1 ){
    for(int i = 0; i < nBins; ++i){
      if(hiHF >= binTableDataUp[i] && hiHF < binTableDataUp[i+1]){
        binPos = i;
        break;
      }
    }
  }
  if( variation == 2 ){
    for(int i = 0; i < nBins; ++i){
      if(hiHF >= binTableDataDown[i] && hiHF < binTableDataDown[i+1]){
        binPos = i;
        break;
      }
    }
  }
  if( variation == 3 ){
    for(int i = 0; i < nBins; ++i){
      if(hiHF >= binTableMC_Drum5F[i] && hiHF < binTableMC_Drum5F[i+1]){
        binPos = i;
        break;
      }
    }
  }

  binPos = nBins - 1 - binPos;

  return (Int_t)(200*((Double_t)binPos)/((Double_t)nBins));
}


#endif
