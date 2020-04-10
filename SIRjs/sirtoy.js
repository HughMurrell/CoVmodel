// -----------------------------------------------------------
// get the clean John Hopkins data
console.log("started")

var csv = `country,population,lat,long,first,last,intervention,beta_before,beta_after,gamma_before,gamma_after,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15,D16,D17,D18,D19,D20,D21,D22,D23,D24,D25,D26,D27,D28,D29,D30,D31,D32,D33,D34,D35,D36,D37,D38,D39,D40,D41,D42,D43,D44,D45,D46,D47,D48,D49,D50,D51,D52,D53,D54,D55,D56,D57,D58,D59,D60,D61,D62,D63,D64,D65,D66,D67,D68,D69,D70,D71,D72,D73,D74,D75,D76,D77,D78
Algeria,42228429,28.0339,1.6596,0020-01-22,0020-04-08,48,0.35,0.18,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,3,5,12,12,17,17,19,20,20,20,24,26,37,48,54,60,74,87,90,139,201,230,264,302,367,409,454,511,584,716,847,986,1171,1251,1320,1423,1468,1572
Argentina,44494502,-38.4161,-63.6167,0020-01-22,0020-04-08,53,0.54,0.16,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,2,8,12,12,17,19,19,31,34,45,56,68,79,97,128,158,266,301,387,387,502,589,690,745,820,1054,1054,1133,1265,1451,1451,1554,1628,1715
Austria,8847037,47.5162,14.5501,0020-01-22,0020-04-08,64,0.4,0.06,0.06,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,3,3,9,14,18,21,29,41,55,79,104,131,182,246,302,504,655,860,1018,1332,1646,2013,2388,2814,3582,4474,5283,5588,6909,7657,8271,8788,9618,10180,10711,11129,11524,11781,12051,12297,12639,12942
Belarus,9485386,53.7098,27.9534,0020-01-22,0020-04-08,53,0.17,0.28,0.06,0.04,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,6,6,6,6,6,6,9,9,12,27,27,27,36,36,51,51,69,76,76,81,81,86,86,94,94,94,152,152,163,304,351,440,562,700,861,1066
Belgium,11422068,50.8333,4.0,0020-01-22,0020-04-08,57,0.24,0.16,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,8,13,23,50,109,169,200,239,267,314,314,559,689,886,1058,1243,1486,1795,2257,2815,3401,3743,4269,4937,6235,7284,9134,10836,11899,12775,13964,15348,16770,18431,19691,20814,22194,23403
Brazil,209469333,-14.235,-51.9253,0020-01-22,0020-04-08,62,0.36,0.18,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,2,2,2,2,4,4,13,13,20,25,31,38,52,151,151,162,200,321,372,621,793,1021,1546,1924,2247,2554,2985,3417,3904,4256,4579,5717,6836,8044,9056,10360,11130,12161,14034,16170
Chile,18729160,-35.6751,-71.543,0020-01-22,0020-04-08,60,0.46,0.15,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,4,4,4,8,8,13,23,23,43,61,74,155,201,238,238,434,537,632,746,922,1142,1306,1610,1909,2139,2449,2738,3031,3404,3737,4161,4471,4815,5116,5546
Colombia,49648685,4.5709,-74.2973,0020-01-22,0020-04-08,55,0.61,0.17,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,3,9,9,13,22,34,54,65,93,102,128,196,231,277,378,470,491,539,608,702,798,906,1065,1161,1267,1406,1485,1579,1780,2054
Croatia,4089400,45.1,15.2,0020-01-22,0020-04-08,64,0.27,0.1,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,3,3,5,6,7,7,9,10,10,11,12,12,12,14,19,19,32,38,49,57,65,81,105,128,206,254,315,382,442,495,586,657,713,790,867,963,1011,1079,1126,1182,1222,1282,1343
Denmark,5797446,56.2639,9.5018,0020-01-22,0020-04-08,52,0.55,0.11,0.06,0.04,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,3,4,4,6,10,10,23,23,35,90,262,442,615,801,827,864,914,977,1057,1151,1255,1326,1395,1450,1591,1724,1877,2046,2201,2395,2577,2860,3107,3386,3757,4077,4369,4681,5071,5402
Dominican Republic,10627165,18.7357,-70.1627,0020-01-22,0020-04-08,65,0.32,0.13,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,2,2,5,5,5,5,5,5,11,11,11,21,21,34,72,112,202,245,312,392,488,581,719,859,901,1109,1284,1380,1488,1488,1745,1828,1956,2111
Ecuador,17084357,-1.8312,-78.1834,0020-01-22,0020-04-08,65,0.38,0.11,0.06,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,6,7,10,13,13,13,14,15,15,17,17,17,28,28,37,58,111,199,367,506,789,981,1082,1173,1403,1595,1823,1924,1962,2240,2748,3163,3368,3465,3646,3747,3747,4450
Estonia,1320884,58.5953,25.0136,0020-01-22,0020-04-08,56,0.35,0.12,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,2,2,3,10,10,10,10,12,16,16,79,115,171,205,225,258,267,283,306,326,352,369,404,538,575,645,679,715,745,779,858,961,1039,1097,1108,1149,1185
Finland,5518050,64.0,26.0,0020-01-22,0020-04-08,55,0.17,0.14,0.04,0.06,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,3,6,6,6,6,12,15,15,23,30,40,59,59,155,225,244,277,321,336,400,450,523,626,700,792,880,958,1041,1167,1240,1352,1418,1446,1518,1615,1882,1927,2176,2308,2487
France,66987244,46.2276,2.2137,0020-01-22,0020-04-08,53,0.23,0.17,0.04,0.06,0,0,2,3,3,3,4,5,5,5,6,6,6,6,6,6,6,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,14,18,38,57,100,130,191,204,285,377,653,949,1126,1209,1784,2281,2281,3661,4469,4499,6633,7652,9043,10871,12612,14282,16018,19856,22304,25233,29155,32964,37575,40174,44550,52128,56989,59105,64338,89953,92839,98010,109069,112950
Germany,82927922,51.0,9.0,0020-01-22,0020-04-08,54,0.25,0.14,0.04,0.06,0,0,0,0,0,1,4,4,4,5,8,10,12,12,12,12,13,13,14,14,16,16,16,16,16,16,16,16,16,16,16,16,16,16,17,27,46,48,79,130,159,196,262,482,670,799,1040,1176,1457,1908,2078,3675,4585,5795,7272,9257,12327,15320,19848,22213,24873,29056,32986,37323,43938,50871,57695,62095,66885,71808,77872,84794,91159,96092,100123,103374,107663,113296
Greece,10727668,39.0742,21.8243,0020-01-22,0020-04-08,49,0.51,0.12,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,3,4,4,7,7,7,9,31,45,46,73,73,89,99,99,190,228,331,331,387,418,418,495,530,624,695,743,821,892,966,1061,1156,1212,1314,1415,1544,1613,1673,1735,1755,1832,1884
Iceland,353574,64.9631,-19.0208,0020-01-22,0020-04-08,50,0.54,0.12,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,3,6,11,26,34,43,50,50,58,69,85,103,134,156,171,180,220,250,330,409,473,568,588,648,737,802,890,963,1020,1086,1135,1220,1319,1364,1417,1486,1562,1586,1616
India,1352617328,21.0,78.0,0020-01-22,0020-04-08,55,0.15,0.24,0.04,0.06,0,0,0,0,0,0,0,0,1,1,1,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5,5,28,30,31,34,39,43,56,62,73,82,102,113,119,142,156,194,244,330,396,499,536,657,727,887,987,1024,1251,1397,1998,2543,2567,3082,3588,4778,5311,5916
Indonesia,267663435,-0.7893,113.9213,0020-01-22,0020-04-08,54,0.53,0.15,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,4,4,6,19,27,34,34,69,96,117,134,172,227,311,369,450,514,579,686,790,893,1046,1155,1285,1414,1528,1677,1790,1986,2092,2273,2491,2738,2956
Iraq,38433600,33.0,44.0,0020-01-22,0020-04-08,48,0.37,0.15,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,5,7,7,13,19,26,32,35,35,40,54,60,60,71,71,71,101,110,116,124,154,164,192,208,214,233,266,316,346,382,458,506,547,630,694,728,772,820,878,961,1031,1122,1202
Ireland,4853506,53.1424,-7.6921,0020-01-22,0020-04-08,61,0.4,0.14,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,2,6,6,18,18,19,21,34,43,43,90,129,129,169,223,292,557,683,785,906,1125,1329,1564,1819,2121,2415,2615,2910,3235,3447,3849,4273,4604,4994,5364,5709,6074
Israel,8883800,31.0,35.0,0020-01-22,0020-04-08,62,0.31,0.14,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,2,3,4,7,10,10,12,15,20,37,43,61,61,75,79,100,126,155,213,218,250,304,427,529,712,883,1071,1238,2369,2693,3035,3619,4247,4695,5358,6092,6857,7428,7851,8430,8904,9248,9404
Italy,60431283,43.0,12.0,0020-01-22,0020-04-08,55,0.3,0.09,0.04,0.06,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,20,62,155,229,322,453,655,888,1128,1694,2036,2502,3089,3858,4636,5883,7375,9172,10149,12462,12462,17660,21157,24747,27980,31506,35713,41035,47021,53578,59138,63927,69176,74386,80589,86498,92472,97689,101739,105792,110574,115242,119827,124632,128948,132547,135586,139422
Japan,126529100,36.0,138.0,0020-01-22,0020-04-08,36,0.2,0.12,0.06,0.04,2,2,2,2,4,4,7,7,11,15,20,20,20,22,22,22,25,25,26,26,26,28,28,29,43,59,66,74,84,94,105,122,147,159,170,189,214,228,241,256,274,293,331,360,420,461,502,511,581,639,639,701,773,839,839,878,889,924,963,1007,1101,1128,1193,1307,1387,1468,1693,1866,1866,1953,2178,2495,2617,3139,3139,3654,3906,4257
Luxembourg,607728,49.8153,6.1296,0020-01-22,0020-04-08,64,0.38,0.07,0.06,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,3,3,5,7,19,34,51,59,77,140,203,335,484,670,798,875,1099,1333,1453,1605,1831,1950,1988,2178,2319,2487,2612,2729,2804,2843,2970,3034
Malaysia,31528585,2.5,112.5,0020-01-22,0020-04-08,53,0.17,0.12,0.04,0.06,0,0,0,3,4,4,4,7,8,8,8,8,8,10,12,12,12,16,16,18,18,18,19,19,22,22,22,22,22,22,22,22,22,22,22,22,23,23,25,29,29,36,50,50,83,93,99,117,129,149,149,197,238,428,566,673,790,900,1030,1183,1306,1518,1624,1796,2031,2161,2320,2470,2626,2766,2908,3116,3333,3483,3662,3793,3963,4119
Mexico,126190788,23.6345,-102.5528,0020-01-22,0020-04-08,65,0.31,0.16,0.06,0.04,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,4,5,5,5,5,5,6,6,7,7,7,8,12,12,26,41,53,82,93,118,164,203,251,316,367,405,475,585,717,848,993,1094,1215,1378,1510,1688,1890,2143,2439,2785
Moldova,3545883,47.4116,28.3699,0020-01-22,0020-04-08,56,0.5,0.22,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,3,3,3,6,12,23,23,30,30,49,66,80,94,109,125,149,177,199,231,263,298,353,423,505,591,752,864,965,1056,1174
Morocco,36029138,31.7917,-7.0926,0020-01-22,0020-04-08,66,0.3,0.14,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,2,2,2,2,2,3,5,6,7,17,28,29,38,49,63,77,96,115,143,170,225,275,345,402,479,556,617,654,708,791,919,1021,1120,1184,1275
Netherlands,17231017,52.1326,5.2913,0020-01-22,0020-04-08,50,0.73,0.14,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,6,10,18,24,38,82,128,188,265,321,382,503,503,804,959,1135,1413,1705,2051,2460,2994,3631,4204,4749,5560,6412,7431,8603,9762,10866,11750,12595,13614,14697,15723,16627,17851,18803,19580,20549
New Zealand,4885500,-40.9006,174.886,0020-01-22,0020-04-08,65,0.27,0.14,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,3,3,4,5,5,5,5,5,5,5,6,8,8,12,20,28,39,52,102,102,155,205,283,368,451,514,589,647,708,797,868,950,1039,1106,1160,1210
Norway,5314336,60.472,8.4689,0020-01-22,0020-04-08,49,0.67,0.11,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,6,15,19,25,32,56,87,108,147,176,205,400,598,702,996,1090,1221,1333,1463,1550,1746,1914,2118,2385,2621,2863,3084,3369,3755,4015,4284,4445,4641,4863,5147,5370,5550,5687,5865,6086,6086
Pakistan,212215030,30.3753,69.3451,0020-01-22,0020-04-08,64,0.33,0.13,0.06,0.04,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,4,4,4,5,5,5,6,6,6,6,16,19,20,28,31,53,136,236,299,454,501,730,776,875,972,1063,1201,1373,1495,1597,1717,1938,2118,2421,2686,2818,3157,3766,4035,4263
Panama,4176873,8.538,-80.7821,0020-01-22,0020-04-08,58,0.78,0.16,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,8,11,27,36,43,55,69,86,109,137,200,313,345,345,443,558,674,786,901,989,1181,1181,1317,1475,1673,1801,1988,2100,2249
Peru,31989256,-9.19,-75.0152,0020-01-22,0020-04-08,55,0.57,0.21,0.06,0.04,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,6,7,11,11,15,28,38,43,86,117,145,234,234,318,363,395,416,480,580,635,671,852,950,1065,1323,1414,1595,1746,2281,2561,2954,4342
Philippines,106651922,13.0,122.0,0020-01-22,0020-04-08,55,0.16,0.18,0.04,0.06,0,0,0,0,0,0,0,0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5,6,10,20,33,49,52,64,111,140,142,187,202,217,230,307,380,462,552,636,707,803,1075,1418,1546,2084,2311,2633,3018,3094,3246,3660,3764,3870
Poland,37978548,51.9194,19.1451,0020-01-22,0020-04-08,54,0.65,0.17,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,5,5,11,16,22,31,49,68,103,119,177,238,251,355,425,536,634,749,901,1051,1221,1389,1638,1862,2055,2311,2554,2946,3383,3627,4102,4413,4848,5205
Portugal,10281762,39.3999,-8.2245,0020-01-22,0020-04-08,66,0.45,0.11,0.06,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,5,8,13,20,30,30,41,59,59,112,169,245,331,448,448,785,1020,1280,1600,2060,2362,2995,3544,4268,5170,5962,6408,7443,8251,9034,9886,10524,11278,11730,12442,13141
Qatar,2781677,25.3548,51.1839,0020-01-22,0020-04-08,53,0.46,0.14,0.06,0.04,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,3,3,7,8,8,8,8,15,18,24,262,262,320,337,401,439,439,452,460,470,481,494,501,526,537,549,562,590,634,693,781,835,949,1075,1325,1604,1832,2057,2210
Romania,19473936,45.9432,24.9668,0020-01-22,0020-04-08,64,0.31,0.15,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,3,3,3,3,3,4,6,9,9,15,15,25,45,49,89,123,131,158,184,260,277,308,367,433,576,794,906,1029,1292,1452,1815,2109,2245,2460,2738,3183,3613,3864,4057,4417,4761
Saudi Arabia,33699947,24.0,45.0,0020-01-22,0020-04-08,53,0.57,0.15,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,5,5,5,11,15,20,21,45,86,103,103,118,171,171,274,344,392,511,562,767,900,1012,1104,1203,1299,1453,1563,1720,1885,2039,2179,2402,2605,2795,2932
Serbia,6982084,44.0165,21.0059,0020-01-22,0020-04-08,55,0.56,0.21,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,5,12,19,35,46,48,55,65,83,103,135,171,222,249,303,384,384,457,659,741,785,900,1060,1171,1476,1624,1908,2200,2447,2666
Singapore,5638676,1.2833,103.8333,0020-01-22,0020-04-08,26,0.2,0.12,0.06,0.04,0,1,3,3,4,5,7,7,10,13,16,18,18,24,28,28,30,33,40,45,47,50,58,67,72,75,77,81,84,84,85,85,89,89,91,93,93,93,102,106,108,110,110,117,130,138,150,150,160,178,178,200,212,226,243,266,313,345,385,432,455,509,558,631,683,732,802,844,879,926,1000,1049,1114,1189,1309,1375,1481,1623
Slovenia,2067372,46.1512,14.9955,0020-01-22,0020-04-08,54,0.7,0.1,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,7,7,16,16,31,57,89,141,181,219,253,275,275,286,341,383,414,442,480,528,562,632,684,730,756,802,841,897,934,977,997,1021,1059,1091
South Africa,57779622,-30.5595,22.9375,0020-01-22,0020-04-08,65,0.43,0.05,0.06,0.04,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,3,3,7,13,17,24,38,51,62,62,116,150,202,240,274,402,554,709,927,1170,1187,1280,1326,1353,1380,1462,1505,1585,1655,1686,1749,1845
Spain,46723749,40.0,-4.0,0020-01-22,0020-04-08,56,0.28,0.14,0.04,0.06,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,6,13,15,32,45,84,120,165,222,259,400,500,673,1073,1695,2277,2277,5232,6391,7798,9942,11748,13910,17963,20410,25374,28768,35136,39885,49515,57786,65719,73235,80110,87956,95923,104118,112065,119199,126168,131646,136675,141942,148220
Sweden,10183175,63.0,16.0,0020-01-22,0020-04-08,55,0.2,0.14,0.04,0.06,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,7,7,12,14,15,21,35,94,101,161,203,248,355,500,599,814,961,1022,1103,1190,1279,1439,1639,1763,1934,2046,2286,2526,2840,3069,3447,3700,4028,4435,4947,5568,6131,6443,6830,7206,7693,8419
Switzerland,8516543,46.8182,8.2275,0020-01-22,0020-04-08,48,0.78,0.13,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,8,8,18,27,42,56,90,114,214,268,337,374,491,652,652,1139,1359,2200,2200,2700,3028,4075,5294,6575,7474,8795,9877,10897,11811,12928,14076,14829,15922,16605,17768,18827,19606,20505,21100,21657,22253,23280
Thailand,69428524,15.0,101.0,0020-01-22,0020-04-08,26,0.16,0.15,0.04,0.06,2,3,5,7,8,8,14,14,14,19,19,19,19,25,25,25,25,32,32,32,33,33,33,33,33,34,35,35,35,35,35,35,35,35,37,40,40,41,42,42,43,43,43,47,48,50,50,50,53,59,70,75,82,114,147,177,212,272,322,411,599,721,827,934,1045,1136,1245,1388,1524,1651,1771,1875,1978,2067,2169,2220,2258,2369
Ukraine,44622516,48.3794,31.1656,0020-01-22,0020-04-08,66,0.31,0.17,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,3,3,3,7,14,14,16,29,47,73,73,97,145,196,310,356,475,548,645,794,897,1072,1225,1308,1319,1462,1668
United Arab Emirates,9630959,24.0,54.0,0020-01-22,0020-04-08,30,0.07,0.21,0.06,0.04,0,0,0,0,0,0,0,4,4,4,4,5,5,5,5,5,5,7,7,8,8,8,8,8,8,9,9,9,9,9,9,13,13,13,13,13,13,19,21,21,21,27,27,29,29,45,45,45,74,74,85,85,85,98,98,98,113,140,140,153,153,198,248,333,333,405,468,570,611,664,814,1024,1264,1505,1799,2076,2359,2659
United Kingdom,66488991,55.3781,-3.4360000000000004,0020-01-22,0020-04-08,55,0.23,0.19,0.04,0.06,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,3,3,3,8,8,9,9,9,9,9,9,9,9,9,9,9,9,13,13,13,15,20,23,36,40,51,85,115,163,206,273,321,382,456,456,798,1140,1140,1543,1950,2626,2689,3983,5018,5683,6650,8077,9529,11658,14543,17089,19522,22141,25150,29474,33718,38168,41903,47806,51608,55242,60733
`;

var jh_data = d3.csvParse(csv);

var latest_data = []

var countries = [];
for (var i = 0; i < jh_data.length; i++){
    countries.push(jh_data[i]["country"]);
}
countries.sort();

var country;

function populateCountries(countryElementId){
    // given the id of the <select> tag as function argument, it inserts <option> tags
    var countryElement = document.getElementById(countryElementId);
    countryElement.length=0;
    // countryElement.options[0] = new Option('Select Country','-1');
    countryElement.selectedIndex = 0;
    for (var i=0; i<countries.length; i++) {
        countryElement.options[countryElement.length] = new Option(countries[i],countries[i]);
    }
    // Assigned all countries. Now assign event listener.

    countryElement.onchange = function(){
        country = countries[countryElement.selectedIndex];
        for (var i=0; i<jh_data.length; i++){
            if (jh_data[i]["country"] == country){
                values = Object.values(jh_data[i]);
                // the 9 below is to drop non-data fields
                int_values = values.slice(11,values.length).map(Number);
                index=int_values.findIndex(function(number) {
                  return number > 0;
                });
                latest_data = int_values.slice(index,int_values.length);
                break_point = jh_data[i]["intervention"] - index;
                p_SI = jh_data[i]["beta_before"];
                p_SI_ld = jh_data[i]["beta_after"];
                p_IR = jh_data[i]["gamma_before"];
                p_IR_ld = jh_data[i]["gamma_after"];
                N = jh_data[i]["population"];
                $("#p_N").prop('disabled', true);
                $('#p_N').val(N);
            }
        }
        reset_all();
        console.log(country)
    }
}

populateCountries("country2");


// latest_data = [1, 1, 1, 3, 3, 7, 13, 17, 24, 38, 51, 62, 62, 116, 150, 202, 240, 274, 402, 554, 709, 927, 1170, 1187, 1280, 1326, 1353, 1380, 1462, 1505, 1585]

var valid_data = 1;
var running = 0;
var p_SI_default = 0.46;
var p_IR_default = 0.1;
var p_SI_ld_default = 0.08;
var p_IR_ld_default = 0.1;

var p_SI = p_SI_default;
var p_IR = p_IR_default;
var p_SI_ld = p_SI_ld_default;
var p_IR_ld = p_IR_ld_default;

var time_interval = 200;
var count = 0;
var N = 60000000;
var timeseries;


var sir_color = {D: "#000000", S: "#00ffff", I: "#f00000", R: "#00f000", C: "#0000f0" }

var epi_state = { S: (N-1)/N, I: 1/N, R: 0 };

function reset_params () {
    
    // p_SI = p_SI_default;
    // p_IR = p_IR_default;
    // p_SI_ld = p_SI_ld_default;
    // p_IR_ld = p_IR_ld_default;
    
    $("#p_SI").val(p_SI)
    $("#p_SI").keyup(update_p_SI);

    $("#p_IR").val(p_IR)
    $("#p_IR").keyup(update_p_IR);

    $("#p_SI_ld").val(p_SI_ld)
    $("#p_SI_ld").keyup(update_p_SI_ld);

    $("#p_IR_ld").val(p_IR_ld)
    $("#p_IR_ld").keyup(update_p_IR_ld);
}

function reset_history () {
    
  // console.log(latest_data[latest_data.length-1]);

  timeseries = {D: {label: "Data", color: sir_color.D, data: []},
                S: {label: "Susceptible", color: sir_color.S, data: []},
                I: {label: "Infective", color: sir_color.I, data: []},
                R: {label: "Recovered (or dead)", color: sir_color.R, data: []},
                C: {label: "Cumulative (Infective + Recovered)", color: sir_color.C, data: []},
                B: {label: "Start of Intervention", color: sir_color.B, data: [] }
  };
    
  for (i = 0; i < latest_data.length; i++) {
        timeseries.D.data.push([i, latest_data[i]]);
  }
    
  timeseries.B.data.push([break_point, 0]);
  timeseries.B.data.push([break_point, latest_data[latest_data.length-1]]);

  timeseries.S.data.push([count, N*epi_state.S]);
  timeseries.I.data.push([count, N*epi_state.I]);
  timeseries.R.data.push([count, N*epi_state.R]);
  timeseries.C.data.push([count, N*(epi_state.R+epi_state.I)]);
}

reset_params();
reset_history();

var plotOptions = {
        lines: { show: true },
	    points: { show: true },
        xaxis: {min: 0},
        series: { shadowSize: 0 },
        grid: {hoverable: true, clickable:true},
         legend:{
                   // backgroundOpacity: 0.5,
                   // noColumns: 0,
                   // backgroundColor: "green",
                   position: "nw"
               }
    };

var placeholder = $("#epicurves");
var plot = $.plot($("#epicurves"), [], plotOptions);

$("#epicurves").on("plotclick",function(event,pos,item){
    if(item){
        break_point = item.series.data[item.dataIndex][0];
        console.log(break_point);
        reset_all();
    }
});

// set default country and default intervention break point
var break_point = 22;
var series_point;
var element = document.getElementById('country2');
element.value = 'South Africa';
var event = new Event('change');
element.dispatchEvent(event);

update_plot();
update_counters();

$("#reset-button").click(reset_all);
$("#start-button").click(start_all);
$("#stop-button").click(stop_all);
$("#continue-button").click(continue_all);

setInterval(run_SIR, time_interval);


function run_SIR() {
    
    if (running == 0)
       return;

    
    var s = epi_state.S;
    var i = epi_state.I;
    var r = epi_state.R;
    
    var beta = p_SI;
    var gamma = p_IR;
    
    if (count >= break_point) {
       beta = p_SI_ld;
       gamma = p_IR_ld;
    }
    
    epi_state.S += ( - beta * s * i );
    epi_state.I += ( + beta * s * i - gamma * i );
    epi_state.R += (gamma * i);
                    
    count++;

    timeseries.S.data.push([count, Math.round(N*epi_state.S)]);
    timeseries.I.data.push([count, Math.round(N*epi_state.I)]);
    timeseries.R.data.push([count, Math.round(N*epi_state.R)]);
    timeseries.C.data.push([count, N*(epi_state.R+epi_state.I)]);

    update_plot();

    update_counters();

    if (Math.round(N*epi_state.I) == 0)
       running = 0;
}

function update_plot () {
 plot.setData([timeseries.D, timeseries.I, timeseries.R,
               timeseries.C, timeseries.B ]);
 plot.setupGrid();
 plot.draw();
}


function update_counters () {
  $("#count_S").html(Math.round(N*epi_state.S));
  $("#count_I").html(Math.round(N*epi_state.I));
  $("#count_R").html(Math.round(N*epi_state.R));
}

function update_country () {
    country = $("#country").val
    console.log("update_country()",country)
}

function update_p_SI () {
  p = Number($("#p_SI").val());
  if (isNaN(p) || p<0.0 || p>1.0) {
     valid_data = 0;
     $("#p_SI").css("background-color", "#f88");
  } else {
     p_SI = p;
     valid_data = 1;
     $("#p_SI").css("background-color", "#fff");
  }
}

function update_p_SI_ld () {
  p = Number($("#p_SI_ld").val());
  if (isNaN(p) || p<0.0 || p>1.0) {
     valid_data = 0;
     $("#p_SI_ld").css("background-color", "#f88");
  } else {
     p_SI_ld = p;
     valid_data = 1;
     $("#p_SI_ld").css("background-color", "#fff");
  }
}

function update_p_IR () {
  p = Number($("#p_IR").val());
  if (isNaN(p) || p<0.0 || p>1.0) {
     valid_data = 0;
     $("#p_IR").css("background-color", "#f88");
  } else {
     p_IR = p;
     valid_data = 1;
     $("#p_IR").css("background-color", "#fff");
  }
}

function update_p_IR_ld () {
  p = Number($("#p_IR_ld").val());
  if (isNaN(p) || p<0.0 || p>1.0) {
     valid_data = 0;
     $("#p_IR_ld").css("background-color", "#f88");
  } else {
     p_IR_ld = p;
     valid_data = 1;
     $("#p_IR_ld").css("background-color", "#fff");
  }
}

function reset_all () {
	running = 0;

    count = 0;
    epi_state = { S: (N-1)/N, I: 1/N, R: 0 };

    reset_params();
    reset_history();
    update_plot();

  update_counters();
    
}

function start_all () {
    running = 1;
    count = 0;
    epi_state = { S: (N-1)/N, I: 1/N, R: 0 };
    reset_history();
    update_plot();
    update_counters();
}

function stop_all () {
    running = 0;
    update_plot();
    update_counters();
}

function continue_all () {
    running = 1;
    update_plot();
    update_counters();
}
