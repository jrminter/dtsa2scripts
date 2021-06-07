# Generated by the NIST EPQ Library
# On Jun 7, 2021
# By johnr
# NOTE: TO SAVE FILES AS ENCAPSULATED POST SCRIPT FILES, OR AS OTHER FILE FORMATS,
#     COMMENT OUT THE "TERMINAL WINDOW" SECTION, SET THE TERMINAL AS
#     THE CHOSEN TYPE, AND SET THE OUTPUT.
# set terminal postscript enhanced "Times-Roman" 12
# set terminal png
set terminal window
unset logscale x
set logscale y
set xrange [0.0:5000.0]
set yrange [1.0:140000.0]
set xlabel "(Energy) eV" font "Arial,16"
set ylabel "" font "Arial,16"
set format x "%g"
set format y "%g"
unset label
# B
set label 1 "B" at first 183.3,3241.0 center
# Al
set label 2 "Al" at first 1554.7,2844.6 center
# Si
set label 3 "Si" at first 1833.8,3810.1 center
# Si
set label 4 "Si" at first 1739.7,18942.3 center
# K
set label 5 "K" at first 3589.6,2817.2 center
# Al
set label 6 "Al" at first 1486.5,4504.0 center
# K
set label 7 "K" at first 3313.8,2945.0 center
# K
set label 8 "K" at first 3311.1,2873.2 center
# K
set label 9 "K" at first 259.7,3646.0 center
# O
set label 10 "O" at first 524.9,49649.0 center
# Na
set label 11 "Na" at first 1041.0,7995.0 center
# Na
set label 12 "Na" at first 1041.0,5397.5 center
# O
set label 13 "O" at first 524.9,96498.0 center
# Si
set label 14 "Si" at first 1739.7,34765.0 center
# Al
set label 15 "Al" at first 1486.5,3663.9 center
# K
set label 16 "K" at first 3589.6,2809.7 center
plot "-" with lines  title "Sim-Pyrex-5kV", "-" with impulses lt -1 notitle
# Sim-Pyrex-5kV
-0.5, 0.0
4.5, 0.0
9.5, 0.0
14.5, 2.0
19.5, 3.0
24.5, 1.0
29.5, 5.0
34.5, 8.0
39.5, 8.0
44.6, 23.0
49.6, 32.0
54.6, 39.0
59.6, 43.0
64.6, 63.0
69.6, 82.0
74.6, 106.0
79.6, 109.0
84.6, 127.0
89.6, 148.0
94.6, 145.0
99.6, 190.0
104.6, 159.0
109.6, 167.0
114.6, 141.0
119.6, 166.0
124.6, 168.0
129.6, 152.0
134.6, 202.0
139.6, 214.0
144.6, 241.0
149.6, 246.0
154.6, 264.0
159.6, 277.0
164.6, 297.0
169.6, 311.0
174.6, 351.0
179.6, 351.0
184.6, 361.0
189.6, 389.0
194.6, 428.0
199.6, 441.0
204.7, 439.0
209.7, 488.0
214.7, 480.0
219.7, 529.0
224.7, 585.0
229.7, 622.0
234.7, 671.0
239.7, 695.0
244.7, 718.0
249.7, 773.0
254.7, 804.0
259.7, 752.0
264.7, 799.0
269.7, 846.0
274.7, 749.0
279.7, 737.0
284.7, 756.0
289.7, 687.0
294.7, 671.0
299.7, 568.0
304.7, 561.0
309.7, 530.0
314.7, 512.0
319.7, 460.0
324.7, 422.0
329.7, 453.0
334.7, 431.0
339.7, 428.0
344.7, 457.0
349.7, 436.0
354.7, 468.0
359.7, 471.0
364.8, 462.0
369.8, 508.0
374.8, 532.0
379.8, 566.0
384.8, 529.0
389.8, 593.0
394.8, 591.0
399.8, 598.0
404.8, 631.0
409.8, 602.0
414.8, 592.0
419.8, 616.0
424.8, 725.0
429.8, 752.0
434.8, 846.0
439.8, 1091.0
444.8, 1395.0
449.8, 1886.0
454.8, 2770.0
459.8, 4152.0
464.8, 6206.0
469.8, 9288.0
474.8, 13385.0
479.8, 18609.0
484.8, 26186.0
489.8, 34448.0
494.8, 44086.0
499.8, 54376.0
504.8, 65702.0
509.8, 75567.0
514.8, 84375.0
519.8, 90465.0
524.9, 93547.0
529.9, 93698.0
534.9, 90444.0
539.9, 84236.0
544.9, 75680.0
549.9, 65457.0
554.9, 54790.0
559.9, 43456.0
564.9, 34248.0
569.9, 25721.0
574.9, 18669.0
579.9, 12953.0
584.9, 9109.0
589.9, 6080.0
594.9, 4085.0
599.9, 2662.0
604.9, 1883.0
609.9, 1261.0
614.9, 931.0
619.9, 763.0
624.9, 679.0
629.9, 572.0
634.9, 588.0
639.9, 602.0
644.9, 535.0
649.9, 581.0
654.9, 591.0
659.9, 579.0
664.9, 595.0
669.9, 584.0
674.9, 624.0
679.9, 621.0
685.0, 586.0
690.0, 595.0
695.0, 622.0
700.0, 617.0
705.0, 612.0
710.0, 654.0
715.0, 642.0
720.0, 626.0
725.0, 675.0
730.0, 651.0
735.0, 649.0
740.0, 659.0
745.0, 699.0
750.0, 671.0
755.0, 627.0
760.0, 699.0
765.0, 661.0
770.0, 647.0
775.0, 630.0
780.0, 658.0
785.0, 662.0
790.0, 655.0
795.0, 671.0
800.0, 629.0
805.0, 671.0
810.0, 670.0
815.0, 708.0
820.0, 692.0
825.0, 680.0
830.0, 652.0
835.0, 678.0
840.0, 651.0
845.1, 643.0
850.1, 679.0
855.1, 671.0
860.1, 628.0
865.1, 678.0
870.1, 643.0
875.1, 652.0
880.1, 670.0
885.1, 629.0
890.1, 680.0
895.1, 662.0
900.1, 692.0
905.1, 718.0
910.1, 684.0
915.1, 676.0
920.1, 645.0
925.1, 638.0
930.1, 647.0
935.1, 694.0
940.1, 686.0
945.1, 665.0
950.1, 718.0
955.1, 694.0
960.1, 779.0
965.1, 838.0
970.1, 927.0
975.1, 1061.0
980.1, 1149.0
985.1, 1401.0
990.1, 1672.0
995.1, 1966.0
1000.1, 2360.0
1005.2, 2691.0
1010.2, 3103.0
1015.2, 3580.0
1020.2, 4082.0
1025.2, 4411.0
1030.2, 4778.0
1035.2, 4997.0
1040.2, 5166.0
1045.2, 5195.0
1050.2, 5105.0
1055.2, 4858.0
1060.2, 4427.0
1065.2, 4059.0
1070.2, 3476.0
1075.2, 3196.0
1080.2, 2736.0
1085.2, 2269.0
1090.2, 1987.0
1095.2, 1653.0
1100.2, 1323.0
1105.2, 1156.0
1110.2, 964.0
1115.2, 900.0
1120.2, 747.0
1125.2, 705.0
1130.2, 673.0
1135.2, 603.0
1140.2, 640.0
1145.2, 585.0
1150.2, 633.0
1155.2, 598.0
1160.2, 608.0
1165.3, 609.0
1170.3, 601.0
1175.3, 592.0
1180.3, 586.0
1185.3, 600.0
1190.3, 544.0
1195.3, 569.0
1200.3, 585.0
1205.3, 560.0
1210.3, 562.0
1215.3, 569.0
1220.3, 584.0
1225.3, 564.0
1230.3, 550.0
1235.3, 577.0
1240.3, 569.0
1245.3, 559.0
1250.3, 524.0
1255.3, 543.0
1260.3, 581.0
1265.3, 551.0
1270.3, 590.0
1275.3, 540.0
1280.3, 512.0
1285.3, 550.0
1290.3, 561.0
1295.3, 503.0
1300.3, 543.0
1305.3, 535.0
1310.3, 542.0
1315.3, 537.0
1320.3, 525.0
1325.4, 517.0
1330.4, 529.0
1335.4, 501.0
1340.4, 548.0
1345.4, 522.0
1350.4, 508.0
1355.4, 508.0
1360.4, 519.0
1365.4, 537.0
1370.4, 506.0
1375.4, 527.0
1380.4, 541.0
1385.4, 518.0
1390.4, 544.0
1395.4, 533.0
1400.4, 499.0
1405.4, 520.0
1410.4, 528.0
1415.4, 672.0
1420.4, 634.0
1425.4, 726.0
1430.4, 751.0
1435.4, 801.0
1440.4, 925.0
1445.4, 1001.0
1450.4, 1089.0
1455.4, 1254.0
1460.4, 1362.0
1465.4, 1406.0
1470.4, 1535.0
1475.4, 1547.0
1480.4, 1649.0
1485.5, 1633.0
1490.5, 1704.0
1495.5, 1672.0
1500.5, 1598.0
1505.5, 1539.0
1510.5, 1419.0
1515.5, 1282.0
1520.5, 1250.0
1525.5, 1109.0
1530.5, 979.0
1535.5, 873.0
1540.5, 841.0
1545.5, 700.0
1550.5, 618.0
1555.5, 588.0
1560.5, 564.0
1565.5, 535.0
1570.5, 463.0
1575.5, 442.0
1580.5, 435.0
1585.5, 454.0
1590.5, 455.0
1595.5, 442.0
1600.5, 418.0
1605.5, 450.0
1610.5, 436.0
1615.5, 469.0
1620.5, 458.0
1625.5, 567.0
1630.5, 622.0
1635.5, 671.0
1640.5, 869.0
1645.6, 1019.0
1650.6, 1323.0
1655.6, 1706.0
1660.6, 2332.0
1665.6, 3112.0
1670.6, 3891.0
1675.6, 5129.0
1680.6, 6722.0
1685.6, 8634.0
1690.6, 10544.0
1695.6, 12660.0
1700.6, 15293.0
1705.6, 17893.0
1710.6, 20543.0
1715.6, 23398.0
1720.6, 26033.0
1725.6, 27975.0
1730.6, 29677.0
1735.6, 30752.0
1740.6, 31627.0
1745.6, 31965.0
1750.6, 31153.0
1755.6, 29491.0
1760.6, 27747.0
1765.6, 25911.0
1770.6, 23286.0
1775.6, 20494.0
1780.6, 18032.0
1785.6, 15357.0
1790.6, 12956.0
1795.6, 10516.0
1800.6, 8479.0
1805.7, 6909.0
1810.7, 5373.0
1815.7, 4215.0
1820.7, 3273.0
1825.7, 2447.0
1830.7, 1862.0
1835.7, 1488.0
1840.7, 1140.0
1845.7, 935.0
1850.7, 784.0
1855.7, 636.0
1860.7, 579.0
1865.7, 534.0
1870.7, 489.0
1875.7, 448.0
1880.7, 407.0
1885.7, 402.0
1890.7, 365.0
1895.7, 358.0
1900.7, 382.0
1905.7, 320.0
1910.7, 319.0
1915.7, 329.0
1920.7, 351.0
1925.7, 325.0
1930.7, 304.0
1935.7, 279.0
1940.7, 279.0
1945.7, 310.0
1950.7, 280.0
1955.7, 296.0
1960.7, 299.0
1965.8, 284.0
1970.8, 275.0
1975.8, 305.0
1980.8, 296.0
1985.8, 314.0
1990.8, 271.0
1995.8, 324.0
2000.8, 286.0
2005.8, 266.0
2010.8, 274.0
2015.8, 288.0
2020.8, 306.0
2025.8, 280.0
2030.8, 270.0
2035.8, 273.0
2040.8, 264.0
2045.8, 269.0
2050.8, 252.0
2055.8, 275.0
2060.8, 273.0
2065.8, 278.0
2070.8, 295.0
2075.8, 268.0
2080.8, 277.0
2085.8, 295.0
2090.8, 256.0
2095.8, 253.0
2100.8, 282.0
2105.8, 245.0
2110.8, 276.0
2115.8, 254.0
2120.8, 288.0
2125.9, 252.0
2130.9, 270.0
2135.9, 282.0
2140.9, 252.0
2145.9, 261.0
2150.9, 239.0
2155.9, 258.0
2160.9, 233.0
2165.9, 262.0
2170.9, 278.0
2175.9, 267.0
2180.9, 241.0
2185.9, 249.0
2190.9, 261.0
2195.9, 249.0
2200.9, 247.0
2205.9, 263.0
2210.9, 250.0
2215.9, 218.0
2220.9, 248.0
2225.9, 267.0
2230.9, 244.0
2235.9, 212.0
2240.9, 240.0
2245.9, 236.0
2250.9, 251.0
2255.9, 242.0
2260.9, 235.0
2265.9, 235.0
2270.9, 249.0
2275.9, 229.0
2280.9, 228.0
2286.0, 214.0
2291.0, 245.0
2296.0, 217.0
2301.0, 233.0
2306.0, 221.0
2311.0, 219.0
2316.0, 232.0
2321.0, 224.0
2326.0, 195.0
2331.0, 231.0
2336.0, 233.0
2341.0, 219.0
2346.0, 252.0
2351.0, 219.0
2356.0, 218.0
2361.0, 240.0
2366.0, 213.0
2371.0, 225.0
2376.0, 218.0
2381.0, 195.0
2386.0, 227.0
2391.0, 233.0
2396.0, 213.0
2401.0, 188.0
2406.0, 208.0
2411.0, 188.0
2416.0, 222.0
2421.0, 210.0
2426.0, 227.0
2431.0, 209.0
2436.0, 192.0
2441.0, 220.0
2446.1, 214.0
2451.1, 214.0
2456.1, 194.0
2461.1, 194.0
2466.1, 196.0
2471.1, 199.0
2476.1, 196.0
2481.1, 217.0
2486.1, 182.0
2491.1, 231.0
2496.1, 223.0
2501.1, 194.0
2506.1, 194.0
2511.1, 190.0
2516.1, 181.0
2521.1, 169.0
2526.1, 206.0
2531.1, 195.0
2536.1, 214.0
2541.1, 205.0
2546.1, 201.0
2551.1, 180.0
2556.1, 195.0
2561.1, 201.0
2566.1, 207.0
2571.1, 171.0
2576.1, 189.0
2581.1, 179.0
2586.1, 181.0
2591.1, 216.0
2596.1, 176.0
2601.1, 192.0
2606.2, 176.0
2611.2, 172.0
2616.2, 226.0
2621.2, 150.0
2626.2, 172.0
2631.2, 179.0
2636.2, 162.0
2641.2, 156.0
2646.2, 184.0
2651.2, 171.0
2656.2, 167.0
2661.2, 187.0
2666.2, 164.0
2671.2, 169.0
2676.2, 146.0
2681.2, 169.0
2686.2, 186.0
2691.2, 180.0
2696.2, 189.0
2701.2, 155.0
2706.2, 174.0
2711.2, 159.0
2716.2, 146.0
2721.2, 169.0
2726.2, 141.0
2731.2, 163.0
2736.2, 140.0
2741.2, 159.0
2746.2, 141.0
2751.2, 145.0
2756.2, 168.0
2761.2, 164.0
2766.3, 137.0
2771.3, 164.0
2776.3, 161.0
2781.3, 186.0
2786.3, 167.0
2791.3, 143.0
2796.3, 138.0
2801.3, 150.0
2806.3, 145.0
2811.3, 161.0
2816.3, 173.0
2821.3, 155.0
2826.3, 159.0
2831.3, 164.0
2836.3, 150.0
2841.3, 161.0
2846.3, 137.0
2851.3, 152.0
2856.3, 129.0
2861.3, 163.0
2866.3, 145.0
2871.3, 144.0
2876.3, 155.0
2881.3, 152.0
2886.3, 152.0
2891.3, 128.0
2896.3, 159.0
2901.3, 127.0
2906.3, 150.0
2911.3, 132.0
2916.3, 147.0
2921.3, 132.0
2926.4, 147.0
2931.4, 139.0
2936.4, 146.0
2941.4, 144.0
2946.4, 140.0
2951.4, 155.0
2956.4, 131.0
2961.4, 135.0
2966.4, 138.0
2971.4, 129.0
2976.4, 138.0
2981.4, 142.0
2986.4, 150.0
2991.4, 142.0
2996.4, 121.0
3001.4, 135.0
3006.4, 145.0
3011.4, 165.0
3016.4, 138.0
3021.4, 104.0
3026.4, 113.0
3031.4, 137.0
3036.4, 115.0
3041.4, 128.0
3046.4, 136.0
3051.4, 124.0
3056.4, 138.0
3061.4, 127.0
3066.4, 127.0
3071.4, 172.0
3076.4, 123.0
3081.4, 119.0
3086.5, 125.0
3091.5, 150.0
3096.5, 121.0
3101.5, 131.0
3106.5, 123.0
3111.5, 130.0
3116.5, 138.0
3121.5, 118.0
3126.5, 117.0
3131.5, 122.0
3136.5, 110.0
3141.5, 122.0
3146.5, 105.0
3151.5, 106.0
3156.5, 116.0
3161.5, 118.0
3166.5, 120.0
3171.5, 130.0
3176.5, 117.0
3181.5, 99.0
3186.5, 110.0
3191.5, 121.0
3196.5, 110.0
3201.5, 116.0
3206.5, 112.0
3211.5, 110.0
3216.5, 94.0
3221.5, 114.0
3226.5, 104.0
3231.5, 118.0
3236.5, 117.0
3241.5, 110.0
3246.6, 114.0
3251.6, 108.0
3256.6, 102.0
3261.6, 139.0
3266.6, 127.0
3271.6, 113.0
3276.6, 125.0
3281.6, 98.0
3286.6, 106.0
3291.6, 145.0
3296.6, 127.0
3301.6, 139.0
3306.6, 140.0
3311.6, 129.0
3316.6, 114.0
3321.6, 127.0
3326.6, 138.0
3331.6, 124.0
3336.6, 122.0
3341.6, 119.0
3346.6, 117.0
3351.6, 113.0
3356.6, 115.0
3361.6, 120.0
3366.6, 102.0
3371.6, 90.0
3376.6, 101.0
3381.6, 114.0
3386.6, 92.0
3391.6, 119.0
3396.6, 91.0
3401.6, 84.0
3406.7, 99.0
3411.7, 97.0
3416.7, 102.0
3421.7, 99.0
3426.7, 87.0
3431.7, 89.0
3436.7, 85.0
3441.7, 76.0
3446.7, 93.0
3451.7, 98.0
3456.7, 76.0
3461.7, 93.0
3466.7, 72.0
3471.7, 83.0
3476.7, 69.0
3481.7, 78.0
3486.7, 91.0
3491.7, 86.0
3496.7, 85.0
3501.7, 74.0
3506.7, 76.0
3511.7, 90.0
3516.7, 77.0
3521.7, 91.0
3526.7, 97.0
3531.7, 96.0
3536.7, 94.0
3541.7, 87.0
3546.7, 77.0
3551.7, 83.0
3556.7, 72.0
3561.7, 84.0
3566.8, 67.0
3571.8, 88.0
3576.8, 75.0
3581.8, 70.0
3586.8, 71.0
3591.8, 80.0
3596.8, 85.0
3601.8, 80.0
3606.8, 81.0
3611.8, 65.0
3616.8, 79.0
3621.8, 85.0
3626.8, 77.0
3631.8, 68.0
3636.8, 75.0
3641.8, 73.0
3646.8, 68.0
3651.8, 64.0
3656.8, 76.0
3661.8, 78.0
3666.8, 73.0
3671.8, 73.0
3676.8, 75.0
3681.8, 65.0
3686.8, 60.0
3691.8, 82.0
3696.8, 51.0
3701.8, 77.0
3706.8, 56.0
3711.8, 58.0
3716.8, 83.0
3721.8, 67.0
3726.9, 66.0
3731.9, 65.0
3736.9, 72.0
3741.9, 61.0
3746.9, 62.0
3751.9, 58.0
3756.9, 62.0
3761.9, 61.0
3766.9, 66.0
3771.9, 68.0
3776.9, 70.0
3781.9, 60.0
3786.9, 64.0
3791.9, 57.0
3796.9, 62.0
3801.9, 56.0
3806.9, 65.0
3811.9, 68.0
3816.9, 47.0
3821.9, 67.0
3826.9, 56.0
3831.9, 56.0
3836.9, 57.0
3841.9, 58.0
3846.9, 60.0
3851.9, 59.0
3856.9, 53.0
3861.9, 56.0
3866.9, 47.0
3871.9, 53.0
3876.9, 52.0
3881.9, 62.0
3887.0, 65.0
3892.0, 68.0
3897.0, 51.0
3902.0, 54.0
3907.0, 70.0
3912.0, 56.0
3917.0, 40.0
3922.0, 59.0
3927.0, 57.0
3932.0, 57.0
3937.0, 52.0
3942.0, 60.0
3947.0, 47.0
3952.0, 44.0
3957.0, 56.0
3962.0, 51.0
3967.0, 42.0
3972.0, 36.0
3977.0, 54.0
3982.0, 42.0
3987.0, 48.0
3992.0, 45.0
3997.0, 51.0
4002.0, 43.0
4007.0, 52.0
4012.0, 44.0
4017.0, 40.0
4022.0, 35.0
4027.0, 48.0
4032.0, 54.0
4037.0, 36.0
4042.0, 39.0
4047.1, 48.0
4052.1, 45.0
4057.1, 37.0
4062.1, 44.0
4067.1, 46.0
4072.1, 47.0
4077.1, 46.0
4082.1, 46.0
4087.1, 47.0
4092.1, 40.0
4097.1, 48.0
4102.1, 52.0
4107.1, 53.0
4112.1, 37.0
4117.1, 40.0
4122.1, 53.0
4127.1, 52.0
4132.1, 63.0
4137.1, 37.0
4142.1, 30.0
4147.1, 38.0
4152.1, 44.0
4157.1, 36.0
4162.1, 36.0
4167.1, 38.0
4172.1, 35.0
4177.1, 40.0
4182.1, 30.0
4187.1, 38.0
4192.1, 29.0
4197.1, 38.0
4202.1, 35.0
4207.2, 37.0
4212.2, 39.0
4217.2, 39.0
4222.2, 36.0
4227.2, 42.0
4232.2, 48.0
4237.2, 40.0
4242.2, 37.0
4247.2, 36.0
4252.2, 36.0
4257.2, 49.0
4262.2, 39.0
4267.2, 40.0
4272.2, 32.0
4277.2, 23.0
4282.2, 38.0
4287.2, 25.0
4292.2, 30.0
4297.2, 35.0
4302.2, 44.0
4307.2, 35.0
4312.2, 32.0
4317.2, 48.0
4322.2, 35.0
4327.2, 32.0
4332.2, 22.0
4337.2, 30.0
4342.2, 28.0
4347.2, 32.0
4352.2, 27.0
4357.2, 24.0
4362.2, 28.0
4367.3, 27.0
4372.3, 35.0
4377.3, 23.0
4382.3, 21.0
4387.3, 27.0
4392.3, 28.0
4397.3, 36.0
4402.3, 27.0
4407.3, 26.0
4412.3, 29.0
4417.3, 20.0
4422.3, 19.0
4427.3, 30.0
4432.3, 26.0
4437.3, 33.0
4442.3, 20.0
4447.3, 32.0
4452.3, 23.0
4457.3, 21.0
4462.3, 23.0
4467.3, 28.0
4472.3, 23.0
4477.3, 21.0
4482.3, 23.0
4487.3, 25.0
4492.3, 23.0
4497.3, 22.0
4502.3, 16.0
4507.3, 24.0
4512.3, 21.0
4517.3, 27.0
4522.3, 18.0
4527.4, 15.0
4532.4, 21.0
4537.4, 15.0
4542.4, 18.0
4547.4, 17.0
4552.4, 23.0
4557.4, 18.0
4562.4, 18.0
4567.4, 16.0
4572.4, 21.0
4577.4, 14.0
4582.4, 12.0
4587.4, 15.0
4592.4, 17.0
4597.4, 14.0
4602.4, 11.0
4607.4, 12.0
4612.4, 13.0
4617.4, 21.0
4622.4, 13.0
4627.4, 9.0
4632.4, 16.0
4637.4, 15.0
4642.4, 7.0
4647.4, 11.0
4652.4, 11.0
4657.4, 10.0
4662.4, 16.0
4667.4, 14.0
4672.4, 16.0
4677.4, 14.0
4682.4, 6.0
4687.4, 10.0
4692.5, 14.0
4697.5, 18.0
4702.5, 10.0
4707.5, 16.0
4712.5, 11.0
4717.5, 11.0
4722.5, 9.0
4727.5, 10.0
4732.5, 18.0
4737.5, 16.0
4742.5, 8.0
4747.5, 10.0
4752.5, 9.0
4757.5, 10.0
4762.5, 13.0
4767.5, 8.0
4772.5, 9.0
4777.5, 6.0
4782.5, 8.0
4787.5, 10.0
4792.5, 9.0
4797.5, 3.0
4802.5, 12.0
4807.5, 8.0
4812.5, 6.0
4817.5, 7.0
4822.5, 3.0
4827.5, 5.0
4832.5, 9.0
4837.5, 5.0
4842.5, 3.0
4847.5, 9.0
4852.6, 4.0
4857.6, 4.0
4862.6, 10.0
4867.6, 3.0
4872.6, 4.0
4877.6, 5.0
4882.6, 5.0
4887.6, 5.0
4892.6, 4.0
4897.6, 3.0
4902.6, 5.0
4907.6, 5.0
4912.6, 2.0
4917.6, 4.0
4922.6, 3.0
4927.6, 1.0
4932.6, 3.0
4937.6, 4.0
4942.6, 3.0
4947.6, 4.0
4952.6, 0.0
4957.6, 1.0
4962.6, 5.0
4967.6, 2.0
4972.6, 3.0
4977.6, 2.0
4982.6, 0.0
4987.6, 0.0
4992.6, 3.0
e
# B
183.3, 441.0
# Al
1554.7, 44.6
# Si
1833.8, 1010.1
# Si
1739.7, 16142.3
# K
3589.6, 17.3
# Al
1486.5, 1704.0
# K
3313.8, 145.0
# K
3311.1, 73.2
# K
259.7, 846.0
# O
524.9, 46849.0
# Na
1041.0, 5195.0
# Na
1041.0, 2597.5
# O
524.9, 93698.0
# Si
1739.7, 31965.0
# Al
1486.5, 863.9
# K
3589.6, 9.7
e
