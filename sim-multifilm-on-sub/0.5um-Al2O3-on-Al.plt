# Generated by the NIST EPQ Library
# On Apr 18, 2016
# By jrminter
# NOTE: TO SAVE FILES AS ENCAPSULATED POST SCRIPT FILES, OR AS OTHER FILE FORMATS,
#     COMMENT OUT THE "TERMINAL WINDOW" SECTION, SET THE TERMINAL AS
#     THE CHOSEN TYPE, AND SET THE OUTPUT.
# set terminal postscript enhanced "Times-Roman" 12
# set terminal png
set terminal window
unset logscale xy
set xrange [0.0:2000.0]
set yrange [0.0:680000.0]
set xlabel "(Energy) eV" font "Arial,16"
set ylabel "" font "Arial,16"
set format x "%g"
set format y "%g"
unset label
# O K-L3
set label 1 "O K-L3" at first 524.9,96050.0 center
# Al K-L3
set label 2 "Al K-L3" at first 1486.5,681650.0 center
plot "-" with lines  title "7 kV", "-" with lines  title "10 kV", "-" with lines  title "15 kV", "-" with lines  title "20 kV", "-" with lines  title "30 kV", "-" with impulses lt -1 notitle
# 7 kV
-2.1, 0.0
2.9, 0.0
7.9, 0.0
12.9, 0.0
17.9, 0.0
22.9, 0.0
27.9, 0.0
32.9, 0.0
37.9, 0.0
42.9, 1.0
47.9, 4.0
52.9, 3.0
57.9, 2.0
62.8, 6.0
67.8, 4.0
72.8, 7.0
77.8, 2.0
82.8, 6.0
87.8, 10.0
92.8, 14.0
97.8, 19.0
102.8, 30.0
107.8, 33.0
112.8, 54.0
117.8, 58.0
122.8, 65.0
127.8, 75.0
132.8, 71.0
137.8, 98.0
142.8, 109.0
147.7, 118.0
152.7, 142.0
157.7, 145.0
162.7, 163.0
167.7, 177.0
172.7, 195.0
177.7, 216.0
182.7, 249.0
187.7, 225.0
192.7, 267.0
197.7, 299.0
202.7, 296.0
207.7, 314.0
212.7, 343.0
217.7, 379.0
222.7, 410.0
227.6, 413.0
232.6, 490.0
237.6, 503.0
242.6, 563.0
247.6, 558.0
252.6, 582.0
257.6, 623.0
262.6, 685.0
267.6, 688.0
272.6, 671.0
277.6, 667.0
282.6, 661.0
287.6, 570.0
292.6, 549.0
297.6, 536.0
302.6, 471.0
307.6, 476.0
312.5, 417.0
317.5, 368.0
322.5, 408.0
327.5, 347.0
332.5, 419.0
337.5, 379.0
342.5, 405.0
347.5, 435.0
352.5, 538.0
357.5, 494.0
362.5, 512.0
367.5, 530.0
372.5, 555.0
377.5, 601.0
382.5, 620.0
387.5, 610.0
392.4, 598.0
397.4, 618.0
402.4, 661.0
407.4, 682.0
412.4, 691.0
417.4, 684.0
422.4, 717.0
427.4, 742.0
432.4, 773.0
437.4, 820.0
442.4, 844.0
447.4, 846.0
452.4, 987.0
457.4, 1189.0
462.4, 1499.0
467.4, 2143.0
472.3, 3205.0
477.3, 4913.0
482.3, 7761.0
487.3, 11682.0
492.3, 17064.0
497.3, 23955.0
502.3, 32633.0
507.3, 41231.0
512.3, 51234.0
517.3, 58457.0
522.3, 64742.0
527.3, 68437.0
532.3, 67866.0
537.3, 64467.0
542.3, 58408.0
547.3, 50371.0
552.3, 41279.0
557.2, 32057.0
562.2, 23640.0
567.2, 16639.0
572.2, 11382.0
577.2, 7515.0
582.2, 4714.0
587.2, 2943.0
592.2, 1778.0
597.2, 1213.0
602.2, 942.0
607.2, 808.0
612.2, 675.0
617.2, 641.0
622.2, 628.0
627.2, 604.0
632.2, 656.0
637.1, 648.0
642.1, 635.0
647.1, 619.0
652.1, 720.0
657.1, 683.0
662.1, 681.0
667.1, 719.0
672.1, 718.0
677.1, 675.0
682.1, 729.0
687.1, 722.0
692.1, 748.0
697.1, 762.0
702.1, 731.0
707.1, 726.0
712.1, 762.0
717.1, 779.0
722.0, 798.0
727.0, 770.0
732.0, 757.0
737.0, 725.0
742.0, 844.0
747.0, 799.0
752.0, 832.0
757.0, 808.0
762.0, 795.0
767.0, 821.0
772.0, 820.0
777.0, 830.0
782.0, 856.0
787.0, 766.0
792.0, 863.0
797.0, 815.0
801.9, 823.0
806.9, 866.0
811.9, 835.0
816.9, 814.0
821.9, 858.0
826.9, 872.0
831.9, 843.0
836.9, 860.0
841.9, 851.0
846.9, 829.0
851.9, 824.0
856.9, 851.0
861.9, 842.0
866.9, 879.0
871.9, 852.0
876.9, 855.0
881.9, 892.0
886.8, 846.0
891.8, 791.0
896.8, 833.0
901.8, 861.0
906.8, 810.0
911.8, 836.0
916.8, 873.0
921.8, 849.0
926.8, 849.0
931.8, 848.0
936.8, 890.0
941.8, 815.0
946.8, 836.0
951.8, 829.0
956.8, 862.0
961.8, 922.0
966.7, 855.0
971.7, 830.0
976.7, 850.0
981.7, 849.0
986.7, 834.0
991.7, 829.0
996.7, 849.0
1001.7, 896.0
1006.7, 847.0
1011.7, 874.0
1016.7, 871.0
1021.7, 866.0
1026.7, 849.0
1031.7, 879.0
1036.7, 843.0
1041.7, 808.0
1046.7, 886.0
1051.6, 823.0
1056.6, 872.0
1061.6, 903.0
1066.6, 856.0
1071.6, 840.0
1076.6, 881.0
1081.6, 826.0
1086.6, 768.0
1091.6, 858.0
1096.6, 870.0
1101.6, 804.0
1106.6, 882.0
1111.6, 831.0
1116.6, 844.0
1121.6, 831.0
1126.6, 826.0
1131.5, 825.0
1136.5, 825.0
1141.5, 855.0
1146.5, 827.0
1151.5, 835.0
1156.5, 809.0
1161.5, 854.0
1166.5, 768.0
1171.5, 821.0
1176.5, 741.0
1181.5, 810.0
1186.5, 811.0
1191.5, 801.0
1196.5, 779.0
1201.5, 813.0
1206.5, 784.0
1211.5, 789.0
1216.4, 784.0
1221.4, 822.0
1226.4, 840.0
1231.4, 802.0
1236.4, 787.0
1241.4, 745.0
1246.4, 771.0
1251.4, 769.0
1256.4, 798.0
1261.4, 779.0
1266.4, 786.0
1271.4, 770.0
1276.4, 702.0
1281.4, 797.0
1286.4, 827.0
1291.4, 773.0
1296.3, 797.0
1301.3, 769.0
1306.3, 758.0
1311.3, 763.0
1316.3, 793.0
1321.3, 764.0
1326.3, 739.0
1331.3, 737.0
1336.3, 771.0
1341.3, 741.0
1346.3, 798.0
1351.3, 750.0
1356.3, 746.0
1361.3, 712.0
1366.3, 731.0
1371.3, 782.0
1376.3, 794.0
1381.2, 794.0
1386.2, 936.0
1391.2, 989.0
1396.2, 1251.0
1401.2, 1651.0
1406.2, 2088.0
1411.2, 2870.0
1416.2, 4178.0
1421.2, 5800.0
1426.2, 7864.0
1431.2, 10565.0
1436.2, 14291.0
1441.2, 18305.0
1446.2, 23315.0
1451.2, 28747.0
1456.2, 34779.0
1461.1, 40950.0
1466.1, 47159.0
1471.1, 52204.0
1476.1, 56180.0
1481.1, 59883.0
1486.1, 60991.0
1491.1, 61224.0
1496.1, 59880.0
1501.1, 56821.0
1506.1, 52112.0
1511.1, 46725.0
1516.1, 41298.0
1521.1, 35211.0
1526.1, 29152.0
1531.1, 23455.0
1536.1, 18702.0
1541.1, 14547.0
1546.0, 11042.0
1551.0, 8078.0
1556.0, 5896.0
1561.0, 4524.0
1566.0, 3198.0
1571.0, 2298.0
1576.0, 1713.0
1581.0, 1330.0
1586.0, 1081.0
1591.0, 947.0
1596.0, 793.0
1601.0, 727.0
1606.0, 616.0
1611.0, 576.0
1616.0, 553.0
1621.0, 520.0
1625.9, 489.0
1630.9, 527.0
1635.9, 500.0
1640.9, 491.0
1645.9, 505.0
1650.9, 499.0
1655.9, 489.0
1660.9, 503.0
1665.9, 518.0
1670.9, 455.0
1675.9, 494.0
1680.9, 493.0
1685.9, 489.0
1690.9, 491.0
1695.9, 512.0
1700.9, 503.0
1705.9, 517.0
1710.8, 512.0
1715.8, 465.0
1720.8, 479.0
1725.8, 471.0
1730.8, 489.0
1735.8, 476.0
1740.8, 497.0
1745.8, 486.0
1750.8, 448.0
1755.8, 463.0
1760.8, 484.0
1765.8, 487.0
1770.8, 489.0
1775.8, 462.0
1780.8, 451.0
1785.8, 456.0
1790.7, 439.0
1795.7, 449.0
1800.7, 475.0
1805.7, 428.0
1810.7, 498.0
1815.7, 463.0
1820.7, 465.0
1825.7, 454.0
1830.7, 428.0
1835.7, 427.0
1840.7, 421.0
1845.7, 459.0
1850.7, 462.0
1855.7, 439.0
1860.7, 433.0
1865.7, 446.0
1870.7, 439.0
1875.6, 427.0
1880.6, 480.0
1885.6, 493.0
1890.6, 468.0
1895.6, 448.0
1900.6, 434.0
1905.6, 451.0
1910.6, 455.0
1915.6, 437.0
1920.6, 420.0
1925.6, 438.0
1930.6, 456.0
1935.6, 431.0
1940.6, 430.0
1945.6, 405.0
1950.6, 420.0
1955.5, 426.0
1960.5, 438.0
1965.5, 401.0
1970.5, 395.0
1975.5, 428.0
1980.5, 455.0
1985.5, 424.0
1990.5, 392.0
1995.5, 426.0
e
# 10 kV
-2.1, 0.0
2.9, 0.0
7.9, 0.0
12.9, 0.0
17.9, 0.0
22.9, 0.0
27.9, 0.0
32.9, 1.0
37.9, 0.0
42.9, 0.0
47.9, 3.0
52.9, 5.0
57.9, 5.0
62.8, 5.0
67.8, 2.0
72.8, 3.0
77.8, 4.0
82.8, 6.0
87.8, 6.0
92.8, 15.0
97.8, 12.0
102.8, 17.0
107.8, 25.0
112.8, 32.0
117.8, 35.0
122.8, 49.0
127.8, 41.0
132.8, 58.0
137.8, 58.0
142.8, 81.0
147.7, 74.0
152.7, 97.0
157.7, 78.0
162.7, 115.0
167.7, 125.0
172.7, 127.0
177.7, 122.0
182.7, 143.0
187.7, 199.0
192.7, 173.0
197.7, 176.0
202.7, 207.0
207.7, 196.0
212.7, 239.0
217.7, 257.0
222.7, 282.0
227.6, 314.0
232.6, 340.0
237.6, 336.0
242.6, 347.0
247.6, 404.0
252.6, 423.0
257.6, 457.0
262.6, 470.0
267.6, 455.0
272.6, 472.0
277.6, 486.0
282.6, 439.0
287.6, 413.0
292.6, 382.0
297.6, 405.0
302.6, 320.0
307.6, 342.0
312.5, 325.0
317.5, 316.0
322.5, 284.0
327.5, 309.0
332.5, 318.0
337.5, 334.0
342.5, 367.0
347.5, 397.0
352.5, 398.0
357.5, 444.0
362.5, 457.0
367.5, 465.0
372.5, 491.0
377.5, 533.0
382.5, 511.0
387.5, 561.0
392.4, 581.0
397.4, 570.0
402.4, 604.0
407.4, 615.0
412.4, 599.0
417.4, 661.0
422.4, 725.0
427.4, 732.0
432.4, 736.0
437.4, 772.0
442.4, 856.0
447.4, 942.0
452.4, 1091.0
457.4, 1297.0
462.4, 1728.0
467.4, 2549.0
472.3, 3880.0
477.3, 5893.0
482.3, 9130.0
487.3, 14154.0
492.3, 20778.0
497.3, 29143.0
502.3, 39099.0
507.3, 50403.0
512.3, 61038.0
517.3, 71486.0
522.3, 78641.0
527.3, 82450.0
532.3, 82417.0
537.3, 78274.0
542.3, 70984.0
547.3, 61092.0
552.3, 49651.0
557.2, 38607.0
562.2, 28850.0
567.2, 20031.0
572.2, 13639.0
577.2, 8780.0
582.2, 5569.0
587.2, 3369.0
592.2, 1958.0
597.2, 1304.0
602.2, 847.0
607.2, 719.0
612.2, 587.0
617.2, 566.0
622.2, 519.0
627.2, 552.0
632.2, 570.0
637.1, 591.0
642.1, 583.0
647.1, 556.0
652.1, 635.0
657.1, 624.0
662.1, 626.0
667.1, 619.0
672.1, 630.0
677.1, 659.0
682.1, 648.0
687.1, 652.0
692.1, 720.0
697.1, 719.0
702.1, 731.0
707.1, 694.0
712.1, 708.0
717.1, 752.0
722.0, 719.0
727.0, 727.0
732.0, 785.0
737.0, 817.0
742.0, 795.0
747.0, 837.0
752.0, 784.0
757.0, 792.0
762.0, 835.0
767.0, 845.0
772.0, 816.0
777.0, 864.0
782.0, 790.0
787.0, 834.0
792.0, 832.0
797.0, 829.0
801.9, 888.0
806.9, 959.0
811.9, 882.0
816.9, 895.0
821.9, 921.0
826.9, 962.0
831.9, 874.0
836.9, 954.0
841.9, 933.0
846.9, 905.0
851.9, 944.0
856.9, 968.0
861.9, 991.0
866.9, 958.0
871.9, 1007.0
876.9, 953.0
881.9, 936.0
886.8, 990.0
891.8, 1000.0
896.8, 978.0
901.8, 955.0
906.8, 966.0
911.8, 958.0
916.8, 1015.0
921.8, 1007.0
926.8, 988.0
931.8, 1007.0
936.8, 1040.0
941.8, 1044.0
946.8, 1038.0
951.8, 1019.0
956.8, 959.0
961.8, 1047.0
966.7, 997.0
971.7, 1019.0
976.7, 1013.0
981.7, 1116.0
986.7, 1042.0
991.7, 1060.0
996.7, 1033.0
1001.7, 1011.0
1006.7, 1017.0
1011.7, 1081.0
1016.7, 1034.0
1021.7, 1060.0
1026.7, 1017.0
1031.7, 1058.0
1036.7, 1056.0
1041.7, 1077.0
1046.7, 1062.0
1051.6, 1050.0
1056.6, 1079.0
1061.6, 1128.0
1066.6, 1016.0
1071.6, 1043.0
1076.6, 1042.0
1081.6, 1055.0
1086.6, 1042.0
1091.6, 993.0
1096.6, 1108.0
1101.6, 1080.0
1106.6, 1121.0
1111.6, 1014.0
1116.6, 1133.0
1121.6, 1069.0
1126.6, 1018.0
1131.5, 1095.0
1136.5, 1066.0
1141.5, 1077.0
1146.5, 1104.0
1151.5, 1078.0
1156.5, 1094.0
1161.5, 1030.0
1166.5, 1078.0
1171.5, 1077.0
1176.5, 1044.0
1181.5, 1060.0
1186.5, 1028.0
1191.5, 1126.0
1196.5, 1058.0
1201.5, 1059.0
1206.5, 1076.0
1211.5, 1183.0
1216.4, 1093.0
1221.4, 1122.0
1226.4, 1135.0
1231.4, 1104.0
1236.4, 1145.0
1241.4, 1088.0
1246.4, 1069.0
1251.4, 1003.0
1256.4, 1096.0
1261.4, 1100.0
1266.4, 1092.0
1271.4, 1089.0
1276.4, 1098.0
1281.4, 1024.0
1286.4, 995.0
1291.4, 1046.0
1296.3, 1021.0
1301.3, 1079.0
1306.3, 1040.0
1311.3, 1085.0
1316.3, 991.0
1321.3, 999.0
1326.3, 1053.0
1331.3, 1058.0
1336.3, 991.0
1341.3, 1056.0
1346.3, 1101.0
1351.3, 1045.0
1356.3, 1070.0
1361.3, 1057.0
1366.3, 990.0
1371.3, 1145.0
1376.3, 1102.0
1381.2, 1267.0
1386.2, 1304.0
1391.2, 1613.0
1396.2, 2048.0
1401.2, 2693.0
1406.2, 3846.0
1411.2, 5288.0
1416.2, 7484.0
1421.2, 10479.0
1426.2, 14387.0
1431.2, 19991.0
1436.2, 26365.0
1441.2, 34306.0
1446.2, 43532.0
1451.2, 53976.0
1456.2, 65009.0
1461.1, 76045.0
1466.1, 86888.0
1471.1, 97075.0
1476.1, 105712.0
1481.1, 111921.0
1486.1, 114855.0
1491.1, 114529.0
1496.1, 112175.0
1501.1, 105145.0
1506.1, 97493.0
1511.1, 87687.0
1516.1, 77034.0
1521.1, 65076.0
1526.1, 54198.0
1531.1, 44008.0
1536.1, 34704.0
1541.1, 26674.0
1546.0, 20279.0
1551.0, 14802.0
1556.0, 11047.0
1561.0, 7916.0
1566.0, 5655.0
1571.0, 3983.0
1576.0, 3074.0
1581.0, 2248.0
1586.0, 1832.0
1591.0, 1407.0
1596.0, 1163.0
1601.0, 1081.0
1606.0, 893.0
1611.0, 822.0
1616.0, 769.0
1621.0, 758.0
1625.9, 672.0
1630.9, 682.0
1635.9, 662.0
1640.9, 606.0
1645.9, 674.0
1650.9, 704.0
1655.9, 647.0
1660.9, 648.0
1665.9, 648.0
1670.9, 639.0
1675.9, 611.0
1680.9, 593.0
1685.9, 653.0
1690.9, 614.0
1695.9, 613.0
1700.9, 645.0
1705.9, 642.0
1710.8, 620.0
1715.8, 589.0
1720.8, 625.0
1725.8, 638.0
1730.8, 631.0
1735.8, 615.0
1740.8, 637.0
1745.8, 599.0
1750.8, 613.0
1755.8, 629.0
1760.8, 637.0
1765.8, 563.0
1770.8, 619.0
1775.8, 637.0
1780.8, 618.0
1785.8, 580.0
1790.7, 614.0
1795.7, 575.0
1800.7, 583.0
1805.7, 619.0
1810.7, 640.0
1815.7, 613.0
1820.7, 614.0
1825.7, 609.0
1830.7, 628.0
1835.7, 651.0
1840.7, 627.0
1845.7, 610.0
1850.7, 601.0
1855.7, 620.0
1860.7, 588.0
1865.7, 611.0
1870.7, 638.0
1875.6, 676.0
1880.6, 600.0
1885.6, 625.0
1890.6, 652.0
1895.6, 582.0
1900.6, 609.0
1905.6, 583.0
1910.6, 580.0
1915.6, 600.0
1920.6, 619.0
1925.6, 622.0
1930.6, 579.0
1935.6, 608.0
1940.6, 604.0
1945.6, 577.0
1950.6, 619.0
1955.5, 595.0
1960.5, 611.0
1965.5, 640.0
1970.5, 559.0
1975.5, 570.0
1980.5, 608.0
1985.5, 600.0
1990.5, 662.0
1995.5, 566.0
e
# 15 kV
-2.1, 0.0
2.9, 0.0
7.9, 0.0
12.9, 0.0
17.9, 0.0
22.9, 0.0
27.9, 0.0
32.9, 0.0
37.9, 1.0
42.9, 0.0
47.9, 0.0
52.9, 5.0
57.9, 1.0
62.8, 5.0
67.8, 1.0
72.8, 3.0
77.8, 4.0
82.8, 3.0
87.8, 4.0
92.8, 8.0
97.8, 10.0
102.8, 6.0
107.8, 15.0
112.8, 12.0
117.8, 20.0
122.8, 26.0
127.8, 31.0
132.8, 29.0
137.8, 45.0
142.8, 43.0
147.7, 48.0
152.7, 45.0
157.7, 50.0
162.7, 58.0
167.7, 71.0
172.7, 89.0
177.7, 85.0
182.7, 94.0
187.7, 81.0
192.7, 115.0
197.7, 122.0
202.7, 114.0
207.7, 142.0
212.7, 150.0
217.7, 155.0
222.7, 178.0
227.6, 179.0
232.6, 214.0
237.6, 211.0
242.6, 217.0
247.6, 259.0
252.6, 247.0
257.6, 281.0
262.6, 308.0
267.6, 298.0
272.6, 280.0
277.6, 281.0
282.6, 280.0
287.6, 264.0
292.6, 279.0
297.6, 224.0
302.6, 262.0
307.6, 230.0
312.5, 218.0
317.5, 184.0
322.5, 189.0
327.5, 202.0
332.5, 223.0
337.5, 202.0
342.5, 227.0
347.5, 251.0
352.5, 259.0
357.5, 283.0
362.5, 321.0
367.5, 319.0
372.5, 340.0
377.5, 375.0
382.5, 356.0
387.5, 373.0
392.4, 391.0
397.4, 446.0
402.4, 445.0
407.4, 427.0
412.4, 498.0
417.4, 486.0
422.4, 513.0
427.4, 529.0
432.4, 592.0
437.4, 584.0
442.4, 631.0
447.4, 748.0
452.4, 843.0
457.4, 991.0
462.4, 1403.0
467.4, 2029.0
472.3, 3075.0
477.3, 4669.0
482.3, 7465.0
487.3, 11389.0
492.3, 16682.0
497.3, 23760.0
502.3, 31848.0
507.3, 40672.0
512.3, 49659.0
517.3, 57247.0
522.3, 63518.0
527.3, 66661.0
532.3, 66979.0
537.3, 63089.0
542.3, 57439.0
547.3, 49179.0
552.3, 40540.0
557.2, 31283.0
562.2, 23272.0
567.2, 16497.0
572.2, 10936.0
577.2, 6981.0
582.2, 4290.0
587.2, 2703.0
592.2, 1666.0
597.2, 1010.0
602.2, 643.0
607.2, 495.0
612.2, 435.0
617.2, 412.0
622.2, 386.0
627.2, 407.0
632.2, 408.0
637.1, 404.0
642.1, 460.0
647.1, 465.0
652.1, 467.0
657.1, 460.0
662.1, 512.0
667.1, 537.0
672.1, 515.0
677.1, 501.0
682.1, 581.0
687.1, 552.0
692.1, 558.0
697.1, 608.0
702.1, 592.0
707.1, 570.0
712.1, 601.0
717.1, 634.0
722.0, 687.0
727.0, 669.0
732.0, 678.0
737.0, 695.0
742.0, 730.0
747.0, 760.0
752.0, 753.0
757.0, 700.0
762.0, 744.0
767.0, 748.0
772.0, 801.0
777.0, 767.0
782.0, 818.0
787.0, 804.0
792.0, 793.0
797.0, 856.0
801.9, 877.0
806.9, 888.0
811.9, 926.0
816.9, 899.0
821.9, 922.0
826.9, 937.0
831.9, 936.0
836.9, 977.0
841.9, 956.0
846.9, 983.0
851.9, 1027.0
856.9, 979.0
861.9, 1044.0
866.9, 1024.0
871.9, 1006.0
876.9, 1042.0
881.9, 1053.0
886.8, 1074.0
891.8, 1050.0
896.8, 1018.0
901.8, 1073.0
906.8, 1075.0
911.8, 1094.0
916.8, 1123.0
921.8, 1085.0
926.8, 1117.0
931.8, 1206.0
936.8, 1154.0
941.8, 1186.0
946.8, 1168.0
951.8, 1150.0
956.8, 1220.0
961.8, 1190.0
966.7, 1228.0
971.7, 1218.0
976.7, 1279.0
981.7, 1232.0
986.7, 1255.0
991.7, 1250.0
996.7, 1299.0
1001.7, 1282.0
1006.7, 1293.0
1011.7, 1346.0
1016.7, 1258.0
1021.7, 1258.0
1026.7, 1270.0
1031.7, 1337.0
1036.7, 1355.0
1041.7, 1356.0
1046.7, 1290.0
1051.6, 1405.0
1056.6, 1336.0
1061.6, 1393.0
1066.6, 1354.0
1071.6, 1361.0
1076.6, 1328.0
1081.6, 1425.0
1086.6, 1415.0
1091.6, 1402.0
1096.6, 1372.0
1101.6, 1361.0
1106.6, 1400.0
1111.6, 1340.0
1116.6, 1402.0
1121.6, 1403.0
1126.6, 1438.0
1131.5, 1427.0
1136.5, 1420.0
1141.5, 1425.0
1146.5, 1429.0
1151.5, 1427.0
1156.5, 1411.0
1161.5, 1463.0
1166.5, 1528.0
1171.5, 1438.0
1176.5, 1465.0
1181.5, 1523.0
1186.5, 1468.0
1191.5, 1391.0
1196.5, 1557.0
1201.5, 1589.0
1206.5, 1543.0
1211.5, 1499.0
1216.4, 1499.0
1221.4, 1475.0
1226.4, 1504.0
1231.4, 1506.0
1236.4, 1485.0
1241.4, 1500.0
1246.4, 1480.0
1251.4, 1508.0
1256.4, 1454.0
1261.4, 1412.0
1266.4, 1539.0
1271.4, 1551.0
1276.4, 1536.0
1281.4, 1473.0
1286.4, 1544.0
1291.4, 1528.0
1296.3, 1506.0
1301.3, 1484.0
1306.3, 1542.0
1311.3, 1566.0
1316.3, 1499.0
1321.3, 1503.0
1326.3, 1486.0
1331.3, 1504.0
1336.3, 1497.0
1341.3, 1540.0
1346.3, 1546.0
1351.3, 1529.0
1356.3, 1496.0
1361.3, 1569.0
1366.3, 1532.0
1371.3, 1594.0
1376.3, 1749.0
1381.2, 1957.0
1386.2, 2288.0
1391.2, 2969.0
1396.2, 4013.0
1401.2, 5534.0
1406.2, 7888.0
1411.2, 11400.0
1416.2, 16580.0
1421.2, 23468.0
1426.2, 32846.0
1431.2, 45272.0
1436.2, 60603.0
1441.2, 79237.0
1446.2, 100791.0
1451.2, 125218.0
1456.2, 152260.0
1461.1, 178526.0
1466.1, 205012.0
1471.1, 228274.0
1476.1, 248659.0
1481.1, 260246.0
1486.1, 268526.0
1491.1, 268942.0
1496.1, 261763.0
1501.1, 247743.0
1506.1, 228741.0
1511.1, 205033.0
1516.1, 178916.0
1521.1, 152925.0
1526.1, 126682.0
1531.1, 101882.0
1536.1, 80511.0
1541.1, 61649.0
1546.0, 46491.0
1551.0, 34283.0
1556.0, 24931.0
1561.0, 17558.0
1566.0, 12545.0
1571.0, 8890.0
1576.0, 6415.0
1581.0, 4466.0
1586.0, 3297.0
1591.0, 2394.0
1596.0, 1961.0
1601.0, 1503.0
1606.0, 1268.0
1611.0, 1124.0
1616.0, 946.0
1621.0, 869.0
1625.9, 889.0
1630.9, 798.0
1635.9, 794.0
1640.9, 716.0
1645.9, 744.0
1650.9, 726.0
1655.9, 704.0
1660.9, 686.0
1665.9, 732.0
1670.9, 737.0
1675.9, 697.0
1680.9, 711.0
1685.9, 680.0
1690.9, 661.0
1695.9, 728.0
1700.9, 707.0
1705.9, 700.0
1710.8, 724.0
1715.8, 703.0
1720.8, 700.0
1725.8, 721.0
1730.8, 719.0
1735.8, 767.0
1740.8, 699.0
1745.8, 717.0
1750.8, 698.0
1755.8, 720.0
1760.8, 695.0
1765.8, 762.0
1770.8, 712.0
1775.8, 715.0
1780.8, 757.0
1785.8, 756.0
1790.7, 761.0
1795.7, 793.0
1800.7, 747.0
1805.7, 735.0
1810.7, 787.0
1815.7, 769.0
1820.7, 727.0
1825.7, 757.0
1830.7, 803.0
1835.7, 678.0
1840.7, 747.0
1845.7, 785.0
1850.7, 774.0
1855.7, 758.0
1860.7, 760.0
1865.7, 778.0
1870.7, 741.0
1875.6, 829.0
1880.6, 767.0
1885.6, 758.0
1890.6, 752.0
1895.6, 758.0
1900.6, 732.0
1905.6, 776.0
1910.6, 783.0
1915.6, 712.0
1920.6, 721.0
1925.6, 746.0
1930.6, 731.0
1935.6, 758.0
1940.6, 782.0
1945.6, 749.0
1950.6, 749.0
1955.5, 794.0
1960.5, 753.0
1965.5, 778.0
1970.5, 728.0
1975.5, 767.0
1980.5, 782.0
1985.5, 760.0
1990.5, 735.0
1995.5, 742.0
e
# 20 kV
-2.1, 0.0
2.9, 0.0
7.9, 0.0
12.9, 0.0
17.9, 0.0
22.9, 0.0
27.9, 0.0
32.9, 0.0
37.9, 1.0
42.9, 2.0
47.9, 1.0
52.9, 5.0
57.9, 2.0
62.8, 2.0
67.8, 1.0
72.8, 5.0
77.8, 1.0
82.8, 6.0
87.8, 1.0
92.8, 6.0
97.8, 10.0
102.8, 6.0
107.8, 15.0
112.8, 12.0
117.8, 20.0
122.8, 25.0
127.8, 31.0
132.8, 30.0
137.8, 47.0
142.8, 39.0
147.7, 53.0
152.7, 54.0
157.7, 51.0
162.7, 55.0
167.7, 66.0
172.7, 79.0
177.7, 68.0
182.7, 91.0
187.7, 69.0
192.7, 82.0
197.7, 99.0
202.7, 103.0
207.7, 102.0
212.7, 122.0
217.7, 131.0
222.7, 146.0
227.6, 160.0
232.6, 178.0
237.6, 183.0
242.6, 173.0
247.6, 184.0
252.6, 212.0
257.6, 242.0
262.6, 234.0
267.6, 229.0
272.6, 243.0
277.6, 209.0
282.6, 251.0
287.6, 205.0
292.6, 207.0
297.6, 195.0
302.6, 180.0
307.6, 134.0
312.5, 166.0
317.5, 152.0
322.5, 144.0
327.5, 150.0
332.5, 154.0
337.5, 178.0
342.5, 165.0
347.5, 215.0
352.5, 216.0
357.5, 202.0
362.5, 204.0
367.5, 241.0
372.5, 265.0
377.5, 264.0
382.5, 279.0
387.5, 310.0
392.4, 333.0
397.4, 310.0
402.4, 377.0
407.4, 359.0
412.4, 352.0
417.4, 406.0
422.4, 398.0
427.4, 400.0
432.4, 445.0
437.4, 503.0
442.4, 552.0
447.4, 539.0
452.4, 627.0
457.4, 781.0
462.4, 1127.0
467.4, 1639.0
472.3, 2563.0
477.3, 3943.0
482.3, 6136.0
487.3, 9310.0
492.3, 13512.0
497.3, 19099.0
502.3, 25421.0
507.3, 32897.0
512.3, 40156.0
517.3, 46403.0
522.3, 51347.0
527.3, 53775.0
532.3, 53968.0
537.3, 51158.0
542.3, 46156.0
547.3, 39670.0
552.3, 32133.0
557.2, 25266.0
562.2, 18530.0
567.2, 13300.0
572.2, 8869.0
577.2, 5674.0
582.2, 3577.0
587.2, 2023.0
592.2, 1289.0
597.2, 787.0
602.2, 578.0
607.2, 418.0
612.2, 352.0
617.2, 356.0
622.2, 301.0
627.2, 345.0
632.2, 335.0
637.1, 366.0
642.1, 361.0
647.1, 412.0
652.1, 395.0
657.1, 385.0
662.1, 418.0
667.1, 430.0
672.1, 441.0
677.1, 478.0
682.1, 476.0
687.1, 516.0
692.1, 529.0
697.1, 487.0
702.1, 510.0
707.1, 523.0
712.1, 551.0
717.1, 528.0
722.0, 541.0
727.0, 573.0
732.0, 580.0
737.0, 597.0
742.0, 606.0
747.0, 599.0
752.0, 655.0
757.0, 655.0
762.0, 665.0
767.0, 673.0
772.0, 692.0
777.0, 677.0
782.0, 730.0
787.0, 729.0
792.0, 780.0
797.0, 762.0
801.9, 786.0
806.9, 766.0
811.9, 805.0
816.9, 856.0
821.9, 851.0
826.9, 865.0
831.9, 823.0
836.9, 865.0
841.9, 930.0
846.9, 915.0
851.9, 902.0
856.9, 991.0
861.9, 978.0
866.9, 1006.0
871.9, 1012.0
876.9, 1024.0
881.9, 1012.0
886.8, 1083.0
891.8, 1058.0
896.8, 1026.0
901.8, 1079.0
906.8, 1157.0
911.8, 1074.0
916.8, 1217.0
921.8, 1109.0
926.8, 1152.0
931.8, 1166.0
936.8, 1189.0
941.8, 1194.0
946.8, 1212.0
951.8, 1211.0
956.8, 1295.0
961.8, 1218.0
966.7, 1357.0
971.7, 1238.0
976.7, 1246.0
981.7, 1337.0
986.7, 1335.0
991.7, 1296.0
996.7, 1390.0
1001.7, 1408.0
1006.7, 1364.0
1011.7, 1347.0
1016.7, 1459.0
1021.7, 1404.0
1026.7, 1431.0
1031.7, 1476.0
1036.7, 1443.0
1041.7, 1459.0
1046.7, 1419.0
1051.6, 1476.0
1056.6, 1497.0
1061.6, 1545.0
1066.6, 1469.0
1071.6, 1439.0
1076.6, 1565.0
1081.6, 1533.0
1086.6, 1503.0
1091.6, 1566.0
1096.6, 1487.0
1101.6, 1600.0
1106.6, 1593.0
1111.6, 1586.0
1116.6, 1632.0
1121.6, 1667.0
1126.6, 1636.0
1131.5, 1642.0
1136.5, 1632.0
1141.5, 1669.0
1146.5, 1648.0
1151.5, 1732.0
1156.5, 1753.0
1161.5, 1748.0
1166.5, 1664.0
1171.5, 1730.0
1176.5, 1693.0
1181.5, 1728.0
1186.5, 1768.0
1191.5, 1737.0
1196.5, 1790.0
1201.5, 1674.0
1206.5, 1769.0
1211.5, 1691.0
1216.4, 1709.0
1221.4, 1743.0
1226.4, 1798.0
1231.4, 1810.0
1236.4, 1853.0
1241.4, 1736.0
1246.4, 1761.0
1251.4, 1762.0
1256.4, 1866.0
1261.4, 1787.0
1266.4, 1771.0
1271.4, 1766.0
1276.4, 1788.0
1281.4, 1791.0
1286.4, 1800.0
1291.4, 1803.0
1296.3, 1766.0
1301.3, 1904.0
1306.3, 1786.0
1311.3, 1807.0
1316.3, 1935.0
1321.3, 1866.0
1326.3, 1861.0
1331.3, 1869.0
1336.3, 1838.0
1341.3, 1825.0
1346.3, 1891.0
1351.3, 1865.0
1356.3, 1983.0
1361.3, 1941.0
1366.3, 1974.0
1371.3, 2128.0
1376.3, 2255.0
1381.2, 2632.0
1386.2, 3232.0
1391.2, 4271.0
1396.2, 5849.0
1401.2, 8178.0
1406.2, 12088.0
1411.2, 17600.0
1416.2, 25688.0
1421.2, 36614.0
1426.2, 51654.0
1431.2, 71169.0
1436.2, 95982.0
1441.2, 124878.0
1446.2, 160231.0
1451.2, 197982.0
1456.2, 240613.0
1461.1, 282193.0
1466.1, 324020.0
1471.1, 360673.0
1476.1, 392822.0
1481.1, 413549.0
1486.1, 426367.0
1491.1, 425221.0
1496.1, 413784.0
1501.1, 391351.0
1506.1, 362171.0
1511.1, 323524.0
1516.1, 283431.0
1521.1, 241444.0
1526.1, 200836.0
1531.1, 161489.0
1536.1, 126835.0
1541.1, 97801.0
1546.0, 73461.0
1551.0, 53744.0
1556.0, 38849.0
1561.0, 27563.0
1566.0, 19363.0
1571.0, 13697.0
1576.0, 9491.0
1581.0, 6593.0
1586.0, 4752.0
1591.0, 3495.0
1596.0, 2459.0
1601.0, 2004.0
1606.0, 1587.0
1611.0, 1252.0
1616.0, 971.0
1621.0, 930.0
1625.9, 821.0
1630.9, 703.0
1635.9, 737.0
1640.9, 679.0
1645.9, 649.0
1650.9, 643.0
1655.9, 612.0
1660.9, 594.0
1665.9, 627.0
1670.9, 630.0
1675.9, 667.0
1680.9, 654.0
1685.9, 623.0
1690.9, 618.0
1695.9, 657.0
1700.9, 641.0
1705.9, 626.0
1710.8, 658.0
1715.8, 628.0
1720.8, 642.0
1725.8, 655.0
1730.8, 637.0
1735.8, 666.0
1740.8, 697.0
1745.8, 624.0
1750.8, 627.0
1755.8, 662.0
1760.8, 646.0
1765.8, 673.0
1770.8, 682.0
1775.8, 643.0
1780.8, 651.0
1785.8, 680.0
1790.7, 679.0
1795.7, 696.0
1800.7, 622.0
1805.7, 717.0
1810.7, 655.0
1815.7, 652.0
1820.7, 698.0
1825.7, 691.0
1830.7, 662.0
1835.7, 644.0
1840.7, 734.0
1845.7, 682.0
1850.7, 705.0
1855.7, 724.0
1860.7, 726.0
1865.7, 705.0
1870.7, 690.0
1875.6, 695.0
1880.6, 725.0
1885.6, 716.0
1890.6, 706.0
1895.6, 762.0
1900.6, 732.0
1905.6, 771.0
1910.6, 722.0
1915.6, 759.0
1920.6, 764.0
1925.6, 771.0
1930.6, 729.0
1935.6, 700.0
1940.6, 719.0
1945.6, 701.0
1950.6, 783.0
1955.5, 777.0
1960.5, 754.0
1965.5, 747.0
1970.5, 751.0
1975.5, 777.0
1980.5, 747.0
1985.5, 766.0
1990.5, 832.0
1995.5, 739.0
e
# 30 kV
-2.1, 0.0
2.9, 0.0
7.9, 0.0
12.9, 0.0
17.9, 0.0
22.9, 0.0
27.9, 0.0
32.9, 0.0
37.9, 0.0
42.9, 1.0
47.9, 1.0
52.9, 0.0
57.9, 1.0
62.8, 1.0
67.8, 2.0
72.8, 2.0
77.8, 0.0
82.8, 2.0
87.8, 4.0
92.8, 3.0
97.8, 3.0
102.8, 11.0
107.8, 7.0
112.8, 16.0
117.8, 8.0
122.8, 14.0
127.8, 12.0
132.8, 13.0
137.8, 17.0
142.8, 24.0
147.7, 30.0
152.7, 35.0
157.7, 27.0
162.7, 40.0
167.7, 38.0
172.7, 40.0
177.7, 48.0
182.7, 44.0
187.7, 41.0
192.7, 61.0
197.7, 78.0
202.7, 75.0
207.7, 70.0
212.7, 73.0
217.7, 101.0
222.7, 93.0
227.6, 86.0
232.6, 95.0
237.6, 107.0
242.6, 125.0
247.6, 122.0
252.6, 132.0
257.6, 135.0
262.6, 143.0
267.6, 130.0
272.6, 136.0
277.6, 143.0
282.6, 133.0
287.6, 118.0
292.6, 127.0
297.6, 95.0
302.6, 126.0
307.6, 111.0
312.5, 86.0
317.5, 120.0
322.5, 100.0
327.5, 96.0
332.5, 97.0
337.5, 106.0
342.5, 104.0
347.5, 123.0
352.5, 107.0
357.5, 134.0
362.5, 131.0
367.5, 146.0
372.5, 156.0
377.5, 164.0
382.5, 181.0
387.5, 159.0
392.4, 181.0
397.4, 186.0
402.4, 188.0
407.4, 204.0
412.4, 195.0
417.4, 209.0
422.4, 242.0
427.4, 209.0
432.4, 271.0
437.4, 267.0
442.4, 286.0
447.4, 319.0
452.4, 387.0
457.4, 487.0
462.4, 702.0
467.4, 952.0
472.3, 1537.0
477.3, 2513.0
482.3, 3895.0
487.3, 5945.0
492.3, 8814.0
497.3, 12228.0
502.3, 16269.0
507.3, 21290.0
512.3, 25696.0
517.3, 29455.0
522.3, 32882.0
527.3, 34939.0
532.3, 34369.0
537.3, 33215.0
542.3, 29992.0
547.3, 25638.0
552.3, 20891.0
557.2, 16175.0
562.2, 11972.0
567.2, 8362.0
572.2, 5579.0
577.2, 3722.0
582.2, 2259.0
587.2, 1414.0
592.2, 817.0
597.2, 541.0
602.2, 364.0
607.2, 281.0
612.2, 224.0
617.2, 212.0
622.2, 193.0
627.2, 207.0
632.2, 216.0
637.1, 241.0
642.1, 231.0
647.1, 251.0
652.1, 252.0
657.1, 232.0
662.1, 225.0
667.1, 269.0
672.1, 294.0
677.1, 280.0
682.1, 277.0
687.1, 317.0
692.1, 339.0
697.1, 328.0
702.1, 314.0
707.1, 356.0
712.1, 337.0
717.1, 334.0
722.0, 357.0
727.0, 407.0
732.0, 384.0
737.0, 427.0
742.0, 408.0
747.0, 456.0
752.0, 455.0
757.0, 453.0
762.0, 458.0
767.0, 488.0
772.0, 498.0
777.0, 534.0
782.0, 516.0
787.0, 529.0
792.0, 539.0
797.0, 567.0
801.9, 533.0
806.9, 550.0
811.9, 661.0
816.9, 591.0
821.9, 625.0
826.9, 596.0
831.9, 642.0
836.9, 680.0
841.9, 701.0
846.9, 734.0
851.9, 713.0
856.9, 732.0
861.9, 757.0
866.9, 755.0
871.9, 739.0
876.9, 790.0
881.9, 781.0
886.8, 807.0
891.8, 824.0
896.8, 865.0
901.8, 817.0
906.8, 899.0
911.8, 876.0
916.8, 1006.0
921.8, 915.0
926.8, 921.0
931.8, 953.0
936.8, 971.0
941.8, 958.0
946.8, 1002.0
951.8, 1066.0
956.8, 1041.0
961.8, 1042.0
966.7, 1146.0
971.7, 1154.0
976.7, 1106.0
981.7, 1105.0
986.7, 1150.0
991.7, 1148.0
996.7, 1188.0
1001.7, 1251.0
1006.7, 1298.0
1011.7, 1249.0
1016.7, 1294.0
1021.7, 1267.0
1026.7, 1350.0
1031.7, 1223.0
1036.7, 1315.0
1041.7, 1388.0
1046.7, 1328.0
1051.6, 1429.0
1056.6, 1438.0
1061.6, 1455.0
1066.6, 1422.0
1071.6, 1524.0
1076.6, 1491.0
1081.6, 1487.0
1086.6, 1541.0
1091.6, 1597.0
1096.6, 1575.0
1101.6, 1564.0
1106.6, 1585.0
1111.6, 1644.0
1116.6, 1589.0
1121.6, 1527.0
1126.6, 1638.0
1131.5, 1656.0
1136.5, 1656.0
1141.5, 1694.0
1146.5, 1664.0
1151.5, 1756.0
1156.5, 1811.0
1161.5, 1735.0
1166.5, 1805.0
1171.5, 1785.0
1176.5, 1835.0
1181.5, 1924.0
1186.5, 1832.0
1191.5, 1887.0
1196.5, 1867.0
1201.5, 1818.0
1206.5, 1963.0
1211.5, 1873.0
1216.4, 1962.0
1221.4, 1908.0
1226.4, 1969.0
1231.4, 1927.0
1236.4, 1859.0
1241.4, 1920.0
1246.4, 1986.0
1251.4, 2107.0
1256.4, 1993.0
1261.4, 1992.0
1266.4, 2037.0
1271.4, 2018.0
1276.4, 2021.0
1281.4, 2050.0
1286.4, 2057.0
1291.4, 2060.0
1296.3, 2159.0
1301.3, 2167.0
1306.3, 2127.0
1311.3, 2173.0
1316.3, 2119.0
1321.3, 2153.0
1326.3, 2224.0
1331.3, 2274.0
1336.3, 2161.0
1341.3, 2191.0
1346.3, 2222.0
1351.3, 2183.0
1356.3, 2180.0
1361.3, 2315.0
1366.3, 2358.0
1371.3, 2570.0
1376.3, 2901.0
1381.2, 3514.0
1386.2, 4496.0
1391.2, 5963.0
1396.2, 8402.0
1401.2, 12375.0
1406.2, 18265.0
1411.2, 26935.0
1416.2, 39754.0
1421.2, 56870.0
1426.2, 80893.0
1431.2, 111037.0
1436.2, 148744.0
1441.2, 196526.0
1446.2, 250443.0
1451.2, 310041.0
1456.2, 374907.0
1461.1, 442356.0
1466.1, 507554.0
1471.1, 565021.0
1476.1, 614363.0
1481.1, 648842.0
1486.1, 668050.0
1491.1, 665903.0
1496.1, 649175.0
1501.1, 615372.0
1506.1, 565464.0
1511.1, 508656.0
1516.1, 445192.0
1521.1, 378410.0
1526.1, 314226.0
1531.1, 252782.0
1536.1, 199023.0
1541.1, 152435.0
1546.0, 114660.0
1551.0, 84530.0
1556.0, 61023.0
1561.0, 43154.0
1566.0, 30274.0
1571.0, 21133.0
1576.0, 14394.0
1581.0, 9895.0
1586.0, 6929.0
1591.0, 4933.0
1596.0, 3493.0
1601.0, 2548.0
1606.0, 1955.0
1611.0, 1378.0
1616.0, 1037.0
1621.0, 817.0
1625.9, 664.0
1630.9, 586.0
1635.9, 524.0
1640.9, 462.0
1645.9, 458.0
1650.9, 405.0
1655.9, 392.0
1660.9, 403.0
1665.9, 403.0
1670.9, 400.0
1675.9, 440.0
1680.9, 441.0
1685.9, 395.0
1690.9, 418.0
1695.9, 452.0
1700.9, 390.0
1705.9, 390.0
1710.8, 455.0
1715.8, 451.0
1720.8, 406.0
1725.8, 444.0
1730.8, 418.0
1735.8, 473.0
1740.8, 437.0
1745.8, 398.0
1750.8, 452.0
1755.8, 419.0
1760.8, 424.0
1765.8, 423.0
1770.8, 435.0
1775.8, 448.0
1780.8, 447.0
1785.8, 435.0
1790.7, 452.0
1795.7, 417.0
1800.7, 472.0
1805.7, 455.0
1810.7, 455.0
1815.7, 496.0
1820.7, 504.0
1825.7, 502.0
1830.7, 478.0
1835.7, 448.0
1840.7, 490.0
1845.7, 517.0
1850.7, 486.0
1855.7, 476.0
1860.7, 482.0
1865.7, 478.0
1870.7, 502.0
1875.6, 489.0
1880.6, 478.0
1885.6, 482.0
1890.6, 502.0
1895.6, 509.0
1900.6, 498.0
1905.6, 523.0
1910.6, 488.0
1915.6, 486.0
1920.6, 513.0
1925.6, 536.0
1930.6, 534.0
1935.6, 525.0
1940.6, 536.0
1945.6, 542.0
1950.6, 561.0
1955.5, 547.0
1960.5, 562.0
1965.5, 557.0
1970.5, 511.0
1975.5, 564.0
1980.5, 514.0
1985.5, 581.0
1990.5, 533.0
1995.5, 500.0
e
# O K-L3
524.9, 82450.0
# Al K-L3
1486.5, 668050.0
e
