# Generated by the NIST EPQ Library
# On Jun 21, 2017
# By l837410
# NOTE: TO SAVE FILES AS ENCAPSULATED POST SCRIPT FILES, OR AS OTHER FILE FORMATS,
#     COMMENT OUT THE "TERMINAL WINDOW" SECTION, SET THE TERMINAL AS
#     THE CHOSEN TYPE, AND SET THE OUTPUT.
# set terminal postscript enhanced "Times-Roman" 12
# set terminal png
set terminal window
unset logscale xy
set xrange [0.0:1000.0]
set yrange [0.0:260.0]
set xlabel "(Energy) eV" font "Arial,16"
set ylabel "" font "Arial,16"
set format x "%g"
set format y ""
unset label
# C K-L2
set label 1 "C K-L2" at first 277.4,3425.7 center
# O K-L2
set label 2 "O K-L2" at first 524.9,98.6 center
# O K-L3
set label 3 "O K-L3" at first 524.9,191.9 center
# C K-L3
set label 4 "C K-L3" at first 281.8,6846.2 center
plot "-" with lines  title "7kV-5mm-S3-HR-EHT-200s", "-" with lines  title "7kV-5mm-S4-HR-EHT-200s", "-" with lines  title "7kV-5mm-S5-HR-EHT-200s", "-" with lines  title "7kV-5mm-S6-HR-EHT-200s", "-" with impulses lt -1 notitle
# 7kV-5mm-S3-HR-EHT-200s
-0.9, 2691.0
4.1, 3056.0
9.1, 3378.0
14.1, 3676.0
19.1, 4066.0
24.1, 4276.0
29.1, 4579.0
34.1, 4847.0
39.1, 4957.0
44.1, 4915.0
49.1, 4618.0
54.1, 4113.0
59.1, 3535.0
64.1, 2867.0
69.1, 2125.0
74.1, 1531.0
79.1, 982.0
84.1, 613.0
89.1, 420.0
94.1, 260.0
99.1, 145.0
104.1, 106.0
109.1, 79.0
114.1, 75.0
119.2, 52.0
124.2, 64.0
129.2, 67.0
134.2, 61.0
139.2, 58.0
144.2, 66.0
149.2, 66.0
154.2, 60.0
159.2, 67.0
164.2, 71.0
169.2, 89.0
174.2, 81.0
179.2, 93.0
184.2, 88.0
189.2, 115.0
194.2, 147.0
199.2, 184.0
204.2, 240.0
209.2, 324.0
214.2, 445.0
219.2, 667.0
224.2, 962.0
229.2, 1308.0
234.2, 1808.0
239.2, 2439.0
244.2, 3145.0
249.2, 3958.0
254.2, 4766.0
259.2, 5680.0
264.2, 6238.0
269.2, 6773.0
274.3, 6841.0
279.3, 6734.0
284.3, 6531.0
289.3, 6046.0
294.3, 5153.0
299.3, 4306.0
304.3, 3445.0
309.3, 2642.0
314.3, 1963.0
319.3, 1374.0
324.3, 960.0
329.3, 616.0
334.3, 385.0
339.3, 229.0
344.3, 134.0
349.3, 77.0
354.3, 41.0
359.3, 24.0
364.3, 15.0
369.3, 9.0
374.3, 7.0
379.3, 3.0
384.3, 7.0
389.3, 12.0
394.3, 9.0
399.3, 9.0
404.3, 13.0
409.3, 8.0
414.3, 12.0
419.3, 10.0
424.3, 7.0
429.4, 10.0
434.4, 8.0
439.4, 13.0
444.4, 16.0
449.4, 19.0
454.4, 18.0
459.4, 21.0
464.4, 23.0
469.4, 21.0
474.4, 32.0
479.4, 42.0
484.4, 36.0
489.4, 47.0
494.4, 58.0
499.4, 60.0
504.4, 68.0
509.4, 63.0
514.4, 59.0
519.4, 76.0
524.4, 83.0
529.4, 69.0
534.4, 72.0
539.4, 67.0
544.4, 57.0
549.4, 46.0
554.4, 41.0
559.4, 38.0
564.4, 32.0
569.4, 28.0
574.4, 21.0
579.4, 30.0
584.4, 21.0
589.5, 13.0
594.5, 20.0
599.5, 21.0
604.5, 23.0
609.5, 22.0
614.5, 14.0
619.5, 19.0
624.5, 13.0
629.5, 12.0
634.5, 19.0
639.5, 19.0
644.5, 19.0
649.5, 23.0
654.5, 28.0
659.5, 18.0
664.5, 18.0
669.5, 19.0
674.5, 20.0
679.5, 21.0
684.5, 27.0
689.5, 20.0
694.5, 26.0
699.5, 28.0
704.5, 20.0
709.5, 21.0
714.5, 18.0
719.5, 20.0
724.5, 26.0
729.5, 29.0
734.5, 24.0
739.5, 23.0
744.6, 24.0
749.6, 27.0
754.6, 25.0
759.6, 30.0
764.6, 23.0
769.6, 20.0
774.6, 19.0
779.6, 24.0
784.6, 23.0
789.6, 22.0
794.6, 27.0
799.6, 28.0
804.6, 33.0
809.6, 26.0
814.6, 29.0
819.6, 29.0
824.6, 26.0
829.6, 21.0
834.6, 24.0
839.6, 35.0
844.6, 30.0
849.6, 24.0
854.6, 26.0
859.6, 31.0
864.6, 26.0
869.6, 23.0
874.6, 30.0
879.6, 32.0
884.6, 26.0
889.6, 26.0
894.6, 28.0
899.7, 29.0
904.7, 27.0
909.7, 23.0
914.7, 21.0
919.7, 25.0
924.7, 24.0
929.7, 23.0
934.7, 31.0
939.7, 32.0
944.7, 29.0
949.7, 29.0
954.7, 18.0
959.7, 27.0
964.7, 26.0
969.7, 31.0
974.7, 38.0
979.7, 25.0
984.7, 29.0
989.7, 28.0
994.7, 30.0
e
# 0.2191995898618988?7kV-5mm-S4-HR-EHT-200s
-0.9, 610.9
4.1, 685.0
9.1, 759.5
14.1, 847.2
19.1, 882.9
24.1, 955.3
29.1, 1017.7
34.1, 1064.4
39.1, 1104.1
44.1, 1085.0
49.1, 1011.8
54.1, 923.5
59.1, 796.8
64.1, 659.6
69.1, 505.9
74.1, 364.5
79.1, 263.9
84.1, 187.4
89.1, 134.2
94.1, 99.7
99.1, 76.5
104.1, 68.4
109.1, 61.6
114.1, 59.6
119.2, 58.7
124.2, 58.3
129.2, 55.5
134.2, 58.5
139.2, 63.3
144.2, 58.7
149.2, 60.1
154.2, 66.0
159.2, 73.7
164.2, 81.8
169.2, 82.0
174.2, 88.1
179.2, 98.0
184.2, 101.1
189.2, 127.8
194.2, 157.8
199.2, 195.3
204.2, 252.5
209.2, 350.5
214.2, 480.3
219.2, 659.8
224.2, 955.5
229.2, 1330.5
234.2, 1820.9
239.2, 2421.9
244.2, 3167.4
249.2, 3946.9
254.2, 4781.4
259.2, 5579.3
264.2, 6249.2
269.2, 6671.1
274.3, 6841.0
279.3, 6734.9
284.3, 6365.1
289.3, 5789.3
294.3, 5058.5
299.3, 4194.4
304.3, 3332.3
309.3, 2555.6
314.3, 1907.9
319.3, 1332.3
324.3, 890.2
329.3, 586.8
334.3, 359.0
339.3, 219.4
344.3, 128.7
349.3, 72.1
354.3, 40.1
359.3, 23.5
364.3, 14.2
369.3, 10.7
374.3, 8.1
379.3, 8.1
384.3, 8.5
389.3, 8.8
394.3, 10.5
399.3, 10.1
404.3, 10.7
409.3, 11.2
414.3, 10.7
419.3, 11.0
424.3, 13.6
429.4, 14.0
434.4, 12.9
439.4, 13.2
444.4, 15.6
449.4, 16.7
454.4, 16.4
459.4, 18.9
464.4, 19.1
469.4, 21.0
474.4, 27.0
479.4, 30.7
484.4, 33.8
489.4, 41.6
494.4, 48.9
499.4, 54.6
504.4, 59.2
509.4, 65.3
514.4, 73.2
519.4, 70.8
524.4, 82.4
529.4, 78.0
534.4, 73.9
539.4, 68.8
544.4, 67.7
549.4, 58.7
554.4, 50.4
559.4, 45.8
564.4, 41.0
569.4, 38.1
574.4, 31.6
579.4, 24.8
584.4, 21.7
589.5, 20.8
594.5, 25.6
599.5, 22.8
604.5, 19.1
609.5, 18.9
614.5, 21.0
619.5, 21.0
624.5, 20.6
629.5, 19.1
634.5, 20.2
639.5, 21.5
644.5, 20.4
649.5, 24.3
654.5, 23.0
659.5, 21.0
664.5, 19.5
669.5, 18.9
674.5, 21.3
679.5, 21.3
684.5, 25.2
689.5, 26.7
694.5, 23.7
699.5, 23.7
704.5, 25.0
709.5, 26.7
714.5, 25.2
719.5, 28.9
724.5, 25.2
729.5, 24.1
734.5, 27.8
739.5, 27.0
744.6, 27.0
749.6, 25.0
754.6, 27.4
759.6, 22.4
764.6, 24.6
769.6, 24.8
774.6, 26.5
779.6, 26.3
784.6, 23.9
789.6, 22.8
794.6, 23.7
799.6, 24.6
804.6, 25.4
809.6, 27.2
814.6, 27.4
819.6, 25.6
824.6, 27.2
829.6, 26.3
834.6, 26.1
839.6, 29.4
844.6, 28.7
849.6, 26.5
854.6, 24.3
859.6, 25.9
864.6, 26.1
869.6, 28.1
874.6, 27.4
879.6, 30.5
884.6, 28.9
889.6, 27.2
894.6, 27.2
899.7, 27.4
904.7, 29.4
909.7, 30.0
914.7, 28.9
919.7, 30.0
924.7, 28.3
929.7, 30.0
934.7, 26.7
939.7, 24.3
944.7, 28.1
949.7, 28.9
954.7, 26.1
959.7, 27.0
964.7, 29.8
969.7, 28.1
974.7, 26.1
979.7, 22.1
984.7, 24.6
989.7, 26.5
994.7, 28.7
e
# 0.05318272280614466?7kV-5mm-S5-HR-EHT-200s
-0.9, 143.3
4.1, 165.5
9.1, 181.9
14.1, 198.7
19.1, 213.9
24.1, 219.6
29.1, 233.4
34.1, 252.9
39.1, 256.9
44.1, 254.8
49.1, 241.4
54.1, 222.6
59.1, 201.5
64.1, 170.8
69.1, 140.3
74.1, 109.6
79.1, 92.4
84.1, 81.8
89.1, 72.8
94.1, 65.3
99.1, 63.1
104.1, 60.5
109.1, 59.4
114.1, 60.3
119.2, 58.7
124.2, 59.7
129.2, 59.2
134.2, 60.3
139.2, 59.7
144.2, 62.9
149.2, 63.4
154.2, 67.4
159.2, 68.7
164.2, 74.6
169.2, 81.5
174.2, 87.8
179.2, 93.6
184.2, 108.0
189.2, 125.5
194.2, 147.8
199.2, 184.5
204.2, 241.5
209.2, 333.7
214.2, 473.8
219.2, 659.7
224.2, 936.4
229.2, 1302.2
234.2, 1775.0
239.2, 2385.5
244.2, 3109.5
249.2, 3891.6
254.2, 4721.8
259.2, 5519.5
264.2, 6179.7
269.2, 6652.4
274.3, 6841.0
279.3, 6775.1
284.3, 6437.9
289.3, 5874.6
294.3, 5124.3
299.3, 4269.9
304.3, 3428.5
309.3, 2637.9
314.3, 1932.3
319.3, 1363.2
324.3, 926.8
329.3, 604.0
334.3, 375.3
339.3, 226.8
344.3, 133.8
349.3, 79.3
354.3, 46.7
359.3, 28.9
364.3, 19.3
369.3, 13.8
374.3, 12.0
379.3, 10.8
384.3, 10.7
389.3, 11.2
394.3, 12.0
399.3, 10.7
404.3, 12.3
409.3, 12.3
414.3, 12.8
419.3, 14.4
424.3, 13.6
429.4, 13.8
434.4, 14.4
439.4, 16.0
444.4, 15.7
449.4, 17.5
454.4, 19.8
459.4, 22.4
464.4, 25.0
469.4, 27.9
474.4, 33.0
479.4, 35.6
484.4, 41.2
489.4, 49.9
494.4, 58.9
499.4, 69.0
504.4, 76.0
509.4, 85.7
514.4, 90.8
519.4, 98.9
524.4, 108.3
529.4, 112.9
534.4, 114.9
539.4, 112.2
544.4, 109.0
549.4, 103.2
554.4, 97.6
559.4, 88.8
564.4, 78.6
569.4, 67.8
574.4, 57.2
579.4, 50.2
584.4, 43.9
589.5, 37.4
594.5, 31.0
599.5, 27.8
604.5, 26.6
609.5, 22.6
614.5, 22.0
619.5, 21.5
624.5, 21.0
629.5, 20.2
634.5, 20.8
639.5, 21.6
644.5, 21.8
649.5, 22.0
654.5, 22.1
659.5, 21.6
664.5, 21.9
669.5, 22.2
674.5, 22.2
679.5, 21.4
684.5, 24.6
689.5, 23.4
694.5, 22.9
699.5, 26.0
704.5, 25.1
709.5, 24.7
714.5, 24.7
719.5, 25.3
724.5, 24.8
729.5, 25.5
734.5, 25.3
739.5, 26.9
744.6, 27.8
749.6, 26.0
754.6, 26.4
759.6, 25.6
764.6, 25.8
769.6, 26.3
774.6, 26.4
779.6, 25.5
784.6, 27.0
789.6, 28.5
794.6, 26.8
799.6, 26.2
804.6, 26.9
809.6, 27.9
814.6, 27.2
819.6, 28.6
824.6, 28.7
829.6, 27.9
834.6, 27.3
839.6, 28.8
844.6, 28.7
849.6, 26.6
854.6, 26.5
859.6, 27.1
864.6, 26.8
869.6, 26.5
874.6, 26.8
879.6, 27.7
884.6, 28.3
889.6, 29.0
894.6, 28.8
899.7, 28.0
904.7, 26.2
909.7, 26.8
914.7, 27.8
919.7, 27.1
924.7, 27.3
929.7, 28.9
934.7, 29.6
939.7, 27.6
944.7, 26.9
949.7, 28.6
954.7, 29.2
959.7, 27.2
964.7, 27.8
969.7, 28.1
974.7, 28.6
979.7, 29.3
984.7, 29.1
989.7, 29.4
994.7, 29.4
e
# 0.023137124922211098?7kV-5mm-S6-HR-EHT-200s
-0.9, 60.8
4.1, 70.9
9.1, 78.0
14.1, 85.3
19.1, 89.4
24.1, 92.9
29.1, 98.7
34.1, 102.7
39.1, 103.8
44.1, 103.5
49.1, 101.0
54.1, 95.1
59.1, 86.3
64.1, 78.4
69.1, 72.5
74.1, 65.6
79.1, 60.9
84.1, 59.5
89.1, 57.7
94.1, 57.1
99.1, 57.8
104.1, 57.0
109.1, 58.1
114.1, 57.6
119.2, 56.6
124.2, 55.3
129.2, 55.6
134.2, 56.9
139.2, 59.9
144.2, 61.8
149.2, 65.1
154.2, 68.0
159.2, 70.3
164.2, 74.6
169.2, 78.6
174.2, 85.0
179.2, 93.9
184.2, 105.6
189.2, 121.6
194.2, 149.0
199.2, 186.5
204.2, 243.0
209.2, 331.1
214.2, 460.5
219.2, 655.7
224.2, 921.0
229.2, 1290.4
234.2, 1783.6
239.2, 2380.4
244.2, 3107.5
249.2, 3886.4
254.2, 4705.1
259.2, 5502.0
264.2, 6158.5
269.2, 6606.9
274.3, 6841.0
279.3, 6772.9
284.3, 6414.2
289.3, 5839.0
294.3, 5098.8
299.3, 4268.6
304.3, 3420.8
309.3, 2621.0
314.3, 1925.9
319.3, 1366.1
324.3, 929.3
329.3, 608.3
334.3, 384.6
339.3, 230.6
344.3, 135.9
349.3, 81.3
354.3, 49.4
359.3, 32.1
364.3, 21.7
369.3, 16.4
374.3, 14.2
379.3, 13.9
384.3, 12.8
389.3, 12.8
394.3, 13.0
399.3, 13.1
404.3, 14.1
409.3, 14.5
414.3, 14.1
419.3, 15.5
424.3, 15.1
429.4, 14.9
434.4, 16.8
439.4, 19.1
444.4, 19.6
449.4, 21.4
454.4, 24.1
459.4, 26.8
464.4, 30.3
469.4, 34.5
474.4, 40.9
479.4, 46.5
484.4, 53.8
489.4, 64.1
494.4, 77.2
499.4, 92.0
504.4, 107.4
509.4, 123.2
514.4, 139.9
519.4, 151.5
524.4, 163.3
529.4, 174.1
534.4, 181.3
539.4, 186.7
544.4, 183.5
549.4, 178.3
554.4, 167.7
559.4, 156.0
564.4, 140.0
569.4, 121.8
574.4, 104.2
579.4, 88.7
584.4, 74.2
589.5, 59.8
594.5, 48.7
599.5, 40.5
604.5, 33.9
609.5, 29.2
614.5, 26.4
619.5, 24.9
624.5, 21.9
629.5, 22.2
634.5, 22.2
639.5, 21.6
644.5, 21.9
649.5, 22.0
654.5, 22.0
659.5, 22.7
664.5, 24.0
669.5, 23.9
674.5, 24.1
679.5, 25.2
684.5, 25.2
689.5, 24.9
694.5, 26.0
699.5, 25.5
704.5, 24.6
709.5, 26.1
714.5, 25.5
719.5, 26.6
724.5, 27.0
729.5, 27.0
734.5, 26.7
739.5, 27.8
744.6, 27.5
749.6, 27.8
754.6, 28.4
759.6, 28.0
764.6, 27.2
769.6, 29.2
774.6, 29.8
779.6, 29.2
784.6, 29.7
789.6, 29.5
794.6, 28.9
799.6, 29.3
804.6, 29.6
809.6, 29.5
814.6, 28.6
819.6, 28.5
824.6, 29.3
829.6, 30.3
834.6, 29.2
839.6, 29.3
844.6, 28.6
849.6, 28.1
854.6, 28.6
859.6, 28.6
864.6, 28.3
869.6, 28.0
874.6, 27.9
879.6, 29.3
884.6, 28.1
889.6, 27.5
894.6, 28.4
899.7, 28.1
904.7, 28.8
909.7, 28.3
914.7, 29.2
919.7, 28.6
924.7, 28.4
929.7, 28.1
934.7, 28.4
939.7, 28.5
944.7, 27.9
949.7, 27.9
954.7, 27.9
959.7, 27.5
964.7, 27.2
969.7, 27.3
974.7, 28.3
979.7, 28.9
984.7, 29.0
989.7, 28.4
994.7, 28.9
e
# C K-L2
277.4, 3420.5
# O K-L2
524.9, 93.4
# O K-L3
524.9, 186.7
# C K-L3
281.8, 6841.0
e
