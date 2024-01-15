%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%   initfilterDaub_DWT_discrete.m                       %
%                                                       %
%        D. Veitch   P.Abry                             %
%                                                       %
%   Melbourne  15/1/99                                  %
%   Modified DV May 99 (trun'd to 0.001 and symmetrized)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    This function provides the stored filters I(m) for Daubechies wavelets with N=1,2,..., 10,
%    calculated using calc_initfilterDWT_discrete.m . They are used when using the DWT 
%    to analyse intrinsically discrete time series X(k).
%
%
%  Input:   regu:  (regularity) =  number of vanishing moments (of the Daubechies wavelet).
%
%  Output:  I: The filter coefficents
%           nI: The index into the Matlab vector I corresponding to  the element m=0.
%
%  The output is the same as that of calc_initfilterDWT_discrete.m 


function [I,nI] = initfilterDaub_DWT_discrete(regu);

if ~(regu<11 & regu>0),
	mess=[' Ondelette non implantee'];
	disp(mess);
	I = []; nI=0; return;
else
% symmetric filters: being truncations to 0.00100 of original output of calc_initfilterDWT_discrete(N,0.0001) with nbvoies3 = 14
if regu==1
  I = [ -0.001006
 0.001011
 -0.001016
 0.001021
 -0.001026
 0.001031
 -0.001037
 0.001042
 -0.001047
 0.001053
 -0.001058
 0.001064
 -0.001069
 0.001075
 -0.001081
 0.001087
 -0.001092
 0.001098
 -0.001104
 0.001110
 -0.001116
 0.001123
 -0.001129
 0.001135
 -0.001142
 0.001148
 -0.001155
 0.001161
 -0.001168
 0.001175
 -0.001182
 0.001189
 -0.001196
 0.001203
 -0.001210
 0.001217
 -0.001224
 0.001232
 -0.001239
 0.001247
 -0.001255
 0.001263
 -0.001270
 0.001279
 -0.001287
 0.001295
 -0.001303
 0.001312
 -0.001320
 0.001329
 -0.001338
 0.001346
 -0.001355
 0.001365
 -0.001374
 0.001383
 -0.001393
 0.001402
 -0.001412
 0.001422
 -0.001432
 0.001442
 -0.001453
 0.001463
 -0.001474
 0.001485
 -0.001496
 0.001507
 -0.001518
 0.001529
 -0.001541
 0.001553
 -0.001565
 0.001577
 -0.001589
 0.001602
 -0.001615
 0.001628
 -0.001641
 0.001654
 -0.001668
 0.001682
 -0.001696
 0.001710
 -0.001725
 0.001739
 -0.001754
 0.001770
 -0.001785
 0.001801
 -0.001817
 0.001834
 -0.001851
 0.001868
 -0.001885
 0.001903
 -0.001921
 0.001939
 -0.001958
 0.001977
 -0.001996
 0.002016
 -0.002037
 0.002057
 -0.002078
 0.002100
 -0.002122
 0.002144
 -0.002167
 0.002191
 -0.002215
 0.002239
 -0.002264
 0.002290
 -0.002316
 0.002343
 -0.002370
 0.002398
 -0.002427
 0.002456
 -0.002486
 0.002517
 -0.002549
 0.002581
 -0.002615
 0.002649
 -0.002684
 0.002720
 -0.002757
 0.002795
 -0.002834
 0.002874
 -0.002916
 0.002958
 -0.003002
 0.003047
 -0.003094
 0.003142
 -0.003191
 0.003242
 -0.003295
 0.003350
 -0.003406
 0.003464
 -0.003524
 0.003587
 -0.003651
 0.003718
 -0.003788
 0.003860
 -0.003935
 0.004013
 -0.004094
 0.004178
 -0.004266
 0.004358
 -0.004454
 0.004554
 -0.004659
 0.004768
 -0.004883
 0.005004
 -0.005130
 0.005264
 -0.005404
 0.005552
 -0.005708
 0.005874
 -0.006049
 0.006235
 -0.006433
 0.006644
 -0.006870
 0.007111
 -0.007369
 0.007647
 -0.007947
 0.008272
 -0.008624
 0.009007
 -0.009426
 0.009886
 -0.010393
 0.010955
 -0.011581
 0.012283
 -0.013076
 0.013978
 -0.015014
 0.016216
 -0.017627
 0.019308
 -0.021342
 0.023856
 -0.027042
 0.031211
 -0.036902
 0.045137
 -0.058124
 0.081682
 -0.138078
 0.589459
 0.589520
 -0.138078
 0.081682
 -0.058124
 0.045137
 -0.036902
 0.031211
 -0.027042
 0.023856
 -0.021342
 0.019308
 -0.017627
 0.016216
 -0.015014
 0.013978
 -0.013076
 0.012283
 -0.011581
 0.010955
 -0.010393
 0.009886
 -0.009426
 0.009007
 -0.008624
 0.008272
 -0.007947
 0.007647
 -0.007369
 0.007111
 -0.006870
 0.006644
 -0.006433
 0.006235
 -0.006049
 0.005874
 -0.005708
 0.005552
 -0.005404
 0.005264
 -0.005130
 0.005004
 -0.004883
 0.004768
 -0.004659
 0.004554
 -0.004454
 0.004358
 -0.004266
 0.004178
 -0.004094
 0.004013
 -0.003935
 0.003860
 -0.003788
 0.003718
 -0.003651
 0.003587
 -0.003524
 0.003464
 -0.003406
 0.003350
 -0.003295
 0.003242
 -0.003191
 0.003142
 -0.003094
 0.003047
 -0.003002
 0.002958
 -0.002916
 0.002874
 -0.002834
 0.002795
 -0.002757
 0.002720
 -0.002684
 0.002649
 -0.002615
 0.002581
 -0.002549
 0.002517
 -0.002486
 0.002456
 -0.002427
 0.002398
 -0.002370
 0.002343
 -0.002316
 0.002290
 -0.002264
 0.002239
 -0.002215
 0.002191
 -0.002167
 0.002144
 -0.002122
 0.002100
 -0.002078
 0.002057
 -0.002037
 0.002016
 -0.001996
 0.001977
 -0.001958
 0.001939
 -0.001921
 0.001903
 -0.001885
 0.001868
 -0.001851
 0.001834
 -0.001817
 0.001801
 -0.001785
 0.001770
 -0.001754
 0.001739
 -0.001725
 0.001710
 -0.001696
 0.001682
 -0.001668
 0.001654
 -0.001641
 0.001628
 -0.001615
 0.001602
 -0.001589
 0.001577
 -0.001565
 0.001553
 -0.001541
 0.001529
 -0.001518
 0.001507
 -0.001496
 0.001485
 -0.001474
 0.001463
 -0.001453
 0.001442
 -0.001432
 0.001422
 -0.001412
 0.001402
 -0.001393
 0.001383
 -0.001374
 0.001365
 -0.001355
 0.001346
 -0.001338
 0.001329
 -0.001320
 0.001312
 -0.001303
 0.001295
 -0.001287
 0.001279
 -0.001270
 0.001263
 -0.001255
 0.001247
 -0.001239
 0.001232
 -0.001224
 0.001217
 -0.001210
 0.001203
 -0.001196
 0.001189
 -0.001182
 0.001175
 -0.001168
 0.001161
 -0.001155
 0.001148
 -0.001142
 0.001135
 -0.001129
 0.001123
 -0.001116
 0.001110
 -0.001104
 0.001098
 -0.001092
 0.001087
 -0.001081
 0.001075
 -0.001069
 0.001064
 -0.001058
 0.001053
 -0.001047
 0.001042
 -0.001037
 0.001031
 -0.001026
 0.001021
 -0.001016
 0.001011
 -0.001006
 0.001001
 ];
  nI = 203;
end
if regu==2
  I = [ 0.001007
 -0.001014
 0.001021
 -0.001028
 0.001035
 -0.001042
 0.001049
 -0.001057
 0.001064
 -0.001072
 0.001079
 -0.001087
 0.001095
 -0.001103
 0.001111
 -0.001119
 0.001127
 -0.001136
 0.001144
 -0.001153
 0.001162
 -0.001171
 0.001180
 -0.001189
 0.001199
 -0.001208
 0.001218
 -0.001228
 0.001238
 -0.001248
 0.001258
 -0.001269
 0.001280
 -0.001290
 0.001301
 -0.001313
 0.001324
 -0.001336
 0.001348
 -0.001360
 0.001372
 -0.001385
 0.001397
 -0.001410
 0.001424
 -0.001437
 0.001451
 -0.001465
 0.001479
 -0.001494
 0.001509
 -0.001524
 0.001539
 -0.001555
 0.001571
 -0.001588
 0.001604
 -0.001622
 0.001639
 -0.001657
 0.001675
 -0.001694
 0.001713
 -0.001733
 0.001753
 -0.001773
 0.001794
 -0.001816
 0.001838
 -0.001860
 0.001883
 -0.001907
 0.001931
 -0.001956
 0.001982
 -0.002008
 0.002035
 -0.002063
 0.002091
 -0.002120
 0.002150
 -0.002181
 0.002213
 -0.002246
 0.002280
 -0.002315
 0.002350
 -0.002387
 0.002426
 -0.002465
 0.002506
 -0.002548
 0.002591
 -0.002637
 0.002683
 -0.002732
 0.002782
 -0.002834
 0.002888
 -0.002944
 0.003002
 -0.003063
 0.003126
 -0.003191
 0.003260
 -0.003332
 0.003406
 -0.003485
 0.003567
 -0.003653
 0.003743
 -0.003838
 0.003937
 -0.004042
 0.004153
 -0.004270
 0.004393
 -0.004525
 0.004664
 -0.004812
 0.004969
 -0.005138
 0.005318
 -0.005511
 0.005719
 -0.005943
 0.006186
 -0.006449
 0.006735
 -0.007048
 0.007392
 -0.007771
 0.008190
 -0.008658
 0.009182
 -0.009774
 0.010448
 -0.011221
 0.012118
 -0.013172
 0.014425
 -0.015943
 0.017818
 -0.020194
 0.023303
 -0.027549
 0.033708
 -0.043510
 0.062551
 -0.161059
 0.812359
 0.377953
 -0.103501
 0.061373
 -0.043660
 0.033889
 -0.027694
 0.023414
 -0.020281
 0.017888
 -0.016000
 0.014472
 -0.013211
 0.012152
 -0.011251
 0.010473
 -0.009797
 0.009202
 -0.008676
 0.008206
 -0.007785
 0.007405
 -0.007060
 0.006746
 -0.006458
 0.006195
 -0.005952
 0.005727
 -0.005518
 0.005325
 -0.005144
 0.004975
 -0.004817
 0.004669
 -0.004529
 0.004398
 -0.004274
 0.004157
 -0.004046
 0.003941
 -0.003841
 0.003746
 -0.003656
 0.003570
 -0.003488
 0.003409
 -0.003334
 0.003263
 -0.003194
 0.003128
 -0.003065
 0.003004
 -0.002946
 0.002890
 -0.002836
 0.002783
 -0.002733
 0.002685
 -0.002638
 0.002593
 -0.002549
 0.002507
 -0.002467
 0.002427
 -0.002389
 0.002352
 -0.002316
 0.002281
 -0.002247
 0.002214
 -0.002183
 0.002152
 -0.002121
 0.002092
 -0.002064
 0.002036
 -0.002009
 0.001983
 -0.001957
 0.001932
 -0.001908
 0.001884
 -0.001861
 0.001839
 -0.001817
 0.001795
 -0.001774
 0.001754
 -0.001733
 0.001714
 -0.001695
 0.001676
 -0.001658
 0.001640
 -0.001622
 0.001605
 -0.001588
 0.001572
 -0.001556
 0.001540
 -0.001524
 0.001509
 -0.001494
 0.001480
 -0.001465
 0.001451
 -0.001438
 0.001424
 -0.001411
 0.001398
 -0.001385
 0.001373
 -0.001360
 0.001348
 -0.001336
 0.001325
 -0.001313
 0.001302
 -0.001291
 0.001280
 -0.001269
 0.001259
 -0.001248
 0.001238
 -0.001228
 0.001218
 -0.001209
 0.001199
 -0.001190
 0.001180
 -0.001171
 0.001162
 -0.001153
 0.001145
 -0.001136
 0.001128
 -0.001119
 0.001111
 -0.001103
 0.001095
 -0.001087
 0.001079
 -0.001072
 0.001064
 -0.001057
 0.001050
 -0.001042
 0.001035
 -0.001028
 0.001021
 -0.001014
 0.001008
 -0.001001
 ];
  nI = 152;
end
if regu==3
  I = [ 0.001003
 -0.001208
 0.001486
 -0.001877
 0.002450
 -0.003342
 0.004846
 -0.007686
 0.013983
 -0.018320
 -0.052992
 0.930205
 0.151315
 -0.027005
 0.011595
 -0.006443
 0.004092
 -0.002823
 0.002061
 -0.001567
 0.001230
 -0.000989
 0.000811
 -0.000677
 0.000572
 ];
  nI = 13;
end
if regu==4
  I = [ 0.001017
 -0.001024
 0.001031
 -0.001038
 0.001045
 -0.001052
 0.001059
 -0.001066
 0.001074
 -0.001081
 0.001089
 -0.001097
 0.001105
 -0.001113
 0.001121
 -0.001129
 0.001137
 -0.001146
 0.001154
 -0.001163
 0.001172
 -0.001181
 0.001190
 -0.001199
 0.001208
 -0.001218
 0.001228
 -0.001237
 0.001247
 -0.001258
 0.001268
 -0.001278
 0.001289
 -0.001300
 0.001311
 -0.001322
 0.001334
 -0.001345
 0.001357
 -0.001369
 0.001381
 -0.001394
 0.001407
 -0.001420
 0.001433
 -0.001446
 0.001460
 -0.001474
 0.001488
 -0.001502
 0.001517
 -0.001532
 0.001548
 -0.001563
 0.001579
 -0.001596
 0.001612
 -0.001629
 0.001647
 -0.001664
 0.001683
 -0.001701
 0.001720
 -0.001739
 0.001759
 -0.001780
 0.001800
 -0.001822
 0.001843
 -0.001866
 0.001888
 -0.001912
 0.001936
 -0.001960
 0.001986
 -0.002011
 0.002038
 -0.002065
 0.002093
 -0.002122
 0.002152
 -0.002182
 0.002213
 -0.002246
 0.002279
 -0.002313
 0.002348
 -0.002384
 0.002422
 -0.002460
 0.002500
 -0.002541
 0.002584
 -0.002628
 0.002673
 -0.002720
 0.002769
 -0.002820
 0.002872
 -0.002927
 0.002983
 -0.003042
 0.003103
 -0.003167
 0.003233
 -0.003302
 0.003375
 -0.003450
 0.003529
 -0.003611
 0.003698
 -0.003789
 0.003884
 -0.003984
 0.004090
 -0.004201
 0.004318
 -0.004443
 0.004574
 -0.004714
 0.004862
 -0.005021
 0.005189
 -0.005370
 0.005563
 -0.005771
 0.005995
 -0.006237
 0.006500
 -0.006786
 0.007097
 -0.007439
 0.007816
 -0.008232
 0.008696
 -0.009215
 0.009799
 -0.010463
 0.011223
 -0.012102
 0.013131
 -0.014350
 0.015819
 -0.017624
 0.019892
 -0.022832
 0.026795
 -0.032432
 0.041115
 -0.056360
 0.085847
 -0.135367
 0.158052
 0.916228
 -0.013877
 0.041053
 -0.034280
 0.028496
 -0.024213
 0.021000
 -0.018520
 0.016556
 -0.014964
 0.013649
 -0.012545
 0.011606
 -0.010797
 0.010093
 -0.009475
 0.008928
 -0.008440
 0.008003
 -0.007609
 0.007252
 -0.006927
 0.006630
 -0.006357
 0.006106
 -0.005874
 0.005658
 -0.005459
 0.005272
 -0.005098
 0.004935
 -0.004782
 0.004639
 -0.004503
 0.004376
 -0.004255
 0.004141
 -0.004033
 0.003930
 -0.003833
 0.003740
 -0.003651
 0.003567
 -0.003487
 0.003410
 -0.003336
 0.003265
 -0.003198
 0.003133
 -0.003070
 0.003011
 -0.002953
 0.002898
 -0.002844
 0.002793
 -0.002743
 0.002695
 -0.002649
 0.002604
 -0.002561
 0.002519
 -0.002479
 0.002440
 -0.002402
 0.002365
 -0.002329
 0.002295
 -0.002261
 0.002228
 -0.002197
 0.002166
 -0.002136
 0.002107
 -0.002078
 0.002051
 -0.002024
 0.001998
 -0.001972
 0.001947
 -0.001923
 0.001899
 -0.001876
 0.001854
 -0.001832
 0.001810
 -0.001789
 0.001769
 -0.001749
 0.001729
 -0.001710
 0.001691
 -0.001673
 0.001655
 -0.001637
 0.001620
 -0.001603
 0.001587
 -0.001571
 0.001555
 -0.001539
 0.001524
 -0.001509
 0.001495
 -0.001480
 0.001466
 -0.001453
 0.001439
 -0.001426
 0.001413
 -0.001400
 0.001387
 -0.001375
 0.001363
 -0.001351
 0.001339
 -0.001328
 0.001316
 -0.001305
 0.001294
 -0.001283
 0.001273
 -0.001262
 0.001252
 -0.001242
 0.001232
 -0.001222
 0.001213
 -0.001203
 0.001194
 -0.001185
 0.001176
 -0.001167
 0.001158
 -0.001150
 0.001141
 -0.001133
 0.001125
 -0.001116
 0.001108
 -0.001100
 0.001093
 -0.001085
 0.001077
 -0.001070
 0.001063
 -0.001055
 0.001048
 -0.001041
 0.001034
 -0.001027
 0.001020
 -0.001014
 0.001007
 -0.001001
 0.000994
 -0.000988
 ];
  nI = 155;
end
if regu==5
  I = [ -0.001014
 0.001019
 -0.001024
 0.001028
 -0.001033
 0.001038
 -0.001043
 0.001047
 -0.001052
 0.001057
 -0.001062
 0.001067
 -0.001073
 0.001078
 -0.001083
 0.001088
 -0.001093
 0.001099
 -0.001104
 0.001110
 -0.001115
 0.001121
 -0.001126
 0.001132
 -0.001138
 0.001144
 -0.001150
 0.001156
 -0.001162
 0.001168
 -0.001174
 0.001180
 -0.001186
 0.001192
 -0.001199
 0.001205
 -0.001212
 0.001218
 -0.001225
 0.001232
 -0.001239
 0.001245
 -0.001252
 0.001259
 -0.001267
 0.001274
 -0.001281
 0.001288
 -0.001296
 0.001303
 -0.001311
 0.001319
 -0.001327
 0.001334
 -0.001342
 0.001351
 -0.001359
 0.001367
 -0.001375
 0.001384
 -0.001393
 0.001401
 -0.001410
 0.001419
 -0.001428
 0.001437
 -0.001446
 0.001456
 -0.001465
 0.001475
 -0.001485
 0.001495
 -0.001505
 0.001515
 -0.001525
 0.001536
 -0.001546
 0.001557
 -0.001568
 0.001579
 -0.001590
 0.001602
 -0.001613
 0.001625
 -0.001637
 0.001649
 -0.001661
 0.001673
 -0.001686
 0.001699
 -0.001712
 0.001725
 -0.001738
 0.001752
 -0.001766
 0.001780
 -0.001794
 0.001808
 -0.001823
 0.001838
 -0.001853
 0.001869
 -0.001884
 0.001900
 -0.001917
 0.001933
 -0.001950
 0.001967
 -0.001985
 0.002002
 -0.002020
 0.002039
 -0.002057
 0.002077
 -0.002096
 0.002116
 -0.002136
 0.002156
 -0.002177
 0.002199
 -0.002221
 0.002243
 -0.002265
 0.002289
 -0.002312
 0.002336
 -0.002361
 0.002386
 -0.002412
 0.002438
 -0.002465
 0.002492
 -0.002520
 0.002549
 -0.002578
 0.002608
 -0.002639
 0.002671
 -0.002703
 0.002736
 -0.002770
 0.002804
 -0.002840
 0.002876
 -0.002914
 0.002952
 -0.002991
 0.003032
 -0.003073
 0.003116
 -0.003160
 0.003205
 -0.003252
 0.003300
 -0.003349
 0.003400
 -0.003452
 0.003506
 -0.003562
 0.003620
 -0.003679
 0.003741
 -0.003804
 0.003870
 -0.003938
 0.004008
 -0.004081
 0.004157
 -0.004236
 0.004317
 -0.004402
 0.004490
 -0.004582
 0.004678
 -0.004778
 0.004882
 -0.004990
 0.005104
 -0.005223
 0.005348
 -0.005479
 0.005616
 -0.005760
 0.005912
 -0.006073
 0.006242
 -0.006421
 0.006610
 -0.006811
 0.007025
 -0.007252
 0.007495
 -0.007755
 0.008033
 -0.008332
 0.008654
 -0.009002
 0.009379
 -0.009789
 0.010236
 -0.010727
 0.011267
 -0.011864
 0.012528
 -0.013272
 0.014109
 -0.015059
 0.016147
 -0.017404
 0.018876
 -0.020620
 0.022722
 -0.025305
 0.028558
 -0.032789
 0.038529
 -0.046820
 0.060114
 -0.084888
 0.136125
 -0.241282
 0.424404
 0.782133
 -0.092379
 0.072524
 -0.055165
 0.044306
 -0.036988
 0.031739
 -0.027793
 0.024721
 -0.022260
 0.020245
 -0.018565
 0.017143
 -0.015923
 0.014865
 -0.013940
 0.013123
 -0.012396
 0.011746
 -0.011160
 0.010630
 -0.010149
 0.009709
 -0.009305
 0.008934
 -0.008591
 0.008274
 -0.007979
 0.007705
 -0.007448
 0.007209
 -0.006984
 0.006773
 -0.006574
 0.006387
 -0.006210
 0.006042
 -0.005883
 0.005733
 -0.005590
 0.005454
 -0.005324
 0.005201
 -0.005083
 0.004970
 -0.004862
 0.004759
 -0.004660
 0.004565
 -0.004474
 0.004386
 -0.004302
 0.004221
 -0.004143
 0.004067
 -0.003995
 0.003925
 -0.003857
 0.003792
 -0.003729
 0.003668
 -0.003609
 0.003552
 -0.003496
 0.003442
 -0.003390
 0.003340
 -0.003291
 0.003243
 -0.003197
 0.003152
 -0.003108
 0.003066
 -0.003024
 0.002984
 -0.002945
 0.002907
 -0.002869
 0.002833
 -0.002798
 0.002763
 -0.002730
 0.002697
 -0.002665
 0.002633
 -0.002603
 0.002573
 -0.002544
 0.002515
 -0.002487
 0.002460
 -0.002433
 0.002407
 -0.002381
 0.002356
 -0.002332
 0.002308
 -0.002284
 0.002261
 -0.002239
 0.002216
 -0.002195
 0.002173
 -0.002153
 0.002132
 -0.002112
 0.002092
 -0.002073
 0.002054
 -0.002035
 0.002017
 -0.001999
 0.001981
 -0.001964
 0.001947
 -0.001930
 0.001914
 -0.001897
 0.001881
 -0.001866
 0.001850
 -0.001835
 0.001820
 -0.001806
 0.001791
 -0.001777
 0.001763
 -0.001749
 0.001736
 -0.001722
 0.001709
 -0.001696
 0.001684
 -0.001671
 0.001659
 -0.001646
 0.001634
 -0.001623
 0.001611
 -0.001599
 0.001588
 -0.001577
 0.001566
 -0.001555
 0.001544
 -0.001534
 0.001523
 -0.001513
 0.001503
 -0.001493
 0.001483
 -0.001473
 0.001464
 -0.001454
 0.001445
 -0.001435
 0.001426
 -0.001417
 0.001408
 -0.001400
 0.001391
 -0.001382
 0.001374
 -0.001365
 0.001357
 -0.001349
 0.001341
 -0.001333
 0.001325
 -0.001317
 0.001310
 -0.001302
 0.001294
 -0.001287
 0.001280
 -0.001272
 0.001265
 -0.001258
 0.001251
 -0.001244
 0.001237
 -0.001231
 0.001224
 -0.001217
 0.001211
 -0.001204
 0.001198
 -0.001191
 0.001185
 -0.001179
 0.001173
 -0.001166
 0.001160
 -0.001154
 0.001149
 -0.001143
 0.001137
 -0.001131
 0.001125
 -0.001120
 0.001114
 -0.001109
 0.001103
 -0.001098
 0.001093
 -0.001087
 0.001082
 -0.001077
 0.001072
 -0.001067
 0.001061
 -0.001056
 0.001051
 -0.001047
 0.001042
 -0.001037
 0.001032
 -0.001027
 0.001023
 -0.001018
 0.001013
 -0.001009
 0.001004
 -0.001000
 0.000995
 ];
  nI = 224;
end
if regu==6
  I = [ -0.001013
 0.001019
 -0.001025
 0.001031
 -0.001037
 0.001043
 -0.001049
 0.001055
 -0.001062
 0.001068
 -0.001075
 0.001081
 -0.001088
 0.001095
 -0.001102
 0.001109
 -0.001116
 0.001123
 -0.001130
 0.001138
 -0.001145
 0.001153
 -0.001160
 0.001168
 -0.001176
 0.001184
 -0.001192
 0.001200
 -0.001208
 0.001217
 -0.001225
 0.001234
 -0.001242
 0.001251
 -0.001260
 0.001269
 -0.001279
 0.001288
 -0.001298
 0.001307
 -0.001317
 0.001327
 -0.001337
 0.001347
 -0.001358
 0.001368
 -0.001379
 0.001390
 -0.001401
 0.001413
 -0.001424
 0.001436
 -0.001448
 0.001460
 -0.001472
 0.001484
 -0.001497
 0.001510
 -0.001523
 0.001536
 -0.001550
 0.001564
 -0.001578
 0.001592
 -0.001607
 0.001622
 -0.001637
 0.001652
 -0.001668
 0.001684
 -0.001700
 0.001717
 -0.001734
 0.001751
 -0.001769
 0.001787
 -0.001805
 0.001824
 -0.001843
 0.001863
 -0.001883
 0.001903
 -0.001924
 0.001946
 -0.001967
 0.001990
 -0.002013
 0.002036
 -0.002060
 0.002084
 -0.002109
 0.002135
 -0.002161
 0.002188
 -0.002216
 0.002244
 -0.002274
 0.002303
 -0.002334
 0.002366
 -0.002398
 0.002431
 -0.002465
 0.002501
 -0.002537
 0.002574
 -0.002612
 0.002652
 -0.002693
 0.002735
 -0.002778
 0.002823
 -0.002869
 0.002917
 -0.002966
 0.003017
 -0.003070
 0.003125
 -0.003181
 0.003240
 -0.003301
 0.003365
 -0.003430
 0.003499
 -0.003570
 0.003644
 -0.003722
 0.003803
 -0.003887
 0.003975
 -0.004067
 0.004164
 -0.004265
 0.004372
 -0.004484
 0.004601
 -0.004725
 0.004856
 -0.004995
 0.005141
 -0.005297
 0.005462
 -0.005638
 0.005825
 -0.006025
 0.006240
 -0.006471
 0.006719
 -0.006987
 0.007277
 -0.007593
 0.007937
 -0.008314
 0.008729
 -0.009187
 0.009696
 -0.010266
 0.010906
 -0.011632
 0.012463
 -0.013422
 0.014542
 -0.015868
 0.017463
 -0.019422
 0.021889
 -0.025104
 0.029505
 -0.036088
 0.047700
 -0.072713
 0.133098
 -0.284929
 0.682266
 0.570610
 -0.092553
 0.062583
 -0.046110
 0.036457
 -0.030145
 0.025700
 -0.022399
 0.019852
 -0.017826
 0.016176
 -0.014806
 0.013650
 -0.012662
 0.011808
 -0.011061
 0.010404
 -0.009820
 0.009299
 -0.008830
 0.008406
 -0.008021
 0.007670
 -0.007348
 0.007052
 -0.006779
 0.006527
 -0.006292
 0.006074
 -0.005871
 0.005680
 -0.005502
 0.005335
 -0.005177
 0.005028
 -0.004888
 0.004756
 -0.004630
 0.004511
 -0.004398
 0.004290
 -0.004187
 0.004090
 -0.003996
 0.003907
 -0.003822
 0.003741
 -0.003662
 0.003587
 -0.003516
 0.003446
 -0.003380
 0.003316
 -0.003254
 0.003195
 -0.003138
 0.003083
 -0.003029
 0.002978
 -0.002928
 0.002880
 -0.002834
 0.002789
 -0.002745
 0.002703
 -0.002661
 0.002622
 -0.002583
 0.002546
 -0.002509
 0.002474
 -0.002439
 0.002406
 -0.002373
 0.002342
 -0.002311
 0.002281
 -0.002251
 0.002223
 -0.002195
 0.002168
 -0.002141
 0.002115
 -0.002090
 0.002066
 -0.002042
 0.002018
 -0.001995
 0.001973
 -0.001951
 0.001929
 -0.001908
 0.001888
 -0.001868
 0.001848
 -0.001829
 0.001810
 -0.001791
 0.001773
 -0.001755
 0.001738
 -0.001721
 0.001704
 -0.001688
 0.001672
 -0.001656
 0.001640
 -0.001625
 0.001610
 -0.001596
 0.001581
 -0.001567
 0.001553
 -0.001540
 0.001526
 -0.001513
 0.001500
 -0.001487
 0.001475
 -0.001462
 0.001450
 -0.001438
 0.001427
 -0.001415
 0.001404
 -0.001393
 0.001382
 -0.001371
 0.001360
 -0.001350
 0.001340
 -0.001329
 0.001319
 -0.001310
 0.001300
 -0.001290
 0.001281
 -0.001272
 0.001262
 -0.001253
 0.001244
 -0.001236
 0.001227
 -0.001219
 0.001210
 -0.001202
 0.001194
 -0.001186
 0.001178
 -0.001170
 0.001162
 -0.001154
 0.001147
 -0.001139
 0.001132
 -0.001125
 0.001118
 -0.001110
 0.001103
 -0.001097
 0.001090
 -0.001083
 0.001076
 -0.001070
 0.001063
 -0.001057
 0.001051
 -0.001044
 0.001038
 -0.001032
 0.001026
 -0.001020
 0.001014
 -0.001008
 0.001003
 -0.000997
 0.000991
 ];
  nI = 176;
end
if regu==7
  I = [ -0.000965
 0.000991
 -0.001018
 0.001047
 -0.001078
 0.001110
 -0.001144
 0.001181
 -0.001219
 0.001260
 -0.001304
 0.001351
 -0.001400
 0.001453
 -0.001510
 0.001571
 -0.001635
 0.001705
 -0.001778
 0.001855
 -0.001936
 0.002018
 -0.002100
 0.002173
 -0.002228
 0.002244
 -0.002183
 0.002017
 -0.002015
 0.004102
 -0.016083
 0.063135
 -0.229502
 0.869556
 0.339201
 -0.044783
 0.025143
 -0.016590
 0.012136
 -0.009467
 0.007712
 -0.006481
 0.005574
 -0.004880
 0.004335
 -0.003895
 0.003534
 -0.003233
 0.002977
 -0.002758
 0.002569
 -0.002403
 0.002257
 -0.002127
 0.002011
 -0.001907
 0.001813
 -0.001728
 0.001650
 -0.001579
 0.001514
 -0.001454
 0.001398
 -0.001346
 0.001298
 -0.001254
 0.001212
 -0.001173
 0.001136
 -0.001102
 0.001069
 ];
  nI = 36;
end
if regu==8
  I = [ -0.001048
 0.001057
 -0.001066
 0.001076
 -0.001086
 0.001095
 -0.001105
 0.001116
 -0.001126
 0.001136
 -0.001147
 0.001158
 -0.001169
 0.001181
 -0.001192
 0.001204
 -0.001216
 0.001229
 -0.001241
 0.001254
 -0.001267
 0.001281
 -0.001294
 0.001308
 -0.001322
 0.001337
 -0.001352
 0.001367
 -0.001383
 0.001399
 -0.001415
 0.001432
 -0.001449
 0.001466
 -0.001484
 0.001503
 -0.001522
 0.001541
 -0.001561
 0.001581
 -0.001602
 0.001624
 -0.001646
 0.001668
 -0.001692
 0.001716
 -0.001740
 0.001766
 -0.001792
 0.001818
 -0.001846
 0.001875
 -0.001904
 0.001935
 -0.001966
 0.001998
 -0.002032
 0.002066
 -0.002102
 0.002139
 -0.002178
 0.002218
 -0.002259
 0.002302
 -0.002346
 0.002393
 -0.002441
 0.002491
 -0.002543
 0.002598
 -0.002655
 0.002714
 -0.002776
 0.002841
 -0.002909
 0.002981
 -0.003056
 0.003135
 -0.003218
 0.003306
 -0.003399
 0.003497
 -0.003600
 0.003710
 -0.003827
 0.003952
 -0.004085
 0.004227
 -0.004380
 0.004544
 -0.004720
 0.004911
 -0.005118
 0.005343
 -0.005589
 0.005859
 -0.006155
 0.006484
 -0.006849
 0.007257
 -0.007717
 0.008239
 -0.008837
 0.009528
 -0.010335
 0.011291
 -0.012441
 0.013851
 -0.015618
 0.017900
 -0.020958
 0.025259
 -0.031620
 0.041205
 -0.054606
 0.068174
 -0.061273
 -0.066417
 0.944190
 0.141458
 0.013597
 -0.017159
 0.016410
 -0.014973
 0.013566
 -0.012326
 0.011259
 -0.010346
 0.009561
 -0.008881
 0.008288
 -0.007767
 0.007307
 -0.006897
 0.006529
 -0.006199
 0.005900
 -0.005628
 0.005380
 -0.005152
 0.004943
 -0.004751
 0.004572
 -0.004407
 0.004253
 -0.004109
 0.003975
 -0.003849
 0.003731
 -0.003619
 0.003515
 -0.003416
 0.003322
 -0.003234
 0.003150
 -0.003070
 0.002995
 -0.002922
 0.002854
 -0.002788
 0.002725
 -0.002666
 0.002608
 -0.002553
 0.002501
 -0.002450
 0.002402
 -0.002355
 0.002310
 -0.002267
 0.002225
 -0.002185
 0.002147
 -0.002109
 0.002073
 -0.002038
 0.002005
 -0.001972
 0.001941
 -0.001910
 0.001880
 -0.001852
 0.001824
 -0.001797
 0.001770
 -0.001745
 0.001720
 -0.001696
 0.001673
 -0.001650
 0.001628
 -0.001606
 0.001585
 -0.001565
 0.001545
 -0.001525
 0.001506
 -0.001488
 0.001470
 -0.001452
 0.001435
 -0.001418
 0.001402
 -0.001386
 0.001370
 -0.001355
 0.001340
 -0.001325
 0.001311
 -0.001297
 0.001283
 -0.001270
 0.001257
 -0.001244
 0.001231
 -0.001219
 0.001207
 -0.001195
 0.001183
 -0.001172
 0.001160
 -0.001149
 0.001139
 -0.001128
 0.001118
 -0.001107
 0.001097
 -0.001087
 0.001078
 -0.001068
 0.001059
 -0.001050
 0.001041
 -0.001032
 0.001023
 -0.001015
 0.001006
 -0.000998
 0.000990
 -0.000982
 0.000974
 ];
  nI = 121;
end
if regu==9
  I = [ -0.001020
 0.001025
 -0.001030
 0.001035
 -0.001040
 0.001045
 -0.001050
 0.001055
 -0.001060
 0.001065
 -0.001071
 0.001076
 -0.001081
 0.001087
 -0.001092
 0.001098
 -0.001103
 0.001109
 -0.001115
 0.001120
 -0.001126
 0.001132
 -0.001138
 0.001144
 -0.001150
 0.001156
 -0.001163
 0.001169
 -0.001175
 0.001182
 -0.001188
 0.001195
 -0.001201
 0.001208
 -0.001215
 0.001222
 -0.001229
 0.001236
 -0.001243
 0.001250
 -0.001257
 0.001264
 -0.001272
 0.001279
 -0.001287
 0.001295
 -0.001302
 0.001310
 -0.001318
 0.001326
 -0.001335
 0.001343
 -0.001351
 0.001360
 -0.001368
 0.001377
 -0.001386
 0.001395
 -0.001404
 0.001413
 -0.001422
 0.001432
 -0.001441
 0.001451
 -0.001461
 0.001471
 -0.001481
 0.001491
 -0.001501
 0.001512
 -0.001522
 0.001533
 -0.001544
 0.001555
 -0.001566
 0.001578
 -0.001589
 0.001601
 -0.001613
 0.001625
 -0.001637
 0.001650
 -0.001663
 0.001675
 -0.001689
 0.001702
 -0.001715
 0.001729
 -0.001743
 0.001757
 -0.001771
 0.001786
 -0.001801
 0.001816
 -0.001831
 0.001847
 -0.001863
 0.001879
 -0.001895
 0.001912
 -0.001929
 0.001947
 -0.001964
 0.001982
 -0.002001
 0.002019
 -0.002038
 0.002058
 -0.002077
 0.002097
 -0.002118
 0.002139
 -0.002160
 0.002182
 -0.002204
 0.002227
 -0.002250
 0.002274
 -0.002298
 0.002323
 -0.002348
 0.002373
 -0.002400
 0.002427
 -0.002454
 0.002482
 -0.002511
 0.002541
 -0.002571
 0.002602
 -0.002633
 0.002666
 -0.002699
 0.002733
 -0.002768
 0.002804
 -0.002841
 0.002878
 -0.002917
 0.002957
 -0.002998
 0.003040
 -0.003083
 0.003128
 -0.003174
 0.003221
 -0.003270
 0.003320
 -0.003371
 0.003425
 -0.003480
 0.003537
 -0.003595
 0.003656
 -0.003719
 0.003784
 -0.003851
 0.003921
 -0.003993
 0.004068
 -0.004146
 0.004227
 -0.004312
 0.004399
 -0.004490
 0.004586
 -0.004685
 0.004788
 -0.004897
 0.005010
 -0.005129
 0.005253
 -0.005384
 0.005521
 -0.005666
 0.005818
 -0.005979
 0.006149
 -0.006329
 0.006519
 -0.006722
 0.006937
 -0.007167
 0.007413
 -0.007676
 0.007958
 -0.008262
 0.008590
 -0.008946
 0.009332
 -0.009753
 0.010213
 -0.010720
 0.011279
 -0.011900
 0.012594
 -0.013374
 0.014258
 -0.015267
 0.016430
 -0.017787
 0.019391
 -0.021316
 0.023672
 -0.026626
 0.030449
 -0.035615
 0.043010
 -0.054327
 0.072551
 -0.102090
 0.147274
 -0.202975
 0.180067
 0.895598
 0.011657
 0.052620
 -0.044142
 0.037106
 -0.031820
 0.027798
 -0.024659
 0.022150
 -0.020100
 0.018395
 -0.016956
 0.015726
 -0.014661
 0.013732
 -0.012913
 0.012186
 -0.011536
 0.010953
 -0.010425
 0.009946
 -0.009509
 0.009109
 -0.008741
 0.008402
 -0.008088
 0.007797
 -0.007526
 0.007273
 -0.007036
 0.006815
 -0.006607
 0.006411
 -0.006227
 0.006053
 -0.005888
 0.005732
 -0.005584
 0.005444
 -0.005310
 0.005183
 -0.005062
 0.004946
 -0.004836
 0.004730
 -0.004629
 0.004532
 -0.004439
 0.004350
 -0.004264
 0.004182
 -0.004103
 0.004026
 -0.003953
 0.003882
 -0.003813
 0.003747
 -0.003684
 0.003622
 -0.003562
 0.003505
 -0.003449
 0.003395
 -0.003342
 0.003292
 -0.003242
 0.003194
 -0.003148
 0.003103
 -0.003059
 0.003016
 -0.002975
 0.002935
 -0.002895
 0.002857
 -0.002820
 0.002784
 -0.002748
 0.002714
 -0.002680
 0.002648
 -0.002616
 0.002584
 -0.002554
 0.002524
 -0.002495
 0.002467
 -0.002439
 0.002412
 -0.002385
 0.002359
 -0.002334
 0.002309
 -0.002284
 0.002261
 -0.002237
 0.002214
 -0.002192
 0.002170
 -0.002148
 0.002127
 -0.002107
 0.002086
 -0.002066
 0.002047
 -0.002028
 0.002009
 -0.001990
 0.001972
 -0.001954
 0.001937
 -0.001920
 0.001903
 -0.001886
 0.001870
 -0.001854
 0.001838
 -0.001823
 0.001808
 -0.001793
 0.001778
 -0.001763
 0.001749
 -0.001735
 0.001721
 -0.001708
 0.001694
 -0.001681
 0.001668
 -0.001655
 0.001643
 -0.001631
 0.001618
 -0.001606
 0.001595
 -0.001583
 0.001571
 -0.001560
 0.001549
 -0.001538
 0.001527
 -0.001516
 0.001506
 -0.001495
 0.001485
 -0.001475
 0.001465
 -0.001455
 0.001445
 -0.001436
 0.001426
 -0.001417
 0.001408
 -0.001399
 0.001390
 -0.001381
 0.001372
 -0.001364
 0.001355
 -0.001347
 0.001338
 -0.001330
 0.001322
 -0.001314
 0.001306
 -0.001298
 0.001290
 -0.001283
 0.001275
 -0.001268
 0.001260
 -0.001253
 0.001246
 -0.001239
 0.001232
 -0.001225
 0.001218
 -0.001211
 0.001204
 -0.001198
 0.001191
 -0.001185
 0.001178
 -0.001172
 0.001165
 -0.001159
 0.001153
 -0.001147
 0.001141
 -0.001135
 0.001129
 -0.001123
 0.001117
 -0.001112
 0.001106
 -0.001100
 0.001095
 -0.001089
 0.001084
 -0.001078
 0.001073
 -0.001068
 0.001062
 -0.001057
 0.001052
 -0.001047
 0.001042
 -0.001037
 0.001032
 -0.001027
 0.001022
 -0.001017
 0.001013
 -0.001008
 0.001003
 -0.000999
 0.000994
 -0.000990
 0.000985
 ];
  nI = 217;
end
if regu==10
  I = [ -0.001023
 0.001028
 -0.001033
 0.001039
 -0.001044
 0.001049
 -0.001055
 0.001060
 -0.001066
 0.001071
 -0.001077
 0.001083
 -0.001089
 0.001095
 -0.001101
 0.001107
 -0.001113
 0.001119
 -0.001125
 0.001131
 -0.001138
 0.001144
 -0.001150
 0.001157
 -0.001164
 0.001170
 -0.001177
 0.001184
 -0.001191
 0.001198
 -0.001205
 0.001212
 -0.001219
 0.001227
 -0.001234
 0.001242
 -0.001250
 0.001257
 -0.001265
 0.001273
 -0.001281
 0.001289
 -0.001297
 0.001306
 -0.001314
 0.001323
 -0.001331
 0.001340
 -0.001349
 0.001358
 -0.001367
 0.001377
 -0.001386
 0.001395
 -0.001405
 0.001415
 -0.001425
 0.001435
 -0.001445
 0.001456
 -0.001466
 0.001477
 -0.001488
 0.001498
 -0.001510
 0.001521
 -0.001532
 0.001544
 -0.001556
 0.001568
 -0.001580
 0.001593
 -0.001605
 0.001618
 -0.001631
 0.001644
 -0.001658
 0.001671
 -0.001685
 0.001699
 -0.001713
 0.001728
 -0.001743
 0.001758
 -0.001773
 0.001789
 -0.001805
 0.001821
 -0.001838
 0.001854
 -0.001871
 0.001889
 -0.001907
 0.001925
 -0.001943
 0.001962
 -0.001981
 0.002000
 -0.002020
 0.002041
 -0.002061
 0.002083
 -0.002104
 0.002126
 -0.002149
 0.002172
 -0.002195
 0.002219
 -0.002244
 0.002269
 -0.002294
 0.002321
 -0.002347
 0.002375
 -0.002403
 0.002432
 -0.002461
 0.002492
 -0.002522
 0.002554
 -0.002587
 0.002620
 -0.002654
 0.002690
 -0.002726
 0.002763
 -0.002801
 0.002840
 -0.002880
 0.002922
 -0.002964
 0.003008
 -0.003054
 0.003100
 -0.003148
 0.003198
 -0.003249
 0.003302
 -0.003357
 0.003413
 -0.003471
 0.003532
 -0.003594
 0.003659
 -0.003726
 0.003796
 -0.003868
 0.003944
 -0.004022
 0.004103
 -0.004188
 0.004276
 -0.004368
 0.004464
 -0.004565
 0.004670
 -0.004780
 0.004895
 -0.005016
 0.005143
 -0.005277
 0.005418
 -0.005566
 0.005723
 -0.005890
 0.006066
 -0.006253
 0.006452
 -0.006664
 0.006890
 -0.007133
 0.007393
 -0.007673
 0.007975
 -0.008302
 0.008656
 -0.009043
 0.009466
 -0.009931
 0.010443
 -0.011012
 0.011647
 -0.012360
 0.013167
 -0.014088
 0.015149
 -0.016385
 0.017846
 -0.019600
 0.021749
 -0.024456
 0.027995
 -0.032887
 0.040207
 -0.052223
 0.073493
 -0.112526
 0.184278
 -0.311121
 0.459272
 0.745705
 -0.042259
 0.058621
 -0.046110
 0.037653
 -0.031759
 0.027447
 -0.024164
 0.021582
 -0.019500
 0.017784
 -0.016347
 0.015125
 -0.014073
 0.013158
 -0.012355
 0.011645
 -0.011012
 0.010444
 -0.009932
 0.009468
 -0.009046
 0.008660
 -0.008305
 0.007978
 -0.007676
 0.007396
 -0.007136
 0.006893
 -0.006667
 0.006455
 -0.006256
 0.006069
 -0.005892
 0.005726
 -0.005569
 0.005420
 -0.005279
 0.005145
 -0.005018
 0.004897
 -0.004782
 0.004672
 -0.004566
 0.004466
 -0.004370
 0.004278
 -0.004189
 0.004105
 -0.004023
 0.003945
 -0.003870
 0.003797
 -0.003728
 0.003661
 -0.003596
 0.003533
 -0.003473
 0.003414
 -0.003358
 0.003303
 -0.003250
 0.003199
 -0.003149
 0.003101
 -0.003055
 0.003009
 -0.002965
 0.002923
 -0.002881
 0.002841
 -0.002802
 0.002764
 -0.002726
 0.002690
 -0.002655
 0.002621
 -0.002587
 0.002555
 -0.002523
 0.002492
 -0.002462
 0.002432
 -0.002404
 0.002376
 -0.002348
 0.002321
 -0.002295
 0.002269
 -0.002244
 0.002220
 -0.002196
 0.002172
 -0.002149
 0.002127
 -0.002105
 0.002083
 -0.002062
 0.002041
 -0.002021
 0.002001
 -0.001981
 0.001962
 -0.001943
 0.001925
 -0.001907
 0.001889
 -0.001872
 0.001855
 -0.001838
 0.001821
 -0.001805
 0.001789
 -0.001774
 0.001758
 -0.001743
 0.001728
 -0.001714
 0.001699
 -0.001685
 0.001671
 -0.001658
 0.001644
 -0.001631
 0.001618
 -0.001605
 0.001593
 -0.001580
 0.001568
 -0.001556
 0.001544
 -0.001533
 0.001521
 -0.001510
 0.001499
 -0.001488
 0.001477
 -0.001466
 0.001456
 -0.001445
 0.001435
 -0.001425
 0.001415
 -0.001405
 0.001396
 -0.001386
 0.001377
 -0.001368
 0.001358
 -0.001349
 0.001340
 -0.001332
 0.001323
 -0.001314
 0.001306
 -0.001298
 0.001289
 -0.001281
 0.001273
 -0.001265
 0.001257
 -0.001250
 0.001242
 -0.001234
 0.001227
 -0.001220
 0.001212
 -0.001205
 0.001198
 -0.001191
 0.001184
 -0.001177
 0.001170
 -0.001164
 0.001157
 -0.001151
 0.001144
 -0.001138
 0.001131
 -0.001125
 0.001119
 -0.001113
 0.001107
 -0.001101
 0.001095
 -0.001089
 0.001083
 -0.001077
 0.001072
 -0.001066
 0.001060
 -0.001055
 0.001049
 -0.001044
 0.001039
 -0.001033
 0.001028
 -0.001023
 0.001018
 -0.001013
 0.001008
 -0.001003
 0.000998
 -0.000993
 0.000988
 ];
  nI = 203;
end
I = I';
end