from scidbbiobank import connect
import lookups
import numpy
import scipy.stats as stats
dcauc = [0.997,
0.996,
0.994,
0.991,
0.988,
0.985,
0.985,
0.985,
0.984,
0.981,
0.979,
0.978,
0.978,
0.976,
0.976,
0.976,
0.976,
0.975,
0.972,
0.972,
0.969,
0.966,
0.962,
0.961,
0.961,
0.961,
0.96,
0.96,
0.955,
0.955,
0.952,
0.949,
0.948,
0.947,
0.943,
0.94,
0.94,
0.939,
0.938,
0.937,
0.937,
0.936,
0.936,
0.935,
0.935,
0.935,
0.934,
0.934,
0.933,
0.932,
0.932,
0.931,
0.93,
0.928,
0.928,
0.927,
0.926,
0.926,
0.924,
0.923,
0.923,
0.922,
0.922,
0.921,
0.92,
0.92,
0.92,
0.919,
0.919,
0.918,
0.918,
0.917,
0.917,
0.916,
0.915,
0.915,
0.914,
0.913,
0.913,
0.912,
0.911,
0.911,
0.909,
0.908,
0.907,
0.906,
0.905,
0.905,
0.904,
0.904,
0.904,
0.904,
0.903,
0.903,
0.903,
0.902,
0.902,
0.901,
0.901,
0.9,
0.9,
0.899,
0.899,
0.899,
0.898,
0.897,
0.896,
0.896,
0.896,
0.895,
0.894,
0.892,
0.892,
0.89,
0.887,
0.887,
0.887,
0.886,
0.886,
0.885,
0.884,
0.883,
0.883,
0.882,
0.881,
0.881,
0.88,
0.88,
0.88,
0.88,
0.88,
0.879,
0.878,
0.878,
0.878,
0.877,
0.877,
0.877,
0.877,
0.877,
0.876,
0.875,
0.874,
0.873,
0.873,
0.873,
0.872,
0.872,
0.872,
0.872,
0.872,
0.871,
0.871,
0.871,
0.87,
0.87,
0.869,
0.869,
0.869,
0.868,
0.867,
0.867,
0.866,
0.866,
0.865,
0.865,
0.865,
0.864,
0.864,
0.863,
0.863,
0.863,
0.863,
0.862,
0.861,
0.861,
0.861,
0.861,
0.861,
0.86,
0.86,
0.86,
0.859,
0.857,
0.857,
0.857,
0.855,
0.855,
0.855,
0.854,
0.854,
0.854,
0.853,
0.852,
0.852,
0.852,
0.852,
0.85,
0.85,
0.85,
0.849,
0.849,
0.849,
0.848,
0.848,
0.848,
0.848,
0.847,
0.846,
0.846,
0.846,
0.845,
0.844,
0.842,
0.842,
0.841,
0.84,
0.84,
0.84,
0.838,
0.838,
0.837,
0.837,
0.836,
0.835,
0.835,
0.835,
0.834,
0.834,
0.833,
0.833,
0.833,
0.833,
0.832,
0.832,
0.832,
0.831,
0.831,
0.831,
0.83,
0.83,
0.83,
0.83,
0.829,
0.829,
0.829,
0.826,
0.826,
0.826,
0.826,
0.824,
0.824,
0.823,
0.823,
0.823,
0.823,
0.822,
0.822,
0.821,
0.821,
0.821,
0.821,
0.821,
0.82,
0.82,
0.82,
0.817,
0.816,
0.815,
0.815,
0.815,
0.815,
0.814,
0.814,
0.814,
0.814,
0.813,
0.813,
0.813,
0.813,
0.812,
0.812,
0.812,
0.811,
0.81,
0.81,
0.81,
0.809,
0.809,
0.809,
0.808,
0.808,
0.807,
0.807,
0.807,
0.806,
0.806,
0.806,
0.805,
0.804,
0.804,
0.803,
0.802,
0.801,
0.801,
0.801,
0.8,
0.8,
0.8,
0.799,
0.798,
0.798,
0.798,
0.798,
0.798,
0.798,
0.797,
0.796,
0.796,
0.795,
0.794,
0.793,
0.792,
0.792,
0.792,
0.791,
0.791,
0.79,
0.79,
0.79,
0.79,
0.79,
0.789,
0.788,
0.787,
0.787,
0.787,
0.787,
0.786,
0.785,
0.785,
0.784,
0.784,
0.783,
0.783,
0.783,
0.783,
0.782,
0.781,
0.781,
0.78,
0.779,
0.779,
0.779,
0.778,
0.778,
0.778,
0.778,
0.776,
0.776,
0.775,
0.775,
0.775,
0.774,
0.773,
0.773,
0.772,
0.772,
0.772,
0.772,
0.772,
0.772,
0.77,
0.77,
0.77,
0.769,
0.769,
0.769,
0.769,
0.769,
0.768,
0.768,
0.767,
0.766,
0.765,
0.765,
0.765,
0.764,
0.764,
0.763,
0.763,
0.762,
0.762,
0.76,
0.76,
0.759,
0.757,
0.757,
0.756,
0.756,
0.756,
0.756,
0.756,
0.756,
0.755,
0.755,
0.755,
0.755,
0.754,
0.754,
0.754,
0.753,
0.753,
0.752,
0.752,
0.75,
0.75,
0.75,
0.749,
0.748,
0.748,
0.748,
0.748,
0.747,
0.746,
0.746,
0.745,
0.745,
0.744,
0.744,
0.744,
0.744,
0.743,
0.743,
0.742,
0.742,
0.742,
0.741,
0.741,
0.74,
0.74,
0.74,
0.74,
0.74,
0.74,
0.739,
0.739,
0.739,
0.739,
0.737,
0.737,
0.736,
0.735,
0.734,
0.732,
0.731,
0.731,
0.73,
0.73,
0.73,
0.728,
0.728,
0.727,
0.727,
0.727,
0.726,
0.726,
0.726,
0.724,
0.724,
0.724,
0.723,
0.723,
0.722,
0.722,
0.721,
0.721,
0.719,
0.719,
0.718,
0.717,
0.716,
0.714,
0.713,
0.713,
0.713,
0.712,
0.712,
0.712,
0.712,
0.712,
0.71,
0.707,
0.707,
0.705,
0.703,
0.7,
0.7,
0.698,
0.698,
0.698,
0.698,
0.698,
0.697,
0.697,
0.696,
0.696,
0.696,
0.694,
0.693,
0.693,
0.692,
0.692,
0.692,
0.692,
0.689,
0.686,
0.684,
0.684,
0.684,
0.682,
0.681,
0.677,
0.677,
0.676,
0.675,
0.675,
0.674,
0.673,
0.67,
0.668,
0.666,
0.665,
0.662,
0.661,
0.659,
0.655,
0.655,
0.652,
0.651,
0.648,
0.648,
0.646,
0.646,
0.638,
0.636,
0.634,
0.633,
0.623,
0.621,
0.619,
0.61,
0.608,
0.603,
0.586,
0.567,
0.562,
0.447]
dcgene = ['BACE2',
'CEACAM16',
'FABP3',
'CLU',
'RBP7',
'MMP15',
'SPARC',
'CRHR1',
'GDF6',
'PLXNB2',
'DBI',
'CTSB',
'HS3ST1',
'PARM1',
'CAR14',
'CYBRD1',
'MANSC4',
'FGFR3',
'ACSL3',
'MGST1',
'STMN1',
'WNT7A',
'SMAGP',
'GM19345',
'PDE6D',
'CD81',
'TM4SF1',
'PDE7B',
'FAM101B',
'DPP4',
'ANXA9',
'TSPO',
'CYB5A',
'KLHL14',
'EPB41',
'TMEM79',
'CD47',
'COPS7B',
'GJB6',
'SCD1',
'CERS2',
'TSPAN8',
'GPRC5C',
'MID1IP1',
'TMOD1',
'PTMS',
'S100A1',
'JAG1',
'KCNJ10',
'HERPUD1',
'1190002N15RIK',
'CAR13',
'LDHB',
'SELM',
'S100A13',
'KCNJ16',
'FAR1',
'TMEM221',
'4930523C07RIK',
'TMEM86A',
'KCNB2',
'ZFP36L1',
'SLC6A6',
'CAV2',
'CEP41',
'CMTM7',
'PTPRF',
'OTOG',
'CD9',
'SLC12A7',
'TMPRSS2',
'GLMP',
'ELOVL5',
'SERINC2',
'IGFBP2',
'LPAR3',
'PTPN5',
'UTP14B',
'MDM1',
'CST3',
'TRPM3',
'GAS2',
'TUBA1A',
'GNG5',
'SALL2',
'TMEM82',
'ASPH',
'TSPAN12',
'MARCKSL1',
'RND2',
'PROM2',
'GM2A',
'ATG3',
'FSCN1',
'ADD1',
'CNN3',
'MAN2B1',
'CLSTN1',
'B4GALNT1',
'KREMEN1',
'SIDT2',
'SEC61A1',
'SEPT9',
'SOX10',
'RUSC2',
'PRSS36',
'CREB3L4',
'RGCC',
'FGF9',
'PLBD1',
'MPZL1',
'TMPRSS6',
'RHOC',
'DGKB',
'SYNE2',
'CEMIP',
'F12',
'OTOGL',
'PTPRK',
'GRB14',
'PLLP',
'CCDC69',
'SYPL',
'CES1D',
'LRP5',
'NDRG2',
'TMEM229B',
'ARHGAP44',
'SPINT2',
'RASSF6',
'RGS6',
'EMID1',
'PAQR6',
'NIPAL2',
'ST3GAL6',
'VIM',
'TUBA1B',
'WLS',
'IFI35',
'KSR1',
'ODF3L2',
'FASN',
'BCL2L14',
'EGFR',
'FSTL3',
'MAL',
'SCARA3',
'CUTAL',
'FARP1',
'CMTM8',
'DRC1',
'MFI2',
'LAPTM4A',
'GULP1',
'PLPP2',
'EPHX4',
'TRPV4',
'RCN2',
'INPP5J',
'RASGRP1',
'PPP1R17',
'FKBP9',
'ACTN4',
'PLD5',
'RPS6KA1',
'TSPAN14',
'NENF',
'GAL3ST4',
'MPZL2',
'CD82',
'VANGL1',
'PLEKHB1',
'ABHD14B',
'GPC2',
'PDZD2',
'TUBB4B',
'PIEZO1',
'MXRA8',
'FBLN5',
'SLC44A2',
'SCPEP1',
'MARVELD3',
'CPQ',
'UGDH',
'NPC2',
'CASKIN2',
'PLBD2',
'ENHO',
'BMP6',
'ARHGAP35',
'CCSAP',
'CD99L2',
'RABGAP1L',
'SEL1L',
'EHD3',
'FAM174B',
'PISD',
'RAB31',
'AEBP1',
'RASAL1',
'ILDR1',
'SLC22A17',
'ITPR2',
'UMODL1',
'P2RX2',
'SLC22A15',
'GCA',
'SKP1A',
'PREX1',
'S1PR2',
'2310030G06RIK',
'MBNL2',
'TSPAN9',
'C530008M17RIK',
'AMIGO1',
'RAMP3',
'FBXO30',
'LAD1',
'TWF1',
'NBEAL2',
'FMN1',
'FAM114A1',
'STARD13',
'ELOVL1',
'LAPTM4B',
'SOCS2',
'CPNE8',
'BDH2',
'MCTP1',
'SLITRK6',
'TAPBP',
'SORL1',
'HEY2',
'GM13889',
'ZFP467',
'RORC',
'RASSF4',
'ERGIC1',
'NDST1',
'CD164',
'ACSS3',
'EPCAM',
'SH3BGRL3',
'ITPR1',
'GPX8',
'CPNE2',
'PRDX6',
'LIPH',
'MLLT3',
'BCAM',
'DENND2D',
'RAC1',
'KIF2A',
'FAM107A',
'G6PDX',
'TSPAN7',
'ATP1B3',
'TGFA',
'EPB41L1',
'PVRL3',
'FDFT1',
'SLC9A3R1',
'ITGB8',
'FAM234A',
'AIM1',
'KCNA5',
'RAP1GAP2',
'CTNND1',
'MED10',
'ADAM10',
'HEY1',
'SLC44A3',
'CYP2J12',
'COL9A3',
'UAP1L1',
'ESRP1',
'SOCS6',
'MAP2K1',
'HERC4',
'SLC35A5',
'FBXO2',
'KIF1C',
'B4GALNT3',
'FUCA1',
'SCAMP4',
'GSTM5',
'2610020H08RIK',
'PON2',
'VEPH1',
'PDE4B',
'PRR15L',
'PTX3',
'DEGS2',
'SARAF',
'SLC46A1',
'CDC42EP1',
'NCSTN',
'MGAT5',
'CAP1',
'CYB561A3',
'SEPP1',
'GIPC2',
'PRDX1',
'MGLL',
'TMEM171',
'ANK2',
'NEU1',
'ATRAID',
'EPHB6',
'PAQR4',
'CRYAB',
'SLC27A1',
'PTPRG',
'FBLIM1',
'ARSB',
'SFXN5',
'UTRN',
'PGK1',
'SLC20A1',
'PCP4L1',
'SMIM3',
'SOX9',
'PTPRZ1',
'MASP1',
'CD151',
'TMEM106C',
'TOR4A',
'ACACA',
'HEYL',
'CAST',
'TMC6',
'E2F5',
'SSR2',
'EFEMP2',
'LAMB2',
'AKAP13',
'RNF13',
'FLRT2',
'ENPP1',
'LGR5',
'IL4RA',
'MAGED1',
'GALNT1',
'ARRB1',
'1700047I17RIK2',
'PLPP3',
'SYT13',
'SMS',
'TUBB2A',
'GAK',
'PDIA4',
'TKT',
'PLCD4',
'SNN',
'IL17RE',
'GNAS',
'ITFG1',
'TECPR2',
'TPD52',
'NUPR1',
'DUSP1',
'TCN2',
'SRGAP2',
'SAPCD2',
'UNC13B',
'FAH',
'LZIC',
'CYB5R3',
'SIN3B',
'PRRT3',
'STX7',
'IGSF8',
'SYMPK',
'BMP7',
'H2AFV',
'NRCAM',
'SPECC1',
'PACS1',
'HDC',
'FZD4',
'CAMTA1',
'TMCO1',
'FUT2',
'LIMA1',
'TNS1',
'SCARB1',
'PPT1',
'CD320',
'PSAP',
'TESC',
'TENM4',
'SHROOM1',
'CABP7',
'ZFP185',
'COL27A1',
'TNPO1',
'HOOK3',
'GNG12',
'FAM206A',
'COPS7A',
'SSX2IP',
'CD24A',
'GLRX3',
'SACS',
'ARHGEF5',
'CSRP1',
'SCD2',
'LLGL2',
'SNRPN',
'SLC16A6',
'PBXIP1',
'PLEKHH1',
'AGTRAP',
'ARHGAP32',
'GNA14',
'LAMP2',
'TESK2',
'DHRS3',
'RNF139',
'SLCO2A1',
'FAM110B',
'MGAT4A',
'CHST1',
'ADAMTS1',
'NUDT4',
'FUT4',
'FAM129B',
'CLDN23',
'ERGIC3',
'CCSER1',
'PLA2G12A',
'USB1',
'KIF21A',
'BSG',
'HHATL',
'ASCC2',
'9530077C05RIK',
'PLPP1',
'MSMO1',
'FBN2',
'FGFR1',
'ANXA5',
'CTSZ',
'RASAL2',
'VKORC1',
'LMTK2',
'PCYT2',
'GSTK1',
'PALM',
'1600002H07RIK',
'MED14',
'KIF19A',
'HS6ST2',
'SCRN1',
'MAT2A',
'RLBP1',
'SCAMP1',
'TTLL3',
'SLC39A8',
'PARVA',
'SUSD3',
'EDEM2',
'RNASE1',
'RPL10',
'DERL1',
'TRIM62',
'LPAR1',
'S100B',
'PDE4DIP',
'MTCH1',
'EBF1',
'PDE4D',
'PPIB',
'ZDHHC12',
'MVB12A',
'RPS14',
'RGMA',
'SHROOM3',
'CDK16',
'ADCY3',
'GDPD5',
'ALDH9A1',
'DOCK7',
'TMEM51',
'GAPDH',
'BZW1',
'TIMP3',
'FKBP4',
'HIST1H2BC',
'MAN2A2',
'TFCP2L1',
'ANKRD9',
'SPG21',
'SLC2A10',
'PIP4K2C',
'MAPK9',
'SULF1',
'MRPL51',
'PRSS8',
'ABHD4',
'POR',
'ERP29',
'APP',
'APC',
'SLC35B2',
'SDF4',
'NCKAP5',
'LRPAP1',
'GM10094',
'ELMO1',
'LONP2',
'1110004E09RIK',
'RHPN2',
'FAAH',
'DGKA',
'TM2D3',
'SNX21',
'LURAP1',
'GAA',
'LAMP1',
'RABAC1',
'MACF1',
'PTGDS',
'GPX4',
'AKR1A1',
'PARK7',
'LAMTOR2',
'F11R',
'STK38',
'GPS1',
'TMEM94',
'GAB1',
'SLC9A9',
'ITIH5',
'INPP5B',
'D630045J12RIK',
'CLMN',
'CLDN9',
'ACY1',
'SPATA2L',
'SEPHS2',
'ATP11A',
'JKAMP',
'PIGS',
'FRY',
'FILIP1',
'CDC42',
'STK33',
'SOD2',
'KCTD10',
'ITPKB',
'BTF3',
'DCT',
'HBP1',
'D930048N14RIK',
'HBB-BS',
'HBA-A1',
'NIP7']
ihcgene = ['OTOF',
'ATP2A3',
'TPBGL',
'DNAJC5B',
'CALB2',
'RPS6KA2',
'GSN',
'RP23-354H24.9',
'SLC7A14',
'BCL2',
'BIN1',
'KIFAP3',
'FAM65C',
'MARCH4',
'ATP1B1',
'CALM2',
'FCRLB',
'SPOCK1',
'CALM1',
'ARVCF',
'PHGDH',
'PPP2R2B',
'KCNAB1',
'IGFBP5',
'PDCD6',
'RELT',
'REEP6',
'PRKD1',
'ARMCX2',
'OSBPL6',
'SLC17A8',
'TCTEX1D1',
'ORMDL3',
'MRPS6',
'ARL5A',
'MIA',
'NEFL',
'FTL1',
'QDPR',
'PLTP',
'NMNAT2',
'EPHX1',
'AP3S1',
'GAREM',
'DNM1',
'WDR7',
'CABP2',
'MYL4',
'OTOL1',
'TSPAN13',
'SLC5A3',
'MYO3A',
'PLCXD2',
'PRKRA',
'CYR61',
'VIMP',
'PRDX2',
'LDLRAD3',
'BMPER',
'SFRP1',
'TBX2',
'SLC1A3',
'KCNA10',
'SYNRG',
'ANXA2',
'ANKRD13D',
'TFG',
'COA3',
'ZYX',
'SLC24A3',
'KIF13B',
'0610011F06RIK',
'PARP6',
'PPP2R2C',
'TRP53INP2',
'ANO1',
'TUBG1',
'RARRES1',
'ISOC1',
'SSFA2',
'CLN6',
'PPM1H',
'OSBPL11',
'EDNRB',
'CDKL4',
'ARPC4',
'TMOD2',
'TLN2',
'PLEKHB2',
'FTH1',
'PEBP1',
'SCG5',
'FAT4',
'HSP90AA1',
'COPZ1',
'ANO4',
'SMDT1',
'TGFB2',
'PGM2L1',
'PCBP3',
'ENDOD1',
'MYL6',
'GSTA4',
'PFN2',
'COL18A1',
'MADD',
'FAM32A',
'TJAP1',
'ATP2B1',
'TMEM132A',
'PVALB',
'ARF5',
'PTPRQ',
'LTBP1',
'FECH',
'FAM213B',
'GM1113',
'PPME1',
'CDO1',
'IDH1',
'CPXM2',
'SDCBP',
'ATP1A1',
'NTRK2',
'GAPDH',
'AMOT',
'OGDHL',
'RLBP1',
'NFASC',
'NOTCH1',
'CST6',
'GCNT4',
'NREP',
'ABHD12',
'HSDL1',
'SPATA6',
'TMSB4X',
'AMPH',
'ATP6V1G1',
'CBX5',
'SRGAP3',
'TPI1',
'GM13394',
'LIFR',
'OSBPL1A',
'SIAE',
'ZCCHC18',
'LYZ2',
'FAT1',
'FAM132A',
'FAM69A',
'ATPAF1',
'LGALS3',
'ANXA5',
'NOL4L',
'FXYD6',
'LRRC52',
'TOX2',
'MT-ND5',
'CERS5',
'RDH10',
'CHCHD2',
'CRIP2',
'ODC1',
'TRIM45',
'GM20594',
'GNB1',
'LUM',
'ILK',
'HDAC3',
'KNCN',
'KIF5C',
'CKAP4',
'SERP1',
'EXTL2',
'TMEM216',
'ANPEP',
'BDNF',
'CACFD1',
'PGAM1',
'DYNC2H1',
'CLTA',
'BCAM',
'DNPEP',
'TOM1L2',
'CLRN2',
'FDPS',
'ALDOA',
'WNT7B',
'ZC3H7B',
'ACTR10',
'MLF1',
'POLR2C',
'COL2A1',
'DIMT1',
'CHRNB2',
'GRAMD1B',
'ME1',
'CCT8',
'DPYSL2',
'YWHAG',
'MORN4',
'PARVA',
'SMARCAL1',
'DGKG',
'APBA1',
'TMEM53',
'GUK1',
'TECTB',
'ZYG11B',
'GPX2',
'MTUS1',
'ARPC1B',
'CDKN2D',
'WNK4',
'S100A8',
'CUBN',
'CD164L2',
'GALK1',
'AP3S2',
'SCN3B',
'UCP2',
'RCAN3',
'LOX',
'VPS39',
'SLC25A25',
'DENND4A',
'YEATS4',
'LARS2',
'TSTD3',
'MT-CO3',
'VPS29',
'NGP',
'COCH',
'ANO3',
'GJB2',
'S100A9',
'SRP72',
'UFC1',
'EPYC',
'TMEM117',
'GOT2',
'GORASP1',
'RPL7',
'VPS41',
'ACADVL',
'GM15662',
'GM10925',
'CAR3',
'NRSN1',
'OBOX3',
'GRWD1',
'GANC',
'SH3GL2',
'ATG12',
'GRIP2',
'CRY2',
'GM28437',
'UBE2W',
'CLTB',
'IQCB1',
'GM14667',
'GM19085',
'GDAP1',
'IPMK',
'HSPBP1',
'ABLIM1',
'EPB41L2',
'CNEP1R1',
'UBE2O',
'D6WSU163E',
'ZFP958',
'TASP1',
'CRTC1',
'AC152063.2',
'OAZ3',
'SASS6',
'MLXIPL',
'GALNTL6',
'RABGAP1',
'ASTN2',
'GM17374',
'MBNL1',
'ANKHD1',
'LONP1']
ihcauc = [0.998,
0.997,
0.99,
0.99,
0.989,
0.986,
0.983,
0.977,
0.965,
0.96,
0.959,
0.947,
0.943,
0.943,
0.936,
0.934,
0.933,
0.923,
0.922,
0.918,
0.917,
0.913,
0.913,
0.912,
0.912,
0.911,
0.908,
0.906,
0.906,
0.903,
0.902,
0.896,
0.895,
0.888,
0.882,
0.881,
0.881,
0.88,
0.877,
0.875,
0.874,
0.873,
0.872,
0.872,
0.868,
0.867,
0.862,
0.859,
0.858,
0.857,
0.855,
0.855,
0.854,
0.851,
0.849,
0.848,
0.846,
0.844,
0.842,
0.84,
0.84,
0.837,
0.835,
0.834,
0.834,
0.833,
0.832,
0.83,
0.83,
0.829,
0.826,
0.826,
0.821,
0.821,
0.819,
0.818,
0.816,
0.816,
0.815,
0.815,
0.814,
0.813,
0.81,
0.81,
0.808,
0.806,
0.806,
0.806,
0.805,
0.804,
0.804,
0.802,
0.801,
0.8,
0.8,
0.799,
0.797,
0.796,
0.795,
0.794,
0.793,
0.793,
0.793,
0.791,
0.79,
0.788,
0.788,
0.788,
0.787,
0.785,
0.785,
0.784,
0.782,
0.781,
0.781,
0.78,
0.779,
0.779,
0.778,
0.777,
0.776,
0.776,
0.776,
0.776,
0.776,
0.775,
0.775,
0.775,
0.774,
0.773,
0.773,
0.771,
0.771,
0.771,
0.771,
0.77,
0.77,
0.769,
0.769,
0.768,
0.767,
0.767,
0.764,
0.762,
0.762,
0.761,
0.761,
0.76,
0.76,
0.76,
0.76,
0.76,
0.759,
0.758,
0.758,
0.755,
0.755,
0.754,
0.754,
0.753,
0.751,
0.75,
0.749,
0.748,
0.746,
0.746,
0.746,
0.745,
0.745,
0.743,
0.743,
0.742,
0.741,
0.741,
0.738,
0.738,
0.738,
0.737,
0.735,
0.735,
0.735,
0.734,
0.733,
0.733,
0.732,
0.731,
0.73,
0.73,
0.729,
0.725,
0.725,
0.724,
0.724,
0.723,
0.722,
0.72,
0.719,
0.719,
0.719,
0.718,
0.718,
0.718,
0.713,
0.712,
0.712,
0.712,
0.711,
0.707,
0.706,
0.704,
0.703,
0.702,
0.699,
0.697,
0.696,
0.695,
0.694,
0.694,
0.692,
0.691,
0.689,
0.688,
0.687,
0.685,
0.685,
0.684,
0.683,
0.682,
0.681,
0.68,
0.677,
0.676,
0.675,
0.674,
0.674,
0.67,
0.668,
0.668,
0.666,
0.665,
0.663,
0.661,
0.65,
0.648,
0.648,
0.648,
0.647,
0.646,
0.645,
0.641,
0.36,
0.639,
0.632,
0.623,
0.621,
0.614,
0.612,
0.602,
0.599,
0.596,
0.406,
0.409,
0.41,
0.589,
0.412,
0.588,
0.587,
0.584,
0.582,
0.578,
0.425,
0.563,
0.439,
0.56,
0.558,
0.444,
0.446,
0.545,
0.459,
0.54,
0.462,
0.482,
0.517,
0.515,
0.513]
ohcgene = ['OCM',
'SRI',
'SLC26A5',
'STRC',
'LBH',
'MMD',
'CHRNA10',
'SDR42E2',
'CAR7',
'LPIN2',
'STRIP2',
'DNM3',
'LMO7',
'OLFM1',
'AK1',
'AI593442',
'DNER',
'CHRNA9',
'CACNA1D',
'CABLES1',
'ELOVL2',
'LRTM2',
'HS3ST3B1',
'PLCE1',
'MKRN2OS',
'ATP8A2',
'ATP2B2',
'SLC45A3',
'FAIM2',
'ABLIM2',
'GM6558',
'TTC7',
'USP46',
'EFCAB14',
'GPR156',
'KCNN2',
'TMC1',
'CACHD1',
'DISC1',
'ESPN',
'HOMER2',
'ST8SIA3',
'HSD17B7',
'EPN3',
'GALNT13',
'CPNE9',
'A530016L24RIK',
'LMOD3',
'LGALSL',
'STRBP',
'ME2',
'TRIOBP',
'RABGEF1',
'GPAM',
'NUP210',
'AQP11',
'LOXHD1',
'4930558C23RIK',
'BMP3',
'MYO15',
'FRMD5',
'CGN',
'B230219D22RIK',
'AP2A2',
'SPNS3',
'EFHD2',
'MYOM1',
'P2RX4',
'IKZF2',
'RAB3IP',
'C1QL1',
'KLHDC8B',
'TRIM16',
'NHLRC2',
'4833439L19RIK',
'RIMS3',
'TMIE',
'CIB2',
'TUB',
'LIMD2',
'ST6GALNAC6',
'MLXIP',
'LHFPL5',
'CTNNAL1',
'PRDX5',
'GM1673',
'SACM1L',
'ZBTB7A',
'CDH23',
'SPTB',
'SH3GL3',
'A930033H14RIK',
'FAM57B',
'LBR',
'ABCA5',
'RBM24',
'PCP4',
'DNAJA4',
'OCIAD2',
'YPEL3',
'PPP1R15A',
'TMEM30B',
'AFTPH',
'MAGI3',
'SLAIN2',
'CDC14A',
'DNAH8',
'MPPED1',
'REEP5',
'ATP8A1',
'TUBA8',
'GFOD1',
'AQP5',
'SEC16B',
'PADI2',
'ACSL1',
'BTBD9',
'MPDZ',
'SYNE4',
'1700025G04RIK',
'CARF',
'MCOLN3',
'FBXL20',
'TMEM38A',
'ELMOD1',
'NF2',
'HERC3',
'MICAL2',
'FGD2',
'SRD5A1',
'TMEM254C',
'STK39',
'IQSEC1',
'PRMT3',
'TMEM145',
'SIPA1L3',
'ZFP365',
'CDV3',
'PPP1CB',
'ZCCHC17',
'CRIP3',
'GRXCR1',
'DST',
'TTC24',
'SEPW1',
'CEP164',
'ABCA2',
'KCP',
'ANKFN1',
'LRRC8B',
'NDUFA8',
'MBOAT7',
'CALML4',
'UBL3',
'TMEM183A',
'CKMT1',
'USP48',
'SLC35E4',
'MEGF6',
'SH3RF3',
'KCNQ4',
'ACER2',
'PLEKHG6',
'PLSCR5',
'ATP8B1',
'PROB1',
'SLC25A42',
'PITPNA',
'TBC1D23',
'PKIG',
'DUSP8',
'SLK',
'AGPAT2',
'EPS8L2',
'TMEM170B',
'PTPN3',
'CHST2',
'MYO1H',
'CDCP2',
'TMEM254B',
'TMEM254A',
'KIF13A',
'POU4F3',
'TMEM255B',
'CCDC186',
'SH2D4A',
'ATF6',
'DYTN',
'NECAB2',
'SLC38A1',
'KALRN',
'MTMR11',
'SH2D4B',
'SENP2',
'SGPP1',
'CPEB3',
'PIK3CD',
'ANKRD24',
'SPOP',
'FSCN2',
'ACSL4',
'DUSP14',
'MYO6',
'GSE1',
'CACNA2D4',
'ERICH3',
'TOMT',
'OSBP2',
'ADARB1',
'SALL1',
'DONSON',
'PROSER2',
'PSMD8',
'PLS1',
'RTN2',
'FAM78A',
'EFHD1',
'ST13',
'MMP23',
'RNF11',
'KLHDC7B',
'DBNDD1',
'SORBS2',
'GRXCR2',
'ATP13A1',
'ABCB9',
'MFSD4',
'PPARGC1B',
'NRP2',
'LHFPL4',
'LMOD1',
'SNAP91',
'GARS',
'OGT',
'CCDC74A',
'USP25',
'CAMK2B',
'CEP85L',
'ZFP618',
'PFKFB2',
'STRA6',
'RHOBTB1',
'CDK5RAP2',
'PKNOX2',
'NBEAL1',
'RGS12',
'SORBS1',
'B3GALT5',
'IRF6',
'PLCH2',
'BCL7A',
'TMTC4',
'MKNK1',
'LASP1',
'DLG2',
'CPE',
'ARHGAP8',
'SEMA4C',
'CLASP1',
'FAM20C',
'DPYSL5',
'CCDC68',
'OSBPL9',
'SLC9A2',
'SMYD2',
'RUNX1T1',
'ARID3B',
'GM1043',
'GIPC3',
'MOB2',
'PDCD11',
'PRKACB',
'IGF2BP3',
'TSKU',
'CALB1',
'ADIPOR2',
'ATE1',
'FBXO9',
'FAM184A',
'ANKRD22',
'FAM104A',
'GM13420',
'CPT1B',
'D7ERTD443E',
'DFFA',
'UACA',
'ARV1',
'KLRK1',
'PYGO1',
'TPM1',
'PUM1',
'2210011C24RIK',
'TMEM191C',
'TMCC2',
'EIF4G2',
'STK24',
'LUC7L2',
'MFF',
'VWA3A',
'IGF1R',
'AMOTL2',
'GSK3B',
'GFI1',
'WDR12',
'RELL1',
'TMF1',
'MARS',
'JUP',
'CERS1',
'NRD1',
'CRELD2',
'TBC1D7',
'RBFOX2',
'SIX4',
'KCNS3',
'GPR152',
'CARD19',
'INPP1',
'D11WSU47E',
'PRKAA2',
'YARS',
'MED1',
'SLC7A6',
'ZFP532',
'HYOU1',
'REEP3',
'TUNAR',
'AP3B2',
'PPP2R5B',
'DUSP27',
'BCAR1',
'PHTF2',
'LHX3',
'ZBTB38',
'WRAP73',
'OTUD7B',
'RNPS1',
'UBE2V1',
'MAML2',
'SUGP2',
'MED25',
'ARHGEF28',
'GOLGA4',
'UBAC1',
'FAM234B',
'PRKAB1',
'INTU',
'CNOT3',
'MTPN',
'ERMP1',
'SIK2',
'PLPP7',
'MYCBP2',
'CEL',
'LRPPRC',
'UPF2',
'E130308A19RIK',
'PLEKHA7',
'PXDC1',
'EFR3A',
'TXNDC15',
'FGGY',
'AIF1L',
'EFTUD1',
'FBXL5',
'ANP32B',
'MAU2',
'IFIT1',
'WDR48',
'PSMD6',
'ACSS2',
'CLTC',
'NUDCD3',
'ZRANB2',
'AKAP9',
'PDZD7',
'BPGM',
'ADAMTSL1',
'PRKCA',
'TMTC1',
'CD164L2',
'TRAK1',
'GM21092',
'TK2',
'HSPA4L',
'NPTN',
'SUN1',
'PREP',
'GORASP2',
'BRAF',
'FAM217B',
'NSF',
'HPS5',
'ZGPAT',
'GABARAPL1',
'LRIG2',
'DCLRE1A',
'RIMKLB',
'DUSP10',
'AK2',
'CDC42BPA',
'GM42669',
'TNRC6B',
'CCZ1',
'ABCG1',
'FAM49A',
'KANK1',
'ANXA6',
'CCNL2',
'ARHGEF16',
'LINGO3',
'APBA1',
'SYNJ2',
'TMEM63B',
'MREG',
'KLHL21',
'TWF2',
'ZDHHC16',
'CCDC6',
'KMT2C',
'FAAP100',
'SH2B2',
'LIMK2',
'APLP1',
'SPOCK2',
'GABRB3',
'DNAJC6',
'CHKA',
'GTF3C5',
'LRRC20',
'OBFC1',
'PTTG1IP',
'ICMT',
'SBNO1',
'AKR1B3',
'ULK3',
'LRRC30',
'UNKL',
'DYM',
'CBFA2T2',
'ATP6V0A1',
'RAB19',
'RYR1',
'FCHO1',
'FAM214A',
'MATR3',
'CLRN1',
'AFAP1',
'ATRX',
'KNCN',
'STRADA',
'DTNB',
'CAR2',
'MYCBPAP',
'ZFP687',
'PLB1',
'WSCD2',
'ST8SIA2',
'ATP6V1B2',
'RNF150',
'TMEM64',
'HDGFRP2',
'MDM2',
'HERC2',
'NPTX1',
'HSPA2',
'CPSF6',
'SOAT1',
'DISP2',
'SPG20',
'FASTK',
'CYC1',
'ZKSCAN3',
'ABCF2',
'CRAMP1L',
'CTBP2',
'OSBPL3',
'KLHDC3',
'AHSA2',
'UBXN2A',
'SHANK2',
'MYO7A',
'IKBKB',
'RBM5',
'DCTN1',
'NSMAF',
'LRP8',
'TMEM201',
'TMEM2',
'BAG2',
'OSER1',
'ABLIM3',
'USP19',
'COPS8',
'LMAN1',
'B230217C12RIK',
'RB1CC1',
'HLCS',
'ZFP266',
'LANCL2',
'ATP6V1C2',
'TLDC1',
'WAPL',
'ROGDI',
'CXXC5',
'WWP2',
'WDR6',
'PPP4R4',
'TXNL4B',
'HMG20A',
'NCOA3',
'PCSK9',
'NOL9',
'NTNG2',
'DLAT',
'INSC',
'RALGAPB',
'ARTN',
'KCNMA1',
'LMBRD1',
'PXK',
'PALM3',
'RGS11',
'CRLF3',
'RAB28',
'ZFP444',
'SCYL1',
'GYS1',
'KRIT1',
'OPA3',
'4930402H24RIK',
'9330151L19RIK',
'UNC45A',
'FLNB',
'ZFP101',
'KCTD20',
'ZDHHC13',
'PAK3',
'DTNA',
'SLC1A5',
'ASIC1',
'GPD1L',
'DHODH',
'FAM124A',
'KSR2',
'GIGYF2',
'USP33',
'IGDCC4',
'EYA4',
'ARCN1',
'PEX26',
'DEPDC5',
'CDC23',
'PDLIM1',
'CTBP1',
'LRRC8C',
'BBS1',
'ANO10',
'UGGT1',
'DHX30',
'TET3',
'LAMA5',
'APPL2',
'RBM33',
'ISL1',
'NEBL',
'MGAT4B',
'GM42927',
'OGDH',
'RAPGEF1',
'PLAT',
'MBOAT2',
'MUC15',
'BTN2A2',
'ACIN1',
'TMEM19',
'UBP1',
'PAKAP',
'DHX38',
'FARP2',
'ASPHD2',
'DNM2',
'OCLN',
'TMEM87A',
'PLOD3',
'SALL3',
'LMAN2L',
'KIT',
'WDR1',
'COG3',
'MUM1',
'MIA3',
'ANK',
'MKRN1',
'UAP1',
'STIP1',
'VPS11',
'WBP1L',
'RALGPS1',
'NR1H2',
'FUBP1',
'ZFYVE19',
'SCAP',
'USP2',
'KLHL11',
'STK35',
'ENTPD5',
'MKRN2',
'SF3B3',
'ARHGAP17',
'PMM2',
'SLC25A12',
'SLC22A23',
'ABCC5',
'AIFM1',
'NHLRC1',
'NOMO1',
'IGHM',
'BORCS7',
'UBTD2',
'PPP6R3',
'NUP153',
'SLC37A3',
'NKTR',
'STAMBP',
'GATSL2',
'RNF180',
'KCNH8',
'EPM2AIP1',
'FAF1',
'DOLPP1',
'NECAP1',
'CMTR1',
'MARVELD2',
'ZFP346',
'SPINT1',
'STAT3',
'INPP5E',
'CCP110',
'PDRG1',
'DDI2',
'MSRB3',
'HIPK3',
'TOP2B',
'IDH3A',
'TTC14',
'PHIP',
'PIK3R2',
'FADS2',
'UBR4',
'FBXO42',
'RNF123',
'EIF4E',
'WAC',
'GART',
'RABEP2',
'JMJD1C',
'ST8SIA5',
'ANAPC5',
'ACVR2A',
'PDS5A',
'URGCP',
'PPP4C',
'CABIN1',
'NUP214',
'HK1',
'UBAP2',
'GPCPD1',
'CYB561D2',
'FLYWCH1',
'EIF3J2',
'GLE1',
'LAPTM5',
'GATC',
'BFAR',
'RAF1',
'USP20',
'GCC1',
'VPS13B',
'PTP4A2',
'ATF7IP',
'FOXK1',
'MLYCD',
'LIMD1',
'ATG4D',
'CRKL',
'SLC37A4',
'PTPN11',
'TRP53BP1',
'GM5464',
'ZMIZ2',
'FLCN',
'CCDC157',
'RP24-199E5.2',
'KDM5C',
'PRDM4',
'PLCB3',
'ANKRD12',
'MAST3',
'TAF1C',
'FBXO28',
'SH3BP1',
'FLII',
'DKK3',
'COG1']
ohcauc = [0.999,
0.994,
0.993,
0.99,
0.989,
0.986,
0.983,
0.98,
0.978,
0.977,
0.976,
0.976,
0.976,
0.976,
0.975,
0.972,
0.97,
0.967,
0.963,
0.958,
0.956,
0.951,
0.95,
0.949,
0.947,
0.944,
0.943,
0.942,
0.942,
0.942,
0.939,
0.937,
0.935,
0.934,
0.932,
0.93,
0.928,
0.926,
0.925,
0.923,
0.923,
0.92,
0.919,
0.917,
0.916,
0.916,
0.915,
0.914,
0.914,
0.913,
0.911,
0.911,
0.91,
0.91,
0.91,
0.91,
0.91,
0.906,
0.905,
0.904,
0.904,
0.902,
0.902,
0.902,
0.901,
0.899,
0.899,
0.898,
0.898,
0.895,
0.895,
0.895,
0.893,
0.893,
0.893,
0.892,
0.892,
0.892,
0.891,
0.89,
0.888,
0.887,
0.887,
0.887,
0.886,
0.885,
0.884,
0.884,
0.882,
0.882,
0.881,
0.881,
0.881,
0.88,
0.879,
0.879,
0.876,
0.874,
0.874,
0.874,
0.873,
0.871,
0.871,
0.87,
0.87,
0.869,
0.868,
0.865,
0.865,
0.865,
0.864,
0.864,
0.863,
0.862,
0.861,
0.861,
0.861,
0.86,
0.86,
0.859,
0.859,
0.858,
0.857,
0.856,
0.855,
0.855,
0.855,
0.854,
0.854,
0.854,
0.854,
0.853,
0.852,
0.852,
0.851,
0.85,
0.85,
0.85,
0.849,
0.849,
0.848,
0.848,
0.847,
0.847,
0.847,
0.845,
0.845,
0.845,
0.844,
0.844,
0.844,
0.842,
0.841,
0.841,
0.84,
0.84,
0.839,
0.839,
0.838,
0.838,
0.838,
0.837,
0.837,
0.837,
0.836,
0.836,
0.836,
0.835,
0.835,
0.835,
0.834,
0.834,
0.833,
0.831,
0.831,
0.831,
0.83,
0.83,
0.83,
0.83,
0.83,
0.829,
0.829,
0.829,
0.829,
0.828,
0.828,
0.828,
0.828,
0.827,
0.827,
0.826,
0.826,
0.826,
0.825,
0.825,
0.825,
0.824,
0.824,
0.824,
0.824,
0.822,
0.821,
0.821,
0.821,
0.82,
0.82,
0.819,
0.819,
0.819,
0.817,
0.817,
0.817,
0.817,
0.817,
0.816,
0.816,
0.816,
0.815,
0.815,
0.815,
0.814,
0.814,
0.814,
0.813,
0.813,
0.812,
0.812,
0.811,
0.811,
0.81,
0.81,
0.809,
0.809,
0.809,
0.808,
0.808,
0.807,
0.807,
0.806,
0.805,
0.805,
0.805,
0.805,
0.805,
0.804,
0.804,
0.804,
0.802,
0.802,
0.802,
0.801,
0.8,
0.799,
0.799,
0.799,
0.799,
0.798,
0.796,
0.796,
0.796,
0.795,
0.794,
0.794,
0.794,
0.794,
0.793,
0.792,
0.792,
0.791,
0.791,
0.791,
0.791,
0.791,
0.79,
0.79,
0.789,
0.789,
0.789,
0.789,
0.789,
0.789,
0.787,
0.787,
0.787,
0.785,
0.785,
0.785,
0.784,
0.784,
0.784,
0.784,
0.783,
0.783,
0.783,
0.783,
0.781,
0.78,
0.78,
0.779,
0.779,
0.778,
0.776,
0.776,
0.775,
0.774,
0.774,
0.773,
0.772,
0.772,
0.772,
0.772,
0.772,
0.772,
0.771,
0.771,
0.771,
0.77,
0.77,
0.77,
0.769,
0.769,
0.767,
0.767,
0.767,
0.767,
0.767,
0.766,
0.766,
0.766,
0.766,
0.765,
0.765,
0.765,
0.764,
0.764,
0.764,
0.764,
0.764,
0.763,
0.763,
0.762,
0.762,
0.762,
0.761,
0.761,
0.761,
0.761,
0.761,
0.76,
0.76,
0.76,
0.759,
0.758,
0.758,
0.757,
0.757,
0.757,
0.754,
0.754,
0.754,
0.753,
0.752,
0.752,
0.752,
0.752,
0.751,
0.751,
0.751,
0.751,
0.751,
0.751,
0.75,
0.75,
0.75,
0.75,
0.749,
0.749,
0.749,
0.749,
0.749,
0.748,
0.748,
0.747,
0.747,
0.747,
0.746,
0.746,
0.746,
0.746,
0.745,
0.745,
0.744,
0.744,
0.743,
0.743,
0.743,
0.741,
0.741,
0.741,
0.741,
0.741,
0.741,
0.74,
0.74,
0.739,
0.739,
0.738,
0.738,
0.737,
0.737,
0.737,
0.736,
0.736,
0.736,
0.735,
0.735,
0.735,
0.733,
0.733,
0.733,
0.733,
0.732,
0.732,
0.732,
0.732,
0.731,
0.731,
0.731,
0.73,
0.73,
0.729,
0.729,
0.729,
0.728,
0.728,
0.728,
0.727,
0.726,
0.726,
0.726,
0.726,
0.725,
0.725,
0.724,
0.724,
0.723,
0.723,
0.723,
0.723,
0.722,
0.722,
0.722,
0.721,
0.721,
0.721,
0.721,
0.72,
0.72,
0.72,
0.719,
0.719,
0.718,
0.718,
0.717,
0.717,
0.717,
0.717,
0.717,
0.716,
0.716,
0.716,
0.714,
0.714,
0.714,
0.713,
0.713,
0.713,
0.713,
0.713,
0.713,
0.712,
0.712,
0.712,
0.712,
0.711,
0.711,
0.71,
0.71,
0.71,
0.71,
0.71,
0.709,
0.709,
0.709,
0.709,
0.709,
0.709,
0.709,
0.708,
0.708,
0.707,
0.706,
0.706,
0.706,
0.705,
0.704,
0.704,
0.704,
0.704,
0.703,
0.702,
0.702,
0.702,
0.7,
0.699,
0.699,
0.699,
0.698,
0.698,
0.698,
0.697,
0.696,
0.696,
0.695,
0.695,
0.695,
0.695,
0.695,
0.694,
0.694,
0.694,
0.694,
0.693,
0.693,
0.692,
0.692,
0.692,
0.691,
0.691,
0.691,
0.69,
0.69,
0.69,
0.69,
0.69,
0.69,
0.69,
0.69,
0.69,
0.689,
0.688,
0.688,
0.687,
0.687,
0.687,
0.687,
0.687,
0.686,
0.686,
0.686,
0.685,
0.685,
0.685,
0.685,
0.685,
0.685,
0.685,
0.683,
0.683,
0.683,
0.683,
0.682,
0.681,
0.681,
0.68,
0.68,
0.679,
0.679,
0.679,
0.679,
0.679,
0.678,
0.678,
0.677,
0.677,
0.677,
0.675,
0.675,
0.675,
0.675,
0.675,
0.674,
0.674,
0.674,
0.673,
0.673,
0.672,
0.671,
0.67,
0.67,
0.669,
0.667,
0.667,
0.667,
0.666,
0.666,
0.666,
0.666,
0.666,
0.666,
0.665,
0.663,
0.663,
0.663,
0.662,
0.662,
0.66,
0.66,
0.656,
0.655,
0.655,
0.655,
0.655,
0.654,
0.654,
0.653,
0.653,
0.652,
0.651,
0.646,
0.646,
0.643,
0.641,
0.641,
0.64,
0.638,
0.637,
0.636,
0.635,
0.635,
0.632,
0.631,
0.63,
0.629,
0.629,
0.628,
0.627,
0.627,
0.623,
0.622,
0.622,
0.622,
0.621,
0.62,
0.619,
0.619,
0.617,
0.616,
0.615,
0.615,
0.614,
0.614,
0.613,
0.612,
0.611,
0.61,
0.61,
0.608,
0.605,
0.605,
0.605,
0.604,
0.603,
0.602,
0.602,
0.598,
0.597,
0.594,
0.594,
0.59,
0.589,
0.578,
0.575,
0.575,
0.572,
0.571,
0.571,
0.57,
0.569,
0.567,
0.567,
0.566,
0.565,
0.564,
0.562,
0.556,
0.554,
0.553,
0.451,
0.546,
0.541,
0.534,
0.526,
0.514]
aucthr = [.5, .6, .7, .8, .9]
pvalarr = [.001, .0001, .00001, .000001, .0000001]
bb = connect(scidb_auth=('scidbadmin', 'Paradigm4'), namespace='RIVAS_HG19')
assocset = str(bb.list_association_sets()['name'][0])

def compute_enrichment_stats(resdf, geneset, aucarr, threshold):
    genesubset = [geneset[idx] for idx in range(0, len(geneset)) if float(aucarr[idx]) >= float(threshold)]
    yescnt = len(set([item[0] for item in filter(len,list(resdf.apply(lambda row: ([gene for gene in row['gene_symbol'].split(',') if gene in genesubset]), axis = 1)))]))
    listyes = list(set([item[0] for item in filter(len,list(resdf.apply(lambda row: ([gene for gene in row['gene_symbol'].split(',') if gene in genesubset]), axis = 1)))]))
    auclist = [aucarr[geneset.index(gene)] for gene in listyes] 
    yesnocnt = len(set(filter(len,resdf['gene_symbol']))) - yescnt
    noyes = len(genesubset) - yescnt
    nono = 20000 - len(genesubset) - yesnocnt - yescnt
    oddsratio, pvalue = stats.fisher_exact([[yescnt, yesnocnt], [noyes, nono]])
    return oddsratio, pvalue, ','.join(listyes), ','.join([str(i) for i in auclist])

for pval in pvalarr:
    res = lookups.get_icd_variant_by_icd_id_pvalue(bb, 'BIN_FC3002247', 2344, pvalue=pval)
    resdf = res.query('(category == "lof_variant" or category == "missense_variant")', inplace = False)
    for aucval in aucthr:
        resorihc, respihc, genesihc, aucarrihc = compute_enrichment_stats(resdf, ihcgene, ihcauc, aucval)
        resorohc, respohc, genesohc, aucarrohc = compute_enrichment_stats(resdf, ohcgene, ohcauc, aucval)
        resordc, respdc, genesdc, aucarrdc = compute_enrichment_stats(resdf, dcgene, dcauc, aucval)
        print -numpy.log10(pval), pval, aucval, resorihc, respihc,"Inner_Hair_Cells", genesihc , aucarrihc
        print -numpy.log10(pval), pval, aucval, resorohc, respohc,"Outer_Hair_Cells", genesohc , aucarrohc
        print -numpy.log10(pval), pval, aucval, resordc, respdc,"Deiters_Cells", genesdc, aucarrdc
















