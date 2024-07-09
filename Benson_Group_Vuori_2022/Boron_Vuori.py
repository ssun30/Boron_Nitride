#!/usr/bin/env python
# encoding: utf-8

name = "Boron_Vuori"
shortDesc = u"Boron compound thermo database calculated by Vuori et al. 2022"
longDesc = u"""
Calculations done by Vuori et al. 2022 using the W1X-1 compound method.

Vuori, H. T., Rautiainen, J. M., Kolehmainen, E. T., & Tuononen, H. M. (2022). Computational thermochemistry: extension of Benson group additivity approach to organoboron compounds and reliable predictions of their thermochemical properties. Dalton Transactions, 51(41), 15816-15829. https://doi.org/10.1039/D2DT02659G
"""
entry(
    index = 0,
    label = 'BH3',
    molecule =
"""
1 B u0 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([27.8664,31.8056,36.2,40.4417,47.7394,53.4,62.5454],'J/(mol*K)'),
        H298 = (103.9,'kJ/mol'),
        S298 = (188.2,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (83.1447,'J/(mol*K)'),
        E0 = (95.5691,'kJ/mol'),
    ),
    shortDesc = u"""BH3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 1,
    label = 'BH2Me',
    molecule =
"""
1 B u0 p0 c0 {2,S} {3,S} {4,S}
2 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([47.0054,57.4614,67.5,76.5681,91.4794,102.7,120.372],'J/(mol*K)'),
        H298 = (32.2,'kJ/mol'),
        S298 = (250,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (153.818,'J/(mol*K)'),
        E0 = (21.3713,'kJ/mol'),
    ),
    shortDesc = u"""CH5B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 2,
    label = 'BH2Et',
    molecule =
"""
1  B u0 p0 c0 {2,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
3  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([69.4302,85.6503,100.7,114.173,136.366,153.2,179.943],'J/(mol*K)'),
        H298 = (21.3,'kJ/mol'),
        S298 = (288,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (224.491,'J/(mol*K)'),
        E0 = (7.16778,'kJ/mol'),
    ),
    shortDesc = u"""C2H7B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 3,
    label = 'BH2Vi',
    molecule =
"""
1 B u0 p0 c0 {2,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {3,D} {6,S}
3 C u0 p0 c0 {2,D} {7,S} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([56.7183,72.25,86.1,97.8422,115.725,128.1,145.871],'J/(mol*K)'),
        H298 = (135.3,'kJ/mol'),
        S298 = (265.2,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (178.761,'J/(mol*K)'),
        E0 = (123.712,'kJ/mol'),
    ),
    shortDesc = u"""C2H5B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 4,
    label = 'BH2Pr',
    molecule =
"""
1  B u0 p0 c0 {2,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {11,S} {12,S} {13,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([91.8768,114.793,135.4,153.483,182.738,204.6,238.905],'J/(mol*K)'),
        H298 = (-0.7,'kJ/mol'),
        S298 = (322.1,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (295.164,'J/(mol*K)'),
        E0 = (-17.8822,'kJ/mol'),
    ),
    shortDesc = u"""C3H9B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 5,
    label = 'BH2iPr',
    molecule =
"""
1  B u0 p0 c0 {2,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
3  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([94.2894,117.592,138.2,156.042,184.5,205.5,238.161],'J/(mol*K)'),
        H298 = (3.2,'kJ/mol'),
        S298 = (312.9,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (295.164,'J/(mol*K)'),
        E0 = (-14.2885,'kJ/mol'),
    ),
    shortDesc = u"""C3H9B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 6,
    label = 'BH2sBu',
    molecule =
"""
1  B u0 p0 c0 {2,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {14,S} {15,S} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([115.665,147.015,174.2,197.299,233.294,259.2,298.439],'J/(mol*K)'),
        H298 = (-16.9,'kJ/mol'),
        S298 = (345.9,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (365.837,'J/(mol*K)'),
        E0 = (-36.9178,'kJ/mol'),
    ),
    shortDesc = u"""C4H11B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 7,
    label = 'BH2tBu',
    molecule =
"""
1  B u0 p0 c0 {2,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([119.163,150.264,177.1,199.872,235.384,261,299.923],'J/(mol*K)'),
        H298 = (-22.6,'kJ/mol'),
        S298 = (339.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (365.837,'J/(mol*K)'),
        E0 = (-43.3225,'kJ/mol'),
    ),
    shortDesc = u"""C4H11B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 8,
    label = 'BH2Ph',
    molecule =
"""
1  B u0 p0 c0 {2,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {7,D}
3  C u0 p0 c0 {2,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {5,S} {11,S}
5  C u0 p0 c0 {4,S} {6,D} {12,S}
6  C u0 p0 c0 {5,D} {7,S} {13,S}
7  C u0 p0 c0 {2,D} {6,S} {14,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([99.3239,133.668,163,187.058,222.293,245.6,277.175],'J/(mol*K)'),
        H298 = (160.9,'kJ/mol'),
        S298 = (313.9,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (328.422,'J/(mol*K)'),
        E0 = (144.822,'kJ/mol'),
    ),
    shortDesc = u"""C6H7B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 9,
    label = 'BHMe2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  B u0 p0 c0 {1,S} {3,S} {7,S}
3  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([68.3325,84.6727,99.8,113.292,135.408,152.1,178.529],'J/(mol*K)'),
        H298 = (-39.7,'kJ/mol'),
        S298 = (291,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (224.491,'J/(mol*K)'),
        E0 = (-53.5826,'kJ/mol'),
    ),
    shortDesc = u"""C2H7B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 10,
    label = 'BHEtMe',
    molecule =
"""
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  B u0 p0 c0 {1,S} {3,S} {8,S}
3  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {11,S} {12,S} {13,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([90.2501,112.079,132.1,150,179.601,202.2,238.345],'J/(mol*K)'),
        H298 = (-50.4,'kJ/mol'),
        S298 = (332.4,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (295.164,'J/(mol*K)'),
        E0 = (-67.5746,'kJ/mol'),
    ),
    shortDesc = u"""C3H9B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 11,
    label = 'BHMeVi',
    molecule =
"""
1  C u0 p0 c0 {2,D} {5,S} {6,S}
2  C u0 p0 c0 {1,D} {3,S} {7,S}
3  B u0 p0 c0 {2,S} {4,S} {8,S}
4  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([80.146,101.586,120.6,136.906,162.37,180.6,207.772],'J/(mol*K)'),
        H298 = (64,'kJ/mol'),
        S298 = (315.6,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (249.434,'J/(mol*K)'),
        E0 = (48.929,'kJ/mol'),
    ),
    shortDesc = u"""C3H7B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 12,
    label = 'BHMePh',
    molecule =
"""
1  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
2  B u0 p0 c0 {1,S} {3,S} {12,S}
3  C u0 p0 c0 {2,S} {4,S} {8,D}
4  C u0 p0 c0 {3,S} {5,D} {13,S}
5  C u0 p0 c0 {4,D} {6,S} {14,S}
6  C u0 p0 c0 {5,S} {7,D} {15,S}
7  C u0 p0 c0 {6,D} {8,S} {16,S}
8  C u0 p0 c0 {3,D} {7,S} {17,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([123.045,162.737,196.4,224.083,265.099,292.8,331.654],'J/(mol*K)'),
        H298 = (91.8,'kJ/mol'),
        S298 = (364.6,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (399.095,'J/(mol*K)'),
        E0 = (72.1756,'kJ/mol'),
    ),
    shortDesc = u"""C7H9B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 13,
    label = 'BHEt2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  B u0 p0 c0 {2,S} {4,S} {11,S}
4  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {14,S} {15,S} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([113.78,141.542,166.6,188.84,225.472,253.4,298.053],'J/(mol*K)'),
        H298 = (-60.4,'kJ/mol'),
        S298 = (354.9,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (365.837,'J/(mol*K)'),
        E0 = (-81.0291,'kJ/mol'),
    ),
    shortDesc = u"""C4H11B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 14,
    label = 'BHEtVi',
    molecule =
"""
1  C u0 p0 c0 {2,D} {6,S} {7,S}
2  C u0 p0 c0 {1,D} {3,S} {8,S}
3  B u0 p0 c0 {2,S} {4,S} {9,S}
4  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {12,S} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([100.668,127.703,151.5,171.931,204.082,227.4,262.845],'J/(mol*K)'),
        H298 = (53.7,'kJ/mol'),
        S298 = (345,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (320.107,'J/(mol*K)'),
        E0 = (35.6963,'kJ/mol'),
    ),
    shortDesc = u"""C4H9B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 15,
    label = 'BHVi2',
    molecule =
"""
1  C u0 p0 c0 {2,D} {6,S} {7,S}
2  C u0 p0 c0 {1,D} {3,S} {8,S}
3  B u0 p0 c0 {2,S} {4,S} {9,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {4,D} {11,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([90.7388,116.243,138.2,156.526,184.202,203.3,230.676],'J/(mol*K)'),
        H298 = (166.9,'kJ/mol'),
        S298 = (317.4,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (274.378,'J/(mol*K)'),
        E0 = (150.696,'kJ/mol'),
    ),
    shortDesc = u"""C4H7B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 16,
    label = 'BHPhVi',
    molecule =
"""
1  C u0 p0 c0 {2,D} {10,S} {11,S}
2  C u0 p0 c0 {1,D} {3,S} {12,S}
3  B u0 p0 c0 {2,S} {4,S} {13,S}
4  C u0 p0 c0 {3,S} {5,S} {9,D}
5  C u0 p0 c0 {4,S} {6,D} {14,S}
6  C u0 p0 c0 {5,D} {7,S} {15,S}
7  C u0 p0 c0 {6,S} {8,D} {16,S}
8  C u0 p0 c0 {7,D} {9,S} {17,S}
9  C u0 p0 c0 {4,D} {8,S} {18,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([133.024,176.417,213.2,243.447,288.232,318.4,360.337],'J/(mol*K)'),
        H298 = (195,'kJ/mol'),
        S298 = (377.3,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (424.038,'J/(mol*K)'),
        E0 = (174.148,'kJ/mol'),
    ),
    shortDesc = u"""C8H9B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 17,
    label = 'BHiPr2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {11,S}
3  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
4  B u0 p0 c0 {2,S} {5,S} {15,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {16,S}
6  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([162.367,203.041,238.6,269.476,319.315,356.7,415.776],'J/(mol*K)'),
        H298 = (-96.7,'kJ/mol'),
        S298 = (417,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (507.183,'J/(mol*K)'),
        E0 = (-124.198,'kJ/mol'),
    ),
    shortDesc = u"""C6H15B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 18,
    label = 'BHPh2',
    molecule =
"""
1  B u0 p0 c0 {2,S} {8,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {7,D}
3  C u0 p0 c0 {2,S} {4,D} {15,S}
4  C u0 p0 c0 {3,D} {5,S} {16,S}
5  C u0 p0 c0 {4,S} {6,D} {17,S}
6  C u0 p0 c0 {5,D} {7,S} {18,S}
7  C u0 p0 c0 {2,D} {6,S} {19,S}
8  C u0 p0 c0 {1,S} {9,S} {13,D}
9  C u0 p0 c0 {8,S} {10,D} {20,S}
10 C u0 p0 c0 {9,D} {11,S} {21,S}
11 C u0 p0 c0 {10,S} {12,D} {22,S}
12 C u0 p0 c0 {11,D} {13,S} {23,S}
13 C u0 p0 c0 {8,D} {12,S} {24,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {12,S}
24 H u0 p0 c0 {13,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([177.743,240.343,292.5,334.634,395.508,435.3,488.645],'J/(mol*K)'),
        H298 = (229.8,'kJ/mol'),
        S298 = (425.1,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (573.699,'J/(mol*K)'),
        E0 = (204.247,'kJ/mol'),
    ),
    shortDesc = u"""C12H11B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 19,
    label = 'BMe3',
    molecule =
"""
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  B u0 p0 c0 {1,S} {3,S} {4,S}
3  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([91.052,112.86,132.7,150.334,179.338,201.4,236.663],'J/(mol*K)'),
        H298 = (-109.6,'kJ/mol'),
        S298 = (357.5,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (295.164,'J/(mol*K)'),
        E0 = (-126.889,'kJ/mol'),
    ),
    shortDesc = u"""C3H9B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 20,
    label = 'BEtMe2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([114.57,141.817,166.4,188.261,224.406,252.1,296.662],'J/(mol*K)'),
        H298 = (-121.3,'kJ/mol'),
        S298 = (375.6,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (365.837,'J/(mol*K)'),
        E0 = (-142.178,'kJ/mol'),
    ),
    shortDesc = u"""C4H11B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 21,
    label = 'BMe2Vi',
    molecule =
"""
1  C u0 p0 c0 {2,D} {6,S} {7,S}
2  C u0 p0 c0 {1,D} {3,S} {8,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {12,S} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([103.535,128.843,151.1,170.344,201.069,223.8,259.271],'J/(mol*K)'),
        H298 = (-2.3,'kJ/mol'),
        S298 = (357.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (320.107,'J/(mol*K)'),
        E0 = (-21.1961,'kJ/mol'),
    ),
    shortDesc = u"""C4H9B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 22,
    label = 'BMe2Ph',
    molecule =
"""
1  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
2  B u0 p0 c0 {1,S} {3,S} {4,S}
3  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {5,S} {9,D}
5  C u0 p0 c0 {4,S} {6,D} {16,S}
6  C u0 p0 c0 {5,D} {7,S} {17,S}
7  C u0 p0 c0 {6,S} {8,D} {18,S}
8  C u0 p0 c0 {7,D} {9,S} {19,S}
9  C u0 p0 c0 {4,D} {8,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([147.337,191.197,228.6,259.825,307.339,340.6,389.361],'J/(mol*K)'),
        H298 = (26.1,'kJ/mol'),
        S298 = (404.4,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (469.768,'J/(mol*K)'),
        E0 = (2.4578,'kJ/mol'),
    ),
    shortDesc = u"""C8H11B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 23,
    label = 'BEtMePh',
    molecule =
"""
1  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {14,S} {15,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {3,S} {6,S} {10,D}
6  C u0 p0 c0 {5,S} {7,D} {19,S}
7  C u0 p0 c0 {6,D} {8,S} {20,S}
8  C u0 p0 c0 {7,S} {9,D} {21,S}
9  C u0 p0 c0 {8,D} {10,S} {22,S}
10 C u0 p0 c0 {5,D} {9,S} {23,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([165.268,215.423,258.5,294.742,350.458,389.9,448.344],'J/(mol*K)'),
        H298 = (16.1,'kJ/mol'),
        S298 = (459.1,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (540.441,'J/(mol*K)'),
        E0 = (-9.90865,'kJ/mol'),
    ),
    shortDesc = u"""C9H13B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 24,
    label = 'BMeVi2',
    molecule =
"""
1  C u0 p0 c0 {2,D} {7,S} {8,S}
2  C u0 p0 c0 {1,D} {3,S} {9,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {6,D} {13,S}
6  C u0 p0 c0 {5,D} {14,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([115.001,143.152,167.5,188.289,221.052,245,281.962],'J/(mol*K)'),
        H298 = (105.4,'kJ/mol'),
        S298 = (378.6,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (345.051,'J/(mol*K)'),
        E0 = (84.8999,'kJ/mol'),
    ),
    shortDesc = u"""C5H9B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 25,
    label = 'BMePhVi',
    molecule =
"""
1  C u0 p0 c0 {2,D} {11,S} {12,S}
2  C u0 p0 c0 {1,D} {3,S} {13,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {6,S} {10,D}
6  C u0 p0 c0 {5,S} {7,D} {17,S}
7  C u0 p0 c0 {6,D} {8,S} {18,S}
8  C u0 p0 c0 {7,S} {9,D} {19,S}
9  C u0 p0 c0 {8,D} {10,S} {20,S}
10 C u0 p0 c0 {5,D} {9,S} {21,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([160.006,206.767,246.2,278.849,328.104,362.3,412.058],'J/(mol*K)'),
        H298 = (133.7,'kJ/mol'),
        S298 = (426.2,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (494.711,'J/(mol*K)'),
        E0 = (108.259,'kJ/mol'),
    ),
    shortDesc = u"""C9H11B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 26,
    label = 'BMePh2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
2  B u0 p0 c0 {1,S} {3,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {8,D}
4  C u0 p0 c0 {3,S} {5,D} {18,S}
5  C u0 p0 c0 {4,D} {6,S} {19,S}
6  C u0 p0 c0 {5,S} {7,D} {20,S}
7  C u0 p0 c0 {6,D} {8,S} {21,S}
8  C u0 p0 c0 {3,D} {7,S} {22,S}
9  C u0 p0 c0 {2,S} {10,S} {14,D}
10 C u0 p0 c0 {9,S} {11,D} {23,S}
11 C u0 p0 c0 {10,D} {12,S} {24,S}
12 C u0 p0 c0 {11,S} {13,D} {25,S}
13 C u0 p0 c0 {12,D} {14,S} {26,S}
14 C u0 p0 c0 {9,D} {13,S} {27,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {12,S}
26 H u0 p0 c0 {13,S}
27 H u0 p0 c0 {14,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([202.254,269.871,326.4,372.481,440.183,485.5,548.212],'J/(mol*K)'),
        H298 = (163.3,'kJ/mol'),
        S298 = (469.9,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (644.372,'J/(mol*K)'),
        E0 = (133.882,'kJ/mol'),
    ),
    shortDesc = u"""C13H13B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 27,
    label = 'BF2Me',
    molecule =
"""
1 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2 B u0 p0 c0 {1,S} {3,S} {4,S}
3 F u0 p3 c0 {2,S}
4 F u0 p3 c0 {2,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([56.8576,69.1795,79.9,88.8736,102.432,111.8,125.389],'J/(mol*K)'),
        H298 = (-809.4,'kJ/mol'),
        S298 = (296.3,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (153.818,'J/(mol*K)'),
        E0 = (-821.551,'kJ/mol'),
    ),
    shortDesc = u"""CH3BF2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 28,
    label = 'BEt3',
    molecule =
"""
1  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  B u0 p0 c0 {2,S} {4,S} {6,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {15,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {7,S} {18,S} {19,S}
7  C u0 p0 c0 {6,S} {20,S} {21,S} {22,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([162.523,201.262,235.6,265.913,315.894,354.2,415.923],'J/(mol*K)'),
        H298 = (-143.3,'kJ/mol'),
        S298 = (423.4,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (507.183,'J/(mol*K)'),
        E0 = (-171.34,'kJ/mol'),
    ),
    shortDesc = u"""C6H15B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 29,
    label = 'BEtVi2',
    molecule =
"""
1  C u0 p0 c0 {2,D} {8,S} {9,S}
2  C u0 p0 c0 {1,D} {3,S} {10,S}
3  B u0 p0 c0 {2,S} {4,S} {6,S}
4  C u0 p0 c0 {3,S} {5,D} {11,S}
5  C u0 p0 c0 {4,D} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {7,S} {14,S} {15,S}
7  C u0 p0 c0 {6,S} {16,S} {17,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([138.858,174.246,204.7,230.617,271.306,300.9,346.16],'J/(mol*K)'),
        H298 = (96.3,'kJ/mol'),
        S298 = (392.7,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (415.724,'J/(mol*K)'),
        E0 = (72.5503,'kJ/mol'),
    ),
    shortDesc = u"""C6H11B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 30,
    label = 'BEtF2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  F u0 p3 c0 {3,S}
5  F u0 p3 c0 {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([79.076,96.9024,112.5,125.847,146.805,162,185.163],'J/(mol*K)'),
        H298 = (-819,'kJ/mol'),
        S298 = (327.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (224.491,'J/(mol*K)'),
        E0 = (-834.491,'kJ/mol'),
    ),
    shortDesc = u"""C2H5BF2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 31,
    label = 'BVi3',
    molecule =
"""
1  C u0 p0 c0 {2,D} {8,S} {9,S}
2  C u0 p0 c0 {1,D} {3,S} {10,S}
3  B u0 p0 c0 {2,S} {4,S} {6,S}
4  C u0 p0 c0 {3,S} {5,D} {11,S}
5  C u0 p0 c0 {4,D} {12,S} {13,S}
6  C u0 p0 c0 {3,S} {7,D} {14,S}
7  C u0 p0 c0 {6,D} {15,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([123.76,154.622,181.2,203.811,239.292,265.1,304.682],'J/(mol*K)'),
        H298 = (213.5,'kJ/mol'),
        S298 = (394.1,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (369.994,'J/(mol*K)'),
        E0 = (191.834,'kJ/mol'),
    ),
    shortDesc = u"""C6H9B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 32,
    label = 'BPhVi2',
    molecule =
"""
1  C u0 p0 c0 {2,D} {12,S} {13,S}
2  C u0 p0 c0 {1,D} {3,S} {14,S}
3  B u0 p0 c0 {2,S} {4,S} {6,S}
4  C u0 p0 c0 {3,S} {5,D} {15,S}
5  C u0 p0 c0 {4,D} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {7,S} {11,D}
7  C u0 p0 c0 {6,S} {8,D} {18,S}
8  C u0 p0 c0 {7,D} {9,S} {19,S}
9  C u0 p0 c0 {8,S} {10,D} {20,S}
10 C u0 p0 c0 {9,D} {11,S} {21,S}
11 C u0 p0 c0 {6,D} {10,S} {22,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {11,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([172.933,225.341,268.9,304.37,356.61,391.8,441.089],'J/(mol*K)'),
        H298 = (242.1,'kJ/mol'),
        S298 = (431.9,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (519.654,'J/(mol*K)'),
        E0 = (215.423,'kJ/mol'),
    ),
    shortDesc = u"""C10H11B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 33,
    label = 'BPh2Vi',
    molecule =
"""
1  C u0 p0 c0 {2,D} {16,S} {17,S}
2  C u0 p0 c0 {1,D} {3,S} {18,S}
3  B u0 p0 c0 {2,S} {4,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {9,D}
5  C u0 p0 c0 {4,S} {6,D} {19,S}
6  C u0 p0 c0 {5,D} {7,S} {20,S}
7  C u0 p0 c0 {6,S} {8,D} {21,S}
8  C u0 p0 c0 {7,D} {9,S} {22,S}
9  C u0 p0 c0 {4,D} {8,S} {23,S}
10 C u0 p0 c0 {3,S} {11,S} {15,D}
11 C u0 p0 c0 {10,S} {12,D} {24,S}
12 C u0 p0 c0 {11,D} {13,S} {25,S}
13 C u0 p0 c0 {12,S} {14,D} {26,S}
14 C u0 p0 c0 {13,D} {15,S} {27,S}
15 C u0 p0 c0 {10,D} {14,S} {28,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {12,S}
26 H u0 p0 c0 {13,S}
27 H u0 p0 c0 {14,S}
28 H u0 p0 c0 {15,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([215.142,286.693,346.3,394.753,465.674,512.9,577.651],'J/(mol*K)'),
        H298 = (268.4,'kJ/mol'),
        S298 = (480.5,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (669.315,'J/(mol*K)'),
        E0 = (237.294,'kJ/mol'),
    ),
    shortDesc = u"""C14H13B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 34,
    label = 'BF2Vi',
    molecule =
"""
1 C u0 p0 c0 {2,D} {6,S} {7,S}
2 C u0 p0 c0 {1,D} {3,S} {8,S}
3 B u0 p0 c0 {2,S} {4,S} {5,S}
4 F u0 p3 c0 {3,S}
5 F u0 p3 c0 {3,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([71.8492,88.0026,101.3,111.943,127.185,137.1,150.522],'J/(mol*K)'),
        H298 = (-695.8,'kJ/mol'),
        S298 = (306.7,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (178.761,'J/(mol*K)'),
        E0 = (-709.955,'kJ/mol'),
    ),
    shortDesc = u"""C2H3BF2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 35,
    label = 'BiPr3',
    molecule =
"""
1  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {14,S}
3  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
4  B u0 p0 c0 {2,S} {5,S} {8,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {18,S}
6  C u0 p0 c0 {5,S} {19,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {4,S} {9,S} {10,S} {25,S}
9  C u0 p0 c0 {8,S} {26,S} {27,S} {28,S}
10 C u0 p0 c0 {8,S} {29,S} {30,S} {31,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {10,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([236.317,297.218,349.2,393.508,463.674,515.4,595.889],'J/(mol*K)'),
        H298 = (-194.1,'kJ/mol'),
        S298 = (503.5,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (719.202,'J/(mol*K)'),
        E0 = (-231.926,'kJ/mol'),
    ),
    shortDesc = u"""C9H21B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 36,
    label = 'BBu3',
    molecule =
"""
1  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {17,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {19,S} {20,S}
4  C u0 p0 c0 {3,S} {5,S} {21,S} {22,S}
5  B u0 p0 c0 {4,S} {6,S} {10,S}
6  C u0 p0 c0 {5,S} {7,S} {23,S} {24,S}
7  C u0 p0 c0 {6,S} {8,S} {25,S} {26,S}
8  C u0 p0 c0 {7,S} {9,S} {27,S} {28,S}
9  C u0 p0 c0 {8,S} {29,S} {30,S} {31,S}
10 C u0 p0 c0 {5,S} {11,S} {32,S} {33,S}
11 C u0 p0 c0 {10,S} {12,S} {34,S} {35,S}
12 C u0 p0 c0 {11,S} {13,S} {36,S} {37,S}
13 C u0 p0 c0 {12,S} {38,S} {39,S} {40,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {4,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {9,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {12,S}
37 H u0 p0 c0 {12,S}
38 H u0 p0 c0 {13,S}
39 H u0 p0 c0 {13,S}
40 H u0 p0 c0 {13,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([294.235,370.438,436.8,494.563,588.39,659.3,771.987],'J/(mol*K)'),
        H298 = (-274.6,'kJ/mol'),
        S298 = (632.4,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (931.221,'J/(mol*K)'),
        E0 = (-321.178,'kJ/mol'),
    ),
    shortDesc = u"""C12H27B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 37,
    label = 'BPh3',
    molecule =
"""
1  C u0 p0 c0 {2,S} {19,D} {20,S}
2  C u0 p0 c0 {1,S} {3,D} {21,S}
3  C u0 p0 c0 {2,D} {4,S} {22,S}
4  C u0 p0 c0 {3,S} {5,S} {18,D}
5  B u0 p0 c0 {4,S} {6,S} {12,S}
6  C u0 p0 c0 {5,S} {7,S} {11,D}
7  C u0 p0 c0 {6,S} {8,D} {23,S}
8  C u0 p0 c0 {7,D} {9,S} {24,S}
9  C u0 p0 c0 {8,S} {10,D} {25,S}
10 C u0 p0 c0 {9,D} {11,S} {26,S}
11 C u0 p0 c0 {6,D} {10,S} {27,S}
12 C u0 p0 c0 {5,S} {13,S} {17,D}
13 C u0 p0 c0 {12,S} {14,D} {28,S}
14 C u0 p0 c0 {13,D} {15,S} {29,S}
15 C u0 p0 c0 {14,S} {16,D} {30,S}
16 C u0 p0 c0 {15,D} {17,S} {31,S}
17 C u0 p0 c0 {12,D} {16,S} {32,S}
18 C u0 p0 c0 {4,D} {19,S} {33,S}
19 C u0 p0 c0 {1,D} {18,S} {34,S}
20 H u0 p0 c0 {1,S}
21 H u0 p0 c0 {2,S}
22 H u0 p0 c0 {3,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {11,S}
28 H u0 p0 c0 {13,S}
29 H u0 p0 c0 {14,S}
30 H u0 p0 c0 {15,S}
31 H u0 p0 c0 {16,S}
32 H u0 p0 c0 {17,S}
33 H u0 p0 c0 {18,S}
34 H u0 p0 c0 {19,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([258.335,348.143,422.9,483.545,571.968,630.5,710.004],'J/(mol*K)'),
        H298 = (293.3,'kJ/mol'),
        S298 = (517.2,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (818.975,'J/(mol*K)'),
        E0 = (257.435,'kJ/mol'),
    ),
    shortDesc = u"""C18H15B""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 38,
    label = 'BF3',
    molecule =
"""
1 F u0 p3 c0 {2,S}
2 B u0 p0 c0 {1,S} {3,S} {4,S}
3 F u0 p3 c0 {2,S}
4 F u0 p3 c0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([42.8452,49.6224,55.1,59.2725,64.6895,67.7,70.9146],'J/(mol*K)'),
        H298 = (-1134.6,'kJ/mol'),
        S298 = (254.9,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (83.1447,'J/(mol*K)'),
        E0 = (-1144.95,'kJ/mol'),
    ),
    shortDesc = u"""BF3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 39,
    label = 'BCl3',
    molecule =
"""
1 Cl u0 p3 c0 {2,S}
2 B  u0 p0 c0 {1,S} {3,S} {4,S}
3 Cl u0 p3 c0 {2,S}
4 Cl u0 p3 c0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([55.1374,60.9263,64.9,67.5603,70.4969,71.8,72.8211],'J/(mol*K)'),
        H298 = (-404.3,'kJ/mol'),
        S298 = (290,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (83.1447,'J/(mol*K)'),
        E0 = (-417.101,'kJ/mol'),
    ),
    shortDesc = u"""BCl3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 40,
    label = 'BH2OH',
    molecule =
"""
1 B u0 p0 c0 {2,S} {3,S} {4,S}
2 O u0 p2 c0 {1,S} {5,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([34.2617,42.5838,50.4,57.1283,67.3495,74.3,83.9753],'J/(mol*K)'),
        H298 = (-274,'kJ/mol'),
        S298 = (230.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (103.931,'J/(mol*K)'),
        E0 = (-282.629,'kJ/mol'),
    ),
    shortDesc = u"""BH3O""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 41,
    label = 'BH2(OMe)',
    molecule =
"""
1 B u0 p0 c0 {2,S} {4,S} {5,S}
2 O u0 p2 c0 {1,S} {3,S}
3 C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([52.1653,65.6267,78.5,90.0888,109.021,123.1,144.702],'J/(mol*K)'),
        H298 = (-245.8,'kJ/mol'),
        S298 = (269.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (174.604,'J/(mol*K)'),
        E0 = (-257.128,'kJ/mol'),
    ),
    shortDesc = u"""CH5BO""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 42,
    label = 'BH2(OEt)',
    molecule =
"""
1  B u0 p0 c0 {2,S} {5,S} {6,S}
2  O u0 p2 c0 {1,S} {3,S}
3  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([74.1223,94.7502,113.5,129.864,155.876,174.8,203.384],'J/(mol*K)'),
        H298 = (-280.2,'kJ/mol'),
        S298 = (303.2,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (245.277,'J/(mol*K)'),
        E0 = (-294.362,'kJ/mol'),
    ),
    shortDesc = u"""C2H7BO""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 43,
    label = 'BH2(OPh)',
    molecule =
"""
1  B u0 p0 c0 {2,S} {9,S} {10,S}
2  O u0 p2 c0 {1,S} {3,S}
3  C u0 p0 c0 {2,S} {4,S} {8,D}
4  C u0 p0 c0 {3,S} {5,D} {11,S}
5  C u0 p0 c0 {4,D} {6,S} {12,S}
6  C u0 p0 c0 {5,S} {7,D} {13,S}
7  C u0 p0 c0 {6,D} {8,S} {14,S}
8  C u0 p0 c0 {3,D} {7,S} {15,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([113.583,150.251,181,205.929,242.045,265.7,297.468],'J/(mol*K)'),
        H298 = (-121.3,'kJ/mol'),
        S298 = (344.5,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (349.208,'J/(mol*K)'),
        E0 = (-139.575,'kJ/mol'),
    ),
    shortDesc = u"""C6H7BO""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 44,
    label = 'BHMe(OH)',
    molecule =
"""
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 B u0 p0 c0 {1,S} {3,S} {7,S}
3 O u0 p2 c0 {2,S} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([56.192,70.4965,83.4,94.4885,111.704,123.9,141.903],'J/(mol*K)'),
        H298 = (-347.3,'kJ/mol'),
        S298 = (274.2,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (174.604,'J/(mol*K)'),
        E0 = (-359.058,'kJ/mol'),
    ),
    shortDesc = u"""CH5BO""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 45,
    label = 'BHEt(OH)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
3  B u0 p0 c0 {2,S} {4,S} {10,S}
4  O u0 p2 c0 {3,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([78.9116,98.7453,116.5,131.91,156.399,174.3,201.665],'J/(mol*K)'),
        H298 = (-358,'kJ/mol'),
        S298 = (307.9,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (245.277,'J/(mol*K)'),
        E0 = (-373.179,'kJ/mol'),
    ),
    shortDesc = u"""C2H7BO""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 46,
    label = 'BHPh(OH)',
    molecule =
"""
1  O u0 p2 c0 {2,S} {9,S}
2  B u0 p0 c0 {1,S} {3,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {8,D}
4  C u0 p0 c0 {3,S} {5,D} {11,S}
5  C u0 p0 c0 {4,D} {6,S} {12,S}
6  C u0 p0 c0 {5,S} {7,D} {13,S}
7  C u0 p0 c0 {6,D} {8,S} {14,S}
8  C u0 p0 c0 {3,D} {7,S} {15,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([112.8,150.226,181.5,206.712,242.9,266.3,297.203],'J/(mol*K)'),
        H298 = (-210.6,'kJ/mol'),
        S298 = (345.1,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (349.208,'J/(mol*K)'),
        E0 = (-228.532,'kJ/mol'),
    ),
    shortDesc = u"""C6H7BO""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 47,
    label = 'BHMe(OMe)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  B u0 p0 c0 {1,S} {3,S} {8,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([74.6851,93.6526,111.3,127.099,153.05,172.6,203.184],'J/(mol*K)'),
        H298 = (-318,'kJ/mol'),
        S298 = (312.7,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (245.277,'J/(mol*K)'),
        E0 = (-332.696,'kJ/mol'),
    ),
    shortDesc = u"""C2H7BO""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 48,
    label = 'BHEt(OMe)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  B u0 p0 c0 {2,S} {4,S} {11,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {12,S} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([97.0507,123.778,148.1,169.542,204.205,229.9,269.308],'J/(mol*K)'),
        H298 = (-329,'kJ/mol'),
        S298 = (342.9,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (315.95,'J/(mol*K)'),
        E0 = (-346.582,'kJ/mol'),
    ),
    shortDesc = u"""C3H9BO""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 49,
    label = 'BHPh(OMe)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {13,S}
4  C u0 p0 c0 {3,S} {5,S} {9,D}
5  C u0 p0 c0 {4,S} {6,D} {14,S}
6  C u0 p0 c0 {5,D} {7,S} {15,S}
7  C u0 p0 c0 {6,S} {8,D} {16,S}
8  C u0 p0 c0 {7,D} {9,S} {17,S}
9  C u0 p0 c0 {4,D} {8,S} {18,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([134.204,176.722,213,243.079,288.139,318.9,362.174],'J/(mol*K)'),
        H298 = (-182,'kJ/mol'),
        S298 = (384,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (419.881,'J/(mol*K)'),
        E0 = (-203.334,'kJ/mol'),
    ),
    shortDesc = u"""C7H9BO""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 50,
    label = 'BMe2(OH)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  B u0 p0 c0 {1,S} {3,S} {4,S}
3  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
4  O u0 p2 c0 {2,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([77.8207,98.0722,116.1,131.63,156.051,173.7,200.395],'J/(mol*K)'),
        H298 = (-415,'kJ/mol'),
        S298 = (317.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (245.277,'J/(mol*K)'),
        E0 = (-429.854,'kJ/mol'),
    ),
    shortDesc = u"""C2H7BO""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 51,
    label = 'BMe2(OMe)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {12,S} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([99.612,124.091,146,165.214,196.354,219.7,256.478],'J/(mol*K)'),
        H298 = (-381.7,'kJ/mol'),
        S298 = (357.7,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (315.95,'J/(mol*K)'),
        E0 = (-400.089,'kJ/mol'),
    ),
    shortDesc = u"""C3H9BO""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 52,
    label = 'BEtMe(OH)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {11,S} {12,S} {13,S}
5  O u0 p2 c0 {3,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([99.5479,125.681,148.8,168.766,200.432,223.6,259.135],'J/(mol*K)'),
        H298 = (-427.2,'kJ/mol'),
        S298 = (348.7,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (315.95,'J/(mol*K)'),
        E0 = (-445.189,'kJ/mol'),
    ),
    shortDesc = u"""C3H9BO""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 53,
    label = 'BMePh(OH)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
2  B u0 p0 c0 {1,S} {3,S} {4,S}
3  O u0 p2 c0 {2,S} {13,S}
4  C u0 p0 c0 {2,S} {5,S} {9,D}
5  C u0 p0 c0 {4,S} {6,D} {14,S}
6  C u0 p0 c0 {5,D} {7,S} {15,S}
7  C u0 p0 c0 {6,S} {8,D} {16,S}
8  C u0 p0 c0 {7,D} {9,S} {17,S}
9  C u0 p0 c0 {4,D} {8,S} {18,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([136.909,179.221,214.6,243.419,285.708,314,353.282],'J/(mol*K)'),
        H298 = (-278.2,'kJ/mol'),
        S298 = (387,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (419.881,'J/(mol*K)'),
        E0 = (-299.899,'kJ/mol'),
    ),
    shortDesc = u"""C7H9BO""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 54,
    label = 'BMePh(OMe)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {6,S} {10,D}
6  C u0 p0 c0 {5,S} {7,D} {17,S}
7  C u0 p0 c0 {6,D} {8,S} {18,S}
8  C u0 p0 c0 {7,S} {9,D} {19,S}
9  C u0 p0 c0 {8,D} {10,S} {20,S}
10 C u0 p0 c0 {5,D} {9,S} {21,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([155.258,200.171,238.7,271.148,321.165,356.7,409.551],'J/(mol*K)'),
        H298 = (-244.6,'kJ/mol'),
        S298 = (430.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (490.554,'J/(mol*K)'),
        E0 = (-269.665,'kJ/mol'),
    ),
    shortDesc = u"""C8H11BO""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 55,
    label = 'BH(OH)2',
    molecule =
"""
1 O u0 p2 c0 {2,S} {4,S}
2 B u0 p0 c0 {1,S} {3,S} {5,S}
3 O u0 p2 c0 {2,S} {6,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([48.35,60.3243,70.5,78.6539,90.0587,97.1,105.778],'J/(mol*K)'),
        H298 = (-646.3,'kJ/mol'),
        S298 = (256.7,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (124.717,'J/(mol*K)'),
        E0 = (-656.773,'kJ/mol'),
    ),
    shortDesc = u"""BH3O2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 56,
    label = 'BH(OMe)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {9,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([84.4594,106.563,126.3,143.392,170.475,190.2,220.186],'J/(mol*K)'),
        H298 = (-579.8,'kJ/mol'),
        S298 = (338.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (266.063,'J/(mol*K)'),
        E0 = (-595.633,'kJ/mol'),
    ),
    shortDesc = u"""C2H7BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 57,
    label = 'BH(iPrO)2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {13,S}
3  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
4  O u0 p2 c0 {2,S} {5,S}
5  B u0 p0 c0 {4,S} {6,S} {17,S}
6  O u0 p2 c0 {5,S} {7,S}
7  C u0 p0 c0 {6,S} {8,S} {9,S} {18,S}
8  C u0 p0 c0 {7,S} {19,S} {20,S} {21,S}
9  C u0 p0 c0 {7,S} {22,S} {23,S} {24,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([180.916,228.068,268.3,302.397,355.824,394.7,454.365],'J/(mol*K)'),
        H298 = (-719.7,'kJ/mol'),
        S298 = (462.1,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (548.755,'J/(mol*K)'),
        E0 = (-749.328,'kJ/mol'),
    ),
    shortDesc = u"""C6H15BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 58,
    label = 'BMe(OH)2',
    molecule =
"""
1 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2 B u0 p0 c0 {1,S} {3,S} {4,S}
3 O u0 p2 c0 {2,S} {8,S}
4 O u0 p2 c0 {2,S} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {4,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([70.579,88.4457,103.6,116.009,134.234,146.4,163.256],'J/(mol*K)'),
        H298 = (-713.6,'kJ/mol'),
        S298 = (304.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (195.39,'J/(mol*K)'),
        E0 = (-727.27,'kJ/mol'),
    ),
    shortDesc = u"""CH5BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 59,
    label = 'BMe(OH)(OMe)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
5  O u0 p2 c0 {3,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([88.9865,111.993,131.9,148.683,174.474,192.7,219.699],'J/(mol*K)'),
        H298 = (-680,'kJ/mol'),
        S298 = (342.3,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (266.063,'J/(mol*K)'),
        E0 = (-696.378,'kJ/mol'),
    ),
    shortDesc = u"""C2H7BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 60,
    label = 'BMe(OMe)2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  C u0 p0 c0 {3,S} {10,S} {11,S} {12,S}
5  O u0 p2 c0 {3,S} {6,S}
6  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([109.743,135.527,158.5,178.69,211.616,236.5,276.009],'J/(mol*K)'),
        H298 = (-646.1,'kJ/mol'),
        S298 = (379.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (336.736,'J/(mol*K)'),
        E0 = (-666.219,'kJ/mol'),
    ),
    shortDesc = u"""C3H9BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 61,
    label = 'BMe(OEt)(OMe)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  O u0 p2 c0 {2,S} {4,S}
4  B u0 p0 c0 {3,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
6  O u0 p2 c0 {4,S} {7,S}
7  C u0 p0 c0 {6,S} {16,S} {17,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([131.227,165.436,195.2,220.718,261.067,290.6,336.04],'J/(mol*K)'),
        H298 = (-680.4,'kJ/mol'),
        S298 = (411.7,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (407.409,'J/(mol*K)'),
        E0 = (-702.963,'kJ/mol'),
    ),
    shortDesc = u"""C4H11BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 62,
    label = 'BMe(OEt)2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
3  O u0 p2 c0 {2,S} {4,S}
4  B u0 p0 c0 {3,S} {5,S} {6,S}
5  C u0 p0 c0 {4,S} {14,S} {15,S} {16,S}
6  O u0 p2 c0 {4,S} {7,S}
7  C u0 p0 c0 {6,S} {8,S} {17,S} {18,S}
8  C u0 p0 c0 {7,S} {19,S} {20,S} {21,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([154.279,195.361,230.7,260.738,307.814,342,394.278],'J/(mol*K)'),
        H298 = (-713.8,'kJ/mol'),
        S298 = (444.5,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (478.082,'J/(mol*K)'),
        E0 = (-739.459,'kJ/mol'),
    ),
    shortDesc = u"""C5H13BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 63,
    label = 'BEt(OH)2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  O u0 p2 c0 {3,S} {11,S}
5  O u0 p2 c0 {3,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([93.0048,116.689,136.9,153.749,179.301,197.1,223.032],'J/(mol*K)'),
        H298 = (-724.8,'kJ/mol'),
        S298 = (331.1,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (266.063,'J/(mol*K)'),
        E0 = (-741.763,'kJ/mol'),
    ),
    shortDesc = u"""C2H7BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 64,
    label = 'BEt(OH)(OMe)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  O u0 p2 c0 {3,S} {12,S}
5  O u0 p2 c0 {3,S} {6,S}
6  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([111.3,139.594,164.2,185.229,218.285,242.3,278.936],'J/(mol*K)'),
        H298 = (-690.9,'kJ/mol'),
        S298 = (372.4,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (336.736,'J/(mol*K)'),
        E0 = (-710.682,'kJ/mol'),
    ),
    shortDesc = u"""C3H9BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 65,
    label = 'BEt(OMe)2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  B u0 p0 c0 {2,S} {4,S} {6,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
6  O u0 p2 c0 {3,S} {7,S}
7  C u0 p0 c0 {6,S} {16,S} {17,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([131.569,163.21,191.3,215.987,256.327,286.9,335.549],'J/(mol*K)'),
        H298 = (-656.5,'kJ/mol'),
        S298 = (412.5,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (407.409,'J/(mol*K)'),
        E0 = (-679.793,'kJ/mol'),
    ),
    shortDesc = u"""C4H11BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 66,
    label = 'BVi(OH)2',
    molecule =
"""
1  C u0 p0 c0 {2,D} {6,S} {7,S}
2  C u0 p0 c0 {1,D} {3,S} {8,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {3,S} {10,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([87.2537,107.898,124.4,137.284,155.201,166.5,181.44],'J/(mol*K)'),
        H298 = (-600.7,'kJ/mol'),
        S298 = (316.5,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (220.334,'J/(mol*K)'),
        E0 = (-616.82,'kJ/mol'),
    ),
    shortDesc = u"""C2H5BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 67,
    label = 'BVi(OH)(OMe)',
    molecule =
"""
1  C u0 p0 c0 {2,D} {7,S} {8,S}
2  C u0 p0 c0 {1,D} {3,S} {9,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  O u0 p2 c0 {3,S} {10,S}
5  O u0 p2 c0 {3,S} {6,S}
6  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([100.952,126.767,148.7,166.957,194.656,214,242.349],'J/(mol*K)'),
        H298 = (-567.5,'kJ/mol'),
        S298 = (362.5,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (291.007,'J/(mol*K)'),
        E0 = (-585.578,'kJ/mol'),
    ),
    shortDesc = u"""C3H7BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 68,
    label = 'BVi(OMe)2',
    molecule =
"""
1  C u0 p0 c0 {2,D} {8,S} {9,S}
2  C u0 p0 c0 {1,D} {3,S} {10,S}
3  B u0 p0 c0 {2,S} {4,S} {6,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
6  O u0 p2 c0 {3,S} {7,S}
7  C u0 p0 c0 {6,S} {14,S} {15,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([122.497,150.555,175.2,196.691,231.547,257.8,299.389],'J/(mol*K)'),
        H298 = (-531.7,'kJ/mol'),
        S298 = (399.9,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (361.68,'J/(mol*K)'),
        E0 = (-553.811,'kJ/mol'),
    ),
    shortDesc = u"""C4H9BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 69,
    label = 'B(All)(OH)2',
    molecule =
"""
1  C u0 p0 c0 {2,D} {7,S} {8,S}
2  C u0 p0 c0 {1,D} {3,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
4  B u0 p0 c0 {3,S} {5,S} {6,S}
5  O u0 p2 c0 {4,S} {12,S}
6  O u0 p2 c0 {4,S} {13,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([106.765,132.932,154.9,173.053,200.413,219.4,246.997],'J/(mol*K)'),
        H298 = (-627.3,'kJ/mol'),
        S298 = (349.5,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (291.007,'J/(mol*K)'),
        E0 = (-646.387,'kJ/mol'),
    ),
    shortDesc = u"""C3H7BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 70,
    label = 'B(All)(OH)(OMe)',
    molecule =
"""
1  C u0 p0 c0 {2,D} {8,S} {9,S}
2  C u0 p0 c0 {1,D} {3,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
4  B u0 p0 c0 {3,S} {5,S} {6,S}
5  O u0 p2 c0 {4,S} {13,S}
6  O u0 p2 c0 {4,S} {7,S}
7  C u0 p0 c0 {6,S} {14,S} {15,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([122.277,153.918,181,203.821,239.104,264.3,302.075],'J/(mol*K)'),
        H298 = (-594.2,'kJ/mol'),
        S298 = (390.3,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (361.68,'J/(mol*K)'),
        E0 = (-615.375,'kJ/mol'),
    ),
    shortDesc = u"""C4H9BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 71,
    label = 'B(All)(OMe)2',
    molecule =
"""
1  C u0 p0 c0 {2,D} {9,S} {10,S}
2  C u0 p0 c0 {1,D} {3,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {12,S} {13,S}
4  B u0 p0 c0 {3,S} {5,S} {7,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
7  O u0 p2 c0 {4,S} {8,S}
8  C u0 p0 c0 {7,S} {17,S} {18,S} {19,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([143.542,178.288,208.6,234.861,277.119,308.7,358.344],'J/(mol*K)'),
        H298 = (-552.7,'kJ/mol'),
        S298 = (436.1,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (432.353,'J/(mol*K)'),
        E0 = (-577.618,'kJ/mol'),
    ),
    shortDesc = u"""C5H11BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 72,
    label = 'BPh(OH)2',
    molecule =
"""
1  O u0 p2 c0 {2,S} {10,S}
2  B u0 p0 c0 {1,S} {3,S} {4,S}
3  O u0 p2 c0 {2,S} {11,S}
4  C u0 p0 c0 {2,S} {5,S} {9,D}
5  C u0 p0 c0 {4,S} {6,D} {12,S}
6  C u0 p0 c0 {5,D} {7,S} {13,S}
7  C u0 p0 c0 {6,S} {8,D} {14,S}
8  C u0 p0 c0 {7,D} {9,S} {15,S}
9  C u0 p0 c0 {4,D} {8,S} {16,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([126.154,165.689,198.2,224.134,260.984,284.6,315.588],'J/(mol*K)'),
        H298 = (-572.6,'kJ/mol'),
        S298 = (382,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (369.994,'J/(mol*K)'),
        E0 = (-592.593,'kJ/mol'),
    ),
    shortDesc = u"""C6H7BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 73,
    label = 'BPh(OH)(OMe)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  O u0 p2 c0 {3,S} {14,S}
5  C u0 p0 c0 {3,S} {6,S} {10,D}
6  C u0 p0 c0 {5,S} {7,D} {15,S}
7  C u0 p0 c0 {6,D} {8,S} {16,S}
8  C u0 p0 c0 {7,S} {9,D} {17,S}
9  C u0 p0 c0 {8,D} {10,S} {18,S}
10 C u0 p0 c0 {5,D} {9,S} {19,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([146.234,189.644,226,255.789,299.986,330,372.432],'J/(mol*K)'),
        H298 = (-539.8,'kJ/mol'),
        S298 = (418.2,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (440.667,'J/(mol*K)'),
        E0 = (-563.136,'kJ/mol'),
    ),
    shortDesc = u"""C7H9BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 74,
    label = 'BPh(OMe)2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {6,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {15,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {7,S} {11,D}
7  C u0 p0 c0 {6,S} {8,D} {18,S}
8  C u0 p0 c0 {7,D} {9,S} {19,S}
9  C u0 p0 c0 {8,S} {10,D} {20,S}
10 C u0 p0 c0 {9,D} {11,S} {21,S}
11 C u0 p0 c0 {6,D} {10,S} {22,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {11,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([170.295,216.624,256.2,289.595,341.39,378.5,434.119],'J/(mol*K)'),
        H298 = (-497.6,'kJ/mol'),
        S298 = (448.1,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (511.34,'J/(mol*K)'),
        E0 = (-525.315,'kJ/mol'),
    ),
    shortDesc = u"""C8H11BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 75,
    label = 'B(p-Tol)(OH)2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {10,D}
3  C u0 p0 c0 {2,S} {4,D} {14,S}
4  C u0 p0 c0 {3,D} {5,S} {15,S}
5  C u0 p0 c0 {4,S} {6,S} {9,D}
6  B u0 p0 c0 {5,S} {7,S} {8,S}
7  O u0 p2 c0 {6,S} {16,S}
8  O u0 p2 c0 {6,S} {17,S}
9  C u0 p0 c0 {5,D} {10,S} {18,S}
10 C u0 p0 c0 {2,D} {9,S} {19,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([148.066,192.852,230.1,260.362,304.689,334.3,375.288],'J/(mol*K)'),
        H298 = (-605.8,'kJ/mol'),
        S298 = (420.3,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (440.667,'J/(mol*K)'),
        E0 = (-629.145,'kJ/mol'),
    ),
    shortDesc = u"""C7H9BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 76,
    label = 'B(p-Tol)(OH)(OMe)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  O u0 p2 c0 {3,S} {15,S}
5  C u0 p0 c0 {3,S} {6,S} {11,D}
6  C u0 p0 c0 {5,S} {7,D} {16,S}
7  C u0 p0 c0 {6,D} {8,S} {17,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  C u0 p0 c0 {8,S} {18,S} {19,S} {20,S}
10 C u0 p0 c0 {8,D} {11,S} {21,S}
11 C u0 p0 c0 {5,D} {10,S} {22,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {11,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([167.843,216.361,257.4,291.527,343.323,379.5,432.29],'J/(mol*K)'),
        H298 = (-572.9,'kJ/mol'),
        S298 = (457.9,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (511.34,'J/(mol*K)'),
        E0 = (-599.575,'kJ/mol'),
    ),
    shortDesc = u"""C8H11BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 77,
    label = 'B(p-Tol)(OMe)2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {6,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {3,S} {7,S} {12,D}
7  C u0 p0 c0 {6,S} {8,D} {19,S}
8  C u0 p0 c0 {7,D} {9,S} {20,S}
9  C u0 p0 c0 {8,S} {10,S} {11,D}
10 C u0 p0 c0 {9,S} {21,S} {22,S} {23,S}
11 C u0 p0 c0 {9,D} {12,S} {24,S}
12 C u0 p0 c0 {6,D} {11,S} {25,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {10,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {11,S}
25 H u0 p0 c0 {12,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([188.53,241.248,286.4,324.59,383.997,426.7,490.968],'J/(mol*K)'),
        H298 = (-530.4,'kJ/mol'),
        S298 = (487.4,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (582.013,'J/(mol*K)'),
        E0 = (-560.471,'kJ/mol'),
    ),
    shortDesc = u"""C9H13BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 78,
    label = 'B(OH)3',
    molecule =
"""
1 O u0 p2 c0 {2,S} {5,S}
2 B u0 p0 c0 {1,S} {3,S} {4,S}
3 O u0 p2 c0 {2,S} {6,S}
4 O u0 p2 c0 {2,S} {7,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([66.1343,81.2991,93.1,101.888,113.052,119.1,125.26],'J/(mol*K)'),
        H298 = (-1002.1,'kJ/mol'),
        S298 = (272.2,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (145.503,'J/(mol*K)'),
        E0 = (-1015.14,'kJ/mol'),
    ),
    shortDesc = u"""BH3O3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 79,
    label = 'B(OH)2(OMe)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {3,S} {10,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([83.8207,103.273,119.4,132.481,151.648,164.5,182.527],'J/(mol*K)'),
        H298 = (-967.6,'kJ/mol'),
        S298 = (323.5,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (216.176,'J/(mol*K)'),
        E0 = (-983.465,'kJ/mol'),
    ),
    shortDesc = u"""CH5BO3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 80,
    label = 'B(OH)2(OPh)',
    molecule =
"""
1  O u0 p2 c0 {2,S} {11,S}
2  B u0 p0 c0 {1,S} {3,S} {4,S}
3  O u0 p2 c0 {2,S} {12,S}
4  O u0 p2 c0 {2,S} {5,S}
5  C u0 p0 c0 {4,S} {6,S} {10,D}
6  C u0 p0 c0 {5,S} {7,D} {13,S}
7  C u0 p0 c0 {6,D} {8,S} {14,S}
8  C u0 p0 c0 {7,S} {9,D} {15,S}
9  C u0 p0 c0 {8,D} {10,S} {16,S}
10 C u0 p0 c0 {5,D} {9,S} {17,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {9,S}
17 H u0 p0 c0 {10,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([140.022,186.534,223.2,251.039,287.667,308.7,332.317],'J/(mol*K)'),
        H298 = (-842.9,'kJ/mol'),
        S298 = (398.7,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (390.78,'J/(mol*K)'),
        E0 = (-863.727,'kJ/mol'),
    ),
    shortDesc = u"""C6H7BO3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 81,
    label = 'B(OH)(OMe)2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  O u0 p2 c0 {3,S} {10,S}
5  O u0 p2 c0 {3,S} {6,S}
6  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([101.803,125.328,145.6,162.845,189.878,209.5,239.463],'J/(mol*K)'),
        H298 = (-933.1,'kJ/mol'),
        S298 = (365.1,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (286.849,'J/(mol*K)'),
        E0 = (-951.89,'kJ/mol'),
    ),
    shortDesc = u"""C2H7BO3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 82,
    label = 'B(OH)(OMe)(OPh)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {5,S}
4  O u0 p2 c0 {3,S} {15,S}
5  O u0 p2 c0 {3,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {11,D}
7  C u0 p0 c0 {6,S} {8,D} {16,S}
8  C u0 p0 c0 {7,D} {9,S} {17,S}
9  C u0 p0 c0 {8,S} {10,D} {18,S}
10 C u0 p0 c0 {9,D} {11,S} {19,S}
11 C u0 p0 c0 {6,D} {10,S} {20,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
19 H u0 p0 c0 {10,S}
20 H u0 p0 c0 {11,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([161.699,207.603,245.4,275.995,320.82,350.9,392.963],'J/(mol*K)'),
        H298 = (-809.2,'kJ/mol'),
        S298 = (442.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (461.453,'J/(mol*K)'),
        E0 = (-834.903,'kJ/mol'),
    ),
    shortDesc = u"""C7H9BO3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 83,
    label = 'B(OMe)3',
    molecule =
"""
1  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {6,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
6  O u0 p2 c0 {3,S} {7,S}
7  C u0 p0 c0 {6,S} {14,S} {15,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([119.702,148.059,173,194.701,229.713,255.9,297.022],'J/(mol*K)'),
        H298 = (-898.6,'kJ/mol'),
        S298 = (396.4,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (357.522,'J/(mol*K)'),
        E0 = (-920.111,'kJ/mol'),
    ),
    shortDesc = u"""C3H9BO3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 84,
    label = 'B(OEt)(OMe)2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
3  O u0 p2 c0 {2,S} {4,S}
4  B u0 p0 c0 {3,S} {5,S} {7,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
7  O u0 p2 c0 {4,S} {8,S}
8  C u0 p0 c0 {7,S} {17,S} {18,S} {19,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([141.864,177.689,208.8,235.534,278.018,309.3,357.652],'J/(mol*K)'),
        H298 = (-932.7,'kJ/mol'),
        S298 = (437.6,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (428.195,'J/(mol*K)'),
        E0 = (-957.038,'kJ/mol'),
    ),
    shortDesc = u"""C4H11BO3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 85,
    label = 'B(OMe)2(OPh)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {6,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {3,S} {7,S}
7  C u0 p0 c0 {6,S} {8,S} {12,D}
8  C u0 p0 c0 {7,S} {9,D} {19,S}
9  C u0 p0 c0 {8,D} {10,S} {20,S}
10 C u0 p0 c0 {9,S} {11,D} {21,S}
11 C u0 p0 c0 {10,D} {12,S} {22,S}
12 C u0 p0 c0 {7,D} {11,S} {23,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {10,S}
22 H u0 p0 c0 {11,S}
23 H u0 p0 c0 {12,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([180.085,230.204,272.2,306.947,359.517,396.2,449.842],'J/(mol*K)'),
        H298 = (-775.5,'kJ/mol'),
        S298 = (484.2,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (532.126,'J/(mol*K)'),
        E0 = (-804.16,'kJ/mol'),
    ),
    shortDesc = u"""C8H11BO3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 86,
    label = 'B(OEt)2(OMe)',
    molecule =
"""
1  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {3,S} {13,S} {14,S}
3  O u0 p2 c0 {2,S} {4,S}
4  B u0 p0 c0 {3,S} {5,S} {7,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
7  O u0 p2 c0 {4,S} {8,S}
8  C u0 p0 c0 {7,S} {9,S} {18,S} {19,S}
9  C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([163.833,207.368,244.8,276.633,326.556,362.8,418.02],'J/(mol*K)'),
        H298 = (-966.8,'kJ/mol'),
        S298 = (469.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (498.868,'J/(mol*K)'),
        E0 = (-993.862,'kJ/mol'),
    ),
    shortDesc = u"""C5H13BO3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 87,
    label = 'B(OMe)(OPh)2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {18,S} {19,S} {20,S}
2  O u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {11,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,S} {10,D}
6  C u0 p0 c0 {5,S} {7,D} {21,S}
7  C u0 p0 c0 {6,D} {8,S} {22,S}
8  C u0 p0 c0 {7,S} {9,D} {23,S}
9  C u0 p0 c0 {8,D} {10,S} {24,S}
10 C u0 p0 c0 {5,D} {9,S} {25,S}
11 O u0 p2 c0 {3,S} {12,S}
12 C u0 p0 c0 {11,S} {13,S} {17,D}
13 C u0 p0 c0 {12,S} {14,D} {26,S}
14 C u0 p0 c0 {13,D} {15,S} {27,S}
15 C u0 p0 c0 {14,S} {16,D} {28,S}
16 C u0 p0 c0 {15,D} {17,S} {29,S}
17 C u0 p0 c0 {12,D} {16,S} {30,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {1,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {10,S}
26 H u0 p0 c0 {13,S}
27 H u0 p0 c0 {14,S}
28 H u0 p0 c0 {15,S}
29 H u0 p0 c0 {16,S}
30 H u0 p0 c0 {17,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([244.59,317.22,376.5,424.187,493.612,539.9,604.179],'J/(mol*K)'),
        H298 = (-652.5,'kJ/mol'),
        S298 = (559.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (706.73,'J/(mol*K)'),
        E0 = (-688.859,'kJ/mol'),
    ),
    shortDesc = u"""C13H13BO3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 88,
    label = 'B(OEt)3',
    molecule =
"""
1  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {14,S} {15,S}
3  O u0 p2 c0 {2,S} {4,S}
4  B u0 p0 c0 {3,S} {5,S} {8,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
7  C u0 p0 c0 {6,S} {18,S} {19,S} {20,S}
8  O u0 p2 c0 {4,S} {9,S}
9  C u0 p0 c0 {8,S} {10,S} {21,S} {22,S}
10 C u0 p0 c0 {9,S} {23,S} {24,S} {25,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
23 H u0 p0 c0 {10,S}
24 H u0 p0 c0 {10,S}
25 H u0 p0 c0 {10,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([186.004,237.35,281.1,317.957,375.067,416,477.575],'J/(mol*K)'),
        H298 = (-1000.8,'kJ/mol'),
        S298 = (492.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (569.541,'J/(mol*K)'),
        E0 = (-1030.58,'kJ/mol'),
    ),
    shortDesc = u"""C6H15BO3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 89,
    label = 'C6H4O2BH',
    molecule =
"""
1  B u0 p0 c0 {2,S} {9,S} {10,S}
2  O u0 p2 c0 {1,S} {3,S}
3  C u0 p0 c0 {2,S} {4,S} {8,D}
4  C u0 p0 c0 {3,S} {5,D} {11,S}
5  C u0 p0 c0 {4,D} {6,S} {12,S}
6  C u0 p0 c0 {5,S} {7,D} {13,S}
7  C u0 p0 c0 {6,D} {8,S} {14,S}
8  C u0 p0 c0 {3,D} {7,S} {9,S}
9  O u0 p2 c0 {1,S} {8,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([105.095,142.317,173.2,197.772,232.239,253.8,281.024],'J/(mol*K)'),
        H298 = (-409.1,'kJ/mol'),
        S298 = (322,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (332.579,'J/(mol*K)'),
        E0 = (-425.486,'kJ/mol'),
    ),
    shortDesc = u"""C6H5BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 90,
    label = 'C6H4O2BMe',
    molecule =
"""
1  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
2  B u0 p0 c0 {1,S} {3,S} {10,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {9,D}
5  C u0 p0 c0 {4,S} {6,D} {14,S}
6  C u0 p0 c0 {5,D} {7,S} {15,S}
7  C u0 p0 c0 {6,S} {8,D} {16,S}
8  C u0 p0 c0 {7,D} {9,S} {17,S}
9  C u0 p0 c0 {4,D} {8,S} {10,S}
10 O u0 p2 c0 {2,S} {9,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([128.704,170.928,206.2,234.745,276.081,303.2,339.846],'J/(mol*K)'),
        H298 = (-482.9,'kJ/mol'),
        S298 = (371.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (403.252,'J/(mol*K)'),
        E0 = (-502.985,'kJ/mol'),
    ),
    shortDesc = u"""C7H7BO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 91,
    label = '3-C6H3FO2BMe',
    molecule =
"""
1  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
2  B u0 p0 c0 {1,S} {3,S} {11,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {10,D}
5  C u0 p0 c0 {4,S} {6,D} {15,S}
6  C u0 p0 c0 {5,D} {7,S} {16,S}
7  C u0 p0 c0 {6,S} {8,D} {17,S}
8  C u0 p0 c0 {7,D} {9,S} {10,S}
9  F u0 p3 c0 {8,S}
10 C u0 p0 c0 {4,D} {8,S} {11,S}
11 O u0 p2 c0 {2,S} {10,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([141.007,182.828,217.2,244.784,284.51,310.5,345.588],'J/(mol*K)'),
        H298 = (-662.3,'kJ/mol'),
        S298 = (390.4,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (403.252,'J/(mol*K)'),
        E0 = (-684.749,'kJ/mol'),
    ),
    shortDesc = u"""C7H6BFO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 92,
    label = '4-C6H3FO2BMe',
    molecule =
"""
1  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
2  B u0 p0 c0 {1,S} {3,S} {11,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {10,D}
5  C u0 p0 c0 {4,S} {6,D} {15,S}
6  C u0 p0 c0 {5,D} {7,S} {16,S}
7  C u0 p0 c0 {6,S} {8,S} {9,D}
8  F u0 p3 c0 {7,S}
9  C u0 p0 c0 {7,D} {10,S} {17,S}
10 C u0 p0 c0 {4,D} {9,S} {11,S}
11 O u0 p2 c0 {2,S} {10,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {9,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([141.312,183.334,217.8,245.394,284.999,310.8,345.459],'J/(mol*K)'),
        H298 = (-673.1,'kJ/mol'),
        S298 = (389.5,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (403.252,'J/(mol*K)'),
        E0 = (-695.547,'kJ/mol'),
    ),
    shortDesc = u"""C7H6BFO2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 93,
    label = '3,4-C6H2F2O2BMe',
    molecule =
"""
1  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
2  B u0 p0 c0 {1,S} {3,S} {12,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {11,D}
5  C u0 p0 c0 {4,S} {6,D} {16,S}
6  C u0 p0 c0 {5,D} {7,S} {17,S}
7  C u0 p0 c0 {6,S} {8,S} {9,D}
8  F u0 p3 c0 {7,S}
9  C u0 p0 c0 {7,D} {10,S} {11,S}
10 F u0 p3 c0 {9,S}
11 C u0 p0 c0 {4,D} {9,S} {12,S}
12 O u0 p2 c0 {2,S} {11,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([153.608,194.916,228.3,254.865,292.933,317.8,351.404],'J/(mol*K)'),
        H298 = (-834.8,'kJ/mol'),
        S298 = (409.4,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (403.252,'J/(mol*K)'),
        E0 = (-859.695,'kJ/mol'),
    ),
    shortDesc = u"""C7H5BF2O2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 94,
    label = '3,6-C6H2F2O2BMe',
    molecule =
"""
1  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
2  B u0 p0 c0 {1,S} {3,S} {12,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {11,D}
5  C u0 p0 c0 {4,S} {6,S} {7,D}
6  F u0 p3 c0 {5,S}
7  C u0 p0 c0 {5,D} {8,S} {16,S}
8  C u0 p0 c0 {7,S} {9,D} {17,S}
9  C u0 p0 c0 {8,D} {10,S} {11,S}
10 F u0 p3 c0 {9,S}
11 C u0 p0 c0 {4,D} {9,S} {12,S}
12 O u0 p2 c0 {2,S} {11,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([153.508,194.808,228.2,254.781,292.892,317.8,351.476],'J/(mol*K)'),
        H298 = (-836.9,'kJ/mol'),
        S298 = (408.5,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (403.252,'J/(mol*K)'),
        E0 = (-861.781,'kJ/mol'),
    ),
    shortDesc = u"""C7H5BF2O2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 95,
    label = '4,5-C6H2F2O2BMe',
    molecule =
"""
1  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
2  B u0 p0 c0 {1,S} {3,S} {12,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {11,D}
5  C u0 p0 c0 {4,S} {6,D} {16,S}
6  C u0 p0 c0 {5,D} {7,S} {8,S}
7  F u0 p3 c0 {6,S}
8  C u0 p0 c0 {6,S} {9,S} {10,D}
9  F u0 p3 c0 {8,S}
10 C u0 p0 c0 {8,D} {11,S} {17,S}
11 C u0 p0 c0 {4,D} {10,S} {12,S}
12 O u0 p2 c0 {2,S} {11,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {10,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([153.911,195.355,228.8,255.37,293.358,318.1,351.42],'J/(mol*K)'),
        H298 = (-845.2,'kJ/mol'),
        S298 = (407.6,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (403.252,'J/(mol*K)'),
        E0 = (-870.113,'kJ/mol'),
    ),
    shortDesc = u"""C7H5BF2O2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 96,
    label = '3,4,5-C6H2F2O2BMe',
    molecule =
"""
1  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
2  B u0 p0 c0 {1,S} {3,S} {13,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {12,D}
5  C u0 p0 c0 {4,S} {6,D} {17,S}
6  C u0 p0 c0 {5,D} {7,S} {8,S}
7  F u0 p3 c0 {6,S}
8  C u0 p0 c0 {6,S} {9,S} {10,D}
9  F u0 p3 c0 {8,S}
10 C u0 p0 c0 {8,D} {11,S} {12,S}
11 F u0 p3 c0 {10,S}
12 C u0 p0 c0 {4,D} {10,S} {13,S}
13 O u0 p2 c0 {2,S} {12,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {5,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([166.306,206.977,239.3,264.82,301.268,325.1,357.439],'J/(mol*K)'),
        H298 = (-1004.3,'kJ/mol'),
        S298 = (428.7,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (403.252,'J/(mol*K)'),
        E0 = (-1031.69,'kJ/mol'),
    ),
    shortDesc = u"""C7H4BF3O2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 97,
    label = '3,4,6-C6H2F2O2BMe',
    molecule =
"""
1  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
2  B u0 p0 c0 {1,S} {3,S} {13,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {12,D}
5  C u0 p0 c0 {4,S} {6,S} {7,D}
6  F u0 p3 c0 {5,S}
7  C u0 p0 c0 {5,D} {8,S} {17,S}
8  C u0 p0 c0 {7,S} {9,S} {10,D}
9  F u0 p3 c0 {8,S}
10 C u0 p0 c0 {8,D} {11,S} {12,S}
11 F u0 p3 c0 {10,S}
12 C u0 p0 c0 {4,D} {10,S} {13,S}
13 O u0 p2 c0 {2,S} {12,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([166.107,206.829,239.2,264.758,301.251,325.1,357.438],'J/(mol*K)'),
        H298 = (-1006.8,'kJ/mol'),
        S298 = (427.6,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (403.252,'J/(mol*K)'),
        E0 = (-1034.15,'kJ/mol'),
    ),
    shortDesc = u"""C7H4BF3O2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 98,
    label = 'C6F4O2BMe',
    molecule =
"""
1  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
2  B u0 p0 c0 {1,S} {3,S} {14,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,S} {13,D}
5  C u0 p0 c0 {4,S} {6,S} {7,D}
6  F u0 p3 c0 {5,S}
7  C u0 p0 c0 {5,D} {8,S} {9,S}
8  F u0 p3 c0 {7,S}
9  C u0 p0 c0 {7,S} {10,S} {11,D}
10 F u0 p3 c0 {9,S}
11 C u0 p0 c0 {9,D} {12,S} {13,S}
12 F u0 p3 c0 {11,S}
13 C u0 p0 c0 {4,D} {11,S} {14,S}
14 O u0 p2 c0 {2,S} {13,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([178.698,218.532,249.7,274.166,309.114,332.1,363.603],'J/(mol*K)'),
        H298 = (-1159.1,'kJ/mol'),
        S298 = (448.5,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (403.252,'J/(mol*K)'),
        E0 = (-1189,'kJ/mol'),
    ),
    shortDesc = u"""C7H3BF4O2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 99,
    label = 'BH2NH2',
    molecule =
"""
1 B u0 p0 c0 {2,S} {3,S} {4,S}
2 N u0 p1 c0 {1,S} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([39.4029,49.6473,59.1,67.2045,79.5712,88.1,100.344],'J/(mol*K)'),
        H298 = (-81.9,'kJ/mol'),
        S298 = (228.4,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (128.874,'J/(mol*K)'),
        E0 = (-91.154,'kJ/mol'),
    ),
    shortDesc = u"""BH4N""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 100,
    label = 'BH2NHMe',
    molecule =
"""
1 B u0 p0 c0 {2,S} {4,S} {5,S}
2 N u0 p1 c0 {1,S} {3,S} {6,S}
3 C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([58.4104,73.9042,88.4,101.28,122.105,137.5,161.177],'J/(mol*K)'),
        H298 = (-72.6,'kJ/mol'),
        S298 = (270.7,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (199.547,'J/(mol*K)'),
        E0 = (-84.7188,'kJ/mol'),
    ),
    shortDesc = u"""CH6BN""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 101,
    label = 'BH2NMe2',
    molecule =
"""
1  B u0 p0 c0 {2,S} {5,S} {6,S}
2  N u0 p1 c0 {1,S} {3,S} {4,S}
3  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([82.0483,103.828,123.6,140.954,168.866,189.5,221.336],'J/(mol*K)'),
        H298 = (-70,'kJ/mol'),
        S298 = (293.4,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (270.22,'J/(mol*K)'),
        E0 = (-85.5075,'kJ/mol'),
    ),
    shortDesc = u"""C2H8BN""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 102,
    label = 'BHMeNH2',
    molecule =
"""
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 B u0 p0 c0 {1,S} {3,S} {7,S}
3 N u0 p1 c0 {2,S} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([62.5243,78.3071,92.5,104.746,123.961,137.8,158.742],'J/(mol*K)'),
        H298 = (-147,'kJ/mol'),
        S298 = (275.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (199.547,'J/(mol*K)'),
        E0 = (-159.717,'kJ/mol'),
    ),
    shortDesc = u"""CH6BN""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 103,
    label = 'BMe2NH2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  B u0 p0 c0 {1,S} {3,S} {4,S}
3  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
4  N u0 p1 c0 {2,S} {11,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([84.9488,106.435,125.5,141.97,168.084,187.2,216.687],'J/(mol*K)'),
        H298 = (-211,'kJ/mol'),
        S298 = (312.3,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (270.22,'J/(mol*K)'),
        E0 = (-227.007,'kJ/mol'),
    ),
    shortDesc = u"""C2H8BN""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 104,
    label = 'BH(NMe2)2',
    molecule =
"""
1  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
2  N u0 p1 c0 {1,S} {3,S} {4,S}
3  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
4  B u0 p0 c0 {2,S} {5,S} {14,S}
5  N u0 p1 c0 {4,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
7  C u0 p0 c0 {5,S} {18,S} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([146.975,183.736,216.5,245.363,292.574,328.3,384.743],'J/(mol*K)'),
        H298 = (-160.2,'kJ/mol'),
        S298 = (390.4,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (457.296,'J/(mol*K)'),
        E0 = (-185.603,'kJ/mol'),
    ),
    shortDesc = u"""C4H13BN2""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 105,
    label = 'B(NMe2)3',
    molecule =
"""
1  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
2  N u0 p1 c0 {1,S} {3,S} {4,S}
3  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
4  B u0 p0 c0 {2,S} {5,S} {8,S}
5  N u0 p1 c0 {4,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
8  N u0 p1 c0 {4,S} {9,S} {10,S}
9  C u0 p0 c0 {8,S} {23,S} {24,S} {25,S}
10 C u0 p0 c0 {8,S} {26,S} {27,S} {28,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {9,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {10,S}
27 H u0 p0 c0 {10,S}
28 H u0 p0 c0 {10,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([216.44,269.495,315.6,355.658,420.551,469.4,546.443],'J/(mol*K)'),
        H298 = (-239.6,'kJ/mol'),
        S298 = (461,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (644.372,'J/(mol*K)'),
        E0 = (-275.315,'kJ/mol'),
    ),
    shortDesc = u"""C6H18BN3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 106,
    label = 'B(SMe)3',
    molecule =
"""
1  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
2  S u0 p2 c0 {1,S} {3,S}
3  B u0 p0 c0 {2,S} {4,S} {6,S}
4  S u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
6  S u0 p2 c0 {3,S} {7,S}
7  C u0 p0 c0 {6,S} {14,S} {15,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([153.426,181.094,203.2,221.366,249.494,270.1,302.55],'J/(mol*K)'),
        H298 = (-199.1,'kJ/mol'),
        S298 = (435.1,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (357.522,'J/(mol*K)'),
        E0 = (-226.794,'kJ/mol'),
    ),
    shortDesc = u"""C3H9BS3""",
    longDesc = u"""Calculated with W1X-1""",
)

entry(
    index = 107,
    label = 'B3N3H6',
    molecule =
"""
1  B u0 p0 c0 {2,S} {6,S} {7,S}
2  N u0 p1 c0 {1,S} {3,S} {8,S}
3  B u0 p0 c0 {2,S} {4,S} {9,S}
4  N u0 p1 c0 {3,S} {5,S} {10,S}
5  B u0 p0 c0 {4,S} {6,S} {11,S}
6  N u0 p1 c0 {1,S} {5,S} {12,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([87.6366,117.583,142.6,162.603,190.846,208.7,231.763],'J/(mol*K)'),
        H298 = (-495.7,'kJ/mol'),
        S298 = (287.8,'J/(mol*K)'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (282.692,'J/(mol*K)'),
        E0 = (-510.185,'kJ/mol'),
    ),
    shortDesc = u"""B3H6N3""",
    longDesc = u"""Calculated with W1X-1""",
)

