# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary','DFT_QCI_thermo','CHOCl_G4','CHON_G4','BoronAllendorf', 'CH','CHN','CHO','CHON','CN','NISTThermoLibrary','thermo_DFT_CCSDTF12_BAC','GRI-Mech3.0-N'],
    reactionLibraries = ['Nitrogen_Dean_and_Bozzelli', 'NOx2018'],
    seedMechanisms = ['BoronAllendorf'],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['halogens',
                        # 'intra_H_migration',
                        # 'intra_NO2_ONO_conversion',
                        # 'intra_OH_migration',
                        # 'intra_substitutionCS_cyclization',
                        # 'intra_substitutionCS_isomerization',
                        # 'intra_substitutionS_cyclization',
                        # 'intra_substitutionS_isomerization',
                        # 'H_Abstraction',
                        # 'intra_H_migration',
                        # '1,2_NH3_elimination',
                        # '1,3_NH3_elimination',
                        # 'Retroene',
                        'default'],
    kineticsEstimator = 'rate rules',
)

# Constraints on generated species
generatedSpeciesConstraints(
    allowed = ['input species', 'seed mechanisms', 'reaction libraries'],
    #maximumCarbonAtoms = 4,
    maximumBoronAtoms = 4,
    maximumNitrogenAtoms = 4,
    # maximumRadicalElectrons = 1,
    allowSingletO2 = False,
)

# List of species
species(
    label='BCl3',
    reactive=True,
    structure=adjacencyList(
        """
        1 B u0 p0 c0 {2,S} {3,S} {4,S}
        2 Cl u0 p3 c0 {1,S}
        3 Cl u0 p3 c0 {1,S}
        4 Cl u0 p3 c0 {1,S}
        """),
)
species(
    label='NH3',
    reactive=True,
    structure=SMILES("N"),
)
species(
    label='H2',
    reactive=True,
    structure=adjacencyList(
        """
        1 H u0 p0 {2,S}
        2 H u0 p0 {1,S}
        """),
)
species(
    label='He',
    reactive=False,
    structure=SMILES("[He]"),
)
species(
    label='Cl3BNH3',
    reactive=True,
    structure=adjacencyList(
        """
        1 B u0 p0 c-1 {2,S} {3,S} {4,S} {5,S}
        2 N u0 p0 c+1 {1,S} {6,S} {7,S} {8,S}
        3 Cl u0 p3 c0 {1,S}
        4 Cl u0 p3 c0 {1,S}
        5 Cl u0 p3 c0 {1,S}
        6 H u0 p0 c0 {2,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {2,S}
        """),
)
forbidden(
    label='Boron1plus',
    structure=adjacencyListGroup(
            """
            1 B ux px c+1
            """
    )
)
forbidden(
    label='Boron1minus',
    structure=adjacencyListGroup(
            """
            1 B ux px c-1
            """
    )
)
forbidden(
    label='Boron2plus',
    structure=adjacencyListGroup(
            """
            1 B ux px c+2
            """
    )
)
forbidden(
    label='Boron2minus',
    structure=adjacencyListGroup(
            """
            1 B ux px c-2
            """
    )
)

# Reaction systems
simpleReactor(
    temperature=(1073,'K'),
    pressure=(0.01,'bar'),
    initialMoleFractions={
        "BCl3": 0.05,
        "NH3": 0.15,
        "H2": 0.8,
    },
    # terminationConversion={
    #     'NH3': 0.9,
    # },
    terminationTime=(1.0e1,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=1e-3,
    toleranceInterruptSimulation=1e3,
    maximumEdgeSpecies=10000
)

# pressureDependence(
#     method='modified strong collision',
#     maximumGrainSize=(0.5,'kcal/mol'),
#     minimumNumberOfGrains=250,
#     temperatures=(500,1300,'K',8),
#     pressures=(0.01,100,'bar',5),
#     interpolation=('Chebyshev', 6, 4),
# )

options(
    units='si',
    generateOutputHTML=True,
    generatePlots=True,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
)
