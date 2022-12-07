# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary','DFT_QCI_thermo','CH','CHN','CHO','CHON','CN','NISTThermoLibrary','thermo_DFT_CCSDTF12_BAC','GRI-Mech3.0-N','Chlorinated_Hydrocarbons', 'CHOCl_G4'],
    reactionLibraries = ['Nitrogen_Dean_and_Bozzelli', 'CF2BrCl'], 
    seedMechanisms = [],
    kineticsDepositories = ['training'], 
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# Constraints on generated species
generatedSpeciesConstraints(
    allowed = ['input species', 'seed mechanisms', 'reaction libraries'],
    maximumCarbonAtoms = 4,
    #maximumOxygenAtoms = 5,
    maximumNitrogenAtoms = 4,
    #maximumSiliconAtoms = 0,
    #maximumSulfurAtoms = 0,
    #maximumHeavyAtoms = 3,
    maximumRadicalElectrons = 2,
    allowSingletO2 = False,
)

# List of species
species(
    label='CCl4',
    reactive=True,
    structure=adjacencyList(
        """
        1 Cl u0 p3 c0 {2,S}
        2 C  u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
        3 Cl u0 p3 c0 {2,S}
        4 Cl u0 p3 c0 {2,S}
        5 Cl u0 p3 c0 {2,S}
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

# Reaction systems
simpleReactor(
    temperature=(1200,'K'),
    pressure=(1.0,'bar'),
    initialMoleFractions={
        "CCl4": 0.05,
        "NH3": 0.15,
        "H2": 0.8,
    },
    terminationConversion={
        'NH3': 0.9,
    },
    terminationTime=(1e0,'s'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.05,
    toleranceInterruptSimulation=0.05,
    maximumEdgeSpecies=100000
)


# quantumMechanics(
#     software='mopac',
#     method='pm3',
#     # fileStore='QMfiles', # relative to where it is run. Defaults within the output folder.
#     scratchDirectory = None, # not currently used
#     onlyCyclics = True,
#     maxRadicalNumber = 0,
#     )

# pressureDependence(
#     method='modified strong collision',
#     maximumGrainSize=(0.5,'kcal/mol'),
#     minimumNumberOfGrains=250,
#     temperatures=(300,2000,'K',8),
#     pressures=(0.01,100,'bar',5),
#     interpolation=('Chebyshev', 6, 4),
# )

options(
    units='si',
    generateOutputHTML=True,
    generatePlots=True,
)
