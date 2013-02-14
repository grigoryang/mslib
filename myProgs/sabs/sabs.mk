MYSOURCE  =  

MYPROGS   = discardSimilarConformers selectCavities \
	buildDunbrackLibrary makeShettyLibrary jobDistributor \
	makePdbFromRotamerLibrary getInteractionGraph hbondRecovery makeHonigLibrary \
	genHeteroUniverse genHomoUniverseForSequence genHomoUniverse CAHTM clusterHeteroCandidates clusterCandidates tmHeteroRulesCreator tmRulesCreator \
	chainizeMissingLoops analyzeTMStructures convertToPdbNames getUniqueChains \
	createPhiPsi conformerEnergies selectBackBones binConformerLibrary \
	minimizeConformers removeBadConformers findDiscardedConformers rotamerEnergies buildCharmmSystemForPDB \
	seedAtoms prepareData prepareData_bbdep rankRotamers rankRotamersOnTestCavities repackSidechains addHydrogensFromPDB \
	plotEnergyWell plotchi1 plotchi2 EnergyWell generateEnergyWellInput collectRotamerStats getDunbrackConformerEnergies testDunbrackConformersOnCavities \
	maximumCover collectStatsFromRotamerLibrary minimiseRotamerLibrary createBackBone extractFragments peelWaters myOptionParser selectTestCavities evaluateConformerLib selectBestCavities \
	evaluateDunbrackLib evaluate5xDunbrackLib evaluate9xDunbrackLib wellEnergyStats collectAllRotlibs conformationEnergies getChi1Chi2MinimisedEnergies analyseHydrogenBonds evaluateOneConformerLibrary collectHBondStats collectHBondStats_Bins \
	collectAllHBondStats maximumCover_fixed_threshold multipleRepack multipleRepackWithClassification makeReducedDunbrack makeSamplingDataSet getCentreOfMass calibrateCavities multipleRepackTuning measureEnergy collectBestEnergies greedyOptimizer optimiseSideChains testDeadEndElimination homoTrimerStats interTerminiDistances interHelicalStats replicateTrimerBiounits getSelfEnergies getSelfEnergiesPerLevel standAloneSelfEnergy computeRMSD estimateGMEC rotamerCavityEnergies makePdbFormRotamerLibrary trimDunbrackLib measureRotamerEnergies singleRepackHbondRecovery computeSearchSize buildSequence buildRotamerLibraryFromDihedrals jansfirstprogram CCBupdate computeSideChainRMSD test reorientHelicesSystematically reorientHeteroHelicesSystematically buildCoiledCoil buildHeteroCoiledCoil clusterStructuresByRMSD getClusterCentroids measureInteractionEnergy\
	mutateSideChains repackBasedOnTemplate seperateHelicesAndMeasureEnergy collectInterHelicalHbondData testGreedyOptimizer iterativeXShiftRepacks clusterSolutions moveToCentreOfMass transformHetero transform getAtomsWithTooMuchDisplacement countAlternateConformations repackWithAlternateConformations \
	build69ALA prepGlyDimer antiCAHTM getDimerUniverseParameters trimerQsy qsy predictHelixDimer predictHelixTrimer visualizePacking stitchFtsB combineCCDTM generateHBondingUniverse generateHBondingTrimerUniverse connectWithFragments mutateAndRepackHelixTrimer mutateAndRepackHelixDimer mutateAndCompareHomoVsHeteroDimers homoUniverseData oldHomoUniverse genHomoMockAlaUniverse CAHTM_test measureDeltaEnergyOfDimer measureDeltaEnergyOfDimerWithRepack bensMeasureEnergyOfDimer dockFtsB screenFtsLFtsB predictFtsLFtsB createFtsLFtsBTrimer screenFtsB_gln_trp screenZipA myDesignSideChains listInteractionEnergies createHomoDimers helixBaseline getInterfacialResidues repackAndPrintTermEnergiesForFitting \


