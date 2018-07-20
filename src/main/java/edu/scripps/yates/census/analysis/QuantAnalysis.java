package edu.scripps.yates.census.analysis;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.clustering.ProteinCluster;
import edu.scripps.yates.census.analysis.clustering.ProteinClusterUtils;
import edu.scripps.yates.census.analysis.util.KeyUtils;
import edu.scripps.yates.census.analysis.wrappers.IntegrationResultWrapper;
import edu.scripps.yates.census.analysis.wrappers.OutStatsLine;
import edu.scripps.yates.census.analysis.wrappers.SanXotAnalysisResult;
import edu.scripps.yates.census.read.model.IonCountRatio;
import edu.scripps.yates.census.read.model.IonSerie.IonSerieType;
import edu.scripps.yates.census.read.model.IsoRatio;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPSM;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPeptide;
import edu.scripps.yates.census.read.model.QuantifiedPSM;
import edu.scripps.yates.census.read.model.QuantifiedPeptide;
import edu.scripps.yates.census.read.model.interfaces.IsobaricQuantParser;
import edu.scripps.yates.census.read.model.interfaces.QuantParser;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.IonExclusion;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.dbindex.DBIndexInterface;
import edu.scripps.yates.dbindex.model.DBIndexSearchParams;
import edu.scripps.yates.utilities.grouping.GroupablePeptide;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.PAnalyzer;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.model.enums.AmountType;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class QuantAnalysis implements PropertyChangeListener {
	private static final Logger log = Logger.getLogger(QuantAnalysis.class);
	private boolean overrideFilesIfExists = true;
	private final List<QuantExperiment> quantExperiments = new ArrayList<QuantExperiment>();
	private final File workingFolder;
	private final String NL = System.getProperty("line.separator");
	private final QuantCondition condition1;
	private final QuantCondition condition2;
	private FileMappingResults fileMappingResults;
	private QuantParameters quantParameters = new QuantParameters();
	private final Map<String, List<String>> replicateAndExperimentNames = new THashMap<String, List<String>>();
	private DBIndexInterface dbIndex;
	private boolean ignorePTMs = true;
	private static final String QUANT_FOLDER = "quant";
	private SanXotAnalysisResult result;
	private final ANALYSIS_LEVEL_OUTCOME analysisOutCome;

	private int minAlignmentScore;

	private double minPercentajeOfsmilirarity;

	private int minConsecutiveLength;
	private Set<Set<String>> proteinAccClusters;
	private final QuantificationType quantType;
	private Boolean keepExperimentsSeparated;

	public static enum ANALYSIS_LEVEL_OUTCOME {
		PEPTIDE, PROTEIN, PROTEINGROUP, PROTEIN_CLUSTER, FORCED_CLUSTERS
	};

	/**
	 * Quantification analysis at protein group level
	 *
	 * @param workingPath
	 * @param condition1
	 * @param condition2
	 */
	public QuantAnalysis(QuantificationType quantType, String workingPath, QuantCondition condition1,
			QuantCondition condition2) {
		this.condition1 = condition1;
		this.condition2 = condition2;
		analysisOutCome = ANALYSIS_LEVEL_OUTCOME.PROTEINGROUP;
		workingFolder = createWorkingFolder(workingPath, analysisOutCome);
		this.quantType = quantType;
	}

	/**
	 * Quantification analysis at protein group level
	 *
	 * @param workingFolder
	 * @param condition1
	 * @param condition2
	 * @param analysisOutCome
	 */
	public QuantAnalysis(QuantificationType quantType, File workingFolder, QuantCondition condition1,
			QuantCondition condition2) {
		this.condition1 = condition1;
		this.condition2 = condition2;
		analysisOutCome = ANALYSIS_LEVEL_OUTCOME.PROTEINGROUP;
		this.workingFolder = createWorkingFolder(workingFolder, analysisOutCome);
		this.quantType = quantType;
	}

	/**
	 * Quantification analysis specifying the outcome level of the analysis
	 *
	 * @param workingPath
	 * @param condition1
	 * @param condition2
	 */
	public QuantAnalysis(QuantificationType quantType, String workingPath, QuantCondition condition1,
			QuantCondition condition2, ANALYSIS_LEVEL_OUTCOME analysisOutCome) {
		this.condition1 = condition1;
		this.condition2 = condition2;
		this.analysisOutCome = analysisOutCome;

		workingFolder = createWorkingFolder(workingPath, analysisOutCome);
		this.quantType = quantType;
	}

	/**
	 * Quantification analysis specifying the outcome level of the analysis
	 *
	 * @param workingFolder
	 * @param condition1
	 * @param condition2
	 * @param analysisOutCome
	 */
	public QuantAnalysis(QuantificationType quantType, File workingFolder, QuantCondition condition1,
			QuantCondition condition2, ANALYSIS_LEVEL_OUTCOME analysisOutCome) {
		this.condition1 = condition1;
		this.condition2 = condition2;
		this.analysisOutCome = analysisOutCome;
		this.workingFolder = createWorkingFolder(workingFolder, analysisOutCome);
		this.quantType = quantType;

	}

	/**
	 * @return the ignorePTMs
	 */
	public boolean isIgnorePTMs() {
		return ignorePTMs;
	}

	/**
	 * @param ignorePTMs
	 *            the ignorePTMs to set
	 */
	public void setIgnorePTMs(boolean ignorePTMs) {
		this.ignorePTMs = ignorePTMs;
	}

	private File createWorkingFolder(String workingPath, ANALYSIS_LEVEL_OUTCOME outcome) {
		final File folder = new File(workingPath + File.separator + QUANT_FOLDER + File.separator + outcome.name());
		if (!folder.exists())
			folder.mkdirs();
		return folder;
	}

	public static File createWorkingFolder(File workingFolder, ANALYSIS_LEVEL_OUTCOME outcome) {
		final File folder = new File(
				workingFolder.getAbsolutePath() + File.separator + QUANT_FOLDER + File.separator + outcome.name());
		if (!folder.exists())
			folder.mkdirs();
		return folder;
	}

	/**
	 * @return the keepExperimentsSeparated
	 */
	public Boolean getKeepExperimentsSeparated() {
		return keepExperimentsSeparated;
	}

	/**
	 * @param keepExperimentsSeparated
	 *            the keepExperimentsSeparated to set
	 */
	public void setKeepExperimentsSeparated(Boolean keepExperimentsSeparated) {
		this.keepExperimentsSeparated = keepExperimentsSeparated;
	}

	public void addQuantExperiment(QuantExperiment exp) {

		quantExperiments.add(exp);
	}

	/**
	 * @return the quantExperiments
	 */
	public List<QuantExperiment> getQuantExperiments() {
		return quantExperiments;
	}

	/**
	 * @return the workingPath
	 */
	public File getWorkingPath() {
		return workingFolder;
	}

	public void runSanXot() throws IOException {

		// clear static data from censusParsers
		// log.info("Clearing static data from parser");
		// QuantStaticMaps.clearInfo();

		if (fileMappingResults == null) {
			writeFiles();
		}

		final SanXotInterfaze sanxot = new SanXotInterfaze(fileMappingResults, quantParameters);

		if (keepExperimentsSeparated != null) {
			sanxot.setKeepExperimentsSeparated(keepExperimentsSeparated);
		}
		// sanxot.addPropertyChangeListener(this);
		// sanxot.execute();
		try {
			sanxot.analyze();
			result = sanxot.getResult();
		} catch (final InterruptedException e) {
			e.printStackTrace();
		} catch (final ExecutionException e) {
			e.printStackTrace();
		}
	}

	public FileMappingResults writeFiles() throws IOException {
		// assign db index to replicates (if exists)
		assignDBIndexToRuns();
		// read experiment and replicate names and stores at
		// replicateAndExperimentNames Map
		readExperimentAndReplicateNames();
		// write data file, which is the one with the ratios of the lower level
		writeDataFile();

		switch (analysisOutCome) {
		case PROTEINGROUP:
			writeRelationshipsFilesForProteinGroupOutcome();
			break;
		case PROTEIN:
			writeRelationshipsFilesForProteinOutcome();
			break;
		case FORCED_CLUSTERS:
			writeRelationshipsFilesForForcedProteinClustersOutcome();
			break;
		case PEPTIDE:
			writeRelationshipsFilesForPeptideOutcome();
			break;
		case PROTEIN_CLUSTER:
			if (minAlignmentScore == 0) {
				throw new IllegalArgumentException(
						"Minimum alignment score is needed. call to setMinAligmentScore(int min)");
			}
			if (Double.compare(minPercentajeOfsmilirarity, 0.0) == 0) {
				throw new IllegalArgumentException(
						"Minimum percentage of similarity is needed. call to setMinPercentajeOfsmilirarity(double min)");
			}
			if (minConsecutiveLength == 0) {
				throw new IllegalArgumentException(
						"Minimum consecutive length is needed. call to setMinConsecutiveLength(int consecutiveLength)");
			}
			writeRelationshipsFilesForProteinClustersOutcome(minAlignmentScore, minPercentajeOfsmilirarity,
					minConsecutiveLength);
		default:
			break;
		}

		fileMappingResults = new FileMappingResults(quantType, workingFolder, analysisOutCome,
				replicateAndExperimentNames);
		return fileMappingResults;
	}

	/**
	 * Iterates over the list of experiments and save the names of the
	 * experiments and replicates for its further use.
	 */
	private void readExperimentAndReplicateNames() {
		for (final QuantExperiment exp : quantExperiments) {
			final String experimentName = exp.getName();
			if (!replicateAndExperimentNames.containsKey(experimentName)) {
				replicateAndExperimentNames.put(experimentName, new ArrayList<String>());
			}
			for (final QuantReplicate rep : exp.getReplicates()) {
				final String replicateName = rep.getName();
				replicateAndExperimentNames.get(experimentName).add(replicateName);
			}
		}
	}

	private void assignDBIndexToRuns() {
		if (dbIndex != null) {
			log.info("Setting index to all replicates in the analysis");
			for (final QuantExperiment quantExperiment : quantExperiments) {
				final List<QuantReplicate> replicates = quantExperiment.getReplicates();
				for (final QuantReplicate quantReplicate : replicates) {
					final QuantParser parser = quantReplicate.getParser();

					if (parser instanceof IsobaricQuantParser) {
						// by default, remove B1 and Y1
						final Set<IonExclusion> ionExclusions = new THashSet<IonExclusion>();
						ionExclusions.add(new IonExclusion(IonSerieType.B, 1));
						ionExclusions.add(new IonExclusion(IonSerieType.Y, 1));
						((IsobaricQuantParser) parser).addIonExclusions(ionExclusions);
					}
					parser.setDbIndex(dbIndex);
				}
			}
		}
	}

	private void writeDataFile() throws IOException {

		final FileWriter dataFileWriter = new FileWriter(
				getWorkingPath().getAbsolutePath() + File.separator + FileMappingResults.DATA_FILE);
		dataFileWriter.write("#id\tX\tVcal" + NL);
		int numPSMsDiscarded = 0;
		try {
			for (final QuantExperiment exp : quantExperiments) {
				String expName = "";
				if (quantExperiments.size() > 1) {
					expName = exp.getName();
				}
				for (final QuantReplicate rep : exp.getReplicates()) {
					String repName = "";
					if (exp.getReplicates().size() > 1) {
						repName = rep.getName();
					}
					String expRepKey = "";
					if (!"".equals(repName)) {
						expRepKey = "_" + repName;
					}
					if (!"".equals(expName)) {
						expRepKey += "_" + expName;
					}
					final QuantParser parser = rep.getParser();
					final Collection<QuantifiedPSMInterface> quantifiedPSMs = parser.getPSMMap().values();
					for (final QuantifiedPSMInterface quantifiedPSM : quantifiedPSMs) {
						if (ignorePTMs && quantifiedPSM.containsPTMs()) {
							continue;
						}
						if (quantifiedPSM.isDiscarded()) {
							numPSMsDiscarded++;
							continue;
						}
						Double ratioValue = null;
						String key = null;
						Double fittingWeight = null;

						// in case of isobaric isotopologues, we will write a
						// data line per isobaric ratio in the PSM
						final Set<QuantRatio> nonInfinityRatios = quantifiedPSM.getNonInfinityRatios();
						if (nonInfinityRatios.isEmpty()) {
							// skip this one
							continue;
						}
						if (quantifiedPSM instanceof IsobaricQuantifiedPSM) {
							for (final QuantRatio ratio : nonInfinityRatios) {
								ratioValue = ratio.getLog2Ratio(condition1, condition2);
								if (ratio instanceof IsoRatio) {

									final IsoRatio isoRatio = (IsoRatio) ratio;
									key = KeyUtils.getIonKey(isoRatio,
											((IsobaricQuantifiedPSM) quantifiedPSM).getPeptide(), true) + expRepKey;

									fittingWeight = null;

									switch (quantType) {
									case ISOTOPOLOGUES:
										fittingWeight = (isoRatio.getMaxIntensity())
												/ Math.sqrt(isoRatio.getMass(QuantificationLabel.LIGHT));
										dataFileWriter.write(key + "\t" + ratioValue + "\t" + fittingWeight + "\n");
										break;

									default:
										throw new IllegalArgumentException("Quant type " + quantType
												+ " is not suitable for Isobaric Isotopologues. Use instead "
												+ QuantificationType.ISOTOPOLOGUES);
									}
								} else if (ratio instanceof IonCountRatio) {
									throw new IllegalArgumentException(
											"Ion count ratios is not suitable for ratio integration since they are peptide node level.");
								}
								// TODO
								// double qualityMeasurement =
								// ratio.getMaxPeak() *
								// ratio.getMaxPeak()
								// / (ratio.getIon1().getMass() *
								// ratio.getIon1().getMass());

							}
						} else {
							switch (quantType) {
							case iTRAQ:
								fittingWeight = quantifiedPSM.getMaxPeak();
								break;

							case SILAC:
								// fittingWeight =
								// QuantUtil.getRegressionFactor(quantifiedPSM.getAmounts());
								fittingWeight = QuantUtils.getMaxAmountValueByAmountType(quantifiedPSM.getAmounts(),
										AmountType.AREA);
								if (fittingWeight == null) {
									log.info("no regression factor");
								}
								break;
							case UNKNOWN:
								// the PSM has to have only one ratio, and it
								// has to have associated a confidence score
								// which is the weight
								Double.valueOf(quantifiedPSM.getRatios().iterator().next()
										.getAssociatedConfidenceScore().getValue());
							default:
								throw new IllegalArgumentException("Quant type " + quantType
										+ " is not supported with this analysis configuration");
							}

							key = KeyUtils.getSpectrumKey(quantifiedPSM, true) + expRepKey;
							// in case of not having isobaric isotopologues, we
							// have one ratio per PSM in the replicate, not
							// matters if it is comming from a TMT, where we
							// have more than one ratio per PSM, because we will
							// write each ratio in different replicates
							if (quantifiedPSM instanceof QuantifiedPSM) {
								QuantRatio validRatio = null;

								if (quantParameters.getRatioName() != null
										&& !"".equals(quantParameters.getRatioName())) {
									validRatio = QuantUtils.getRatioByName(quantifiedPSM,
											quantParameters.getRatioName());
								} else {
									validRatio = QuantUtils.getRatioValidForIntegrationAnalysis(quantifiedPSM);
								}
								if (validRatio != null) {
									ratioValue = validRatio.getLog2Ratio(condition1, condition2);
									if (ratioValue == null || Double.isInfinite(ratioValue)
											|| Double.isNaN(ratioValue)) {
										// do not print
										continue;
									}
								} else {
									// dont print
									continue;
								}
							} else {
								ratioValue = nonInfinityRatios.iterator().next().getLog2Ratio(condition1, condition2);
							}

							dataFileWriter.write(key + "\t" + ratioValue + "\t" + fittingWeight + "\n");
						}

					}

				}
			}

		} finally {
			log.info(numPSMsDiscarded + " PSMs were tagged as discarded and will not be considered in the analysis");
			if (dataFileWriter != null) {
				dataFileWriter.close();
			}
		}
	}

	private void writeRelationshipsFilesForProteinGroupOutcome() throws IOException {
		writeIonToSpectrumMap();
		writeSpectrumToPeptideExperimentReplicateMap();
		writePeptideExperimentReplicateToProteinGroupExperimentReplicateMap();
		try {
			// writeProteinExperimentReplicateToProteinExperimentMap();
			writeProteinGroupExperimentReplicateToProteinGroupExperimentMap();
			try {
				// writeProteinExperimentToProteinMap();
				writeProteinGroupExperimentToProteinGroupMap();
			} catch (final IllegalArgumentException e) {
				log.info(e.getMessage());
				// "There is not more than 1 experiment"
				// do nothing and perform the last step: proteinToAll
			}

		} catch (final IllegalArgumentException e) {
			log.info(e.getMessage());
			// "There is not any experiment with some replicates"
			try {
				// writeProteinExperimentToProteinMap();
				writeProteinGroupExperimentReplicateToProteinGroupExperimentMap();
				writeProteinGroupExperimentToProteinGroupMap();
			} catch (final IllegalArgumentException e2) {
				log.info(e2.getMessage());
				// "There is not more than 1 experiment"
				// do nothing and perform the last step: proteinToAll
			}
		}

		// writeProteinToAllMap();

		writeProteinGroupToAllMap();

	}

	private void writeRelationshipsFilesForForcedProteinClustersOutcome() throws IOException {
		writeIonToSpectrumMap();

		writeSpectrumToPeptideExperimentReplicateMap();
		try {
			writePeptideExperimentReplicateToPeptideExperimentMap();
			try {
				writePeptideExperimentToPeptideMap();
			} catch (final IllegalArgumentException e) {
				log.info(e.getMessage());
				// "There is not more than 1 experiment"
				// do nothing and perform the last step: proteinToAll
			}

		} catch (final IllegalArgumentException e) {
			log.info(e.getMessage());
			// "There is not any experiment with some replicates"
			try {
				writePeptideExperimentToPeptideMap();
			} catch (final IllegalArgumentException e2) {
				log.info(e2.getMessage());
				// "There is not more than 1 experiment"
				// do nothing and perform the last step: proteinToAll
			}
		}
		final Set<ProteinCluster> proteinClusterMap = writePeptideToForcedProteinClusterMap(proteinAccClusters);
		writeProteinClusterToAllMap(minAlignmentScore, minPercentajeOfsmilirarity, minConsecutiveLength,
				proteinClusterMap);

	}

	private void writeRelationshipsFilesForProteinClustersOutcome(int minAlignmentScore,
			double minPercentajeOfsmilirarity, int minConsecutiveLength) throws IOException {

		writeIonToSpectrumMap();
		writeSpectrumToPeptideExperimentReplicateMap();
		try {
			writePeptideExperimentReplicateToPeptideExperimentMap();
			try {
				writePeptideExperimentToPeptideMap();
			} catch (final IllegalArgumentException e) {
				log.info(e.getMessage());
				// "There is not more than 1 experiment"
				// do nothing and perform the last step: proteinToAll
			}

		} catch (final IllegalArgumentException e) {
			log.info(e.getMessage());
			// "There is not any experiment with some replicates"
			try {
				writePeptideExperimentToPeptideMap();
			} catch (final IllegalArgumentException e2) {
				log.info(e2.getMessage());
				// "There is not more than 1 experiment"
				// do nothing and perform the last step: proteinToAll
			}
		}
		final Set<ProteinCluster> proteinClusterMap = writePeptideToProteinClusterMap(minAlignmentScore,
				minPercentajeOfsmilirarity, minConsecutiveLength);
		writeProteinClusterToAllMap(minAlignmentScore, minPercentajeOfsmilirarity, minConsecutiveLength,
				proteinClusterMap);

	}

	private void writeRelationshipsFilesForProteinOutcome() throws IOException {
		if (quantType == QuantificationType.ISOTOPOLOGUES) {
			writeIonToSpectrumMap();
		}
		writeSpectrumToPeptideExperimentReplicateMap();
		writePeptideToProteinMap();
		try {
			writeProteinExperimentReplicateToProteinExperimentMap();
			try {
				writeProteinExperimentToProteinMap();
			} catch (final IllegalArgumentException e) {
				log.info(e.getMessage());
				// "There is not more than 1 experiment"
				// do nothing and perform the last step: proteinToAll
			}

		} catch (final IllegalArgumentException e) {
			log.info(e.getMessage());
			// "There is not any experiment with some replicates"
			try {
				writeProteinExperimentToProteinMap();
			} catch (final IllegalArgumentException e2) {
				log.info(e2.getMessage());
				// "There is not more than 1 experiment"
				// do nothing and perform the last step: proteinToAll
			}
		}
		writeProteinToAllMap();
	}

	private void writeRelationshipsFilesForPeptideOutcome() throws IOException {
		writeIonToSpectrumMap();
		writeSpectrumToPeptideExperimentReplicateMap();

		try {
			writePeptideExperimentReplicateToPeptideExperimentMap();
			try {
				writePeptideExperimentToPeptideMap();
			} catch (final IllegalArgumentException e) {
				log.info(e.getMessage());
				// "There is not more than 1 experiment"
				// do nothing and perform the last step: proteinToAll
			}

		} catch (final IllegalArgumentException e) {
			log.info(e.getMessage());
			// "There is not any experiment with some replicates"
			try {
				writePeptideExperimentToPeptideMap();
			} catch (final IllegalArgumentException e2) {
				log.info(e2.getMessage());
				// "There is not more than 1 experiment"
				// do nothing and perform the last step: proteinToAll
			}
		}
		writePeptideToAllMap();
	}

	/**
	 * ALL --> ACC1<br>
	 * ALL --> ACC2<br>
	 * ALL --> ACC3<br>
	 *
	 * @throws IOException
	 */
	private void writeProteinToAllMap() throws IOException {
		final FileWriter writer = new FileWriter(
				getWorkingPath().getAbsolutePath() + File.separator + FileMappingResults.PROTEIN_TO_ALL_6);
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();
		final String all = "all";
		for (final QuantExperiment exp : quantExperiments) {
			for (final QuantReplicate rep : exp.getReplicates()) {
				final QuantParser parser = rep.getParser();
				final Map<String, QuantifiedProteinInterface> quantifiedProteinMap = parser.getProteinMap();
				for (final String proteinKey : quantifiedProteinMap.keySet()) {
					if (map.containsKey(all)) {
						map.get(all).add(proteinKey);
					} else {
						final Set<String> set = new THashSet<String>();
						set.add(proteinKey);
						map.put(all, set);
					}
				}
			}
		}
		final String header = "all" + "\t" + "acc" + "\t" + "all --> protein";
		writeMapToFile(header, map, writer);
	}

	/**
	 * PEP1 --> CLUSTER1<br>
	 * PEP2 --> CLUSTER1<br>
	 * PEP3 --> CLUSTER2<br>
	 * PEP3 --> CLUSTER2<br>
	 * <br>
	 *
	 * @return
	 * @throws IOException
	 */
	private Set<ProteinCluster> writePeptideToForcedProteinClusterMap(Set<Set<String>> forcedProteinClusters)
			throws IOException {

		final File file = new File(
				getWorkingPath().getAbsolutePath() + File.separator + FileMappingResults.PEPTIDE_TO_PROTEIN_CLUSTER_5);

		// create the clusters according to the parameter
		final Set<ProteinCluster> proteinClusters = new THashSet<ProteinCluster>();
		for (final Set<String> proteinAccs : forcedProteinClusters) {
			final List<String> proteinAccList = new ArrayList<String>();
			proteinAccList.addAll(proteinAccs);
			Collections.sort(proteinAccList);
			final StringBuilder proteinClusterKey = new StringBuilder();
			for (final String proteinAcc : proteinAccList) {
				proteinClusterKey.append(proteinAcc + ":");
			}
			final ProteinCluster cluster = new ProteinCluster();
			cluster.setProteinClusterKey(proteinClusterKey.toString());
			proteinClusters.add(cluster);
			for (final String proteinACC : proteinAccs) {
				for (final QuantExperiment exp : quantExperiments) {
					for (final QuantReplicate rep : exp.getReplicates()) {
						final QuantParser parser = rep.getParser();
						final Map<String, QuantifiedProteinInterface> proteinMap = parser.getProteinMap();
						if (proteinMap.containsKey(proteinACC)) {
							final QuantifiedProteinInterface quantifiedProtein = proteinMap.get(proteinACC);
							cluster.addProtein(quantifiedProtein);
						}
					}
				}
			}
		}

		if (!overrideFilesIfExists && file.exists()) {
			return proteinClusters;
		}
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();

		for (final ProteinCluster proteinCluster : proteinClusters) {
			final String proteinClusterKey = proteinCluster.getProteinClusterKey();
			final Set<QuantifiedPeptideInterface> quantifiedPeptides = proteinCluster.getPeptideSet();
			for (final QuantifiedPeptideInterface quantifiedPeptide : quantifiedPeptides) {
				final String peptideKey = quantifiedPeptide.getKey();
				if (map.containsKey(proteinClusterKey)) {
					map.get(proteinClusterKey).add(peptideKey);
				} else {
					final Set<String> set = new THashSet<String>();
					set.add(peptideKey);
					map.put(proteinClusterKey, set);
				}
			}
		}
		final FileWriter writer = new FileWriter(file);
		final String header = "pep" + "\t" + "proteinCluster" + "\t" + "proteinCluster --> peptide";
		writeMapToFile(header, map, writer);
		return proteinClusters;
	}

	/**
	 * PEP1 --> CLUSTER1<br>
	 * PEP2 --> CLUSTER1<br>
	 * PEP3 --> CLUSTER2<br>
	 * PEP3 --> CLUSTER2<br>
	 * <br>
	 * Parameters for the clustering building:<br>
	 *
	 * @param minAlignmentScore
	 * @param minPercentajeOfsmilirarity
	 * @param minConsecutiveLength
	 * @return
	 * @throws IOException
	 */
	private Set<ProteinCluster> writePeptideToProteinClusterMap(int minAlignmentScore,
			double minPercentajeOfsmilirarity, int minConsecutiveLength) throws IOException {

		final File file = new File(
				getWorkingPath().getAbsolutePath() + File.separator + FileMappingResults.PEPTIDE_TO_PROTEIN_CLUSTER_5);

		final List<QuantifiedPSMInterface> quantPSMs = new ArrayList<QuantifiedPSMInterface>();
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();
		for (final QuantExperiment exp : quantExperiments) {

			for (final QuantReplicate rep : exp.getReplicates()) {

				final QuantParser parser = rep.getParser();
				// in this case, get all proteins and construct protein groups.
				// Then, asign in the map to peptides
				final Map<String, QuantifiedPSMInterface> psmMap = parser.getPSMMap();

				for (final QuantifiedPSMInterface quantPSM : psmMap.values()) {
					quantPSMs.add(quantPSM);
				}

			}
		}

		// create the peptides from all PSMs.
		final Map<String, QuantifiedPeptideInterface> peptideMap = QuantifiedPeptide.getQuantifiedPeptides(quantPSMs);

		final Set<ProteinCluster> proteinClusters = ProteinClusterUtils.getProteinClusters(peptideMap,
				minAlignmentScore, minPercentajeOfsmilirarity, minConsecutiveLength);
		if (!overrideFilesIfExists && file.exists()) {
			return proteinClusters;
		}
		final FileWriter writer = new FileWriter(file);
		for (final ProteinCluster proteinCluster : proteinClusters) {
			final String proteinClusterKey = proteinCluster.getProteinClusterKey();

			final Set<QuantifiedPeptideInterface> quantifiedPeptides = proteinCluster.getPeptideSet();
			for (final QuantifiedPeptideInterface quantifiedPeptide : quantifiedPeptides) {
				final String peptideKey = quantifiedPeptide.getKey();
				if (map.containsKey(proteinClusterKey)) {
					map.get(proteinClusterKey).add(peptideKey);
				} else {
					final Set<String> set = new THashSet<String>();
					set.add(peptideKey);
					map.put(proteinClusterKey, set);
				}
			}
		}

		final String header = "pep" + "\t" + "proteinCluster" + "\t" + "proteinCluster --> peptide";
		writeMapToFile(header, map, writer);
		return proteinClusters;
	}

	/**
	 * ALL --> PEP1<br>
	 * ALL --> PEP2<br>
	 * ALL --> PEP3<br>
	 *
	 * @throws IOException
	 */
	private void writePeptideToAllMap() throws IOException {
		final File file = new File(
				getWorkingPath().getAbsolutePath() + File.separator + FileMappingResults.PEPTIDE_TO_ALL_5);
		if (!overrideFilesIfExists && file.exists()) {
			return;
		}
		final FileWriter writer = new FileWriter(file);
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();
		final String all = "all";
		for (final QuantExperiment exp : quantExperiments) {
			for (final QuantReplicate rep : exp.getReplicates()) {
				final QuantParser parser = rep.getParser();
				final Map<String, Set<String>> peptideToSpectraMap2 = parser.getPeptideToSpectraMap();
				for (final String peptideKey : peptideToSpectraMap2.keySet()) {
					if (map.containsKey(all)) {
						map.get(all).add(peptideKey);
					} else {
						final Set<String> set = new THashSet<String>();
						set.add(peptideKey);
						map.put(all, set);
					}
				}
			}
		}
		final String header = "all" + "\t" + "sequence+charge" + "\t" + "all --> peptide";
		writeMapToFile(header, map, writer);
	}

	/**
	 * ALL --> ProteinCluster1<br>
	 * ALL --> ProteinCluster2<br>
	 * ALL --> ProteinCluster3<br>
	 *
	 * @param proteinClusterMap
	 *
	 * @return
	 *
	 * @throws IOException
	 */
	private void writeProteinClusterToAllMap(int minAlignmentScore, double minPercentajeOfsmilirarity,
			int minConsecutiveLength, Set<ProteinCluster> proteinClusters) throws IOException {
		final File file = new File(
				getWorkingPath().getAbsolutePath() + File.separator + FileMappingResults.PROTEIN_CLUSTER_TO_ALL_6);
		if (!overrideFilesIfExists && file.exists()) {
			return;
		}
		final FileWriter writer = new FileWriter(file);
		final String all = "all";
		final List<GroupableProtein> groupableProteins = new ArrayList<GroupableProtein>();
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();
		for (final QuantExperiment exp : quantExperiments) {

			for (final QuantReplicate rep : exp.getReplicates()) {

				final QuantParser parser = rep.getParser();
				// in this case, get all proteins and construct protein groups.
				// Then, asign in the map to peptides
				final Map<String, QuantifiedProteinInterface> proteinMap = parser.getProteinMap();

				for (final QuantifiedProteinInterface quantProtein : proteinMap.values()) {
					groupableProteins.add(quantProtein);
				}

			}
		}

		final Set<String> set = new THashSet<String>();
		for (final ProteinCluster proteinCluster : proteinClusters) {
			final String proteinClusterKey = proteinCluster.getProteinClusterKey();
			set.add(proteinClusterKey);
		}

		map.put(all, set);
		final String header = "all" + "\t" + "proteinCluster" + "\t" + "all --> proteinCluster";
		writeMapToFile(header, map, writer);

	}

	private void writeProteinGroupToAllMap() throws IOException {
		final File file = new File(
				getWorkingPath().getAbsolutePath() + File.separator + FileMappingResults.PROTEINGROUP_TO_ALL_6);
		if (!overrideFilesIfExists && file.exists()) {
			return;
		}
		final FileWriter writer = new FileWriter(file);
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();
		final String all = "all";
		for (final QuantExperiment exp : quantExperiments) {
			for (final QuantReplicate rep : exp.getReplicates()) {
				final QuantParser parser = rep.getParser();
				final List<GroupableProtein> groupableProteins = new ArrayList<GroupableProtein>();
				groupableProteins.addAll(parser.getProteinMap().values());
				final List<ProteinGroup> proteinGroups = getProteinGroups(groupableProteins);
				for (final ProteinGroup proteinGroup : proteinGroups) {
					final String proteinKey = KeyUtils.getGroupKey(proteinGroup);
					if (map.containsKey(all)) {
						map.get(all).add(proteinKey);
					} else {
						final Set<String> set = new THashSet<String>();
						set.add(proteinKey);
						map.put(all, set);
					}
				}
			}
		}
		final String header = "all" + "\t" + "acc" + "\t" + "all --> protein";
		writeMapToFile(header, map, writer);
	}

	/**
	 * ACC1 --> ACC1_EXP1<br>
	 * ACC1 --> ACC1_EXP2<br>
	 * ACC2 --> ACC2_EXP1<br>
	 * ACC3 --> ACC3_EXP2<br>
	 *
	 * @throws IOException
	 */
	private void writeProteinExperimentToProteinMap() throws IOException {
		if (quantExperiments.size() == 1)
			throw new IllegalArgumentException("There is not more than 1 experiment");

		final File file = new File(getWorkingPath().getAbsolutePath() + File.separator
				+ FileMappingResults.PROTEIN_EXPERIMENT_TO_PROTEIN_5);
		if (!overrideFilesIfExists && file.exists()) {
			return;
		}
		final FileWriter writer = new FileWriter(file);
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();

		for (final QuantExperiment exp : quantExperiments) {
			String expName = "";
			if (quantExperiments.size() > 1) {
				expName = exp.getName();
			}
			for (final QuantReplicate rep : exp.getReplicates()) {
				final QuantParser parser = rep.getParser();
				final Map<String, QuantifiedProteinInterface> quantifiedProteinMap = parser.getProteinMap();
				for (final String proteinKey : quantifiedProteinMap.keySet()) {
					if (map.containsKey(proteinKey)) {
						map.get(proteinKey).add(proteinKey + "_" + expName);
					} else {
						final Set<String> set = new THashSet<String>();
						set.add(proteinKey + "_" + expName);
						map.put(proteinKey, set);
					}
				}
			}
		}
		final String header = "acc" + "\t" + "acc+experiment" + "\t" + "protein --> protein-experiment";
		writeMapToFile(header, map, writer);

	}

	/**
	 * PEP1 --> PEP1_EXP1<br>
	 * PEP1 --> PEP1_EXP2<br>
	 * PEP2 --> PEP2_EXP1<br>
	 * PEP3 --> PEP3_EXP2<br>
	 *
	 * @throws IOException
	 */
	private void writePeptideExperimentToPeptideMap() throws IOException {
		if (quantExperiments.size() == 1)
			throw new IllegalArgumentException("There is not more than 1 experiment");

		final File file = new File(getWorkingPath().getAbsolutePath() + File.separator
				+ FileMappingResults.PEPTIDE_EXPERIMENT_TO_PEPTIDE_4);
		if (!overrideFilesIfExists && file.exists()) {
			return;
		}
		final FileWriter writer = new FileWriter(file);
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();

		for (final QuantExperiment exp : quantExperiments) {
			String expName = "";
			if (quantExperiments.size() > 1) {
				expName = exp.getName();
			}
			for (final QuantReplicate rep : exp.getReplicates()) {
				final QuantParser parser = rep.getParser();
				for (final String peptideKey : parser.getPeptideToSpectraMap().keySet()) {
					if (map.containsKey(peptideKey)) {
						map.get(peptideKey).add(peptideKey + "_" + expName);
					} else {
						final Set<String> set = new THashSet<String>();
						set.add(peptideKey + "_" + expName);
						map.put(peptideKey, set);
					}
				}
			}
		}
		final String header = "sequence+charge" + "\t" + "sequence+charge+experiment" + "\t"
				+ "peptide --> peptide-experiment";
		writeMapToFile(header, map, writer);

	}

	private void writeProteinGroupExperimentToProteinGroupMap() throws IOException {
		if (quantExperiments.size() == 1)
			throw new IllegalArgumentException("There is not more than 1 experiment");

		final File file = new File(getWorkingPath().getAbsolutePath() + File.separator
				+ FileMappingResults.PROTEINGROUP_EXPERIMENT_TO_PROTEINGROUP_5);
		if (!overrideFilesIfExists && file.exists()) {
			return;
		}
		final FileWriter writer = new FileWriter(file);
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();

		for (final QuantExperiment exp : quantExperiments) {
			String expName = "";
			if (quantExperiments.size() > 1) {
				expName = exp.getName();
			}
			for (final QuantReplicate rep : exp.getReplicates()) {
				final QuantParser parser = rep.getParser();
				final List<GroupableProtein> groupableProteins = new ArrayList<GroupableProtein>();
				groupableProteins.addAll(parser.getProteinMap().values());
				final List<ProteinGroup> proteinGroups = getProteinGroups(groupableProteins);
				for (final ProteinGroup proteinGroup : proteinGroups) {
					final String proteinGroupKey = KeyUtils.getGroupKey(proteinGroup);
					if (map.containsKey(proteinGroupKey)) {
						map.get(proteinGroupKey).add(proteinGroupKey + "_" + expName);
					} else {
						final Set<String> set = new THashSet<String>();
						set.add(proteinGroupKey + "_" + expName);
						map.put(proteinGroupKey, set);
					}
				}
			}
		}
		final String header = "acc" + "\t" + "acc+experiment" + "\t" + "proteinGroup --> proteinGroup-experiment";
		writeMapToFile(header, map, writer);

	}

	/**
	 * This step is only valid if some experiment has replicates. Otherwise, an
	 * exception will be thrown.<br>
	 * ACC_EXP1 --> ACC_REP1_EXP1<br>
	 * ACC_EXP1 --> ACC_REP2_EXP1<br>
	 * ACC_EXP2 --> ACC_REP1_EXP2<br>
	 * ACC_EXP2 --> ACC_REP2_EXP2<br>
	 *
	 * @throws IOException
	 */
	private void writeProteinExperimentReplicateToProteinExperimentMap() throws IOException {

		boolean someReplicates = false;
		for (final QuantExperiment exp : quantExperiments) {
			if (exp.getReplicates().size() > 1)
				someReplicates = true;
		}
		if (!someReplicates)
			throw new IllegalArgumentException("There is not any experiment with some replicates");

		final File file = new File(getWorkingPath().getAbsolutePath() + File.separator
				+ FileMappingResults.PROTEIN_REPLICATE_EXPERIMENT_TO_PROTEIN_EXPERIMENT_4);
		if (!overrideFilesIfExists && file.exists()) {
			return;
		}
		final FileWriter writer = new FileWriter(file);
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();

		for (final QuantExperiment exp : quantExperiments) {
			String expName = "";
			if (quantExperiments.size() > 1) {
				expName = exp.getName();
			}
			for (final QuantReplicate rep : exp.getReplicates()) {
				String repName = "";
				if (exp.getReplicates().size() > 1) {
					repName = rep.getName();
				}
				String expRepKey = "";
				if (!"".equals(repName)) {
					expRepKey = "_" + repName;
				}
				if (!"".equals(expName)) {
					expRepKey += "_" + expName;
				}
				final QuantParser parser = rep.getParser();
				final Map<String, QuantifiedProteinInterface> quantifiedProteinMap = parser.getProteinMap();
				for (final String proteinKey : quantifiedProteinMap.keySet()) {
					if (map.containsKey(proteinKey + "_" + expName)) {
						map.get(proteinKey + "_" + expName).add(proteinKey + expRepKey);
					} else {
						final Set<String> set = new THashSet<String>();
						set.add(proteinKey + expRepKey);
						map.put(proteinKey + "_" + expName, set);
					}
				}

			}
		}
		final String header = "acc+experiment" + "\t" + "acc+replicate+experiment" + "\t"
				+ "protein-experiment --> protein-replicate-experiment";
		writeMapToFile(header, map, writer);

	}

	/**
	 * This step is only valid if some experiment has replicates. Otherwise, an
	 * exception will be thrown.<br>
	 * PEP_EXP1 --> PEP_REP1_EXP1<br>
	 * PEP_EXP1 --> PEP_REP2_EXP1<br>
	 * PEP_EXP2 --> PEP_REP1_EXP2<br>
	 * PEP_EXP2 --> PEP_REP2_EXP2<br>
	 *
	 * @throws IOException
	 */
	private void writePeptideExperimentReplicateToPeptideExperimentMap() throws IOException {

		boolean someReplicates = false;
		for (final QuantExperiment exp : quantExperiments) {
			if (exp.getReplicates().size() > 1)
				someReplicates = true;
		}
		if (!someReplicates)
			throw new IllegalArgumentException("There is not any experiment with some replicates");

		final File file = new File(getWorkingPath().getAbsolutePath() + File.separator
				+ FileMappingResults.PEPTIDE_REPLICATE_EXPERIMENT_TO_PEPTIDE_EXPERIMENT_3);
		if (!overrideFilesIfExists && file.exists()) {
			return;
		}
		final FileWriter writer = new FileWriter(file);
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();

		for (final QuantExperiment exp : quantExperiments) {
			String expName = "";
			if (quantExperiments.size() > 1) {
				expName = exp.getName();
			}
			for (final QuantReplicate rep : exp.getReplicates()) {
				String repName = "";
				if (exp.getReplicates().size() > 1) {
					repName = rep.getName();
				}
				String expRepKey = "";
				if (!"".equals(repName)) {
					expRepKey = "_" + repName;
				}
				if (!"".equals(expName)) {
					expRepKey += "_" + expName;
				}
				final QuantParser parser = rep.getParser();
				for (final String peptideKey : parser.getPeptideToSpectraMap().keySet()) {
					if (map.containsKey(peptideKey + "_" + expName)) {
						map.get(peptideKey + "_" + expName).add(peptideKey + expRepKey);
					} else {
						final Set<String> set = new THashSet<String>();
						set.add(peptideKey + expRepKey);
						map.put(peptideKey + "_" + expName, set);
					}
				}
			}
		}
		final String header = "sequence+charge+experiment" + "\t" + "sequence+charge+replicate+experiment" + "\t"
				+ "peptide-replicate-experiment --> peptide-experiment";
		writeMapToFile(header, map, writer);

	}

	private void writeProteinGroupExperimentReplicateToProteinGroupExperimentMap() throws IOException {

		boolean someReplicates = false;
		for (final QuantExperiment exp : quantExperiments) {
			if (exp.getReplicates().size() > 1)
				someReplicates = true;
		}
		if (!someReplicates)
			throw new IllegalArgumentException("There is not any experiment with some replicates");

		final File file = new File(getWorkingPath().getAbsolutePath() + File.separator
				+ FileMappingResults.PROTEINGROUP_REPLICATE_EXPERIMENT_TO_PROTEINGROUP_EXPERIMENT_4);
		if (!overrideFilesIfExists && file.exists()) {
			return;
		}
		final FileWriter writer = new FileWriter(file);
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();

		for (final QuantExperiment exp : quantExperiments) {
			String expName = "";
			if (quantExperiments.size() > 1) {
				expName = exp.getName();
			}
			for (final QuantReplicate rep : exp.getReplicates()) {
				String repName = "";
				if (exp.getReplicates().size() > 1) {
					repName = rep.getName();
				}
				String expRepKey = "";
				if (!"".equals(repName)) {
					expRepKey = "_" + repName;
				}
				if (!"".equals(expName)) {
					expRepKey += "_" + expName;
				}
				final QuantParser parser = rep.getParser();
				final List<GroupableProtein> groupableProteins = new ArrayList<GroupableProtein>();
				groupableProteins.addAll(parser.getProteinMap().values());
				final List<ProteinGroup> proteinGroups = getProteinGroups(groupableProteins);
				for (final ProteinGroup proteinGroup : proteinGroups) {
					final String proteinGroupKey = KeyUtils.getGroupKey(proteinGroup);
					if (map.containsKey(proteinGroupKey + "_" + expName)) {
						map.get(proteinGroupKey + "_" + expName).add(proteinGroupKey + expRepKey);
					} else {
						final Set<String> set = new THashSet<String>();
						set.add(proteinGroupKey + expRepKey);
						map.put(proteinGroupKey + "_" + expName, set);
					}
				}

			}
		}
		final String header = "acc+experiment" + "\t" + "acc+replicate+experiment" + "\t"
				+ "proteinGroup-experiment --> proteinGroup-replicate-experiment";
		writeMapToFile(header, map, writer);

	}

	/**
	 * ACC1_REP1_EXP1 --> PEPTIDEA_3_REP1_EXP1<br>
	 * ACC1_REP1_EXP1 --> PEPTIDEB_2_REP1_EXP1<br>
	 * ACC1_REP2_EXP1 --> PEPTIDEA_3_REP2_EXP1<br>
	 * ACC1_REP1_EXP2 --> PEPTIDEA_3_REP1_EXP2<br>
	 * ACC1_REP1_EXP2 --> PEPTIDEB_2_REP1_EXP2<br>
	 * ACC1_REP2_EXP2 --> PEPTIDEA_3_REP2_EXP2<br>
	 *
	 * @throws IOException
	 */
	private void writePeptideExperimentReplicateToProteinGroupExperimentReplicateMap() throws IOException {
		final File file = new File(getWorkingPath().getAbsolutePath() + File.separator
				+ FileMappingResults.PEPTIDE_REPLICATE_EXPERIMENT_TO_PROTEINGROUP_REPLICATE_EXPERIMENT_3);
		if (!overrideFilesIfExists && file.exists()) {
			return;
		}
		final FileWriter writer = new FileWriter(file);
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();

		for (final QuantExperiment exp : quantExperiments) {
			String expName = "";
			if (quantExperiments.size() > 1) {
				expName = exp.getName();
			}
			for (final QuantReplicate rep : exp.getReplicates()) {
				String repName = "";
				if (exp.getReplicates().size() > 1) {
					repName = rep.getName();
				}
				String expRepKey = "";
				if (!"".equals(repName)) {
					expRepKey = "_" + repName;
				}
				if (!"".equals(expName)) {
					expRepKey += "_" + expName;
				}
				final QuantParser parser = rep.getParser();
				// in this case, get all proteins and construct protein groups.
				// Then, asign in the map to peptides
				final List<GroupableProtein> groupableProteins = new ArrayList<GroupableProtein>();
				groupableProteins.addAll(parser.getProteinMap().values());
				final List<ProteinGroup> proteinGroups = getProteinGroups(groupableProteins);
				final Map<String, Set<String>> proteinGroupToPeptideMap2 = getProteinGroupToPeptideMap(proteinGroups);

				mergeMaps(map, proteinGroupToPeptideMap2, expRepKey, expRepKey);
			}
		}
		final String header = "acc+replicate+experiment" + "\t" + "sequence+charge+replicate+experiment" + "\t"
				+ "proteinGroup-replicate-experiment --> peptide-replicate-experiment";
		writeMapToFile(header, map, writer);

	}

	/**
	 * ACC1_REP1_EXP1 --> PEPTIDEA_3_REP1_EXP1<br>
	 * ACC1_REP1_EXP1 --> PEPTIDEB_2_REP1_EXP1<br>
	 * ACC1_REP2_EXP1 --> PEPTIDEA_3_REP2_EXP1<br>
	 * ACC1_REP1_EXP2 --> PEPTIDEA_3_REP1_EXP2<br>
	 * ACC1_REP1_EXP2 --> PEPTIDEB_2_REP1_EXP2<br>
	 * ACC1_REP2_EXP2 --> PEPTIDEA_3_REP2_EXP2<br>
	 *
	 * @throws IOException
	 */
	private void writePeptideToProteinMap() throws IOException {
		final File file = new File(getWorkingPath().getAbsolutePath() + File.separator
				+ FileMappingResults.PEPTIDE_REPLICATE_EXPERIMENT_TO_PROTEIN_REPLICATE_EXPERIMENT_3);
		if (!overrideFilesIfExists && file.exists()) {
			return;
		}
		final FileWriter writer = new FileWriter(file);
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();

		for (final QuantExperiment exp : quantExperiments) {
			String expName = "";
			if (quantExperiments.size() > 1) {
				expName = exp.getName();
			}
			for (final QuantReplicate rep : exp.getReplicates()) {
				String repName = "";
				if (exp.getReplicates().size() > 1) {
					repName = rep.getName();
				}
				String expRepKey = "";
				if (!"".equals(repName)) {
					expRepKey = "_" + repName;
				}
				if (!"".equals(expName)) {
					expRepKey += "_" + expName;
				}
				final QuantParser parser = rep.getParser();
				// in this case, get all proteins and construct protein groups.
				// Then, asign in the map to peptides
				final Map<String, QuantifiedProteinInterface> proteinMap = parser.getProteinMap();
				final Map<String, Set<String>> map2 = new THashMap<String, Set<String>>();
				for (final String proteinKey : proteinMap.keySet()) {
					final QuantifiedProteinInterface quantifiedProtein = proteinMap.get(proteinKey);
					if (quantifiedProtein.isDiscarded()) {
						continue;
					}
					final Set<QuantifiedPSMInterface> quantifiedPSMs = quantifiedProtein.getQuantifiedPSMs();
					for (final QuantifiedPSMInterface quantifiedPSM : quantifiedPSMs) {
						if (quantifiedPSM.isDiscarded()) {
							continue;
						}
						final String peptideKey = KeyUtils.getSequenceKey(quantifiedPSM, true);
						if (map2.containsKey(proteinKey)) {
							map2.get(proteinKey).add(peptideKey);
						} else {
							final Set<String> set = new THashSet<String>();
							set.add(peptideKey);
							map2.put(proteinKey, set);
						}
					}
					// }
				}
				mergeMaps(map, map2, expRepKey, expRepKey);
			}
		}
		final String header = "acc+replicate+experiment" + "\t" + "sequence+charge+replicate+experiment" + "\t"
				+ "protein-replicate-experiment --> peptide-replicate-experiment";
		writeMapToFile(header, map, writer);

	}

	private Map<String, Set<String>> getProteinGroupToPeptideMap(List<ProteinGroup> proteinGroups) {
		final Map<String, Set<String>> ret = new THashMap<String, Set<String>>();
		for (final ProteinGroup proteinGroup : proteinGroups) {
			final String groupKey = KeyUtils.getGroupKey(proteinGroup);
			final Set<String> set = new THashSet<String>();
			final List<GroupablePeptide> psMs = proteinGroup.getPSMs();
			for (final GroupablePeptide groupablePSM : psMs) {
				if (groupablePSM instanceof QuantifiedPSMInterface) {
					final QuantifiedPSMInterface psm = (QuantifiedPSMInterface) groupablePSM;
					final String peptideKey = KeyUtils.getSequenceKey(psm, true);
					set.add(peptideKey);
				}
			}
			ret.put(groupKey, set);
		}
		return ret;
	}

	private List<ProteinGroup> getProteinGroups(List<GroupableProtein> groupableProteins) {
		final PAnalyzer pa = new PAnalyzer(false);
		final List<ProteinGroup> proteinGroups = pa.run(groupableProteins);
		return proteinGroups;
	}

	/**
	 * PEPTIDEA_3_REP1_EXP1 --> SCAN1_RAWFILE1<br>
	 * PEPTIDEA_3_REP1_EXP1 --> SCAN2_RAWFILE1<br>
	 * PEPTIDEA_3_REP2_EXP1 --> SCAN1_RAWFILE2<br>
	 * PEPTIDEB_2_REP1_EXP1 --> SCAN3_RAWFILE1<br>
	 *
	 * @throws IOException
	 */
	private void writeSpectrumToPeptideExperimentReplicateMap() throws IOException {
		final FileWriter writer = new FileWriter(getWorkingPath().getAbsolutePath() + File.separator
				+ FileMappingResults.SPECTRUM_TO_PEPTIDE_REPLICATE_EXPERIMENT_2);
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();

		for (final QuantExperiment exp : quantExperiments) {
			String expName = "";
			if (quantExperiments.size() > 1) {
				expName = exp.getName();
			}
			for (final QuantReplicate rep : exp.getReplicates()) {
				String repName = "";
				if (exp.getReplicates().size() > 1) {
					repName = rep.getName();
				}
				String expRepKey = "";
				if (!"".equals(repName)) {
					expRepKey = "_" + repName;
				}
				if (!"".equals(expName)) {
					expRepKey += "_" + expName;
				}
				final QuantParser parser = rep.getParser();
				final Map<String, Set<String>> tmpMap = parser.getPeptideToSpectraMap();
				final Map<String, Set<String>> peptideToSpectraMap2 = addRepNameToMap(tmpMap, expRepKey);
				mergeMaps(map, peptideToSpectraMap2, expRepKey, "");
			}
		}
		final String header = "sequence+charge+replicate+experiment" + "\t" + "scan+raw_file" + "\t"
				+ "spectrum --> peptide-replicate-experiment";
		writeMapToFile(header, map, writer);

	}

	/**
	 * SCAN1_RAWFILE1 --> IONY1_SCAN1_RAWFILE1<br>
	 * SCAN1_RAWFILE1 --> IONB2_SCAN1_RAWFILE1<br>
	 * SCAN2_RAWFILE1 --> IONY1_SCAN2_RAWFILE1<br>
	 *
	 * @throws IOException
	 */
	private void writeIonToSpectrumMap() throws IOException {
		if (quantType != QuantificationType.ISOTOPOLOGUES) {
			return;
		}
		final File file = new File(
				getWorkingPath().getAbsolutePath() + File.separator + FileMappingResults.ION_TO_SPECTRUM_1);
		if (!overrideFilesIfExists && file.exists()) {
			return;
		}
		final FileWriter writer = new FileWriter(file);
		final Map<String, Set<String>> map = new THashMap<String, Set<String>>();

		for (final QuantExperiment exp : quantExperiments) {

			String expName = "";
			if (quantExperiments.size() > 1) {
				expName = exp.getName();
			}
			for (final QuantReplicate rep : exp.getReplicates()) {

				String repName = "";
				if (exp.getReplicates().size() > 1) {
					repName = rep.getName();
				}
				String expRepKey = "";
				if (!"".equals(repName)) {
					expRepKey = "_" + repName;
				}
				if (!"".equals(expName)) {
					expRepKey += "_" + expName;
				}
				final IsobaricQuantParser parser = (IsobaricQuantParser) rep.getParser();
				final Map<String, Set<String>> tmpMap = parser.getSpectrumToIonsMap();

				// add repName to the elements of the map
				final Map<String, Set<String>> spectrumToIonsMap2 = addRepNameToMap(tmpMap, expRepKey);

				mergeMaps(map, spectrumToIonsMap2, "", "");
			}
		}

		final String header = "scan+raw_file" + "\t" + "ion_type+scan+raw_file" + "\t" + "ion --> spectrum";
		writeMapToFile(header, map, writer);
	}

	private Map<String, Set<String>> addRepNameToMap(Map<String, Set<String>> map, String repName) {
		final Map<String, Set<String>> ret = new THashMap<String, Set<String>>();
		for (final String key : map.keySet()) {
			final Set<String> keys = map.get(key);
			final Set<String> correctedKeys = new THashSet<String>();
			for (final String key2 : keys) {
				String newKey = key2;
				if (!key2.endsWith(repName))
					newKey += repName;
				correctedKeys.add(newKey);
			}
			String newKey2 = key;
			if (!newKey2.endsWith(repName))
				newKey2 += repName;
			ret.put(newKey2, correctedKeys);
		}
		return ret;
	}

	private void writeMapToFile(String header, Map<String, Set<String>> map, FileWriter writer) throws IOException {
		try {
			writer.write("#" + header + NL);
			final Set<String> keySet = map.keySet();
			final List<String> list = new ArrayList<String>();
			list.addAll(keySet);
			Collections.sort(list);
			for (final String key : list) {
				final Set<String> set = map.get(key);
				for (final String value : set) {
					writer.write(key + "\t" + value + NL);
				}
			}
		} finally {
			if (writer != null) {
				writer.close();
			}
		}
	}

	private void mergeMaps(Map<String, Set<String>> receiver, Map<String, Set<String>> donor, String suffixForKey,
			String suffixForValue) {
		for (final String originalkey : donor.keySet()) {
			String key = originalkey;
			if (!key.endsWith(suffixForKey))
				key += suffixForKey;
			final Set<String> donorValues = donor.get(originalkey);
			if (receiver.containsKey(key)) {
				for (final String donorValue : donorValues) {
					receiver.get(key).add(donorValue + suffixForValue);
				}
			} else {
				final Set<String> set = new THashSet<String>();
				for (final String donorValue : donorValues) {
					String key2 = donorValue;
					if (!key2.endsWith(suffixForKey))
						key2 += suffixForKey;
					set.add(key2);
				}
				receiver.put(key, set);
			}
		}
	}

	@Override
	public void propertyChange(PropertyChangeEvent prop) {
		if (prop.getPropertyName().equals(SanXotInterfaze.CALIBRATING)) {
			log.info("Waiting for calibrating data...");
		} else if (prop.getPropertyName().equals(SanXotInterfaze.STARTING_COMMAND)) {
			final String commandString = (String) prop.getNewValue();
			log.info("Waiting for running command: " + commandString);
		} else if (prop.getPropertyName().equals(SanXotInterfaze.END_COMMAND)) {
			final Long command = (Long) prop.getNewValue();
			log.info("Command finished with ext code: " + command);
			if (command == -999)
				System.exit(-999);
		} else if (prop.getPropertyName().equals(SanXotInterfaze.END_ANALYSIS)) {
			result = (SanXotAnalysisResult) prop.getNewValue();
			log.info(result);
			final IntegrationResultWrapper lastIntegrationResults = result.getLastIntegrationResults();
			final List<OutStatsLine> resultData = lastIntegrationResults.getResultData();
			for (final OutStatsLine outStatsLine : resultData) {
				System.out.println(outStatsLine.getIdsup() + "\t" + outStatsLine.getFDR());
			}
		}
	}

	public void setCalibration(boolean b) {
		quantParameters.setPerformCalibration(b);
	}

	public void setOutlierRemovalFDR(Double fdr) {
		quantParameters.setOutlierRemovalFDR(fdr);
	}

	public void setSanxotScriptsLocationFolder(File folder) {
		quantParameters.setSanxotScriptsFolder(folder);
	}

	public void setFastaFile(File fastaFile) {
		log.info("Constructing index from fasta file: " + fastaFile + " using default parameters");
		setFastaFile(DBIndexInterface.getDefaultDBIndexParams(fastaFile));
	}

	public void setFastaFile(DBIndexSearchParams dbIndexSearchParams) {
		log.info("Constructing index using provided parameters");
		dbIndex = DBIndexInterface.getByParam(dbIndexSearchParams);
	}

	/**
	 * @return the result
	 */
	public SanXotAnalysisResult getResult() {
		return result;
	}

	/**
	 * @return the minAlignmentScore
	 */
	public int getMinAlignmentScore() {
		return minAlignmentScore;
	}

	/**
	 * @param minAlignmentScore
	 *            the minAlignmentScore to set
	 */
	public void setMinAlignmentScore(int minAlignmentScore) {
		this.minAlignmentScore = minAlignmentScore;
	}

	/**
	 * @return the minPercentajeOfsmilirarity
	 */
	public double getMinPercentajeOfsmilirarity() {
		return minPercentajeOfsmilirarity;
	}

	/**
	 * @param minPercentajeOfsmilirarity
	 *            the minPercentajeOfsmilirarity to set
	 */
	public void setMinPercentajeOfsmilirarity(double minPercentajeOfsmilirarity) {
		this.minPercentajeOfsmilirarity = minPercentajeOfsmilirarity;
	}

	/**
	 * @return the minConsecutiveLength
	 */
	public int getMinConsecutiveLength() {
		return minConsecutiveLength;
	}

	/**
	 * @param minConsecutiveLength
	 *            the minConsecutiveLength to set
	 */
	public void setMinConsecutiveLength(int minConsecutiveLength) {
		this.minConsecutiveLength = minConsecutiveLength;
	}

	/**
	 * @return the overrideFilesIfExists
	 */
	public boolean isOverrideFilesIfExists() {
		return overrideFilesIfExists;
	}

	/**
	 * Determine if override or not the relationaship and data files if they
	 * already exist.
	 *
	 * @param overrideFilesIfExists
	 *            the overrideFilesIfExists to set
	 */
	public void setOverrideFilesIfExists(boolean overrideFilesIfExists) {
		this.overrideFilesIfExists = overrideFilesIfExists;
	}

	/**
	 * force the analysis to retrieve the protein level values by grouping
	 * proteins according to the proteinAccClusters provided in the parameter
	 *
	 * @param proteinAccClusters
	 * @throws IllegalArgumentException
	 *             if the {@link ANALYSIS_LEVEL_OUTCOME} is not set to
	 *             ANALYSIS_LEVEL_OUTCOME.FORCED_CLUSTERS
	 */
	public void forceProteinClusters(Set<Set<String>> proteinAccClusters) throws IllegalArgumentException {
		if (analysisOutCome != ANALYSIS_LEVEL_OUTCOME.FORCED_CLUSTERS) {
			throw new IllegalArgumentException(
					"Analysis outcome has to be set to " + ANALYSIS_LEVEL_OUTCOME.FORCED_CLUSTERS);
		}

		this.proteinAccClusters = proteinAccClusters;
	}

	public void setQuantParameters(QuantParameters quantParameters2) {
		quantParameters = quantParameters2;
	}

	public void setTimeout(long timeout) {
		quantParameters.setTimeout(timeout);

	}

}
