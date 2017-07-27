package edu.scripps.yates.census.read;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.analysis.util.KeyUtils;
import edu.scripps.yates.census.read.model.CensusRatio;
import edu.scripps.yates.census.read.model.QuantAmount;
import edu.scripps.yates.census.read.model.QuantifiedPSM;
import edu.scripps.yates.census.read.model.QuantifiedPeptide;
import edu.scripps.yates.census.read.model.QuantifiedProtein;
import edu.scripps.yates.census.read.model.QuantifiedProteinFromDBIndexEntry;
import edu.scripps.yates.census.read.model.RatioDescriptor;
import edu.scripps.yates.census.read.model.RatioScore;
import edu.scripps.yates.census.read.model.StaticQuantMaps;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.MyHashMap;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.dbindex.IndexedProtein;
import edu.scripps.yates.dbindex.util.PeptideNotFoundInDBIndexException;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.model.enums.AmountType;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import edu.scripps.yates.utilities.strings.StringUtils;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class CensusOutParser extends AbstractQuantParser {
	private static final Logger log = Logger.getLogger(CensusOutParser.class);

	private static final String H = "H";
	private static final String P = "P";
	private static final String S = "S";
	private static final String SINGLETON_S = "&S";
	private static final String PLINE = "PLINE";
	private static final String SLINE = "SLINE";
	private static final String SINGLETON_SLINE = "&SLINE";
	private static final String LOCUS = "LOCUS";
	private static final String UNIQUE = "UNIQUE";
	private static final String SEQUENCE = "SEQUENCE";
	public static final String FILENAME = "FILE_NAME";
	// synonym from previous versions
	public static final String FILE_name = "Filename";
	private static final String DESCRIPTION = "DESCRIPTION";

	// PLine columns
	private static final String AVERAGE_RATIO = "AVERAGE_RATIO";
	private static final String STANDARD_DEVIATION = "STANDARD_DEVIATION";
	private static final String NORM_STDEV = "NORM_STDEV";
	private static final String COMPOSITE_RATIO = "COMPOSITE_RATIO";
	private static final String NORM_COMPOSITE_RATIO = "NORM_COMPOSITE_RATIO";
	private static final String COMPOSITE_RATIO_STANDARD_DEVIATION = "COMPOSITE_RATIO_STANDARD_DEVIATION";
	private static final String NORM_COMPOSITE_RATIO_STDEV = "NORM_COMPOSITE_RATIO_STDEV";
	private static final String COMPOSITE_RATIO_STDEV = "COMPOSITE_RATIO_STDEV";
	private static final String MEDIAN_NORM_RATIO = "MEDIAN_NORM_RATIO";
	private static final String MEDIAN_AREA_RATIO = "MEDIAN_AREA_RATIO";
	// SLine columns
	public static final String RATIO = "RATIO";
	public static final String NORM_RATIO = "NORM_RATIO";
	private static final String PVALUE = "PVALUE";
	private static final String PROBABILITY_SCORE = "PROBABILITY_SCORE";
	private static final String PROFILE_SCORE = "PROFILE_SCORE";
	public static final String AREA_RATIO = "AREA_RATIO";
	private static final String SAM_INT = "SAM_INT";
	private static final String PEAK_AREA = "PEAK_AREA";
	private static final String PEAK_AREA_L = "PEAK_AREA_L";
	private static final String PEAK_AREA_M = "PEAK_AREA_M";
	private static final String PEAK_AREA_H = "PEAK_AREA_H";
	private static final String REF_INT = "REF_INT";
	private static final String PEAK_INT = "PEAK_INT";
	private static final String REGRESSION_FACTOR = "REGRESSION_FACTOR";

	// &SLine columns
	private static final String SINGLETON_SCORE = "SINGLETON_SCORE";

	public static final String CS = "CS";
	// synonym of CS, from older versions of census out files:
	public static final String CState = "CState";
	public static final String SCAN = "SCAN";
	// synonym of SCAN, from older versions of census out files
	public static final String SCAN_NUM = "ScanNum";

	private static final String XCORR = "XCorr";
	private static final String DELTACN = "deltaCN";

	private static final String DET_FACTOR = "DET_FACTOR";
	private boolean onlyOneSpectrumPerChromatographicPeakAndPerSaltStep = false;
	private boolean skipSingletons = false; // by default
	private boolean skipNonResolvedPeaks = true; // by default

	public CensusOutParser() {
		super();
	}

	public CensusOutParser(List<RemoteSSHFileReference> remoteSSHServers,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditions, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) {
		super(remoteSSHServers, labelsByConditions, labelNumerator, labelDenominator);
	}

	public CensusOutParser(Map<QuantCondition, QuantificationLabel> labelsByConditions,
			Collection<RemoteSSHFileReference> remoteSSHServers, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) {
		super(labelsByConditions, remoteSSHServers, labelNumerator, labelDenominator);
	}

	public CensusOutParser(RemoteSSHFileReference remoteSSHServer,
			Map<QuantCondition, QuantificationLabel> labelsByConditions, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) throws FileNotFoundException {
		super(remoteSSHServer, labelsByConditions, labelNumerator, labelDenominator);
	}

	public CensusOutParser(File xmlFile, Map<QuantCondition, QuantificationLabel> labelsByConditions,
			QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) throws FileNotFoundException {
		super(xmlFile, labelsByConditions, labelNumerator, labelDenominator);
	}

	public CensusOutParser(File[] xmlFiles, Map<QuantCondition, QuantificationLabel> labelsByConditions,
			QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) throws FileNotFoundException {
		super(xmlFiles, labelsByConditions, labelNumerator, labelDenominator);
	}

	public CensusOutParser(File[] xmlFiles, Map<QuantCondition, QuantificationLabel>[] labelsByConditions,
			QuantificationLabel[] labelNumerator, QuantificationLabel[] labelDenominator) throws FileNotFoundException {
		super(xmlFiles, labelsByConditions, labelNumerator, labelDenominator);
	}

	public CensusOutParser(Collection<File> xmlFiles, Map<QuantCondition, QuantificationLabel> labelsByConditions,
			QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) throws FileNotFoundException {
		super(xmlFiles, labelsByConditions, labelNumerator, labelDenominator);
	}

	public CensusOutParser(RemoteSSHFileReference remoteServer, QuantificationLabel label1, QuantCondition cond1,
			QuantificationLabel label2, QuantCondition cond2) {
		super(remoteServer, label1, cond1, label2, cond2);
	}

	public CensusOutParser(File inputFile, QuantificationLabel label1, QuantCondition cond1, QuantificationLabel label2,
			QuantCondition cond2) throws FileNotFoundException {
		super(inputFile, label1, cond1, label2, cond2);
	}

	public CensusOutParser(File inputFile, Map<QuantCondition, QuantificationLabel> map, QuantificationLabel light,
			QuantificationLabel medium, QuantificationLabel heavy) throws FileNotFoundException {
		super(inputFile, map, light, medium, heavy);
	}

	/**
	 *
	 * @param writeFiles
	 *            whether to write output files necessary to run SanXot program
	 */
	@Override
	protected void process() {
		processed = false;
		log.info("Processing quant file...");

		try {
			int numDecoy = 0;
			boolean someValidFile = false;
			for (RemoteSSHFileReference remoteFileRetriever : remoteFileRetrievers) {
				final Map<QuantCondition, QuantificationLabel> labelsByConditions = labelsByConditionsByFile
						.get(remoteFileRetriever);
				final Map<QuantificationLabel, QuantCondition> conditionsByLabels = new THashMap<QuantificationLabel, QuantCondition>();
				for (QuantCondition cond : labelsByConditions.keySet()) {
					conditionsByLabels.put(labelsByConditions.get(cond), cond);
				}
				List<RatioDescriptor> ratioDescriptors = ratioDescriptorsByFile.get(remoteFileRetriever);
				// QuantificationLabel labelNumerator =
				// numeratorLabelByFile.get(remoteFileRetriever);
				// QuantificationLabel labelDenominator =
				// denominatorLabelByFile.get(remoteFileRetriever);

				String experimentKey = FilenameUtils.getBaseName(remoteFileRetriever.getOutputFile().getAbsolutePath());
				String fileName = FilenameUtils.getName(remoteFileRetriever.getOutputFile().getAbsolutePath());
				log.info(experimentKey);
				// get all the Quantified PSMs first
				// Set<QuantifiedPSMInterface> psms = new
				// HashSet<QuantifiedPSMInterface>();

				log.info("Reading " + remoteFileRetriever.getRemoteFileName() + " from "
						+ remoteFileRetriever.getRemotePath());
				someValidFile = true;
				BufferedReader br = null;
				try {
					br = new BufferedReader(
							new InputStreamReader(new BufferedInputStream(remoteFileRetriever.getRemoteInputStream())));

					String line;
					List<String> pLineHeaderList = new ArrayList<String>();
					List<String> sLineHeaderList = new ArrayList<String>();
					List<String> singletonSLineHeaderList = new ArrayList<String>();
					Set<QuantifiedProteinInterface> proteinGroup = new THashSet<QuantifiedProteinInterface>();
					boolean itWasPeptides = false;
					while ((line = br.readLine()) != null) {
						if (line.startsWith(H)) {
							if (line.contains("\t")) {
								final String[] split = line.split("\t");
								// if second element is PLINE
								if (split[1].equals(PLINE) && split[2].equals(LOCUS)) {
									for (int i = 1; i < split.length; i++) {
										pLineHeaderList.add(split[i]);
									}
								} else // if second element is SLINE
								if (split[1].equals(SLINE) && split[2].equals(UNIQUE)) {
									for (int i = 1; i < split.length; i++) {
										sLineHeaderList.add(split[i]);
									}
								} else // if second element is &SLINE
								if (split[1].equals(SINGLETON_SLINE) && split[2].equals(UNIQUE)) {
									for (int i = 1; i < split.length; i++) {
										if (split[i].equals(PVALUE)) {
											continue; // SINGLETONG SLINE DOESNT
														// HAVE PVALUES!
										}
										singletonSLineHeaderList.add(split[i]);
									}
								}
							}
						} else if (line.startsWith(P)) {
							try {
								if (itWasPeptides) {
									proteinGroup.clear();
								}
								itWasPeptides = false;
								QuantifiedProteinInterface quantifiedProtein = processProteinLine(line, pLineHeaderList,
										conditionsByLabels, ratioDescriptors, numDecoy);
								proteinGroup.add(quantifiedProtein);
							} catch (DiscardProteinException e) {
								continue;
							} catch (IllegalArgumentException e) {
								log.error(e);
							}
						} else if (line.startsWith(S)) {
							itWasPeptides = true;
							// if there is not protein is because it was
							// discarded by the decoy pattern, so ignore any psm
							if (proteinGroup.isEmpty()) {
								continue;
							}

							processPSMLine(line, sLineHeaderList, proteinGroup, conditionsByLabels, labelsByConditions,
									ratioDescriptors, experimentKey, remoteFileRetriever, false);

						} else if (line.startsWith(SINGLETON_S)) {
							if (skipSingletons) {
								continue;
							}
							itWasPeptides = true;
							// if there is not protein is because it was
							// discarded by the decoy pattern, so ignore any psm
							if (proteinGroup.isEmpty()) {
								continue;
							}
							if (singletonSLineHeaderList.isEmpty()) {
								// it is because sometimes there is not header
								// for singletons becase is the same as no
								// singletons
								singletonSLineHeaderList.addAll(sLineHeaderList);
							}
							processPSMLine(line, singletonSLineHeaderList, proteinGroup, conditionsByLabels,
									labelsByConditions, ratioDescriptors, experimentKey, remoteFileRetriever, true);

						}

					}

					br.close();

				} catch (PeptideNotFoundInDBIndexException e) {
					if (!super.ignoreNotFoundPeptidesInDB) {
						throw e;
					}
				} catch (Exception e) {

					e.printStackTrace();
					throw e;
				} finally {
					if (br != null) {
						br.close();
					}
				}

				log.info("(" + experimentKey + ") " + localPsmMap.size() + " PSMs from this parser. "
						+ StaticQuantMaps.psmMap.size() + " PSMs in the system");
				log.info("(" + experimentKey + ") " + localProteinMap.size() + " Proteins from this parser. "
						+ StaticQuantMaps.proteinMap.size() + " Proteins in the system");
				log.info("(" + experimentKey + ") " + localPeptideMap.size() + " Peptides from this parser. "
						+ StaticQuantMaps.peptideMap.size() + " Peptides in the system");
				if (decoyPattern != null) {
					log.info(numDecoy + " decoy Proteins were discarded  in " + experimentKey);
				}

				if (onlyOneSpectrumPerChromatographicPeakAndPerSaltStep) {
					log.info(
							"Reviewing data in order to remove redundant measurements of the same chromatographic peak in the same salt step");
					// create a map in which the key is the peptideSequence +
					// rawFile (removing the H) + chargeState
					// and the values are Sets of psms
					Map<String, Set<QuantifiedPSMInterface>> map = new THashMap<String, Set<QuantifiedPSMInterface>>();
					for (QuantifiedPSMInterface psm : localPsmMap.values()) {
						String key = getSpectrumPerChromatographicPeakAndPerSaltStepKey(psm);
						if (map.containsKey(key)) {
							map.get(key).add(psm);
						} else {
							Set<QuantifiedPSMInterface> set = new THashSet<QuantifiedPSMInterface>();
							set.add(psm);
							map.put(key, set);
						}
					}
					// once the map is populated,
					// look for each key, if we have more than one psm
					// in that case, select the best one
					int numRemoved = 0;
					for (String key : map.keySet()) {
						final Set<QuantifiedPSMInterface> psmSet = map.get(key);
						if (psmSet.size() > 1) {
							QuantifiedPSMInterface bestPSM = getBestPSM(psmSet);
							Set<QuantifiedPSMInterface> toIgnore = new THashSet<QuantifiedPSMInterface>();
							for (QuantifiedPSMInterface psm : psmSet) {
								if (!psm.equals(bestPSM)) {
									toIgnore.add(psm);
								}
							}
							// remove the psms in toIgnore Set
							if (!toIgnore.isEmpty()) {
								for (QuantifiedPSMInterface psmToIgnore : toIgnore) {
									numRemoved++;
									localPsmMap.remove(psmToIgnore.getKey());
									StaticQuantMaps.psmMap.remove(psmToIgnore);
									// remove it from its peptide
									QuantifiedPeptideInterface quantifiedPeptide = psmToIgnore.getQuantifiedPeptide();
									if (quantifiedPeptide != null) {
										quantifiedPeptide.getQuantifiedPSMs().remove(psmToIgnore);
									}
									// remove it from its proteins
									Set<QuantifiedProteinInterface> quantifiedProteins = psmToIgnore
											.getQuantifiedProteins();
									for (QuantifiedProteinInterface protein : quantifiedProteins) {
										protein.getQuantifiedPSMs().remove(psmToIgnore);
									}
								}
							}
						}
					}
					log.info(numRemoved + " PSMs where redundant and removed.");

				}
			}
			if (!someValidFile)
				throw new IllegalArgumentException("some error occurred while reading the files");

			processed = true;
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			// if (processed) {
			// // to create the peptides at the end
			// peptideMap.clear();
			// peptideMap.putAll(
			// QuantifiedPeptide.getQuantifiedPeptides(getPSMMap().values(),
			// distinguishModifiedPeptides));
			// }
		}
	}

	/**
	 * Get the best PSM from the set of PSMs, looking into (in this order) the
	 * REGRESSION_FACTOR or the XCorr
	 *
	 * @param psmSet
	 * @return
	 */
	private QuantifiedPSMInterface getBestPSM(Set<QuantifiedPSMInterface> psmSet) {
		List<QuantifiedPSMInterface> list = new ArrayList<QuantifiedPSMInterface>();
		list.addAll(psmSet);
		Collections.sort(list, new Comparator<QuantifiedPSMInterface>() {

			@Override
			public int compare(QuantifiedPSMInterface o1, QuantifiedPSMInterface o2) {
				// regression_factor
				Double regressionFactor1 = getRatioScore(REGRESSION_FACTOR, o1);
				Double regressionFactor2 = getRatioScore(REGRESSION_FACTOR, o2);
				if (regressionFactor1 != null && regressionFactor2 != null) {
					final int compare = Double.compare(regressionFactor2, regressionFactor1);
					if (compare != 0) {
						return compare;
					}
				}
				// xcorr
				Float xcorr1 = o1.getXcorr();
				Float xcorr2 = o2.getXcorr();
				if (xcorr1 != null && xcorr2 != null) {
					final int compare = Float.compare(xcorr2, xcorr1);
					if (compare != 0) {
						return compare;
					}
				}
				return 0;
			}

			private Double getRatioScore(String scoreName, QuantifiedPSMInterface o1) {
				if (o1.getRatios() != null) {
					for (QuantRatio quantRatio : o1.getRatios()) {
						if (quantRatio.getAssociatedConfidenceScore() != null) {
							if (quantRatio.getAssociatedConfidenceScore().getScoreName().equals(scoreName)) {
								try {
									return Double.valueOf(quantRatio.getAssociatedConfidenceScore().getValue());
								} catch (NumberFormatException e) {

								}
							}
						}
					}
				}
				if (o1.getAmounts() != null) {
					for (Amount amount : o1.getAmounts()) {
						if (amount.getAmountType().name().equals(scoreName)) {
							return amount.getValue();
						}

					}
				}
				return null;
			}
		});
		// return the first element
		return list.get(0);
	}

	/**
	 * peptideSequence + rawFile (removing the H) + chargeState
	 *
	 * @param psm
	 * @return
	 */
	private String getSpectrumPerChromatographicPeakAndPerSaltStepKey(QuantifiedPSMInterface psm) {
		StringBuilder sb = new StringBuilder();
		String rawFileName = psm.getRawFileNames().iterator().next();

		if (rawFileName.startsWith("H") || rawFileName.startsWith("M")) {
			String substring = rawFileName.substring(1);
			// only take it if there is another one like that in the static map
			// this will avoid to consider an H or an M in the file name that
			// were part of the name and not meaning anything
			if (StaticQuantMaps.rawFileNames.contains(substring)) {
				rawFileName = substring;
			}
		}

		sb.append(psm.getSequence()).append("_").append(rawFileName).append("_").append(psm.getCharge());
		return sb.toString();
	}

	private void processPSMLine(String line, List<String> sLineHeaderList,
			Set<QuantifiedProteinInterface> quantifiedProteins,
			Map<QuantificationLabel, QuantCondition> conditionsByLabels,
			Map<QuantCondition, QuantificationLabel> labelsByConditions, List<RatioDescriptor> ratioDescriptors,
			String experimentKey, RemoteSSHFileReference remoteFileRetriever, boolean singleton) throws IOException {

		// new psm
		try {
			MyHashMap<String, String> mapValues = getMapFromSLine(sLineHeaderList, line);

			String sequence = mapValues.get(SEQUENCE);

			// dont look into the QuantifiedPSM.map because each
			// line is always a new PSM
			String inputFileName = FilenameUtils.getName(remoteFileRetriever.getOutputFile().getAbsolutePath());
			String rawFileName = null;
			if (mapValues.containsKey(FILENAME)) {
				rawFileName = mapValues.get(FILENAME);
			} else {
				rawFileName = inputFileName;
			}
			StaticQuantMaps.rawFileNames.add(rawFileName);
			// scan number
			int scanNumber = 0;
			if (mapValues.containsKey(SCAN)) {
				scanNumber = Double.valueOf(mapValues.get(SCAN)).intValue();
			}
			QuantifiedPSMInterface quantifiedPSM = new QuantifiedPSM(sequence, labelsByConditions, peptideToSpectraMap,
					scanNumber, Double.valueOf(mapValues.get(CS)).intValue(), rawFileName, singleton);

			// xcorr
			Float xcorr = null;
			if (mapValues.containsKey(XCORR)) {
				try {
					xcorr = Float.valueOf(mapValues.get(XCORR));
					((QuantifiedPSM) quantifiedPSM).setXcorr(xcorr);
				} catch (NumberFormatException e) {

				}
			}
			// deltacn
			Float deltaCn = null;
			if (mapValues.containsKey(DELTACN)) {
				try {
					deltaCn = Float.valueOf(mapValues.get(DELTACN));
					((QuantifiedPSM) quantifiedPSM).setDeltaCN(deltaCn);
				} catch (NumberFormatException e) {

				}
			}

			quantifiedPSM.getFileNames().add(inputFileName);
			final String psmKey = KeyUtils.getSpectrumKey(quantifiedPSM, true);
			// in case of TMT, the psm may have been created before
			if (StaticQuantMaps.psmMap.containsKey(psmKey)) {
				quantifiedPSM = StaticQuantMaps.psmMap.getItem(psmKey);
			}
			StaticQuantMaps.psmMap.addItem(quantifiedPSM);

			// psms.add(quantifiedPSM);
			// add to map
			if (!localPsmMap.containsKey(quantifiedPSM.getKey())) {
				localPsmMap.put(quantifiedPSM.getKey(), quantifiedPSM);
			}

			if (ratioDescriptors.size() == 1) {
				QuantificationLabel labelNumerator = ratioDescriptors.get(0).getLabel1();
				QuantificationLabel labelDenominator = ratioDescriptors.get(0).getLabel2();
				String ratioSuffix = ratioDescriptors.get(0).getRatioSuffix();
				// PSM regular ratio
				// add PSM ratios from census out
				if (mapValues.containsKey(RATIO)) {
					ratioSuffix = "";
				} else if (mapValues.containsKey(RATIO + ratioSuffix)) {

				}
				if (mapValues.containsKey(RATIO + ratioSuffix)) {
					try {

						final double ratioValue = QuantUtils.parseCensusRatioValue(mapValues.get(RATIO + ratioSuffix));
						CensusRatio ratio = new CensusRatio(ratioValue, false, conditionsByLabels, labelNumerator,
								labelDenominator, AggregationLevel.PSM, RATIO);
						RatioScore ratioScore = null;
						// profile score
						String scoreValue = null;
						// check in any case the regressionFactor. If that is
						// -1,
						// then
						// convert the ratio from 0 to inf
						String regressionFactor = null;
						if (mapValues.containsKey(PROFILE_SCORE) && // PROFILE
																	// SCORE IS
																	// NOT
																	// SPECIFIC
																	// OF A PAIR
																	// OF LABELS
						// because profile score will be
						// assigned to area ratio in case of
						// exist
								!mapValues.containsKey(AREA_RATIO + ratioSuffix)) {
							scoreValue = mapValues.get(PROFILE_SCORE);
							if (!"NA".equals(scoreValue)) {
								ratioScore = new RatioScore(scoreValue, PROFILE_SCORE,
										"PSM-level quantification confidence metric",
										"fitting score comparing peak and gaussian distribution");
								ratio.setRatioScore(ratioScore);
							}
							// score regression p-value N15
						} else if (mapValues.containsKey(PVALUE + ratioSuffix)) {
							scoreValue = mapValues.get(PVALUE + ratioSuffix);
							if (!"NA".equals(scoreValue)) {
								ratioScore = new RatioScore(scoreValue, PVALUE, "PSM-level p-value",
										"probability score based on LR");
								ratio.setRatioScore(ratioScore);
							}
							// score regression p-value SILAC
						} else if (mapValues.containsKey(PROBABILITY_SCORE + ratioSuffix)) {
							scoreValue = mapValues.get(PROBABILITY_SCORE + ratioSuffix);
							if (!"NA".equals(scoreValue)) {
								ratioScore = new RatioScore(scoreValue, PROBABILITY_SCORE, "PSM-level p-value",
										"probability score based on LR");
								ratio.setRatioScore(ratioScore);
							}
						} else if (mapValues.containsKey(DET_FACTOR + ratioSuffix)) {
							scoreValue = mapValues.get(DET_FACTOR + ratioSuffix);
							if (!"NA".equals(scoreValue)) {
								ratioScore = new RatioScore(scoreValue, DET_FACTOR,
										"PSM-level quantification confidence metric", "PSM-level Determining factor");
								ratio.setRatioScore(ratioScore);
							}
						}
						// get the regression_factor anyway, but add it only in
						// case
						// there is not any other score
						if (mapValues.containsKey(REGRESSION_FACTOR + ratioSuffix)) {
							regressionFactor = mapValues.get(REGRESSION_FACTOR + ratioSuffix);
							if (!"NA".equals(regressionFactor)) {
								// just in case there was not any other
								// ratioScore:
								if (ratioScore == null) {
									ratioScore = new RatioScore(regressionFactor, REGRESSION_FACTOR,
											"PSM-level quantification confidence metric",
											"Regression factor or linear regression");
									ratio.setRatioScore(ratioScore);
								}
							}
						}

						try {
							// if ratio is 0 and regression factor is -1
							if (Double.compare(ratio.getValue(), 0.0) == 0) {
								if (regressionFactor != null && ("NA".equals(regressionFactor)
										|| Double.valueOf(-1.0).equals(Double.valueOf(regressionFactor)))) {
									// check area_ratio value. If the are_ratio
									// is < 1, leave it as 0. If the area_ratio
									// is > 1,
									// convert it to +INF.
									// note that all numbers are not log
									// numbers.
									if (mapValues.containsKey(AREA_RATIO + ratioSuffix)) {
										double areaRatioValue = QuantUtils
												.parseCensusRatioValue(mapValues.get(AREA_RATIO + ratioSuffix));
										if (Double.isInfinite(areaRatioValue) || areaRatioValue > 1) {
											ratio = new CensusRatio(Double.POSITIVE_INFINITY, false, conditionsByLabels,
													labelNumerator, labelDenominator, AggregationLevel.PSM, RATIO);
											if (ratioScore != null) {
												ratio.setRatioScore(ratioScore);
											}
										}
									}
								}
							}
						} catch (NumberFormatException e) {
							// do nothing
						}
						// add ratio to PSM
						quantifiedPSM.addRatio(ratio);

					} catch (NumberFormatException e) {
						// skip this
					}
				}

				// PSM AREA RATIO
				// add PSM ratios from census out
				if (mapValues.containsKey(AREA_RATIO + ratioSuffix)) {
					try {
						double ratioValue = QuantUtils.parseCensusRatioValue(mapValues.get(AREA_RATIO + ratioSuffix));

						CensusRatio ratio = new CensusRatio(ratioValue, false, conditionsByLabels, labelNumerator,
								labelDenominator, AggregationLevel.PSM, AREA_RATIO);
						// profile score
						// PROFILE SCORE IS NOT SPECIFIC FOR A PAIR OF LABELS
						if (mapValues.containsKey(PROFILE_SCORE)) {
							String scoreValue = mapValues.get(PROFILE_SCORE);
							if (!"NA".equals(scoreValue)) {
								RatioScore ratioScore = new RatioScore(scoreValue, PROFILE_SCORE,
										"PSM-level quantification confidence metric",
										"fitting score comparing peak and gaussian distribution");
								ratio.setRatioScore(ratioScore);
							}
						} else if (mapValues.containsKey(SINGLETON_SCORE)) {
							String scoreValue = mapValues.get(SINGLETON_SCORE);
							if (!"NA".equals(scoreValue)) {
								RatioScore ratioScore = new RatioScore(scoreValue, SINGLETON_SCORE,
										"PSM-level quantification confidence metric", "Singleton score");
								ratio.setRatioScore(ratioScore);
							}
						}
						// add ratio to PSM
						quantifiedPSM.addRatio(ratio);
					} catch (NumberFormatException e) {
						// skip this
					}
				}
				// PSM AREA RATIO
				// add PSM ratios from census out
				if (mapValues.containsKey(NORM_RATIO + ratioSuffix)) {
					try {
						double ratioValue = QuantUtils.parseCensusRatioValue(mapValues.get(NORM_RATIO + ratioSuffix));

						CensusRatio ratio = new CensusRatio(ratioValue, false, conditionsByLabels, labelNumerator,
								labelDenominator, AggregationLevel.PSM, NORM_RATIO);
						// profile score
						// PROFILE SCORE IS NOT SPECIFIC FOR A PAIR OF LABELS
						if (mapValues.containsKey(PROFILE_SCORE)) {
							String scoreValue = mapValues.get(PROFILE_SCORE);
							if (!"NA".equals(scoreValue)) {
								RatioScore ratioScore = new RatioScore(scoreValue, PROFILE_SCORE,
										"PSM-level quantification confidence metric",
										"fitting score comparing peak and gaussian distribution");
								ratio.setRatioScore(ratioScore);
							}
						} else if (mapValues.containsKey(SINGLETON_SCORE)) {
							String scoreValue = mapValues.get(SINGLETON_SCORE);
							if (!"NA".equals(scoreValue)) {
								RatioScore ratioScore = new RatioScore(scoreValue, SINGLETON_SCORE,
										"PSM-level quantification confidence metric", "Singleton score");
								ratio.setRatioScore(ratioScore);
							}
						}
						// add ratio to PSM
						quantifiedPSM.addRatio(ratio);
					} catch (NumberFormatException e) {
						// skip this
					}
				}
				// TMT
				if (isTMT(labelNumerator) && isTMT(labelDenominator)) {
					// numerator
					Double numeratorIntensity = null;
					final String headerNumerator = getHeaderForTMTLabel(labelNumerator);
					if (mapValues.containsKey(headerNumerator)) {
						numeratorIntensity = Double.valueOf(mapValues.get(headerNumerator));
						QuantAmount amount = new QuantAmount(numeratorIntensity, AmountType.NORMALIZED_INTENSITY,
								conditionsByLabels.get(labelDenominator));
						quantifiedPSM.addAmount(amount);
					}
					// denominator
					Double denominatorIntensity = null;
					final String headerDenominator = getHeaderForTMTLabel(labelDenominator);
					if (mapValues.containsKey(headerDenominator)) {
						denominatorIntensity = Double.valueOf(mapValues.get(headerDenominator));
						QuantAmount amount = new QuantAmount(denominatorIntensity, AmountType.NORMALIZED_INTENSITY,
								conditionsByLabels.get(labelDenominator));
						quantifiedPSM.addAmount(amount);
					}
					// build the ratio
					Double ratioValue = numeratorIntensity / denominatorIntensity;
					CensusRatio ratio = new CensusRatio(ratioValue, false, conditionsByLabels, labelNumerator,
							labelDenominator, AggregationLevel.PSM, labelNumerator + "/" + labelDenominator);
					quantifiedPSM.addRatio(ratio);
				}
			} else {
				// having more than one ratioDescriptor, means that we have L,
				// M, H
				for (RatioDescriptor ratioDescriptor : ratioDescriptors) {
					QuantificationLabel labelNumerator = ratioDescriptor.getLabel1();
					QuantificationLabel labelDenominator = ratioDescriptor.getLabel2();
					String ratioSuffix = ratioDescriptor.getRatioSuffix();
					if (mapValues.containsKey(RATIO + ratioSuffix)) {
						final double ratioValue = QuantUtils.parseCensusRatioValue(mapValues.get(RATIO + ratioSuffix));
						CensusRatio ratio = new CensusRatio(ratioValue, false, conditionsByLabels, labelNumerator,
								labelDenominator, AggregationLevel.PSM, RATIO);
						quantifiedPSM.addRatio(ratio);
					}
					if (mapValues.containsKey(NORM_RATIO + ratioSuffix)) {
						final double ratioValue = QuantUtils
								.parseCensusRatioValue(mapValues.get(NORM_RATIO + ratioSuffix));
						CensusRatio ratio = new CensusRatio(ratioValue, false, conditionsByLabels, labelNumerator,
								labelDenominator, AggregationLevel.PSM, NORM_RATIO);
						quantifiedPSM.addRatio(ratio);
					}
					if (mapValues.containsKey(AREA_RATIO + ratioSuffix)) {
						try {
							double ratioValue = QuantUtils
									.parseCensusRatioValue(mapValues.get(AREA_RATIO + ratioSuffix));
							CensusRatio ratio = new CensusRatio(ratioValue, false, conditionsByLabels, labelNumerator,
									labelDenominator, AggregationLevel.PSM, AREA_RATIO);
							// add ratio to PSM
							quantifiedPSM.addRatio(ratio);
						} catch (NumberFormatException e) {
							// skip this
						}
					}
				}
			}
			// PSM amounts

			// SAM_INT
			// light peptide peak area from reconstructed
			// chromatogram
			if (mapValues.containsKey(SAM_INT)) {
				try {
					final double value = Double.valueOf(mapValues.get(SAM_INT));
					QuantAmount amount = new QuantAmount(value, AmountType.AREA, getLightCondition(conditionsByLabels));
					if (singleton && amount.getValue() != 0.0) {
						amount.setSingleton(true);
					}
					// add amount to PSM
					quantifiedPSM.addAmount(amount);
				} catch (NumberFormatException e) {
					// skip this
				}
			}

			// get the peak areas of all the conditions (all the labels)
			Set<String> usedSuffixes = new HashSet<String>();
			for (RatioDescriptor ratioDescriptor : ratioDescriptors) {
				Map<String, QuantCondition> conditionsByIndividualRatioSuffixes = ratioDescriptor
						.getConditionsByIndividualRatioSuffixes();
				Set<String> differentValuesOfPeakArea = new HashSet<String>();
				for (String suffix : conditionsByIndividualRatioSuffixes.keySet()) {
					if (usedSuffixes.contains(suffix)) {
						continue;
					}
					usedSuffixes.add(suffix);
					QuantCondition quantCondition = conditionsByIndividualRatioSuffixes.get(suffix);

					if (mapValues.containsKey(PEAK_AREA + suffix)) {
						try {
							String stringValue = mapValues.get(PEAK_AREA + suffix);
							if (isSkipNonResolvedPeaks() && differentValuesOfPeakArea.contains(stringValue)) {
								log.warn("PSM '" + quantifiedPSM.getPSMIdentifier()
										+ "' contains not resolved quantitation values (Repeated peak area '"
										+ stringValue + "'). Skipping it...");
								// removing it from local and static maps
								StaticQuantMaps.psmMap.remove(quantifiedPSM);
								localPsmMap.remove(quantifiedPSM.getKey());
								return;
							}
							differentValuesOfPeakArea.add(stringValue);
							final double value = Double.valueOf(stringValue);
							QuantAmount amount = new QuantAmount(value, AmountType.AREA, quantCondition);
							if (singleton && amount.getValue() != 0.0) {
								amount.setSingleton(true);
							}
							// add amount to PSM
							quantifiedPSM.addAmount(amount);
						} catch (NumberFormatException e) {
							// skip this
						}
					}
				}
			}
			// REF_INT
			// heavy peptide peak area from reconstructed
			// chromatogram
			if (mapValues.containsKey(REF_INT)) {
				try {
					final double value = Double.valueOf(mapValues.get(REF_INT));
					QuantAmount amount = new QuantAmount(value, AmountType.AREA, getHeavyCondition(conditionsByLabels));
					if (singleton && amount.getValue() != 0.0) {
						amount.setSingleton(true);
					}
					// add amount to PSM
					quantifiedPSM.addAmount(amount);
				} catch (NumberFormatException e) {
					// skip this
				}
			}

			// REGRESSION_FACTOR
			// regression score (r)
			if (mapValues.containsKey(AmountType.REGRESSION_FACTOR.name())) {
				try {
					final double value = Double.valueOf(mapValues.get(AmountType.REGRESSION_FACTOR.name()));
					QuantAmount amount = new QuantAmount(value, AmountType.REGRESSION_FACTOR,
							getLightCondition(conditionsByLabels));
					// add amount to PSM
					quantifiedPSM.addAmount(amount);
				} catch (NumberFormatException e) {
					// skip this
				}
			}

			// create the peptide
			QuantifiedPeptideInterface quantifiedPeptide = null;
			final String peptideKey = KeyUtils.getSequenceKey(quantifiedPSM, true);
			if (StaticQuantMaps.peptideMap.containsKey(peptideKey)) {
				quantifiedPeptide = StaticQuantMaps.peptideMap.getItem(peptideKey);
			} else {
				quantifiedPeptide = new QuantifiedPeptide(quantifiedPSM);
			}
			StaticQuantMaps.peptideMap.addItem(quantifiedPeptide);

			quantifiedPSM.setQuantifiedPeptide(quantifiedPeptide, true);
			// add peptide to map
			if (!localPeptideMap.containsKey(peptideKey)) {
				localPeptideMap.put(peptideKey, quantifiedPeptide);
			}

			if (dbIndex != null) {
				String cleanSeq = quantifiedPSM.getSequence();
				final Set<IndexedProtein> indexedProteins = dbIndex.getProteins(cleanSeq);
				if (indexedProteins.isEmpty()) {
					if (!ignoreNotFoundPeptidesInDB) {
						throw new PeptideNotFoundInDBIndexException("The peptide " + cleanSeq
								+ " is not found in Fasta DB.\nReview the default indexing parameters such as the number of allowed misscleavages.");
					}
					// log.warn("The peptide " + cleanSeq +
					// " is not found in Fasta DB.");
					// continue;
				}
				// create a new Quantified Protein for each
				// indexedProtein
				for (IndexedProtein indexedProtein : indexedProteins) {
					String proteinKey = KeyUtils.getProteinKey(indexedProtein);
					QuantifiedProteinInterface newQuantifiedProtein = null;
					if (StaticQuantMaps.proteinMap.containsKey(proteinKey)) {
						newQuantifiedProtein = StaticQuantMaps.proteinMap.getItem(proteinKey);

					} else {
						newQuantifiedProtein = new QuantifiedProteinFromDBIndexEntry(indexedProtein);

					}
					StaticQuantMaps.proteinMap.addItem(newQuantifiedProtein);
					// add protein to protein map
					taxonomies.addAll(newQuantifiedProtein.getTaxonomies());
					localProteinMap.put(proteinKey, newQuantifiedProtein);
					// add to protein-experiment map
					addToMap(experimentKey, experimentToProteinsMap, proteinKey);
					// add psm to the protein
					newQuantifiedProtein.addPSM(quantifiedPSM, true);
					// add peptide to the protein
					newQuantifiedProtein.addPeptide(quantifiedPeptide, true);
					// add protein to the psm
					quantifiedPSM.addQuantifiedProtein(newQuantifiedProtein, true);
					// add to the map (if it was already
					// there is not a problem, it will be
					// only once)
					addToMap(proteinKey, proteinToPeptidesMap, KeyUtils.getSequenceKey(quantifiedPSM, true));

				}
			}
			// use the already created quantified
			// protein

			for (QuantifiedProteinInterface quantifiedProtein : quantifiedProteins) {
				// add psm to the proteins
				quantifiedProtein.addPSM(quantifiedPSM, true);
				// add protein to the psm
				quantifiedPSM.addQuantifiedProtein(quantifiedProtein, true);
				// add peptide to the protein
				quantifiedProtein.addPeptide(quantifiedPeptide, true);
				// add to the map (if it was already there
				// is not a problem, it will be only once)
				String proteinKey = quantifiedProtein.getKey();
				addToMap(proteinKey, proteinToPeptidesMap, KeyUtils.getSequenceKey(quantifiedPSM, true));
				// add protein to protein map
				localProteinMap.put(proteinKey, quantifiedProtein);
				// add to protein-experiment map
				addToMap(experimentKey, experimentToProteinsMap, proteinKey);
			}

			// in case of quantifing sites, set the sites to the ratios in case
			// of no ambiguities
			if (!getQuantifiedAAs().isEmpty()) {
				for (QuantRatio ratio : quantifiedPSM.getRatios()) {
					for (Character c : getQuantifiedAAs()) {
						if (quantifiedPSM.getSequence().contains(String.valueOf(c))) {
							ratio.setQuantifiedAA(c);
						}
					}
					// check for ambiguity on the quantified site
					int numSites = 0;
					int quantifiedSitePositionInPeptide = -1;
					for (Character c : getQuantifiedAAs()) {
						List<Integer> allPositionsOf = StringUtils.allPositionsOf(quantifiedPSM.getSequence(), c);
						numSites = +allPositionsOf.size();
						if (allPositionsOf.size() == 1) {
							quantifiedSitePositionInPeptide = allPositionsOf.get(0);
						}
					}
					// if no ambiguities
					if (numSites == 1) {
						ratio.setQuantifiedSitePositionInPeptide(quantifiedSitePositionInPeptide);
					}
				}
			}
		} catch (IllegalArgumentException e) {
			e.printStackTrace();
			log.warn(e);
			log.info("Error reading line '" + line + "' from file. Skipping it...");

		}

	}

	private boolean isSkipNonResolvedPeaks() {
		return skipNonResolvedPeaks;
	}

	public void setSkipNonResolvedPeaks(boolean skipNonResolvedPeaks) {
		this.skipNonResolvedPeaks = skipNonResolvedPeaks;
	}

	private QuantCondition getLightCondition(Map<QuantificationLabel, QuantCondition> conditionsByLabels) {
		for (QuantificationLabel label : conditionsByLabels.keySet()) {
			if (label.isLight()) {
				return conditionsByLabels.get(label);
			}
		}
		return null;
	}

	private QuantCondition getHeavyCondition(Map<QuantificationLabel, QuantCondition> conditionsByLabels) {
		for (QuantificationLabel label : conditionsByLabels.keySet()) {
			if (label.isHeavy()) {
				return conditionsByLabels.get(label);
			}
		}
		return null;
	}

	private QuantCondition getMediumCondition(Map<QuantificationLabel, QuantCondition> conditionsByLabels) {
		for (QuantificationLabel label : conditionsByLabels.keySet()) {
			if (label.isMedium()) {
				return conditionsByLabels.get(label);
			}
		}
		return null;
	}

	private boolean isTMT(QuantificationLabel label) {
		if (QuantificationLabel.TMT_6PLEX_126 == label || QuantificationLabel.TMT_6PLEX_127 == label
				|| QuantificationLabel.TMT_6PLEX_128 == label || QuantificationLabel.TMT_6PLEX_129 == label
				|| QuantificationLabel.TMT_6PLEX_130 == label || QuantificationLabel.TMT_6PLEX_131 == label) {
			return true;
		}
		return false;
	}

	/**
	 * Gets how the header of the normalized intensity for each channel should
	 * start
	 *
	 * @param label
	 * @return
	 */
	private String getHeaderForTMTLabel(QuantificationLabel label) {
		switch (label) {
		case TMT_6PLEX_126:
			return "norm_ratio(126.127725)";
		case TMT_6PLEX_127:
			return "norm_ratio(127.12476)";
		case TMT_6PLEX_128:
			return "norm_ratio(128.134433)";
		case TMT_6PLEX_129:
			return "norm_ratio(129.131468)";
		case TMT_6PLEX_130:
			return "norm_ratio(130.141141)";
		case TMT_6PLEX_131:
			return "norm_ratio(131.138176)";
		default:
			return null;
		}
	}

	private QuantifiedProteinInterface processProteinLine(String line, List<String> pLineHeaderList,
			Map<QuantificationLabel, QuantCondition> conditionsByLabels, List<RatioDescriptor> ratioDescriptors,
			int numDecoy) throws DiscardProteinException {
		// new protein
		MyHashMap<String, String> mapValues = getMapFromPLine(pLineHeaderList, line);
		String proteinACC = mapValues.get(LOCUS);

		QuantifiedProteinInterface quantifiedProtein = null;
		if (StaticQuantMaps.proteinMap.containsKey(proteinACC)) {
			quantifiedProtein = StaticQuantMaps.proteinMap.getItem(proteinACC);
		} else {
			quantifiedProtein = new QuantifiedProtein(proteinACC);
			String description = mapValues.get(DESCRIPTION);
			quantifiedProtein.setDescription(description);
		}
		StaticQuantMaps.proteinMap.addItem(quantifiedProtein);

		// apply the pattern if available
		if (decoyPattern != null) {
			final Matcher matcher = decoyPattern.matcher(proteinACC);
			if (matcher.find()) {
				quantifiedProtein = null;
				numDecoy++;
				throw new DiscardProteinException("Protein " + proteinACC + " is DECOY.");
			}
		}

		// if we only have one ratio descriptor, it is L over H
		if (ratioDescriptors.size() == 1) {
			QuantificationLabel labelNumerator = ratioDescriptors.get(0).getLabel1();
			QuantificationLabel labelDenominator = ratioDescriptors.get(0).getLabel2();
			String ratioSuffix = ratioDescriptors.get(0).getRatioSuffix();
			// add protein ratio
			// first look if the composite ratio is calculated
			if (mapValues.containsKey(COMPOSITE_RATIO)) {
				try {
					final double ratioValue = Double.valueOf(mapValues.get(COMPOSITE_RATIO));
					String stdValue = null;
					if (mapValues.containsKey(COMPOSITE_RATIO_STANDARD_DEVIATION)) {
						stdValue = mapValues.get(COMPOSITE_RATIO_STANDARD_DEVIATION);
						if ("".equals(stdValue)) {
							stdValue = "0.0";
						}
					}
					QuantRatio ratio = new CensusRatio(ratioValue, stdValue, false, conditionsByLabels, labelNumerator,
							labelDenominator, AggregationLevel.PROTEIN, COMPOSITE_RATIO);
					quantifiedProtein.addRatio(ratio);
				} catch (NumberFormatException e) {
					// skip this
				}
				// if there is not composite ratio, use the
				// regular ratio
			} else if (mapValues.containsKey(COMPOSITE_RATIO + ratioSuffix)) {
				try {
					final double ratioValue = Double.valueOf(mapValues.get(COMPOSITE_RATIO + ratioSuffix));
					String stdValue = null;
					if (mapValues.containsKey(COMPOSITE_RATIO_STANDARD_DEVIATION + ratioSuffix)) {
						stdValue = mapValues.get(COMPOSITE_RATIO_STANDARD_DEVIATION + ratioSuffix);
						if ("".equals(stdValue)) {
							stdValue = "0.0";
						}
					}
					QuantRatio ratio = new CensusRatio(ratioValue, stdValue, false, conditionsByLabels, labelNumerator,
							labelDenominator, AggregationLevel.PROTEIN, COMPOSITE_RATIO);
					quantifiedProtein.addRatio(ratio);
				} catch (NumberFormatException e) {
					// skip this
				}
				// if there is not composite ratio, use the
				// regular ratio
			} else if (mapValues.containsKey(AVERAGE_RATIO) || mapValues.containsKey(AREA_RATIO)) {
				try {
					if (mapValues.containsKey(AVERAGE_RATIO)) {
						final double ratioValue = Double.valueOf(mapValues.get(AVERAGE_RATIO));
						String stdValue = null;
						if (mapValues.containsKey(STANDARD_DEVIATION)) {
							stdValue = mapValues.get(STANDARD_DEVIATION);
							if ("".equals(stdValue)) {
								stdValue = "0.0";
							}
						}
						QuantRatio ratio = new CensusRatio(ratioValue, stdValue, false, conditionsByLabels,
								labelNumerator, labelDenominator, AggregationLevel.PROTEIN, AVERAGE_RATIO);
						quantifiedProtein.addRatio(ratio);
					}
				} catch (NumberFormatException e) {
					// skip this
				}
				try {
					if (mapValues.containsKey(AREA_RATIO)) {
						final double ratioValue = Double.valueOf(mapValues.get(AREA_RATIO));
						QuantRatio ratio = new CensusRatio(ratioValue, null, false, conditionsByLabels, labelNumerator,
								labelDenominator, AggregationLevel.PROTEIN, AREA_RATIO);
						quantifiedProtein.addRatio(ratio);
					}
				} catch (NumberFormatException e) {
					// skip this
				}
			} else if (mapValues.containsKey(AVERAGE_RATIO + ratioSuffix)
					|| mapValues.containsKey(AREA_RATIO + ratioSuffix)) {
				try {
					if (mapValues.containsKey(AVERAGE_RATIO + ratioSuffix)) {
						final double ratioValue = Double.valueOf(mapValues.get(AVERAGE_RATIO + ratioSuffix));
						String stdValue = null;
						if (mapValues.containsKey(STANDARD_DEVIATION)) {
							stdValue = mapValues.get(STANDARD_DEVIATION);
							if ("".equals(stdValue)) {
								stdValue = "0.0";
							}
						}
						QuantRatio ratio = new CensusRatio(ratioValue, stdValue, false, conditionsByLabels,
								labelNumerator, labelDenominator, AggregationLevel.PROTEIN, AVERAGE_RATIO);
						quantifiedProtein.addRatio(ratio);
					}
				} catch (NumberFormatException e) {
					// skip this
				}
				try {
					if (mapValues.containsKey(AREA_RATIO + ratioSuffix)) {
						final double ratioValue = Double.valueOf(mapValues.get(AREA_RATIO + ratioSuffix));
						QuantRatio ratio = new CensusRatio(ratioValue, null, false, conditionsByLabels, labelNumerator,
								labelDenominator, AggregationLevel.PROTEIN, AREA_RATIO);
						quantifiedProtein.addRatio(ratio);
					}
				} catch (NumberFormatException e) {
					// skip this
				}
			}
		} else {
			// having more than one ratio descriptor, it is L, M, H
			for (RatioDescriptor ratioDescriptor : ratioDescriptors) {
				QuantificationLabel labelNumerator = ratioDescriptor.getLabel1();
				QuantificationLabel labelDenominator = ratioDescriptor.getLabel2();
				String ratioSuffix = ratioDescriptor.getRatioSuffix();
				// add protein ratio
				// first look if the normalized composite ratio is calculated
				if (mapValues.containsKey(NORM_COMPOSITE_RATIO + ratioSuffix)
						|| mapValues.containsKey(COMPOSITE_RATIO + ratioSuffix)) {
					if (mapValues.containsKey(NORM_COMPOSITE_RATIO + ratioSuffix)) {
						try {
							final double ratioValue = Double.valueOf(mapValues.get(NORM_COMPOSITE_RATIO + ratioSuffix));
							String stdValue = null;
							if (mapValues.containsKey(NORM_COMPOSITE_RATIO_STDEV + ratioSuffix)) {
								stdValue = mapValues.get(NORM_COMPOSITE_RATIO_STDEV + ratioSuffix);
								if ("".equals(stdValue)) {
									stdValue = "0.0";
								}
							}
							QuantRatio ratio = new CensusRatio(ratioValue, stdValue, false, conditionsByLabels,
									labelNumerator, labelDenominator, AggregationLevel.PROTEIN, NORM_COMPOSITE_RATIO);
							quantifiedProtein.addRatio(ratio);
						} catch (NumberFormatException e) {
							// skip this
						}
						// if there is not composite ratio, use the
						// regular composite ratio
					}
					if (mapValues.containsKey(COMPOSITE_RATIO + ratioSuffix)) {
						try {
							final double ratioValue = Double.valueOf(mapValues.get(COMPOSITE_RATIO + ratioSuffix));
							String stdValue = null;
							if (mapValues.containsKey(COMPOSITE_RATIO_STDEV + ratioSuffix)) {
								stdValue = mapValues.get(COMPOSITE_RATIO_STDEV + ratioSuffix);
								if ("".equals(stdValue)) {
									stdValue = "0.0";
								}
							}
							QuantRatio ratio = new CensusRatio(ratioValue, stdValue, false, conditionsByLabels,
									labelNumerator, labelDenominator, AggregationLevel.PROTEIN, COMPOSITE_RATIO);
							quantifiedProtein.addRatio(ratio);
						} catch (NumberFormatException e) {
							// skip this
						}
						// if there is not composite ratio, use the
						// regular composite ratio
					}
				} else if (mapValues.containsKey(MEDIAN_NORM_RATIO + ratioSuffix)
						|| mapValues.containsKey(MEDIAN_AREA_RATIO + ratioSuffix)) {
					try {
						if (mapValues.containsKey(MEDIAN_NORM_RATIO + ratioSuffix)) {
							final double ratioValue = Double.valueOf(mapValues.get(MEDIAN_NORM_RATIO + ratioSuffix));
							String stdValue = null;
							if (mapValues.containsKey(NORM_STDEV + ratioSuffix)) {
								stdValue = mapValues.get(NORM_STDEV + ratioSuffix);
								if ("".equals(stdValue)) {
									stdValue = "0.0";
								}
							}
							QuantRatio ratio = new CensusRatio(ratioValue, stdValue, false, conditionsByLabels,
									labelNumerator, labelDenominator, AggregationLevel.PROTEIN, NORM_STDEV);
							quantifiedProtein.addRatio(ratio);
						}
					} catch (NumberFormatException e) {
						// skip this
					}
					try {
						if (mapValues.containsKey(MEDIAN_AREA_RATIO + ratioSuffix)) {
							final double ratioValue = Double.valueOf(mapValues.get(MEDIAN_AREA_RATIO + ratioSuffix));
							QuantRatio ratio = new CensusRatio(ratioValue, null, false, conditionsByLabels,
									labelNumerator, labelDenominator, AggregationLevel.PROTEIN, MEDIAN_AREA_RATIO);
							quantifiedProtein.addRatio(ratio);
						}
					} catch (NumberFormatException e) {
						// skip this
					}
				}
			}
		}
		return quantifiedProtein;

	}

	private MyHashMap<String, String> getMapFromPLine(List<String> pLineHeaderList, String line) {
		MyHashMap<String, String> map = new MyHashMap<String, String>();
		String[] split = line.split("\t");
		if (split.length != pLineHeaderList.size()) {
			// try to see if there is one that is empty
			String[] splitTMP = new String[split.length - 1];
			if (split[5].equals("")) {

				for (int i = 0; i < split.length - 1; i++) {
					if (i < 5) {
						splitTMP[i] = split[i];
					} else {
						splitTMP[i] = split[i + 1];
					}
				}
				split = splitTMP;
			}
		}
		if (split.length == pLineHeaderList.size()) {
			int i = 1;
			for (String header : pLineHeaderList) {
				if (header.equals(PLINE)) {
					continue;
				}
				map.put(header, split[i]);
				i++;
			}
		} else {
			throw new IllegalArgumentException(
					line + " has different number of columns than the header which have " + pLineHeaderList.size());
		}
		return map;
	}

	private MyHashMap<String, String> getMapFromSLine(List<String> sLineHeaderList, String line) {
		MyHashMap<String, String> map = new MyHashMap<String, String>();
		String[] split = line.split("\t");
		// not remove the last element:
		final int upToThisIndex = split.length - 1;
		split = removeElements(split, "N/A", upToThisIndex);

		// trying to recover some lines in which the peptide information is
		// shifted to the right one column from the sequence
		if ("".equals(split[1]) && "".equals(split[2]) && split.length == sLineHeaderList.size() + 1) {
			// copy into a new array
			String[] newSplit = new String[sLineHeaderList.size()];
			for (int index = 0; index < newSplit.length; index++) {
				if (index <= 1) {
					newSplit[index] = split[index];
				} else {
					newSplit[index] = split[index + 1];
				}
			}
			split = newSplit;
		}
		if (split.length == sLineHeaderList.size() || split.length + 1 == sLineHeaderList.size()) {
			int i = 0;
			for (String header : sLineHeaderList) {
				if (i < split.length) {
					map.put(header, split[i]);
				}
				i++;
			}
		} else {
			throw new IllegalArgumentException(line + " has different number of columns than the header which have "
					+ sLineHeaderList.size() + ". LINE HAS " + split.length);
		}
		return map;
	}

	private String[] removeElements(String[] split, String elementToRemove, int upToThisIndex) {
		List<String> list = new ArrayList<String>();
		int index = -1;
		for (String splitElement : split) {
			index++;
			if (splitElement.equals(elementToRemove) && index < upToThisIndex) {
				continue;
			}
			list.add(splitElement);

		}
		return list.toArray(new String[0]);
	}

	/**
	 * @param onlyOneSpectrumPerChromatographicPeakAndPerSaltStep
	 *            the onlyOneSpectrumPerChromatographicPeakAndPerSaltStep to set
	 */
	public void setOnlyOneSpectrumPerChromatographicPeakAndPerSaltStep(
			boolean onlyOneSpectrumPerChromatographicPeakAndPerSaltStep) {
		this.onlyOneSpectrumPerChromatographicPeakAndPerSaltStep = onlyOneSpectrumPerChromatographicPeakAndPerSaltStep;
	}

	/**
	 * @return the skipSingletons
	 */
	public boolean isSkipSingletons() {
		return skipSingletons;
	}

	/**
	 * @param skipSingletons
	 *            the skipSingletons to set
	 */
	public void setSkipSingletons(boolean skipSingletons) {
		this.skipSingletons = skipSingletons;
	}

}
