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
import java.util.stream.Collectors;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantCondition;
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
import edu.scripps.yates.census.read.util.QuantKeyUtils;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.dbindex.util.PeptideNotFoundInDBIndexException;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.proteomicsmodel.PTMSite;
import edu.scripps.yates.utilities.proteomicsmodel.Ratio;
import edu.scripps.yates.utilities.proteomicsmodel.Score;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AggregationLevel;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AmountType;
import edu.scripps.yates.utilities.proteomicsmodel.factories.PTMSiteEx;
import edu.scripps.yates.utilities.proteomicsmodel.factories.ScoreEx;
import edu.scripps.yates.utilities.proteomicsmodel.utils.KeyUtils;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import edu.scripps.yates.utilities.strings.StringUtils;
import gnu.trove.list.array.TIntArrayList;
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
	public static final String AVERAGE_RATIO = "AVERAGE_RATIO";
	private static final String STANDARD_DEVIATION = "STANDARD_DEVIATION";
	private static final String NORM_STDEV = "NORM_STDEV";
	public static final String COMPOSITE_RATIO = "COMPOSITE_RATIO";
	public static final String NORM_COMPOSITE_RATIO = "NORM_COMPOSITE_RATIO";
	private static final String COMPOSITE_RATIO_STANDARD_DEVIATION = "COMPOSITE_RATIO_STANDARD_DEVIATION";
	private static final String NORM_COMPOSITE_RATIO_STDEV = "NORM_COMPOSITE_RATIO_STDEV";
	private static final String COMPOSITE_RATIO_STDEV = "COMPOSITE_RATIO_STDEV";
	public static final String MEDIAN_NORM_RATIO = "MEDIAN_NORM_RATIO";
	public static final String MEDIAN_AREA_RATIO = "MEDIAN_AREA_RATIO";
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

	public static final String XCORR = "XCorr";
	public static final String DELTACN = "deltaCN";

	public static final String DET_FACTOR = "DET_FACTOR";

	public static final String TMT_PURITY = "TMT_purity";
	public static final String SIGNAL_TO_NOISE = "Signal-noise";
	public static final String LOCALIZATION_SCORE = "Localization_Score";
	public static final String ION_COUNT = "Ion Count";
	private boolean onlyOneSpectrumPerChromatographicPeakAndPerSaltStep = false;
	private boolean skipSingletons = false; // by default
	private boolean skipNonResolvedPeaks = true; // by default

	/**
	 * TMT plex. Only one per parser, and it will throw exception if different
	 * plexes are found in different files of the parser
	 */
	private Integer tmtPlex = null;

	public CensusOutParser() {
		super();
	}

	public CensusOutParser(List<RemoteSSHFileReference> remoteSSHServers,
			List<Map<QuantificationLabel, QuantCondition>> conditionsByLabels, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) {
		super(remoteSSHServers, conditionsByLabels, labelNumerator, labelDenominator);
	}

	public CensusOutParser(Map<QuantificationLabel, QuantCondition> conditionsByLabels,
			Collection<RemoteSSHFileReference> remoteSSHServers, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) {
		super(conditionsByLabels, remoteSSHServers, labelNumerator, labelDenominator);
	}

	public CensusOutParser(RemoteSSHFileReference remoteSSHServer,
			Map<QuantificationLabel, QuantCondition> conditionsByLabels, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) throws FileNotFoundException {
		super(remoteSSHServer, conditionsByLabels, labelNumerator, labelDenominator);
	}

	/**
	 * This is suitable for TMT data, since the map between labels and conditions
	 * allow multiple labels for the same condition
	 * 
	 * @param inputFile
	 * @param conditionsByLabels
	 * @throws FileNotFoundException
	 */
	public CensusOutParser(File inputFile, Map<QuantificationLabel, QuantCondition> conditionsByLabels)
			throws FileNotFoundException {
		super(inputFile, conditionsByLabels);
	}

	public CensusOutParser(File xmlFile, Map<QuantificationLabel, QuantCondition> conditionsByLabels,
			QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) throws FileNotFoundException {
		super(xmlFile, conditionsByLabels, labelNumerator, labelDenominator);
	}

	public CensusOutParser(File[] xmlFiles, Map<QuantificationLabel, QuantCondition> conditionsByLabels,
			QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) throws FileNotFoundException {
		super(xmlFiles, conditionsByLabels, labelNumerator, labelDenominator);
	}

	public CensusOutParser(File[] xmlFiles, Map<QuantificationLabel, QuantCondition>[] conditionsByLabels,
			QuantificationLabel[] labelNumerator, QuantificationLabel[] labelDenominator) throws FileNotFoundException {
		super(xmlFiles, conditionsByLabels, labelNumerator, labelDenominator);
	}

	public CensusOutParser(Collection<File> xmlFiles, Map<QuantificationLabel, QuantCondition> conditionsByLabels,
			QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) throws FileNotFoundException {
		super(xmlFiles, conditionsByLabels, labelNumerator, labelDenominator);
	}

	public CensusOutParser(RemoteSSHFileReference remoteServer, QuantificationLabel label1, QuantCondition cond1,
			QuantificationLabel label2, QuantCondition cond2) {
		super(remoteServer, label1, cond1, label2, cond2);
	}

	public CensusOutParser(File inputFile, QuantificationLabel label1, QuantCondition cond1, QuantificationLabel label2,
			QuantCondition cond2) throws FileNotFoundException {
		super(inputFile, label1, cond1, label2, cond2);
	}

	public CensusOutParser(File inputFile, Map<QuantificationLabel, QuantCondition> conditionsByLabels,
			QuantificationLabel light, QuantificationLabel medium, QuantificationLabel heavy)
			throws FileNotFoundException {
		super(inputFile, conditionsByLabels, light, medium, heavy);
	}

	public Integer getTmtPlex() throws IOException {
		if (tmtPlex == null) {
			for (final RemoteSSHFileReference remoteFileRetriever : super.remoteFileRetrievers) {
				BufferedReader br = null;
				try {
					br = new BufferedReader(
							new InputStreamReader(new BufferedInputStream(remoteFileRetriever.getRemoteInputStream())));
					String line;
					final List<String> pLineHeaderList = new ArrayList<String>();
					final List<String> sLineHeaderList = new ArrayList<String>();
					final List<String> singletonSLineHeaderList = new ArrayList<String>();
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
						} else {
							final List<String> tmtHeaders = getTMTHeaders(sLineHeaderList);
							final int plex = tmtHeaders.size();
							// check if we find different plexes in other file, in that case, we throw an
							// error (this complicates things a lot)
							if (tmtPlex != null && Integer.compare(tmtPlex, plex) != 0) {
								throw new IllegalArgumentException(
										"Multiple TMT plexes have been found (" + tmtPlex + " and " + plex
												+ ". Only one type of TMT plex is allowed in a single parser.");
							}
							tmtPlex = plex;
							break;
						}
					}

				} catch (final IOException e) {
					e.printStackTrace();
					throw e;
				} finally {
					if (br != null) {
						br.close();
					}
				}
			}
		}
		if (tmtPlex == 0) {
			tmtPlex = null;
		}
		return tmtPlex;
	}

	/**
	 * Returns a list of headers that start by "m/z_"
	 * 
	 * @param sLineHeaderList
	 * @return
	 */
	private List<String> getTMTHeaders(List<String> sLineHeaderList) {
		final List<String> ret = new ArrayList<String>();
		for (final String header : sLineHeaderList) {
			if (header.startsWith("m/z_")) {
				ret.add(header);
			}
		}
		return ret;
	}

	/**
	 *
	 * @param writeFiles whether to write output files necessary to run SanXot
	 *                   program
	 * @throws WrongTMTLabels, IOException, QuantParserException
	 */
	@Override
	protected void process() throws QuantParserException {
		processed = false;
		log.info("Processing quant file...");

		try {
			checkLabels();
			int numDecoy = 0;
			boolean someValidFile = false;
			for (final RemoteSSHFileReference remoteFileRetriever : remoteFileRetrievers) {
				final Map<QuantCondition, Set<QuantificationLabel>> labelsByConditions = QuantUtils
						.getLabelsByConditions(conditionsByLabelsByFile.get(remoteFileRetriever));
				final Map<QuantificationLabel, QuantCondition> conditionsByLabels = conditionsByLabelsByFile
						.get(remoteFileRetriever);

				final List<RatioDescriptor> ratioDescriptors = ratioDescriptorsByFile.get(remoteFileRetriever);
				// QuantificationLabel labelNumerator =
				// numeratorLabelByFile.get(remoteFileRetriever);
				// QuantificationLabel labelDenominator =
				// denominatorLabelByFile.get(remoteFileRetriever);

				final String experimentKey = FilenameUtils
						.getBaseName(remoteFileRetriever.getOutputFile().getAbsolutePath());
				final String fileName = FilenameUtils.getName(remoteFileRetriever.getOutputFile().getAbsolutePath());
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
					final List<String> pLineHeaderList = new ArrayList<String>();
					final List<String> sLineHeaderList = new ArrayList<String>();
					final List<String> singletonSLineHeaderList = new ArrayList<String>();
					final Set<QuantifiedProteinInterface> proteinGroup = new THashSet<QuantifiedProteinInterface>();
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
								final QuantifiedProteinInterface quantifiedProtein = processProteinLine(line,
										pLineHeaderList, conditionsByLabels, ratioDescriptors, experimentKey);
								proteinGroup.add(quantifiedProtein);
							} catch (final DiscardProteinException e) {
								numDecoy++;
								continue;
							} catch (final IllegalArgumentException e) {
								e.printStackTrace();
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

				} catch (final PeptideNotFoundInDBIndexException e) {
					e.printStackTrace();
					if (!super.ignoreNotFoundPeptidesInDB) {
						throw e;
					}
				} catch (final Exception e) {
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
					final Map<String, Set<QuantifiedPSMInterface>> map = new THashMap<String, Set<QuantifiedPSMInterface>>();
					for (final QuantifiedPSMInterface psm : localPsmMap.values()) {
						final String key = getSpectrumPerChromatographicPeakAndPerSaltStepKey(psm);
						if (map.containsKey(key)) {
							map.get(key).add(psm);
						} else {
							final Set<QuantifiedPSMInterface> set = new THashSet<QuantifiedPSMInterface>();
							set.add(psm);
							map.put(key, set);
						}
					}
					// once the map is populated,
					// look for each key, if we have more than one psm
					// in that case, select the best one
					int numRemoved = 0;
					for (final String key : map.keySet()) {
						final Set<QuantifiedPSMInterface> psmSet = map.get(key);
						if (psmSet.size() > 1) {
							final QuantifiedPSMInterface bestPSM = getBestPSM(psmSet);
							final Set<QuantifiedPSMInterface> toIgnore = new THashSet<QuantifiedPSMInterface>();
							for (final QuantifiedPSMInterface psm : psmSet) {
								if (!psm.equals(bestPSM)) {
									toIgnore.add(psm);
								}
							}
							// remove the psms in toIgnore Set
							if (!toIgnore.isEmpty()) {
								for (final QuantifiedPSMInterface psmToIgnore : toIgnore) {
									numRemoved++;
									localPsmMap.remove(psmToIgnore.getKey());
									StaticQuantMaps.psmMap.remove(psmToIgnore);
									// remove it from its peptide
									final QuantifiedPeptideInterface quantifiedPeptide = psmToIgnore
											.getQuantifiedPeptide();
									if (quantifiedPeptide != null) {
										quantifiedPeptide.getQuantifiedPSMs().remove(psmToIgnore);
									}
									if (quantifiedPeptide.getQuantifiedPSMs().isEmpty()) {
										localPeptideMap.remove(quantifiedPeptide.getKey());
										StaticQuantMaps.peptideMap.remove(quantifiedPeptide);
									}
									// remove it from its proteins
									final Set<QuantifiedProteinInterface> quantifiedProteins = psmToIgnore
											.getQuantifiedProteins();
									for (final QuantifiedProteinInterface protein : quantifiedProteins) {
										protein.getQuantifiedPSMs().remove(psmToIgnore);
										if (protein.getQuantifiedPSMs().isEmpty()) {
											localProteinMap.remove(protein.getKey());
											StaticQuantMaps.proteinMap.remove(protein);
										}
									}

								}
							}
						}
					}
					log.info(numRemoved + " PSMs were redundant and removed.");

				}
			}
			if (!someValidFile)
				throw new IllegalArgumentException("some error occurred while reading the files");

			processed = true;
		} catch (final IOException e) {
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
	 * This will check whether the labels provided in the constructor match the ones
	 * detected in case of having a TMT file
	 * 
	 * @throws IOException
	 * @throws WrongTMTLabels
	 */
	private void checkLabels() throws IOException {
		if (conditionsByLabelsByFile == null || conditionsByLabelsByFile.isEmpty()) {
			return;
		}
		if (getTmtPlex() != null & getTmtPlex() > 0) {
			final int tmtPlex = getTmtPlex();
			final List<QuantificationLabel> possibleLabels = QuantificationLabel.getTMTPlexLabels(tmtPlex);

			for (final RemoteSSHFileReference remoteFile : remoteFileRetrievers) {
				if (conditionsByLabelsByFile.containsKey(remoteFile)) {
					final Set<QuantificationLabel> labelsFromConstructor = conditionsByLabelsByFile.get(remoteFile)
							.keySet();
					for (final QuantificationLabel possibleLabel : possibleLabels) {
						if (!labelsFromConstructor.contains(possibleLabel)) {
							// now we check whether the labels from constructor are the same plex
							for (final QuantificationLabel labelFromConstructor : labelsFromConstructor) {
								if (QuantificationLabel.getTMTPlex(labelFromConstructor) != tmtPlex) {
									throw new IllegalArgumentException(
											"Labels in the constructor must be TMT " + tmtPlex + " labels");
								}
							}
							// if all are from the same plex, we just write a warning saying that the label
							// is not used
							log.warn("Label from TMT plex " + tmtPlex + " (" + possibleLabel + ") is not used");
						}
					}
				}
			}
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
		final List<QuantifiedPSMInterface> list = new ArrayList<QuantifiedPSMInterface>();
		list.addAll(psmSet);
		Collections.sort(list, new Comparator<QuantifiedPSMInterface>() {

			@Override
			public int compare(QuantifiedPSMInterface o1, QuantifiedPSMInterface o2) {
				// regression_factor
				final Double regressionFactor1 = getRatioScore(REGRESSION_FACTOR, o1);
				final Double regressionFactor2 = getRatioScore(REGRESSION_FACTOR, o2);
				if (regressionFactor1 != null && regressionFactor2 != null) {
					final int compare = Double.compare(regressionFactor2, regressionFactor1);
					if (compare != 0) {
						return compare;
					}
				}
				// xcorr
				final Float xcorr1 = o1.getXCorr();
				final Float xcorr2 = o2.getXCorr();
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
					for (final Ratio quantRatio : o1.getRatios()) {
						if (quantRatio.getAssociatedConfidenceScore() != null) {
							if (quantRatio.getAssociatedConfidenceScore().getScoreName().equals(scoreName)) {
								try {
									return Double.valueOf(quantRatio.getAssociatedConfidenceScore().getValue());
								} catch (final NumberFormatException e) {

								}
							}
						}
					}
				}
				if (o1.getAmounts() != null) {
					for (final Amount amount : o1.getAmounts()) {
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
		final StringBuilder sb = new StringBuilder();
		String rawFileName = psm.getRawFileNames().iterator().next();

		if (rawFileName.startsWith("H") || rawFileName.startsWith("M")) {
			final String substring = rawFileName.substring(1);
			// only take it if there is another one like that in the static map
			// this will avoid to consider an H or an M in the file name that
			// were part of the name and not meaning anything
			if (StaticQuantMaps.rawFileNamesContains(substring)) {
				rawFileName = substring;
			}
		}

		sb.append(psm.getSequence()).append("_").append(rawFileName).append("_").append(psm.getChargeState());
		return sb.toString();
	}

	private void processPSMLine(String line, List<String> sLineHeaderList,
			Set<QuantifiedProteinInterface> quantifiedProteins,
			Map<QuantificationLabel, QuantCondition> conditionsByLabels,
			Map<QuantCondition, Set<QuantificationLabel>> labelsByConditions, List<RatioDescriptor> ratioDescriptors,
			String experimentKey, RemoteSSHFileReference remoteFileRetriever, boolean singleton)
			throws IOException, WrongTMTLabels {

		// new psm
		try {
			final MyHashMap<String, String> mapValues = getMapFromSLine(sLineHeaderList, line);

			final String sequence = mapValues.get(SEQUENCE);

			// dont look into the QuantifiedPSM.map because each
			// line is always a new PSM
			final String inputFileName = FilenameUtils.getName(remoteFileRetriever.getOutputFile().getAbsolutePath());
			String rawFileName = null;
			if (mapValues.containsKey(FILENAME)) {
				rawFileName = mapValues.get(FILENAME);
			} else {
				rawFileName = inputFileName;
			}
			StaticQuantMaps.addRawFileName(rawFileName);
			// scan number
			String scanNumber = "0";
			if (mapValues.containsKey(SCAN)) {
				scanNumber = String.valueOf(Double.valueOf(mapValues.get(SCAN)).intValue());
			}
			QuantifiedPSMInterface quantifiedPSM = null;
			// if (!isGetPTMInProteinMap()) {
			// quantifiedPSM = new QuantifiedPSM(sequence, labelsByConditions,
			// peptideToSpectraMap, scanNumber,
			// Double.valueOf(mapValues.get(CS)).intValue(), rawFileName,
			// singleton);
			// } else {
			quantifiedPSM = new QuantifiedPSM(sequence, peptideToSpectraMap, scanNumber,
					Double.valueOf(mapValues.get(CS)).intValue(), rawFileName, singleton,
					isDistinguishModifiedSequences(), isChargeSensible());
			// }

			// xcorr
			Float xcorr = null;
			if (mapValues.containsKey(XCORR)) {
				try {
					xcorr = Float.valueOf(mapValues.get(XCORR));
					quantifiedPSM.setXCorr(xcorr);
				} catch (final NumberFormatException e) {

				}
			}
			// deltacn
			Float deltaCn = null;
			if (mapValues.containsKey(DELTACN)) {
				try {
					deltaCn = Float.valueOf(mapValues.get(DELTACN));
					quantifiedPSM.setDeltaCn(deltaCn);
				} catch (final NumberFormatException e) {

				}
			}
			// tmt purity
			if (mapValues.containsKey(TMT_PURITY)) {
				try {
					final Float tmtPurity = Float.valueOf(mapValues.get(TMT_PURITY));
					final Score score = new ScoreEx(String.valueOf(tmtPurity), TMT_PURITY, TMT_PURITY, TMT_PURITY);
					quantifiedPSM.addScore(score);
				} catch (final NumberFormatException e) {

				}
			}
			// signal to noise
			if (mapValues.containsKey(SIGNAL_TO_NOISE)) {
				try {
					final Float signalToNoise = Float.valueOf(mapValues.get(SIGNAL_TO_NOISE));
					final Score score = new ScoreEx(String.valueOf(signalToNoise), SIGNAL_TO_NOISE, SIGNAL_TO_NOISE,
							SIGNAL_TO_NOISE);
					quantifiedPSM.addScore(score);
				} catch (final NumberFormatException e) {

				}
			}
			// ion count
			if (mapValues.containsKey(ION_COUNT)) {
				try {
					final Float ionCount = Float.valueOf(mapValues.get(ION_COUNT));
					final Score score = new ScoreEx(String.valueOf(ionCount), ION_COUNT, ION_COUNT, ION_COUNT);
					quantifiedPSM.addScore(score);
				} catch (final NumberFormatException e) {

				}
			}
			// localization score
			String localizationScore = null;
			if (mapValues.containsKey(LOCALIZATION_SCORE)) {

				localizationScore = mapValues.get(LOCALIZATION_SCORE);

			}

			quantifiedPSM.getFileNames().add(inputFileName);
			final String psmKey = KeyUtils.getInstance().getSpectrumKey(quantifiedPSM, isDistinguishModifiedSequences(),
					isChargeSensible());
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

			// if we have ratios
			if (ratioDescriptors != null && !ratioDescriptors.isEmpty()) {
				// if there is only one ratio
				if (ratioDescriptors.size() == 1) {
					final QuantificationLabel labelNumerator = ratioDescriptors.get(0).getLabel1();
					final QuantificationLabel labelDenominator = ratioDescriptors.get(0).getLabel2();
					String ratioSuffix = ratioDescriptors.get(0).getRatioSuffix();
					// PSM regular ratio
					// add PSM ratios from census out
					if (mapValues.containsKey(RATIO)) {
						ratioSuffix = "";
					} else if (mapValues.containsKey(RATIO + ratioSuffix)) {

					}
					if (mapValues.containsKey(RATIO + ratioSuffix)) {
						try {

							final double ratioValue = QuantUtils
									.parseCensusRatioValue(mapValues.get(RATIO + ratioSuffix));
							CensusRatio ratio = new CensusRatio(ratioValue, false, conditionsByLabels, labelNumerator,
									labelDenominator, AggregationLevel.PSM, RATIO);
							RatioScore ratioScore = null;
							// profile score
							String scoreValue = null;
							// check in any case the regressionFactor. If that
							// is
							// -1,
							// then
							// convert the ratio from 0 to inf
							String regressionFactor = null;
							if (mapValues.containsKey(PROFILE_SCORE) && // PROFILE
																		// SCORE
																		// IS
																		// NOT
																		// SPECIFIC
																		// OF A
																		// PAIR
																		// OF
																		// LABELS
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
											"PSM-level quantification confidence metric",
											"PSM-level Determining factor");
									ratio.setRatioScore(ratioScore);
								}
							}
							// get the regression_factor anyway, but add it only
							// in
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
										// check area_ratio value. If the
										// are_ratio
										// is < 1, leave it as 0. If the
										// area_ratio
										// is > 1,
										// convert it to +INF.
										// note that all numbers are not log
										// numbers.
										if (mapValues.containsKey(AREA_RATIO + ratioSuffix)) {
											final double areaRatioValue = QuantUtils
													.parseCensusRatioValue(mapValues.get(AREA_RATIO + ratioSuffix));
											if (Double.isInfinite(areaRatioValue) || areaRatioValue > 1) {
												ratio = new CensusRatio(Double.POSITIVE_INFINITY, false,
														conditionsByLabels, labelNumerator, labelDenominator,
														AggregationLevel.PSM, RATIO);
												if (ratioScore != null) {
													ratio.setRatioScore(ratioScore);
												}
											}
										}
									}
								}
							} catch (final NumberFormatException e) {
								// do nothing
							}
							// add ratio to PSM
							quantifiedPSM.addRatio(ratio);

						} catch (final NumberFormatException e) {
							// skip this
						}
					}

					// PSM AREA RATIO
					// add PSM ratios from census out
					if (mapValues.containsKey(AREA_RATIO + ratioSuffix)) {
						try {
							final double ratioValue = QuantUtils
									.parseCensusRatioValue(mapValues.get(AREA_RATIO + ratioSuffix));

							final CensusRatio ratio = new CensusRatio(ratioValue, false, conditionsByLabels,
									labelNumerator, labelDenominator, AggregationLevel.PSM, AREA_RATIO);
							// profile score
							// PROFILE SCORE IS NOT SPECIFIC FOR A PAIR OF
							// LABELS
							if (mapValues.containsKey(PROFILE_SCORE)) {
								final String scoreValue = mapValues.get(PROFILE_SCORE);
								if (!"NA".equals(scoreValue)) {
									final RatioScore ratioScore = new RatioScore(scoreValue, PROFILE_SCORE,
											"PSM-level quantification confidence metric",
											"fitting score comparing peak and gaussian distribution");
									ratio.setRatioScore(ratioScore);
								}
							} else if (mapValues.containsKey(SINGLETON_SCORE)) {
								final String scoreValue = mapValues.get(SINGLETON_SCORE);
								if (!"NA".equals(scoreValue)) {
									final RatioScore ratioScore = new RatioScore(scoreValue, SINGLETON_SCORE,
											"PSM-level quantification confidence metric", "Singleton score");
									ratio.setRatioScore(ratioScore);
								}
							}
							// add ratio to PSM
							quantifiedPSM.addRatio(ratio);
						} catch (final NumberFormatException e) {
							// skip this
						}
					}
					// PSM AREA RATIO
					// add PSM ratios from census out
					if (mapValues.containsKey(NORM_RATIO + ratioSuffix)) {
						try {
							final double ratioValue = QuantUtils
									.parseCensusRatioValue(mapValues.get(NORM_RATIO + ratioSuffix));

							final CensusRatio ratio = new CensusRatio(ratioValue, false, conditionsByLabels,
									labelNumerator, labelDenominator, AggregationLevel.PSM, NORM_RATIO);
							// profile score
							// PROFILE SCORE IS NOT SPECIFIC FOR A PAIR OF
							// LABELS
							if (mapValues.containsKey(PROFILE_SCORE)) {
								final String scoreValue = mapValues.get(PROFILE_SCORE);
								if (!"NA".equals(scoreValue)) {
									final RatioScore ratioScore = new RatioScore(scoreValue, PROFILE_SCORE,
											"PSM-level quantification confidence metric",
											"fitting score comparing peak and gaussian distribution");
									ratio.setRatioScore(ratioScore);
								}
							} else if (mapValues.containsKey(SINGLETON_SCORE)) {
								final String scoreValue = mapValues.get(SINGLETON_SCORE);
								if (!"NA".equals(scoreValue)) {
									final RatioScore ratioScore = new RatioScore(scoreValue, SINGLETON_SCORE,
											"PSM-level quantification confidence metric", "Singleton score");
									ratio.setRatioScore(ratioScore);
								}
							}
							// add ratio to PSM
							quantifiedPSM.addRatio(ratio);
						} catch (final NumberFormatException e) {
							// skip this
						}
					}
					// TMT
					if (getTmtPlex() != null) {
						final int plex = getTmtPlex();
						final List<QuantificationLabel> labels = QuantificationLabel.getTMTPlexLabels(plex);
						if (labels.contains(labelNumerator) && labels.contains(labelDenominator)) {
							final int channelNumerator = labels.indexOf(labelNumerator) + 1;
							final int channelDenominator = labels.indexOf(labelDenominator) + 1;
							if (channelNumerator <= 0 || channelDenominator <= 0) {
								log.warn("Intensities for labels " + labelNumerator + " and " + labelDenominator
										+ " are not found in file with TMT plex " + plex);
							}
							// numerator
							Double numeratorIntensity = null;
							final String headerNumerator = getHeaderForPeptideNormalizedIntensityInTMT(channelNumerator,
									sLineHeaderList);
							if (mapValues.containsKey(headerNumerator)) {
								numeratorIntensity = Double.valueOf(mapValues.get(headerNumerator));
							}
							// denominator
							Double denominatorIntensity = null;
							final String headerDenominator = getHeaderForPeptideNormalizedIntensityInTMT(
									channelDenominator, sLineHeaderList);
							if (mapValues.containsKey(headerDenominator)) {
								denominatorIntensity = Double.valueOf(mapValues.get(headerDenominator));
							}
							// build the ratio
							if (numeratorIntensity != null && denominatorIntensity != null) {
								final Double ratioValue = numeratorIntensity / denominatorIntensity;
								final CensusRatio ratio = new CensusRatio(ratioValue, false, conditionsByLabels,
										labelNumerator, labelDenominator, AggregationLevel.PSM,
										labelNumerator + "/" + labelDenominator);
								quantifiedPSM.addRatio(ratio);
							}
						}
					}

				} else {
					// having more than one ratioDescriptor, means that we have
					// L,
					// M, H
					for (final RatioDescriptor ratioDescriptor : ratioDescriptors) {
						final QuantificationLabel labelNumerator = ratioDescriptor.getLabel1();
						final QuantificationLabel labelDenominator = ratioDescriptor.getLabel2();
						final String ratioSuffix = ratioDescriptor.getRatioSuffix();
						if (mapValues.containsKey(RATIO + ratioSuffix)) {
							final double ratioValue = QuantUtils
									.parseCensusRatioValue(mapValues.get(RATIO + ratioSuffix));
							final CensusRatio ratio = new CensusRatio(ratioValue, false, conditionsByLabels,
									labelNumerator, labelDenominator, AggregationLevel.PSM, RATIO);
							quantifiedPSM.addRatio(ratio);
						}
						if (mapValues.containsKey(NORM_RATIO + ratioSuffix)) {
							final double ratioValue = QuantUtils
									.parseCensusRatioValue(mapValues.get(NORM_RATIO + ratioSuffix));
							final CensusRatio ratio = new CensusRatio(ratioValue, false, conditionsByLabels,
									labelNumerator, labelDenominator, AggregationLevel.PSM, NORM_RATIO);
							quantifiedPSM.addRatio(ratio);
						}
						if (mapValues.containsKey(AREA_RATIO + ratioSuffix)) {
							try {
								final double ratioValue = QuantUtils
										.parseCensusRatioValue(mapValues.get(AREA_RATIO + ratioSuffix));
								final CensusRatio ratio = new CensusRatio(ratioValue, false, conditionsByLabels,
										labelNumerator, labelDenominator, AggregationLevel.PSM, AREA_RATIO);
								// add ratio to PSM
								quantifiedPSM.addRatio(ratio);
							} catch (final NumberFormatException e) {
								// skip this
							}
						}
					}
				}
			} else {
				// no ratios
			}
			int plex = 0;
			if (getTmtPlex() != null) {
				plex = getTmtPlex();
				final List<QuantificationLabel> labels = QuantificationLabel.getTMTPlexLabels(plex);
				for (int channel = 1; channel <= plex; channel++) {
					// normalized intensity
					final QuantificationLabel label = labels.get(channel - 1);
					QuantCondition condition = null;
					if (conditionsByLabels != null) {
						condition = conditionsByLabels.get(label);
					}
					String header = getHeaderForPeptideNormalizedIntensityInTMT(channel, sLineHeaderList);
					if (mapValues.containsKey(header)) {
						final Double normalizedIntensity = Double.valueOf(mapValues.get(header));

						final QuantAmount amount = new QuantAmount(normalizedIntensity, AmountType.NORMALIZED_INTENSITY,
								condition, label);
						quantifiedPSM.addAmount(amount);
					}
					// raw intensity
					header = getHeaderForPeptideRawIntensityInTMT(channel, sLineHeaderList);
					if (mapValues.containsKey(header)) {
						final Double rawIntensity = Double.valueOf(mapValues.get(header));
						final QuantAmount amount = new QuantAmount(rawIntensity, AmountType.INTENSITY, condition,
								label);
						quantifiedPSM.addAmount(amount);
					}
				}
			}
			// PSM amounts

			if (conditionsByLabels != null) {
				// SAM_INT
				// light peptide peak area from reconstructed
				// chromatogram
				if (mapValues.containsKey(SAM_INT)) {
					try {
						final double value = Double.valueOf(mapValues.get(SAM_INT));
						final QuantificationLabel lightLabel = conditionsByLabels.keySet().stream()
								.filter(label -> label.isLight()).findAny().get();
						final QuantAmount amount = new QuantAmount(value, AmountType.AREA,
								getLightCondition(conditionsByLabels), lightLabel);
						if (singleton && amount.getValue() != 0.0) {
							amount.setSingleton(true);
						}
						// add amount to PSM
						quantifiedPSM.addAmount(amount);
					} catch (final NumberFormatException e) {
						// skip this
					}
				}
			}
			// get the peak areas of all the conditions (all the labels)
			if (ratioDescriptors != null) {
				final Set<String> usedSuffixes = new HashSet<String>();
				for (final RatioDescriptor ratioDescriptor : ratioDescriptors) {
					final Map<String, QuantCondition> conditionsByIndividualRatioSuffixes = ratioDescriptor
							.getConditionsByIndividualRatioSuffixes();
					final Set<String> differentValuesOfPeakArea = new HashSet<String>();
					for (final String suffix : conditionsByIndividualRatioSuffixes.keySet()) {
						if (usedSuffixes.contains(suffix)) {
							continue;
						}
						usedSuffixes.add(suffix);
						QuantificationLabel label = null;
						switch (suffix) {
						case "_L":
							label = QuantificationLabel.LIGHT;
							break;
						case "_M":
							label = QuantificationLabel.MEDIUM;
							break;
						case "_H":
							label = QuantificationLabel.HEAVY;
							break;
						default:
							break;
						}
						final QuantCondition quantCondition = conditionsByIndividualRatioSuffixes.get(suffix);

						if (mapValues.containsKey(PEAK_AREA + suffix)) {
							try {
								final String stringValue = mapValues.get(PEAK_AREA + suffix);
								if (isSkipNonResolvedPeaks() && differentValuesOfPeakArea.contains(stringValue)) {
									log.warn("PSM '" + quantifiedPSM.getIdentifier()
											+ "' contains not resolved quantitation values (Repeated peak area '"
											+ stringValue + "'). Skipping it...");
									// removing it from local and static maps
									StaticQuantMaps.psmMap.remove(quantifiedPSM);
									localPsmMap.remove(quantifiedPSM.getKey());
									return;
								}
								if (!stringValue.equals("0.0")) {
									differentValuesOfPeakArea.add(stringValue);
								}
								final double value = Double.valueOf(stringValue);
								final QuantAmount amount = new QuantAmount(value, AmountType.AREA, quantCondition,
										label);
								if (singleton && amount.getValue() != 0.0) {
									amount.setSingleton(true);
								}
								// add amount to PSM
								quantifiedPSM.addAmount(amount);
							} catch (final NumberFormatException e) {
								// skip this
							}
						}
					}
				}
			}
			if (conditionsByLabels != null) {
				// REF_INT
				// heavy peptide peak area from reconstructed
				// chromatogram
				if (mapValues.containsKey(REF_INT)) {
					try {
						final double value = Double.valueOf(mapValues.get(REF_INT));
						final QuantAmount amount = new QuantAmount(value, AmountType.AREA,
								getHeavyCondition(conditionsByLabels), QuantificationLabel.HEAVY);
						if (singleton && amount.getValue() != 0.0) {
							amount.setSingleton(true);
						}
						// add amount to PSM
						quantifiedPSM.addAmount(amount);
					} catch (final NumberFormatException e) {
						// skip this
					}
				}
			}
			if (conditionsByLabels != null) {
				// REGRESSION_FACTOR
				// regression score (r)
				if (mapValues.containsKey(AmountType.REGRESSION_FACTOR.name())) {
					try {
						final double value = Double.valueOf(mapValues.get(AmountType.REGRESSION_FACTOR.name()));
						final QuantAmount amount = new QuantAmount(value, AmountType.REGRESSION_FACTOR,
								getLightCondition(conditionsByLabels), QuantificationLabel.LIGHT);
						// add amount to PSM
						quantifiedPSM.addAmount(amount);
					} catch (final NumberFormatException e) {
						// skip this
					}
				}
			}
			// PTM localization score
			if (localizationScore != null) {
				parseLocalizationScore(localizationScore, quantifiedPSM);
			}

			// create the peptide
			QuantifiedPeptideInterface quantifiedPeptide = null;
			final String peptideKey = KeyUtils.getInstance().getSequenceChargeKey(quantifiedPSM,
					isDistinguishModifiedSequences(), isChargeSensible());
			if (StaticQuantMaps.peptideMap.containsKey(peptideKey)) {
				quantifiedPeptide = StaticQuantMaps.peptideMap.getItem(peptideKey);
			} else {
				quantifiedPeptide = new QuantifiedPeptide(quantifiedPSM, isIgnoreTaxonomies(),
						isDistinguishModifiedSequences(), isChargeSensible());
			}
			StaticQuantMaps.peptideMap.addItem(quantifiedPeptide);

			quantifiedPSM.setQuantifiedPeptide(quantifiedPeptide, true);
			// add peptide to map
			if (!localPeptideMap.containsKey(peptideKey)) {
				localPeptideMap.put(peptideKey, quantifiedPeptide);
			}

			if (dbIndex != null) {
				final String cleanSeq = quantifiedPSM.getSequence();
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
				for (final IndexedProtein indexedProtein : indexedProteins) {
					final String proteinKey = QuantKeyUtils.getInstance().getProteinKey(indexedProtein,
							isIgnoreACCFormat());
					QuantifiedProteinInterface newQuantifiedProtein = null;
					if (StaticQuantMaps.proteinMap.containsKey(proteinKey)) {
						newQuantifiedProtein = StaticQuantMaps.proteinMap.getItem(proteinKey);

					} else {
						newQuantifiedProtein = new QuantifiedProteinFromDBIndexEntry(indexedProtein,
								isIgnoreTaxonomies(), isIgnoreACCFormat());

					}
					StaticQuantMaps.proteinMap.addItem(newQuantifiedProtein);
					// add protein to protein map
					if (newQuantifiedProtein.getTaxonomies() != null) {
						taxonomies.addAll(newQuantifiedProtein.getTaxonomies());
					}
					final QuantifiedProteinInterface tmp = localProteinMap.put(proteinKey, newQuantifiedProtein);
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
					addToMap(proteinKey, proteinToPeptidesMap, KeyUtils.getInstance()
							.getSequenceChargeKey(quantifiedPSM, isDistinguishModifiedSequences(), isChargeSensible()));

				}
			}
			// use the already created quantified
			// protein

			for (final QuantifiedProteinInterface quantifiedProtein : quantifiedProteins) {
				// add psm to the proteins
				quantifiedProtein.addPSM(quantifiedPSM, true);
				// add protein to the psm
				quantifiedPSM.addQuantifiedProtein(quantifiedProtein, true);
				// add peptide to the protein
				quantifiedProtein.addPeptide(quantifiedPeptide, true);
				// add to the map (if it was already there
				// is not a problem, it will be only once)
				final String proteinKey = quantifiedProtein.getKey();
				addToMap(proteinKey, proteinToPeptidesMap, KeyUtils.getInstance().getSequenceChargeKey(quantifiedPSM,
						isDistinguishModifiedSequences(), isChargeSensible()));
				// add protein to protein map
				localProteinMap.put(proteinKey, quantifiedProtein);
				// add to protein-experiment map
				addToMap(experimentKey, experimentToProteinsMap, proteinKey);
			}

			// in case of quantifying sites, set the sites to the ratios in case
			// of no ambiguities
			if (!getQuantifiedAAs().isEmpty()) {
				for (final QuantRatio ratio : quantifiedPSM.getQuantRatios()) {
					for (final Character c : getQuantifiedAAs()) {
						if (quantifiedPSM.getSequence().contains(String.valueOf(c))) {
							ratio.setQuantifiedAA(c);
						}
					}
					// check for ambiguity on the quantified site
					int numSites = 0;
					PositionInPeptide quantifiedSitePositionInPeptide = null;
					for (final Character c : getQuantifiedAAs()) {
						final TIntArrayList allPositionsOf = StringUtils.allPositionsOf(quantifiedPSM.getSequence(), c);
						numSites += allPositionsOf.size();
						if (allPositionsOf.size() == 1) {
							quantifiedSitePositionInPeptide = new PositionInPeptide(allPositionsOf.get(0), c,
									quantifiedPSM.getSequence());
						}
						ratio.setQuantifiedAA(c);
					}
					// if no ambiguities
					if (numSites == 1) {
						ratio.addQuantifiedSitePositionInPeptide(quantifiedSitePositionInPeptide);
					}
				}
			}
		} catch (final IllegalArgumentException e) {
			e.printStackTrace();
			log.warn(e);
			log.info("Error reading line '" + line + "' from file. Skipping it...");

		} catch (final NullPointerException e) {
			e.printStackTrace();
			log.warn(e);
			log.warn("Error reading line '" + line + "' from file. Skipping it...");
		} catch (final DBIndexStoreException e) {

			e.printStackTrace();
			log.warn(e);
			log.warn("Error reading line '" + line + "' from file. Skipping it...");
		}

	}

	private void parseLocalizationScore(String localizationScoreString, QuantifiedPSMInterface quantifiedPSM) {
		final List<Float> scores = new ArrayList<Float>();

		if (!"NA".equals(localizationScoreString)) {

			if (localizationScoreString.startsWith("[")) {
				localizationScoreString = localizationScoreString.substring(1);
			}
			if (localizationScoreString.endsWith("]")) {
				localizationScoreString = localizationScoreString.substring(0, localizationScoreString.length() - 1);
			}
			if (localizationScoreString.contains(",")) {
				final String[] split = localizationScoreString.split(",");
				for (final String string : split) {
					if ("".equals(string.trim())) {
						continue;
					}
					scores.add(Float.valueOf(string));
				}
			} else {
				if (!"".equals(localizationScoreString.trim())) {
					scores.add(Float.valueOf(localizationScoreString));
				}
			}
		} else {
			scores.add(Float.NaN);
		}

		// add them to the ptms

		final List<PTMSite> ptms = new ArrayList<PTMSite>();
		quantifiedPSM.getPTMs().stream().forEach(ptm -> ptms.addAll(ptm.getPTMSites()));
		Collections.sort(ptms, new Comparator<PTMSite>() {

			@Override
			public int compare(PTMSite o1, PTMSite o2) {
				return Integer.compare(o1.getPosition(), o2.getPosition());
			}
		});
//		if (scores.size() != ptms.size()) {
//			throw new IllegalArgumentException("psm with sequence " + quantifiedPSM.getFullSequence() + " contains "
//					+ ptms.size() + " PTMs but we have " + scores.size() + " localization scores!");
//		}

		int scoresIndexes = 0;
		for (int i = 0; i < ptms.size(); i++) {

			final Float score = scores.get(scoresIndexes);
			final PTMSiteEx ptm = (PTMSiteEx) ptms.get(i);
			ptm.setScore(
					new ScoreEx(String.valueOf(score), LOCALIZATION_SCORE, LOCALIZATION_SCORE, LOCALIZATION_SCORE));
			if (scoresIndexes + 1 < scores.size()) {
				scoresIndexes++;
			}
		}
	}

	private boolean isTMT4Plex(Collection<QuantificationLabel> labels) {
		for (final QuantificationLabel label : labels) {
			if (QuantificationLabel.isTMT4PLEX(label)) {
				return true;
			}
		}
		return false;
	}

	private boolean isTMT6Plex(Collection<QuantificationLabel> labels) {
		for (final QuantificationLabel label : labels) {
			if (QuantificationLabel.isTMT6PLEX(label)) {
				return true;
			}
		}
		return false;
	}

	private boolean isTMT10Plex(Collection<QuantificationLabel> labels) {
		for (final QuantificationLabel label : labels) {
			if (QuantificationLabel.isTMT10PLEX(label)) {
				return true;
			}
		}

		return false;
	}

	private boolean isTMT11Plex(Collection<QuantificationLabel> labels) {
		for (final QuantificationLabel label : labels) {
			if (QuantificationLabel.isTMT11PLEX(label)) {
				return true;
			}
		}

		return false;
	}

	private boolean isSkipNonResolvedPeaks() {
		return skipNonResolvedPeaks;
	}

	public void setSkipNonResolvedPeaks(boolean skipNonResolvedPeaks) {
		this.skipNonResolvedPeaks = skipNonResolvedPeaks;
	}

	private QuantCondition getLightCondition(Map<QuantificationLabel, QuantCondition> conditionsByLabels) {
		for (final QuantificationLabel label : conditionsByLabels.keySet()) {
			if (label.isLight()) {
				final QuantCondition quantCondition = conditionsByLabels.get(label);
				if (quantCondition == null) {
					log.warn("condition is null");
				}
				return quantCondition;
			}
		}
		return null;
	}

	private QuantCondition getHeavyCondition(Map<QuantificationLabel, QuantCondition> conditionsByLabels) {
		for (final QuantificationLabel label : conditionsByLabels.keySet()) {
			if (label.isHeavy()) {
				return conditionsByLabels.get(label);
			}
		}
		return null;
	}

	private QuantCondition getMediumCondition(Map<QuantificationLabel, QuantCondition> conditionsByLabels) {
		for (final QuantificationLabel label : conditionsByLabels.keySet()) {
			if (label.isMedium()) {
				return conditionsByLabels.get(label);
			}
		}
		return null;
	}

	/**
	 * Gets how the header of the ratio for each channel should start
	 *
	 * @param label
	 * @return
	 */
	private String getHeaderForNormalizedRatioInTMT6Plex(QuantificationLabel label, List<String> headers) {
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

	/**
	 * Gets how the header of the peptide normalized intensity for each channel
	 * should start
	 *
	 * @param label
	 * @return
	 */
	private String getHeaderForPeptideNormalizedIntensityInTMT4Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {

		case TMT_4PLEX_127:
			return findStartingBy(headers, "norm_m/z_127.");
		case TMT_4PLEX_128:
			return findStartingBy(headers, "norm_m/z_128.");
		case TMT_4PLEX_129:
			return findStartingBy(headers, "norm_m/z_129.");
		case TMT_4PLEX_130:
			return findStartingBy(headers, "norm_m/z_130.");

		default:
			return null;
		}
	}

	/**
	 * Gets how the header of the peptide normalized intensity for each channel
	 * should start
	 *
	 * @param label
	 * @return
	 */
	private String getHeaderForPeptideNormalizedIntensityInTMT6Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {
		case TMT_6PLEX_126:
			return findStartingBy(headers, "norm_m/z_126.");
		case TMT_6PLEX_127:
			return findStartingBy(headers, "norm_m/z_127.");
		case TMT_6PLEX_128:
			return findStartingBy(headers, "norm_m/z_128.");
		case TMT_6PLEX_129:
			return findStartingBy(headers, "norm_m/z_129.");
		case TMT_6PLEX_130:
			return findStartingBy(headers, "norm_m/z_130.");
		case TMT_6PLEX_131:
			return findStartingBy(headers, "norm_m/z_131.");
		default:
			return null;
		}
	}

	private String getHeaderForPeptideNormalizedIntensityInTMT(int tmtChannel, List<String> pLineHeaders) {
		final List<String> normalizedIntensityHeaders = pLineHeaders.stream()
				.filter(header -> header.startsWith("norm_m/z_")).collect(Collectors.toList());
		return normalizedIntensityHeaders.get(tmtChannel - 1);
	}

	private String getHeaderForPeptideRawIntensityInTMT(int tmtChannel, List<String> pLineHeaders) {
		final List<String> rawIntensityHeaders = pLineHeaders.stream().filter(header -> header.startsWith("m/z_"))
				.collect(Collectors.toList());
		return rawIntensityHeaders.get(tmtChannel - 1);
	}

	/**
	 * Gets how the header of the peptide raw intensity for each channel should
	 * start
	 *
	 * @param label
	 * @return
	 */
	private String getHeaderForPeptideRawIntensityInTMT4Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {

		case TMT_4PLEX_127:
			return findStartingBy(headers, "m/z_127.");
		case TMT_4PLEX_128:
			return findStartingBy(headers, "m/z_128.");
		case TMT_4PLEX_129:
			return findStartingBy(headers, "m/z_129.");
		case TMT_4PLEX_130:
			return findStartingBy(headers, "m/z_130.");

		default:
			return null;
		}
	}

	/**
	 * Gets how the header of the peptide raw intensity for each channel should
	 * start
	 *
	 * @param label
	 * @return
	 */
	private String getHeaderForPeptideRawIntensityInTMT6Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {
		case TMT_6PLEX_126:
			return findStartingBy(headers, "m/z_126.");
		case TMT_6PLEX_127:
			return findStartingBy(headers, "m/z_127.");
		case TMT_6PLEX_128:
			return findStartingBy(headers, "m/z_128.");
		case TMT_6PLEX_129:
			return findStartingBy(headers, "m/z_129.");
		case TMT_6PLEX_130:
			return findStartingBy(headers, "m/z_130.");
		case TMT_6PLEX_131:
			return findStartingBy(headers, "m/z_131.");
		default:
			return null;
		}
	}

	/**
	 * Gets how the header of the peptide normalized intensity for each channel
	 * should start
	 *
	 * @param label
	 * @return
	 */
	private String getHeaderForPeptideNormalizedIntensityInTMT10Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {
		case TMT_10PLEX_126:
			return findStartingBy(headers, "norm_m/z_126.12");
		case TMT_10PLEX_127N:
			return findStartingBy(headers, "norm_m/z_127.124");
		case TMT_10PLEX_127C:
			return findStartingBy(headers, "norm_m/z_127.131");
		case TMT_10PLEX_128N:
			return findStartingBy(headers, "norm_m/z_128.128");
		case TMT_10PLEX_128C:
			return findStartingBy(headers, "norm_m/z_128.134");
		case TMT_10PLEX_129N:
			return findStartingBy(headers, "norm_m/z_129.131");
		case TMT_10PLEX_129C:
			return findStartingBy(headers, "norm_m/z_129.137");
		case TMT_10PLEX_130N:
			return findStartingBy(headers, "norm_m/z_130.134");
		case TMT_10PLEX_130C:
			return findStartingBy(headers, "norm_m/z_130.141");
		case TMT_10PLEX_131:
			return findStartingBy(headers, "norm_m/z_131.138");
		default:
			return null;
		}
	}

	/**
	 * Gets how the header of the peptide normalized intensity for each channel
	 * should start
	 *
	 * @param label
	 * @return
	 */
	private String getHeaderForPeptideNormalizedIntensityInTMT11Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {
		case TMT_11PLEX_126:
			return findStartingBy(headers, "norm_m/z_126.12");
		case TMT_11PLEX_127N:
			return findStartingBy(headers, "norm_m/z_127.124");
		case TMT_11PLEX_127C:
			return findStartingBy(headers, "norm_m/z_127.131");
		case TMT_11PLEX_128N:
			return findStartingBy(headers, "norm_m/z_128.128");
		case TMT_11PLEX_128C:
			return findStartingBy(headers, "norm_m/z_128.134");
		case TMT_11PLEX_129N:
			return findStartingBy(headers, "norm_m/z_129.131");
		case TMT_11PLEX_129C:
			return findStartingBy(headers, "norm_m/z_129.137");
		case TMT_11PLEX_130N:
			return findStartingBy(headers, "norm_m/z_130.134");
		case TMT_11PLEX_130C:
			return findStartingBy(headers, "norm_m/z_130.141");
		case TMT_11PLEX_131:
			return findStartingBy(headers, "norm_m/z_131.138");
		case TMT_11PLEX_131C:
			return findStartingBy(headers, "norm_m/z_131.144");
		default:
			return null;
		}
	}

	/**
	 * Gets how the header of the peptide normalized intensity for each channel
	 * should start
	 *
	 * @param label
	 * @return
	 */
	private String getHeaderForPeptideNormalizedIntensityInTMT16Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {
		case TMT_16PLEX_126:
			return findStartingBy(headers, "norm_m/z_126.12");
		case TMT_16PLEX_127N:
			return findStartingBy(headers, "norm_m/z_127.124");
		case TMT_16PLEX_127C:
			return findStartingBy(headers, "norm_m/z_127.131");
		case TMT_16PLEX_128N:
			return findStartingBy(headers, "norm_m/z_128.128");
		case TMT_16PLEX_128C:
			return findStartingBy(headers, "norm_m/z_128.134");
		case TMT_16PLEX_129N:
			return findStartingBy(headers, "norm_m/z_129.131");
		case TMT_16PLEX_129C:
			return findStartingBy(headers, "norm_m/z_129.137");
		case TMT_16PLEX_130N:
			return findStartingBy(headers, "norm_m/z_130.134");
		case TMT_16PLEX_130C:
			return findStartingBy(headers, "norm_m/z_130.141");
		case TMT_16PLEX_131N:
			return findStartingBy(headers, "norm_m/z_131.138");
		case TMT_16PLEX_131C:
			return findStartingBy(headers, "norm_m/z_131.144");
		case TMT_16PLEX_132N:
			return findStartingBy(headers, "norm_m/z_132.141");
		case TMT_16PLEX_132C:
			return findStartingBy(headers, "norm_m/z_132.147");
		case TMT_16PLEX_133N:
			return findStartingBy(headers, "norm_m/z_133.144");
		case TMT_16PLEX_133C:
			return findStartingBy(headers, "norm_m/z_133.151");
		case TMT_16PLEX_134:
			return findStartingBy(headers, "norm_m/z_134.14");
		default:
			return null;
		}
	}

	/**
	 * Gets how the header of the peptide raw intensity for each channel should
	 * start
	 *
	 * @param label
	 * @return
	 */
	private String getHeaderForPeptideRawIntensityInTMT16Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {
		case TMT_16PLEX_126:
			return findStartingBy(headers, "m/z_126.127");
		case TMT_16PLEX_127N:
			return findStartingBy(headers, "m/z_127.124");
		case TMT_16PLEX_127C:
			return findStartingBy(headers, "m/z_127.131");
		case TMT_16PLEX_128N:
			return findStartingBy(headers, "m/z_128.128");
		case TMT_16PLEX_128C:
			return findStartingBy(headers, "m/z_128.134");
		case TMT_16PLEX_129N:
			return findStartingBy(headers, "m/z_129.131");
		case TMT_16PLEX_129C:
			return findStartingBy(headers, "m/z_129.137");
		case TMT_16PLEX_130N:
			return findStartingBy(headers, "m/z_130.134");
		case TMT_16PLEX_130C:
			return findStartingBy(headers, "m/z_130.141");
		case TMT_16PLEX_131N:
			return findStartingBy(headers, "m/z_131.138");
		case TMT_16PLEX_131C:
			return findStartingBy(headers, "m/z_131.144");
		case TMT_16PLEX_132N:
			return findStartingBy(headers, "m/z_132.141");
		case TMT_16PLEX_132C:
			return findStartingBy(headers, "m/z_132.147");
		case TMT_16PLEX_133N:
			return findStartingBy(headers, "m/z_133.144");
		case TMT_16PLEX_133C:
			return findStartingBy(headers, "m/z_133.151");
		case TMT_16PLEX_134:
			return findStartingBy(headers, "m/z_134.14");
		default:
			return null;
		}
	}

	private String findStartingBy(List<String> headers, String toFind) {
		for (final String header : headers) {
			if (header.toLowerCase().startsWith(toFind.toLowerCase())) {
				return header;
			}
		}
		throw new IllegalArgumentException(toFind + " is not found in the headers of the line");
	}

	/**
	 * Gets how the header of the peptide raw intensity for each channel should
	 * start
	 *
	 * @param label
	 * @return
	 */
	private String getHeaderForPeptideRawIntensityInTMT10Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {
		case TMT_10PLEX_126:
			return findStartingBy(headers, "m/z_126.127");
		case TMT_10PLEX_127N:
			return findStartingBy(headers, "m/z_127.124");
		case TMT_10PLEX_127C:
			return findStartingBy(headers, "m/z_127.131");
		case TMT_10PLEX_128N:
			return findStartingBy(headers, "m/z_128.128");
		case TMT_10PLEX_128C:
			return findStartingBy(headers, "m/z_128.134");
		case TMT_10PLEX_129N:
			return findStartingBy(headers, "m/z_129.131");
		case TMT_10PLEX_129C:
			return findStartingBy(headers, "m/z_129.137");
		case TMT_10PLEX_130N:
			return findStartingBy(headers, "m/z_130.134");
		case TMT_10PLEX_130C:
			return findStartingBy(headers, "m/z_130.141");
		case TMT_10PLEX_131:
			return findStartingBy(headers, "m/z_131.138");
		default:
			return null;
		}
	}

	/**
	 * Gets how the header of the peptide raw intensity for each channel should
	 * start
	 *
	 * @param label
	 * @return
	 */
	private String getHeaderForPeptideRawIntensityInTMT11Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {
		case TMT_11PLEX_126:
			return findStartingBy(headers, "m/z_126.127");
		case TMT_11PLEX_127N:
			return findStartingBy(headers, "m/z_127.124");
		case TMT_11PLEX_127C:
			return findStartingBy(headers, "m/z_127.131");
		case TMT_11PLEX_128N:
			return findStartingBy(headers, "m/z_128.128");
		case TMT_11PLEX_128C:
			return findStartingBy(headers, "m/z_128.134");
		case TMT_11PLEX_129N:
			return findStartingBy(headers, "m/z_129.131");
		case TMT_11PLEX_129C:
			return findStartingBy(headers, "m/z_129.137");
		case TMT_11PLEX_130N:
			return findStartingBy(headers, "m/z_130.134");
		case TMT_11PLEX_130C:
			return findStartingBy(headers, "m/z_130.141");
		case TMT_11PLEX_131:
			return findStartingBy(headers, "m/z_131.138");
		case TMT_11PLEX_131C:
			return findStartingBy(headers, "m/z_131.144");
		default:
			return null;
		}
	}

	/**
	 * Gets how the header of the protein normalized intensity for each channel
	 * should start
	 *
	 * @param label
	 * @param mapValues
	 * @return
	 */
	private String getHeaderForProteinNormalizedIntensityInTMT4Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {

		case TMT_4PLEX_127:
			return findStartingBy(headers, "norm_total m/z_127.");
		case TMT_4PLEX_128:
			return findStartingBy(headers, "norm_total m/z_128.");
		case TMT_4PLEX_129:
			return findStartingBy(headers, "norm_total m/z_129.");
		case TMT_4PLEX_130:
			return findStartingBy(headers, "norm_total m/z_130.");

		default:
			return null;
		}
	}

	/**
	 * Gets how the header of the protein normalized intensity for each channel
	 * should start
	 *
	 * @param label
	 * @param mapValues
	 * @return
	 */
	private String getHeaderForProteinNormalizedIntensityInTMT6Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {
		case TMT_6PLEX_126:
			return findStartingBy(headers, "norm_total m/z_126.");
		case TMT_6PLEX_127:
			return findStartingBy(headers, "norm_total m/z_127.");
		case TMT_6PLEX_128:
			return findStartingBy(headers, "norm_total m/z_128.");
		case TMT_6PLEX_129:
			return findStartingBy(headers, "norm_total m/z_129.");
		case TMT_6PLEX_130:
			return findStartingBy(headers, "norm_total m/z_130.");
		case TMT_6PLEX_131:
			return findStartingBy(headers, "norm_total m/z_131.");
		default:
			return null;
		}
	}

	/**
	 * Gets how the header of the protein raw intensity for each channel should
	 * start
	 *
	 * @param label
	 * @param mapValues
	 * @return
	 */
	private String getHeaderForProteinRawIntensityInTMT4Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {

		case TMT_4PLEX_127:
			return findStartingBy(headers, "total m/z_127.");
		case TMT_4PLEX_128:
			return findStartingBy(headers, "total m/z_128.");
		case TMT_4PLEX_129:
			return findStartingBy(headers, "total m/z_129.");
		case TMT_4PLEX_130:
			return findStartingBy(headers, "total m/z_130.");

		default:
			return null;
		}
	}

	/**
	 * Gets how the header of the protein raw intensity for each channel should
	 * start
	 *
	 * @param label
	 * @param mapValues
	 * @return
	 */
	private String getHeaderForProteinRawIntensityInTMT6Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {
		case TMT_6PLEX_126:
			return findStartingBy(headers, "total m/z_126.");
		case TMT_6PLEX_127:
			return findStartingBy(headers, "total m/z_127.");
		case TMT_6PLEX_128:
			return findStartingBy(headers, "total m/z_128.");
		case TMT_6PLEX_129:
			return findStartingBy(headers, "total m/z_129.");
		case TMT_6PLEX_130:
			return findStartingBy(headers, "total m/z_130.");
		case TMT_6PLEX_131:
			return findStartingBy(headers, "total m/z_131.");
		default:
			return null;
		}
	}

	/**
	 * Gets how the header of the normalized intensity for each channel should start
	 *
	 * @param label
	 * @return
	 */
	private String getHeaderForProteinNormalizedIntensityInTMT10Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {
		case TMT_10PLEX_126:
			return findStartingBy(headers, "norm_total m/z_126.127");
		case TMT_10PLEX_127N:
			return findStartingBy(headers, "norm_total m/z_127.124");
		case TMT_10PLEX_127C:
			return findStartingBy(headers, "norm_total m/z_127.131");
		case TMT_10PLEX_128N:
			return findStartingBy(headers, "norm_total m/z_128.128");
		case TMT_10PLEX_128C:
			return findStartingBy(headers, "norm_total m/z_128.134");
		case TMT_10PLEX_129N:
			return findStartingBy(headers, "norm_total m/z_129.131");
		case TMT_10PLEX_129C:
			return findStartingBy(headers, "norm_total m/z_129.137");
		case TMT_10PLEX_130N:
			return findStartingBy(headers, "norm_total m/z_130.134");
		case TMT_10PLEX_130C:
			return findStartingBy(headers, "norm_total m/z_130.141");
		case TMT_10PLEX_131:
			return findStartingBy(headers, "norm_total m/z_131.138");
		default:
			return null;
		}
	}

	/**
	 * Gets how the header of the protein raw intensity for each channel should
	 * start
	 *
	 * @param label
	 * @return
	 */
	private String getHeaderForProteinRawIntensityInTMT10Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {
		case TMT_10PLEX_126:
			return findStartingBy(headers, "total m/z_126.127");
		case TMT_10PLEX_127N:
			return findStartingBy(headers, "total m/z_127.124");
		case TMT_10PLEX_127C:
			return findStartingBy(headers, "total m/z_127.131");
		case TMT_10PLEX_128N:
			return findStartingBy(headers, "total m/z_128.128");
		case TMT_10PLEX_128C:
			return findStartingBy(headers, "total m/z_128.134");
		case TMT_10PLEX_129N:
			return findStartingBy(headers, "total m/z_129.131");
		case TMT_10PLEX_129C:
			return findStartingBy(headers, "total m/z_129.137");
		case TMT_10PLEX_130N:
			return findStartingBy(headers, "total m/z_130.134");
		case TMT_10PLEX_130C:
			return findStartingBy(headers, "total m/z_130.141");
		case TMT_10PLEX_131:
			return findStartingBy(headers, "total m/z_131.138");
		default:
			return null;
		}
	}

	/**
	 * Gets how the header of the protein raw intensity for each channel should
	 * start
	 *
	 * @param label
	 * @return
	 */
	private String getHeaderForProteinRawIntensityInTMT11Plex(QuantificationLabel label, List<String> headers) {
		switch (label) {
		case TMT_11PLEX_126:
			return findStartingBy(headers, "total m/z_126.127");
		case TMT_11PLEX_127N:
			return findStartingBy(headers, "total m/z_127.124");
		case TMT_11PLEX_127C:
			return findStartingBy(headers, "total m/z_127.131");
		case TMT_11PLEX_128N:
			return findStartingBy(headers, "total m/z_128.128");
		case TMT_11PLEX_128C:
			return findStartingBy(headers, "total m/z_128.134");
		case TMT_11PLEX_129N:
			return findStartingBy(headers, "total m/z_129.131");
		case TMT_11PLEX_129C:
			return findStartingBy(headers, "total m/z_129.137");
		case TMT_11PLEX_130N:
			return findStartingBy(headers, "total m/z_130.134");
		case TMT_11PLEX_130C:
			return findStartingBy(headers, "total m/z_130.141");
		case TMT_11PLEX_131:
			return findStartingBy(headers, "total m/z_131.138");
		case TMT_11PLEX_131C:
			return findStartingBy(headers, "total m/z_131.144");
		default:
			return null;
		}
	}

	private QuantifiedProteinInterface processProteinLine(String line, List<String> pLineHeaderList,
			Map<QuantificationLabel, QuantCondition> conditionsByLabels, List<RatioDescriptor> ratioDescriptors,
			String experimentKey) throws DiscardProteinException {
		// new protein
		final MyHashMap<String, String> mapValues = getMapFromPLine(pLineHeaderList, line);
		final String proteinACC = mapValues.get(LOCUS);

		QuantifiedProteinInterface quantifiedProtein = null;
		// apply the pattern if available
		if (decoyPattern != null) {
			final Matcher matcher = decoyPattern.matcher(proteinACC);
			if (matcher.find()) {
				quantifiedProtein = null;
				throw new DiscardProteinException("Protein " + proteinACC + " is DECOY.");
			}
		}

		if (StaticQuantMaps.proteinMap.containsKey(proteinACC)) {
			quantifiedProtein = StaticQuantMaps.proteinMap.getItem(proteinACC);
		} else {
			quantifiedProtein = new QuantifiedProtein(proteinACC, proteinACC, isIgnoreTaxonomies(), true);
			final String description = mapValues.get(DESCRIPTION);
			quantifiedProtein.setDescription(description);
		}
		StaticQuantMaps.proteinMap.addItem(quantifiedProtein);

		final QuantifiedProteinInterface tmp = localProteinMap.put(proteinACC, quantifiedProtein);
		// add to protein-experiment map
		addToMap(experimentKey, experimentToProteinsMap, proteinACC);

		// if we have ratios
		if (ratioDescriptors != null && !ratioDescriptors.isEmpty()) {
			// if we only have one ratio descriptor, it is L over H
			if (ratioDescriptors.size() == 1) {
				final QuantificationLabel labelNumerator = ratioDescriptors.get(0).getLabel1();
				final QuantificationLabel labelDenominator = ratioDescriptors.get(0).getLabel2();
				final String ratioSuffix = ratioDescriptors.get(0).getRatioSuffix();
				// add protein ratio
				// first look if the composite ratio is calculated
				boolean hasCompositeRatio = false;
				if (mapValues.containsKey(COMPOSITE_RATIO)) {
					try {
						hasCompositeRatio = true;
						final double ratioValue = Double.valueOf(mapValues.get(COMPOSITE_RATIO));
						String stdValue = null;
						if (mapValues.containsKey(COMPOSITE_RATIO_STANDARD_DEVIATION)) {
							stdValue = mapValues.get(COMPOSITE_RATIO_STANDARD_DEVIATION);
							if ("".equals(stdValue)) {
								stdValue = "0.0";
							}
						}
						final QuantRatio ratio = new CensusRatio(ratioValue, stdValue, false, conditionsByLabels,
								labelNumerator, labelDenominator, AggregationLevel.PROTEIN, COMPOSITE_RATIO);
						quantifiedProtein.addRatio(ratio);
					} catch (final NumberFormatException e) {
						// skip this
					}
					// if there is not composite ratio, use the
					// regular ratio
				}
				boolean hasCompositeRatioSuffix = false;
				if (mapValues.containsKey(COMPOSITE_RATIO + ratioSuffix)
						&& (isCapturingRatioName(COMPOSITE_RATIO) || !hasCompositeRatio)) {
					try {
						hasCompositeRatioSuffix = true;
						final double ratioValue = Double.valueOf(mapValues.get(COMPOSITE_RATIO + ratioSuffix));
						String stdValue = null;
						if (mapValues.containsKey(COMPOSITE_RATIO_STANDARD_DEVIATION + ratioSuffix)) {
							stdValue = mapValues.get(COMPOSITE_RATIO_STANDARD_DEVIATION + ratioSuffix);
							if ("".equals(stdValue)) {
								stdValue = "0.0";
							}
						}
						final QuantRatio ratio = new CensusRatio(ratioValue, stdValue, false, conditionsByLabels,
								labelNumerator, labelDenominator, AggregationLevel.PROTEIN, COMPOSITE_RATIO);
						quantifiedProtein.addRatio(ratio);
					} catch (final NumberFormatException e) {
						// skip this
					}
				}
				// if there is not composite ratio, use the
				// regular ratio
				boolean hasOneOfTheseRatios = false;
				if (mapValues.containsKey(AVERAGE_RATIO) || mapValues.containsKey(AREA_RATIO)) {
					if (isCapturingRatioName(AVERAGE_RATIO) || isCapturingRatioName(AREA_RATIO)
							|| (!hasCompositeRatioSuffix && !hasCompositeRatio)) {
						hasOneOfTheseRatios = true;
						try {
							if (mapValues.containsKey(AVERAGE_RATIO)
									&& ((!hasCompositeRatioSuffix && !hasCompositeRatio)
											|| isCapturingRatioName(AVERAGE_RATIO))) {
								final double ratioValue = Double.valueOf(mapValues.get(AVERAGE_RATIO));
								String stdValue = null;
								if (mapValues.containsKey(STANDARD_DEVIATION)) {
									stdValue = mapValues.get(STANDARD_DEVIATION);
									if ("".equals(stdValue)) {
										stdValue = "0.0";
									}
								}
								final QuantRatio ratio = new CensusRatio(ratioValue, stdValue, false,
										conditionsByLabels, labelNumerator, labelDenominator, AggregationLevel.PROTEIN,
										AVERAGE_RATIO);
								quantifiedProtein.addRatio(ratio);
							}
						} catch (final NumberFormatException e) {
							// skip this
						}
						try {
							if (mapValues.containsKey(AREA_RATIO) && ((!hasCompositeRatioSuffix && !hasCompositeRatio)
									|| isCapturingRatioName(AREA_RATIO))) {
								final double ratioValue = Double.valueOf(mapValues.get(AREA_RATIO));
								final QuantRatio ratio = new CensusRatio(ratioValue, null, false, conditionsByLabels,
										labelNumerator, labelDenominator, AggregationLevel.PROTEIN, AREA_RATIO);
								quantifiedProtein.addRatio(ratio);
							}
						} catch (final NumberFormatException e) {
							// skip this
						}
					}
				}

				if (mapValues.containsKey(AVERAGE_RATIO + ratioSuffix)
						|| mapValues.containsKey(AREA_RATIO + ratioSuffix)) {
					if (isCapturingRatioName(AVERAGE_RATIO) || isCapturingRatioName(AREA_RATIO)
							|| (!hasOneOfTheseRatios && !hasCompositeRatio && !hasCompositeRatioSuffix)) {
						try {
							if (mapValues.containsKey(AVERAGE_RATIO + ratioSuffix)
									&& ((!hasOneOfTheseRatios && !hasCompositeRatio && !hasCompositeRatioSuffix)
											|| isCapturingRatioName(AVERAGE_RATIO))) {
								final double ratioValue = Double.valueOf(mapValues.get(AVERAGE_RATIO + ratioSuffix));
								String stdValue = null;
								if (mapValues.containsKey(STANDARD_DEVIATION)) {
									stdValue = mapValues.get(STANDARD_DEVIATION);
									if ("".equals(stdValue)) {
										stdValue = "0.0";
									}
								}
								final QuantRatio ratio = new CensusRatio(ratioValue, stdValue, false,
										conditionsByLabels, labelNumerator, labelDenominator, AggregationLevel.PROTEIN,
										AVERAGE_RATIO);
								quantifiedProtein.addRatio(ratio);
							}
						} catch (final NumberFormatException e) {
							// skip this
						}
						try {
							if (mapValues.containsKey(AREA_RATIO + ratioSuffix)
									&& ((!hasOneOfTheseRatios && !hasCompositeRatio && !hasCompositeRatioSuffix)
											|| isCapturingRatioName(AREA_RATIO))) {
								final double ratioValue = Double.valueOf(mapValues.get(AREA_RATIO + ratioSuffix));
								final QuantRatio ratio = new CensusRatio(ratioValue, null, false, conditionsByLabels,
										labelNumerator, labelDenominator, AggregationLevel.PROTEIN, AREA_RATIO);
								quantifiedProtein.addRatio(ratio);
							}
						} catch (final NumberFormatException e) {
							// skip this
						}
					}
				}
			} else {
				// having more than one ratio descriptor, it is L, M, H
				for (final RatioDescriptor ratioDescriptor : ratioDescriptors) {
					final QuantificationLabel labelNumerator = ratioDescriptor.getLabel1();
					final QuantificationLabel labelDenominator = ratioDescriptor.getLabel2();
					final String ratioSuffix = ratioDescriptor.getRatioSuffix();
					// add protein ratio
					// first look if the normalized composite ratio is
					// calculated
					boolean hasCompositeRatio = false;
					if (mapValues.containsKey(NORM_COMPOSITE_RATIO + ratioSuffix)
							|| mapValues.containsKey(COMPOSITE_RATIO + ratioSuffix)) {
						hasCompositeRatio = true;
						if (mapValues.containsKey(NORM_COMPOSITE_RATIO + ratioSuffix)) {
							try {
								final double ratioValue = Double
										.valueOf(mapValues.get(NORM_COMPOSITE_RATIO + ratioSuffix));
								String stdValue = null;
								if (mapValues.containsKey(NORM_COMPOSITE_RATIO_STDEV + ratioSuffix)) {
									stdValue = mapValues.get(NORM_COMPOSITE_RATIO_STDEV + ratioSuffix);
									if ("".equals(stdValue)) {
										stdValue = "0.0";
									}
								}
								final QuantRatio ratio = new CensusRatio(ratioValue, stdValue, false,
										conditionsByLabels, labelNumerator, labelDenominator, AggregationLevel.PROTEIN,
										NORM_COMPOSITE_RATIO);
								quantifiedProtein.addRatio(ratio);
							} catch (final NumberFormatException e) {
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
								final QuantRatio ratio = new CensusRatio(ratioValue, stdValue, false,
										conditionsByLabels, labelNumerator, labelDenominator, AggregationLevel.PROTEIN,
										COMPOSITE_RATIO);
								quantifiedProtein.addRatio(ratio);
							} catch (final NumberFormatException e) {
								// skip this
							}
							// if there is not composite ratio, use the
							// regular composite ratio
						}
					}

					if (mapValues.containsKey(MEDIAN_NORM_RATIO + ratioSuffix)
							|| mapValues.containsKey(MEDIAN_AREA_RATIO + ratioSuffix)) {
						if (isCapturingRatioName(MEDIAN_NORM_RATIO) || isCapturingRatioName(MEDIAN_NORM_RATIO)
								|| !hasCompositeRatio) {
							try {
								if (mapValues.containsKey(MEDIAN_NORM_RATIO + ratioSuffix)
										&& (!hasCompositeRatio || isCapturingRatioName(MEDIAN_NORM_RATIO))) {
									final double ratioValue = Double
											.valueOf(mapValues.get(MEDIAN_NORM_RATIO + ratioSuffix));
									String stdValue = null;
									if (mapValues.containsKey(NORM_STDEV + ratioSuffix)) {
										stdValue = mapValues.get(NORM_STDEV + ratioSuffix);
										if ("".equals(stdValue)) {
											stdValue = "0.0";
										}
									}
									final QuantRatio ratio = new CensusRatio(ratioValue, stdValue, false,
											conditionsByLabels, labelNumerator, labelDenominator,
											AggregationLevel.PROTEIN, NORM_STDEV);
									quantifiedProtein.addRatio(ratio);
								}
							} catch (final NumberFormatException e) {
								// skip this
							}
							try {
								if (mapValues.containsKey(MEDIAN_AREA_RATIO + ratioSuffix)
										&& (!hasCompositeRatio || isCapturingRatioName(MEDIAN_AREA_RATIO))) {
									final double ratioValue = Double
											.valueOf(mapValues.get(MEDIAN_AREA_RATIO + ratioSuffix));
									final QuantRatio ratio = new CensusRatio(ratioValue, null, false,
											conditionsByLabels, labelNumerator, labelDenominator,
											AggregationLevel.PROTEIN, MEDIAN_AREA_RATIO);
									quantifiedProtein.addRatio(ratio);
								}
							} catch (final NumberFormatException e) {
								// skip this
							}
						}
					}
				}
			}
		} else {
			// no ratios
		}
		// TMT4PLEX
		boolean isTMT4Plex = false;
		if (conditionsByLabels != null && !conditionsByLabels.isEmpty()) {
			isTMT4Plex = isTMT4Plex(conditionsByLabels.keySet());
		}
		if (isTMT4Plex) {
			for (final QuantificationLabel label : QuantificationLabel.getTMT4PlexLabels()) {

				String header = getHeaderForProteinNormalizedIntensityInTMT4Plex(label, pLineHeaderList);
				if (mapValues.containsKey(header)) {
					final Double normalizedIntensity = Double.valueOf(mapValues.get(header));
					final QuantAmount amount = new QuantAmount(normalizedIntensity, AmountType.NORMALIZED_INTENSITY,
							conditionsByLabels.get(label), label);
					quantifiedProtein.addAmount(amount);
				}
				header = getHeaderForProteinRawIntensityInTMT4Plex(label, pLineHeaderList);
				if (mapValues.containsKey(header)) {
					final Double rawIntensity = Double.valueOf(mapValues.get(header));
					final QuantAmount amount = new QuantAmount(rawIntensity, AmountType.INTENSITY,
							conditionsByLabels.get(label), label);
					quantifiedProtein.addAmount(amount);
				}
			}
		}
		// TMT6PLEX
		boolean isTMT6Plex = false;
		if (conditionsByLabels != null && !conditionsByLabels.isEmpty()) {
			isTMT6Plex = isTMT6Plex(conditionsByLabels.keySet());
		}
		if (isTMT6Plex) {
			for (final QuantificationLabel label : QuantificationLabel.getTMT6PlexLabels()) {

				String header = getHeaderForProteinNormalizedIntensityInTMT6Plex(label, pLineHeaderList);
				if (mapValues.containsKey(header)) {
					final Double normalizedIntensity = Double.valueOf(mapValues.get(header));
					final QuantAmount amount = new QuantAmount(normalizedIntensity, AmountType.NORMALIZED_INTENSITY,
							conditionsByLabels.get(label), label);
					quantifiedProtein.addAmount(amount);
				}
				header = getHeaderForProteinRawIntensityInTMT6Plex(label, pLineHeaderList);
				if (mapValues.containsKey(header)) {
					final Double rawIntensity = Double.valueOf(mapValues.get(header));
					final QuantAmount amount = new QuantAmount(rawIntensity, AmountType.INTENSITY,
							conditionsByLabels.get(label), label);
					quantifiedProtein.addAmount(amount);
				}
			}
		}
		// TMT10PLEX
		boolean isTMT10Plex = false;
		if (conditionsByLabels != null && !conditionsByLabels.isEmpty()) {
			isTMT10Plex = isTMT10Plex(conditionsByLabels.keySet());
		}
		if (isTMT10Plex) {
			for (final QuantificationLabel label : QuantificationLabel.getTMT10PlexLabels()) {
				String header = getHeaderForProteinNormalizedIntensityInTMT10Plex(label, pLineHeaderList);
				if (mapValues.containsKey(header)) {
					final Double normalizedIntensity = Double.valueOf(mapValues.get(header));
					final QuantAmount amount = new QuantAmount(normalizedIntensity, AmountType.NORMALIZED_INTENSITY,
							conditionsByLabels.get(label), label);
					quantifiedProtein.addAmount(amount);
				}
				header = getHeaderForProteinRawIntensityInTMT10Plex(label, pLineHeaderList);
				if (mapValues.containsKey(header)) {
					final Double rawIntensity = Double.valueOf(mapValues.get(header));
					final QuantAmount amount = new QuantAmount(rawIntensity, AmountType.INTENSITY,
							conditionsByLabels.get(label), label);
					quantifiedProtein.addAmount(amount);
				}
			}
		}
		// TMT11PLEX
		boolean isTMT11Plex = false;
		if (conditionsByLabels != null && !conditionsByLabels.isEmpty()) {
			isTMT11Plex = isTMT11Plex(conditionsByLabels.keySet());
		}
		if (isTMT11Plex) {
			for (final QuantificationLabel label : QuantificationLabel.getTMT11PlexLabels()) {
				String header = getHeaderForProteinNormalizedIntensityInTMT10Plex(label, pLineHeaderList);
				if (mapValues.containsKey(header)) {
					final Double normalizedIntensity = Double.valueOf(mapValues.get(header));
					final QuantAmount amount = new QuantAmount(normalizedIntensity, AmountType.NORMALIZED_INTENSITY,
							conditionsByLabels.get(label), label);
					quantifiedProtein.addAmount(amount);
				}
				header = getHeaderForProteinRawIntensityInTMT11Plex(label, pLineHeaderList);
				if (mapValues.containsKey(header)) {
					final Double rawIntensity = Double.valueOf(mapValues.get(header));
					final QuantAmount amount = new QuantAmount(rawIntensity, AmountType.INTENSITY,
							conditionsByLabels.get(label), label);
					quantifiedProtein.addAmount(amount);
				}
			}
		}

		return quantifiedProtein;

	}

	private MyHashMap<String, String> getMapFromPLine(List<String> pLineHeaderList, String line) {
		final MyHashMap<String, String> map = new MyHashMap<String, String>();
		String[] split = line.split("\t");
		if (split.length != pLineHeaderList.size()) {
			// try to see if there is one that is empty
			final String[] splitTMP = new String[split.length - 1];
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
			for (final String header : pLineHeaderList) {
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
		final MyHashMap<String, String> map = new MyHashMap<String, String>();
		String[] split = line.split("\t");
		// not remove the last element:
		final int upToThisIndex = split.length - 1;
		split = removeElements(split, "N/A", upToThisIndex);

		// trying to recover some lines in which the peptide information is
		// shifted to the right one column from the sequence
		if ("".equals(split[1]) && "".equals(split[2]) && split.length == sLineHeaderList.size() + 1) {
			// copy into a new array
			final String[] newSplit = new String[sLineHeaderList.size()];
			for (int index = 0; index < newSplit.length; index++) {
				if (index <= 1) {
					newSplit[index] = split[index];
				} else {
					newSplit[index] = split[index + 1];
				}
			}
			split = newSplit;
		}
		if (split.length == sLineHeaderList.size() || split.length + 1 == sLineHeaderList.size()
				|| split.length > sLineHeaderList.size()) {
			int i = 0;
			for (final String header : sLineHeaderList) {
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
		final List<String> list = new ArrayList<String>();
		int index = -1;
		for (final String splitElement : split) {
			index++;
			if (splitElement.equals(elementToRemove) && index < upToThisIndex) {
				continue;
			}
			list.add(splitElement);

		}
		return list.toArray(new String[0]);
	}

	/**
	 * @param onlyOneSpectrumPerChromatographicPeakAndPerSaltStep the
	 *                                                            onlyOneSpectrumPerChromatographicPeakAndPerSaltStep
	 *                                                            to set
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
	 * @param skipSingletons the skipSingletons to set
	 */
	public void setSkipSingletons(boolean skipSingletons) {
		this.skipSingletons = skipSingletons;
	}

	@Override
	public boolean canRead() {
		try {

			for (final RemoteSSHFileReference remoteFileRetriever : this.remoteFileRetrievers) {
				final File file = remoteFileRetriever.getOutputFile();

				List<String> lines = null;

				// check whether it is an excel file
				if (FileUtils.isExcelFile(file)) {
					lines = FileUtils.readLinesFromXLSX(file, "\t", 0);
				} else {
					lines = FileUtils.readFirstLines(file, 1);
				}

				final String line = lines.get(0);

				final String[] split = line.split("\t");

				if (!line.startsWith(H) || !split[1].startsWith("Census version")) {
					return false;
				}

			}

		} catch (final Exception e) {
			return false;
		}
		return true;
	}
}
