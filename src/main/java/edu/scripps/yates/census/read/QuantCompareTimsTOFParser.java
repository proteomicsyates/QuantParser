package edu.scripps.yates.census.read;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.QuantifiedPSM;
import edu.scripps.yates.census.read.model.QuantifiedPeptide;
import edu.scripps.yates.census.read.model.QuantifiedProtein;
import edu.scripps.yates.census.read.model.StaticQuantMaps;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.proteomicsmodel.AbstractPSM;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.proteomicsmodel.Score;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AmountType;
import edu.scripps.yates.utilities.proteomicsmodel.factories.AmountEx;
import edu.scripps.yates.utilities.proteomicsmodel.factories.ScoreEx;
import edu.scripps.yates.utilities.proteomicsmodel.utils.KeyUtils;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.THashSet;
import gnu.trove.set.hash.TIntHashSet;

public class QuantCompareTimsTOFParser extends QuantCompareParser {
	private final static Logger log = Logger.getLogger(QuantCompareTimsTOFParser.class);

	private static final String CCS = "CCS";
	private static final String XIC = "XIC";
	private static final String ESTIMATED_XIC = "ESTIMATED_XIC";
	private static final String SCAN_NUM = "SCAN_NUM";
	private static final String CHARGE_STATE = "CHARGE_STATE";
	private static final Object LOCUS = "LOCUS";
	private static final Object DESCRIPTION = "DESCRIPTION";

	/**
	 * 
	 * @param psmLevelCensusQuantCompareFile input file that is a PSM level census
	 *                                       quant compare file
	 */
	public QuantCompareTimsTOFParser(File psmLevelCensusQuantCompareFile) throws FileNotFoundException {
		super(psmLevelCensusQuantCompareFile);

	}

	@Override
	protected void process() throws QuantParserException {
		try {
			List<String> lines = null;
			final Map<String, Set<String>> peptideToSpectraMap = new THashMap<String, Set<String>>();
			// check whether it is an excel file
			boolean excelFile;
			if (isExcel.containsKey(file.getAbsolutePath())) {
				excelFile = isExcel.get(file.getAbsolutePath());
			} else {
				excelFile = FileUtils.isExcelFile(file);
				isExcel.put(file.getAbsolutePath(), excelFile);
			}
			if (excelFile) {
				lines = FileUtils.readLinesFromXLSX(file, "\t", 0);
			} else {
				lines = Files.readAllLines(this.file.toPath());
			}
			final Set<String> uniqueLineStrings = new THashSet<String>();
			int numLine = 0;
			for (final String line : lines) {
				numLine++;
				final String[] split = line.split("\t");

				if (numLine == 1) {
					// for now, we need to fix it because there are columns missing
					processColumns(split);
				} else {

					// first of all, check whether the line has been already seen before except for
					// the protein columns
					final String uniqueLineString = getUniqueLineStringWithNoProtein(split);
					// take the sequence
					final String rawSequence = split[indexByColumn.get(SEQUENCE)];
					int chargeState = -1;
					if (indexByColumn.containsKey(CHARGE_STATE)) {
						final String chargeStateString = split[indexByColumn.get(CHARGE_STATE)];
						try {
							chargeState = Double.valueOf(chargeStateString).intValue();
						} catch (final NumberFormatException e) {

						}
					}

					QuantifiedPeptideInterface quantPeptide = null;
					int repWithPeptidePresent = 0;
					for (int rep = 1; rep <= columnsByExperiments.size(); rep++) {

						// create a PSM per experiment

						// SCAN
						int scanNumber = -1;
						if (getIndexByColumnAndExperiment(rep, SCAN_NUM) != -1) {
							final String scanNumberString = split[getIndexByColumnAndExperiment(rep, SCAN_NUM)];
							// if scanNumberString is NA, this psm has not been detected in this replicate,
							// and therefore we dont create it
							if ("NA".equals(scanNumberString)) {
								continue;
							}
							repWithPeptidePresent++;

							try {
								scanNumber = Double.valueOf(scanNumberString).intValue();
							} catch (final NumberFormatException e) {

							}
						}

						// in principle, there is no redundancy column for now:
						int redundancy = 1;
						if (getIndexByColumnAndExperiment(rep, REDUNDANCY) != -1) {
							final String redundancyString = split[getIndexByColumnAndExperiment(rep, REDUNDANCY)];
							try {
								redundancy = Double.valueOf(redundancyString).intValue();
							} catch (final NumberFormatException e) {

							}
						}
						String rawFileName = split[getIndexByColumnAndExperiment(rep, FILENAME)];
						if ("NA".equals(rawFileName)) {
							rawFileName = "NA_" + rep;
						}

						final boolean singleton = false;// not used here
						// create a PSM per peptide
						final List<QuantifiedPSMInterface> quantPSMs = new ArrayList<QuantifiedPSMInterface>();

						// the scan number we get is from one PSM only
						// we will create fake scan numbers for all PSMs that are behind this peptide by
						// adding up 1 to each scan number
						for (int i = 1; i <= redundancy; i++) {
							String scanNumberStringForPSM = String.valueOf(scanNumber);
							if (i != 1) {
								scanNumberStringForPSM += "_" + i;
							}
							QuantifiedPSMInterface quantPSM = new QuantifiedPSM(rawSequence, peptideToSpectraMap,
									scanNumberStringForPSM, chargeState, rawFileName, singleton,
									isDistinguishModifiedSequences(), isChargeSensible());
							// make this PSM to be from replicate rep
							quantPSM.addCondition(conditionByExp.get(rep));

							quantPSMs.add(quantPSM);
							final String key = quantPSM.getKey();
							if (!localPsmMap.containsKey(key)) {
								localPsmMap.put(key, quantPSM);
							}
							if (StaticQuantMaps.psmMap.containsKey(key)) {
								quantPSM = StaticQuantMaps.psmMap.getItem(key);
							}
							StaticQuantMaps.psmMap.addItem(quantPSM);
						}
						if (quantPeptide == null) {
							final String peptideKey = KeyUtils.getInstance().getSequenceChargeKey(quantPSMs.get(0),
									isDistinguishModifiedSequences(), isChargeSensible());
							if (StaticQuantMaps.peptideMap.containsKey(peptideKey)) {
								quantPeptide = StaticQuantMaps.peptideMap.getItem(peptideKey);
							} else {
								quantPeptide = new QuantifiedPeptide(quantPSMs.get(0), isIgnoreTaxonomies(),
										isDistinguishModifiedSequences(), isChargeSensible());
							}
							StaticQuantMaps.peptideMap.addItem(quantPeptide);
							if (!localPeptideMap.containsKey(peptideKey)) {
								localPeptideMap.put(quantPeptide.getKey(), quantPeptide);
							}
						}
						// add the rest of PSMs
						for (final QuantifiedPSMInterface psm : quantPSMs) {
							quantPeptide.addPSM(psm, true);
						}

						// if the line already appear is because it was the same peptide with different
						// protein
						if (repWithPeptidePresent == 1 && uniqueLineStrings.contains(uniqueLineString)) {
							break;
						}
						// add the intensities to the PEPTIDE
						// INTENSITY
						final String intensityString = split[getIndexByColumnAndExperiment(rep, INTENSITY)];
						if (intensityString != null && !"NA".equals(intensityString)
								&& !"null*".equals(intensityString)) {
							final double intensity = Double.valueOf(intensityString);
							final Amount intensityAmount = new AmountEx(intensity, AmountType.INTENSITY,
									conditionByExp.get(rep));
							quantPeptide.addAmount(intensityAmount);
						}
						// NORM_INTENSITY is not present, I have to normalize myself, I will at the end
						// of the loop here
						if (normIntensityColumnPerExperiment.containsKey(rep)) {
							final String normIntensityString = split[normIntensityColumnPerExperiment.get(rep)];
							if (normIntensityString != null && !"NA".equals(normIntensityString)) {
								final double normIntensity = Double.valueOf(normIntensityString);
								final Amount normIntensityAmount = new AmountEx(normIntensity,
										AmountType.NORMALIZED_INTENSITY, conditionByExp.get(rep));
								quantPeptide.addAmount(normIntensityAmount);
							}
						}
						// CORRIONINJECTION_INTENSITY
						if (getIndexByColumnAndExperiment(rep, CORRIONINJECTION_INTENSITY) != -1) {
							final String corrInjectionIntensityString = split[getIndexByColumnAndExperiment(rep,
									CORRIONINJECTION_INTENSITY)];
							if (corrInjectionIntensityString != null && !"NA".equals(corrInjectionIntensityString)) {
								final double corrInjectionIntensity = Double.valueOf(corrInjectionIntensityString);
								final Amount corrInjectionIntensityAmount = new AmountEx(corrInjectionIntensity,
										AmountType.CORRIONINJECTION_INTENSITY, conditionByExp.get(rep));
								quantPeptide.addAmount(corrInjectionIntensityAmount);
							}
						}
						// add the scores to the PSM
						// XCorr
						if (getIndexByColumnAndExperiment(rep, XCORR) != -1) {
							final String xcorrString = split[getIndexByColumnAndExperiment(rep, XCORR)];
							if (xcorrString != null && !"NA".equals(xcorrString)) {
								final Score score = new ScoreEx(xcorrString, XCORR, scoreType, null);
								quantPeptide.addScore(score);
							}
						}
						// DCN
						if (getIndexByColumnAndExperiment(rep, DCN) != -1) {
							final String dcnString = split[getIndexByColumnAndExperiment(rep, DCN)];
							if (dcnString != null && !"NA".equals(dcnString)) {
								final Score score = new ScoreEx(dcnString, DCN, scoreType, null);
								quantPeptide.addScore(score);
							}
						}
						// Retention time to the first PSM (this is not entirely correct)
						if (getIndexByColumnAndExperiment(rep, RETENTIONTIME) != -1) {
							final String rtString = split[getIndexByColumnAndExperiment(rep, RETENTIONTIME)];
							if (rtString != null && !"NA".equals(rtString)) {
								try {
									final float rt = Float.valueOf(rtString);
									if (quantPSMs.get(0) instanceof AbstractPSM) {
										((AbstractPSM) quantPSMs.get(0)).setRtInMinutes(rt);
									}
								} catch (final NumberFormatException e) {
									// do nothing
								}
							}
						}
						// CCS
						final String ccsString = split[getIndexByColumnAndExperiment(rep, CCS)];
						if (ccsString != null && !"0".equals(ccsString)) {
							final double ccs = Double.valueOf(ccsString);
							final Amount ccsAmount = new AmountEx(ccs, AmountType.CCS, conditionByExp.get(rep));
							quantPeptide.addAmount(ccsAmount);
						}
						// XIC
						final String xicString = split[getIndexByColumnAndExperiment(rep, XIC)];
						if (xicString != null && !"-1".equals(xicString)) {
							final double xic = Double.valueOf(xicString);
							final Amount xicAmount = new AmountEx(xic, AmountType.XIC, conditionByExp.get(rep));
							quantPeptide.addAmount(xicAmount);
						}
						// ESTIMATED_XIC (I dont grab it for now...unless they tell me is important)
						if (false) {
							final String estimatedXicString = split[getIndexByColumnAndExperiment(rep, ESTIMATED_XIC)];
							if (estimatedXicString != null && !"-1".equals(estimatedXicString)) {
								final double xic = Double.valueOf(estimatedXicString);
								final Amount xicAmount = new AmountEx(xic, AmountType.XIC, conditionByExp.get(rep));
								quantPeptide.addAmount(xicAmount);
							}
						}

						if (repWithPeptidePresent == 1) {
							uniqueLineStrings.add(uniqueLineString);
						}
					}
					// asign the pvalue and qvalue to the peptide
					if (indexByColumn.containsKey(PVALUE)) {
						final String pvalue = split[indexByColumn.get(PVALUE)];
						quantPeptide.addScore(new ScoreEx(pvalue, PVALUE,
								"p-value at peptide level between " + columnsByExperiments.size() + " experiments",
								null));
					}
					if (indexByColumn.containsKey(QVALUE)) {
						final String qvalue = split[indexByColumn.get(QVALUE)];
						quantPeptide.addScore(new ScoreEx(qvalue, QVALUE,
								"q-value at peptide level between " + columnsByExperiments.size() + " experiments",
								null));
					}
					// create protein(s)
					final List<String> accs = new ArrayList<String>();

					if (indexByColumn.containsKey(LOCUS)) {
						String proteinAccs = split[indexByColumn.get(LOCUS)];
						proteinAccs = removeQuotes(proteinAccs);

						if (proteinAccs.contains(",")) {
							final String[] split2 = proteinAccs.split(",");
							for (final String acc : split2) {
								accs.add(acc);
							}
						} else {
							accs.add(proteinAccs);
						}
					}
					final List<String> descriptions = new ArrayList<String>();
					if (indexByColumn.containsKey(DESCRIPTION)) {
						String proteinDescriptions = split[indexByColumn.get(DESCRIPTION)];
						proteinDescriptions = removeQuotes(proteinDescriptions);

						if (proteinDescriptions.contains(",")) {
							final String[] split2 = proteinDescriptions.split(",");
							for (final String description : split2) {
								descriptions.add(description);
							}
						} else {
							descriptions.add(proteinDescriptions);
						}
					}
					for (int i = 0; i < accs.size(); i++) {
						final String proteinACC = accs.get(i);
						QuantifiedProteinInterface protein = null;
						if (StaticQuantMaps.proteinMap.containsKey(proteinACC)) {
							protein = StaticQuantMaps.proteinMap.getItem(proteinACC);
						} else {
							protein = new QuantifiedProtein(proteinACC, isIgnoreTaxonomies());

						}
						StaticQuantMaps.proteinMap.addItem(protein);

						if (!localProteinMap.containsKey(protein.getKey())) {
							localProteinMap.put(protein.getKey(), protein);
						}
						if (descriptions.size() > i) {
							protein.setDescription(descriptions.get(i));
						}
						protein.addPeptide(quantPeptide, true);

					}
				}
			}
		} catch (final IOException e) {
			throw new QuantParserException(e);
		}
	}

	private String getUniqueLineStringWithNoProtein(String[] split) {
		// here the columns with information that may change between same peptide rows
		// with different proteins
		final TIntSet indexToAvoid = new TIntHashSet();
		indexToAvoid.add(indexByColumn.get(PROTEIN));
		indexToAvoid.add(indexByColumn.get(PROTEIN_DESCRIPTION));
		indexToAvoid.add(indexByColumn.get(PVALUE));
		indexToAvoid.add(indexByColumn.get(QVALUE));
		// we include the peptide sequences because they include extra aminoacids that
		// may be different in different proteins
		for (int rep = 1; rep <= columnsByExperiments.size(); rep++) {
			indexToAvoid.add(getIndexByColumnAndExperiment(rep, SEQUENCE));
		}
		final StringBuilder sb = new StringBuilder();
		for (int index = 0; index < split.length; index++) {
			if (!indexToAvoid.contains(index)) {
				sb.append(split[index] + ",");
			}
		}
		return sb.toString();
	}

	private String removeQuotes(String string) {
		if (string.startsWith("\"") && string.endsWith("\"")) {
			final String ret = string.substring(1, string.length() - 1);
			return ret;
		}
		return string;
	}

	private void processColumns(String[] split) {
		int exp = -1;
		int index = 0;
		final Pattern pattern = Pattern.compile("^(.+)_(\\d+)$");

		for (final String header : split) {
			final Matcher matcher = pattern.matcher(header);
			final boolean matchFound = matcher.find();
			if (matchFound) {
				exp = Integer.valueOf(matcher.group(2));
				final String conditionName = REPLICATE_PREFIX + exp;
				conditionByExp.put(exp, new QuantCondition(conditionName));
				if (!columnsByExperiments.containsKey(exp)) {
					columnsByExperiments.put(exp, new TObjectIntHashMap<String>());
				}
				final String column = matcher.group(1);
				columnsByExperiments.get(exp).put(column, index);
			} else {
				indexByColumn.put(header, index);
			}
			index++;
		}

	}

	private int getIndexByColumnAndExperiment(int experiment, String columnName) {
		final TObjectIntMap<String> tObjectIntMap = columnsByExperiments.get(experiment);
		if (tObjectIntMap.containsKey(columnName)) {
			return tObjectIntMap.get(columnName);
		}
		return -1;
	}

	@Override
	public boolean canRead() {
		try {
			List<String> lines = null;
			// check whether it is an excel file
			if (FileUtils.isExcelFile(file)) {
				lines = FileUtils.readLinesFromXLSX(file, "\t", 0);
			} else {
				lines = Files.readAllLines(this.file.toPath());
			}

			final String line = lines.get(0);

			if (line.startsWith("SEQUENCE\tCHARGE_STATE\tFILENAME_")) {

				return true;

			}

			return false;

		} catch (final Exception e) {
			return false;
		}
	}

}
