package edu.scripps.yates.census.read;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.THashSet;
import gnu.trove.set.hash.TIntHashSet;

public class QuantCompareParser extends AbstractQuantParser {
	private final static Logger log = Logger.getLogger(QuantCompareParser.class);
	private static final String H = "H";
	private static final String GROUP_SAMPLE = "GROUP_SAMPLE";
	private static final String SLINE = "SLINE";
	private static final String NORM_INTENSITY = "NORM_INTENSITY";
	private static final String PVALUE = "PVALUE";
	private static final String QVALUE = "QVALUE";
	private static final String PROTEIN = "PROTEIN";
	private static final String PROTEIN_DESCRIPTION = "PROTEIN DESCRIPTION";
	private static final String S = "S";
	private static final String SEQUENCE = "SEQUENCE";
	private static final String SCAN = "SCAN";
	private static final String CSTATE = "CSTATE";
	private static final String FILENAME = "FILENAME";
	private static final String INTENSITY = "INTENSITY";
	private static final String XCORR = "XCORR";
	private static final String scoreType = "PSM-level search engine specific statistic";
	private static final String DCN = "DCN";
	private static final String RETENTIONTIME = "RETENTIONTIME";
	private static final String CORRIONINJECTION_INTENSITY = "CORRIONINJECTION_INTENSITY";
	private static final String REDUNDANCY = "REDUNDANCY";
	public static final String REPLICATE_PREFIX = "Rep_";
	private static final String EXP_HEADER_PREFIX = "EXP_";

	private final File file;
	private final TIntObjectMap<TObjectIntMap<String>> columnsByExperiments = new TIntObjectHashMap<TObjectIntMap<String>>();
	private final TObjectIntMap<String> indexByColumn = new TObjectIntHashMap<String>();

	private final TIntIntMap normIntensityColumnPerExperiment = new TIntIntHashMap();
	private final TIntObjectMap<QuantCondition> conditionByExp = new TIntObjectHashMap<QuantCondition>();

	/**
	 * 
	 * @param psmLevelCensusQuantCompareFile input file that is a PSM level census
	 *                                       quant compare file
	 */
	public QuantCompareParser(File psmLevelCensusQuantCompareFile) {
		this.file = psmLevelCensusQuantCompareFile;
	}

	@Override
	protected void process() throws IOException {
		final Map<String, Set<String>> peptideToSpectraMap = new THashMap<String, Set<String>>();
		List<String> lines = null;
		// check whether it is an excel file
		if (FileUtils.isExcelFile(file)) {
			lines = FileUtils.readLinesFromXLSX(file, "\t", 0);
		} else {
			lines = Files.readAllLines(this.file.toPath());
		}
		final Set<String> uniqueLineStrings = new THashSet<String>();
		int numLine = 0;
		for (final String line : lines) {
			numLine++;
			final String[] split = line.split("\t");

			if (line.startsWith(H)) {
				if (split[1].equals(GROUP_SAMPLE)) {
					// do nothing yet. I have to figure out what is that
				}
			} else if (line.startsWith(SLINE)) {
				processColumns(split);
			} else if (split[0].equals(S)) {
				// first of all, check whether the line has been already seen before except for
				// the protein columns
				final String uniqueLineString = getUniqueLineStringWithNoProtein(split);

				// this is a peptide line
				// take the sequence
				final String rawSequence = split[getIndexByColumnAndExperiment(1, SEQUENCE)];
				QuantifiedPeptideInterface quantPeptide = null;
				for (int rep = 1; rep <= columnsByExperiments.size(); rep++) {
					// create a PSM per experiment
					final String scanNumberString = split[getIndexByColumnAndExperiment(rep, SCAN)];
					int scanNumber = -1;
					try {
						scanNumber = Double.valueOf(scanNumberString).intValue();
					} catch (final NumberFormatException e) {

					}
					int chargeState = -1;
					final String chargeStateString = split[getIndexByColumnAndExperiment(rep, CSTATE)];
					try {
						chargeState = Double.valueOf(chargeStateString).intValue();
					} catch (final NumberFormatException e) {

					}
					int redundancy = 1;
					final String redundancyString = split[getIndexByColumnAndExperiment(rep, REDUNDANCY)];
					try {
						redundancy = Double.valueOf(redundancyString).intValue();
					} catch (final NumberFormatException e) {

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
					int scanNumberForPSM = scanNumber;
					for (int i = 1; i <= redundancy; i++) {
						QuantifiedPSMInterface quantPSM = new QuantifiedPSM(rawSequence, null, peptideToSpectraMap,
								scanNumberForPSM++, chargeState, rawFileName, singleton,
								isDistinguishModifiedSequences(), isChargeSensible());
						// make this PSM to be from replicate rep
						quantPSM.addCondition(conditionByExp.get(rep));

						quantPSMs.add(quantPSM);
						if (!localPsmMap.containsKey(quantPSM.getKey())) {
							localPsmMap.put(quantPSM.getKey(), quantPSM);
						}
						if (StaticQuantMaps.psmMap.containsKey(quantPSM.getKey())) {
							quantPSM = StaticQuantMaps.psmMap.getItem(quantPSM.getKey());
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
					if (rep == 1 && uniqueLineStrings.contains(uniqueLineString)) {
						continue;
					}
					// add the intensities to the PEPTIDE
					// INTENSITY
					final String intensityString = split[getIndexByColumnAndExperiment(rep, INTENSITY)];
					if (intensityString != null && !"NA".equals(intensityString)) {
						final double intensity = Double.valueOf(intensityString);
						final Amount intensityAmount = new AmountEx(intensity, AmountType.INTENSITY,
								conditionByExp.get(rep));
						quantPeptide.addAmount(intensityAmount);
					}
					// NORM_INTENSITY
					final String normIntensityString = split[normIntensityColumnPerExperiment.get(rep)];
					if (normIntensityString != null && !"NA".equals(normIntensityString)) {
						final double normIntensity = Double.valueOf(normIntensityString);
						final Amount normIntensityAmount = new AmountEx(normIntensity, AmountType.NORMALIZED_INTENSITY,
								conditionByExp.get(rep));
						quantPeptide.addAmount(normIntensityAmount);
					}
					// CORRIONINJECTION_INTENSITY
					final String corrInjectionIntensityString = split[getIndexByColumnAndExperiment(rep,
							CORRIONINJECTION_INTENSITY)];
					if (corrInjectionIntensityString != null && !"NA".equals(corrInjectionIntensityString)) {
						final double corrInjectionIntensity = Double.valueOf(corrInjectionIntensityString);
						final Amount corrInjectionIntensityAmount = new AmountEx(corrInjectionIntensity,
								AmountType.CORRIONINJECTION_INTENSITY, conditionByExp.get(rep));
						quantPeptide.addAmount(corrInjectionIntensityAmount);
					}
					// add the scores to the PSM
					// XCorr
					final String xcorrString = split[getIndexByColumnAndExperiment(rep, XCORR)];
					if (xcorrString != null && !"NA".equals(xcorrString)) {
						final Score score = new ScoreEx(xcorrString, XCORR, scoreType, null);
						quantPeptide.addScore(score);
					}
					// XCorr
					final String dcnString = split[getIndexByColumnAndExperiment(rep, DCN)];
					if (dcnString != null && !"NA".equals(dcnString)) {
						final Score score = new ScoreEx(dcnString, DCN, scoreType, null);
						quantPeptide.addScore(score);
					}
					// Retention time to the first PSM (this is not entirely correct)
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
					if (rep == 1) {
						uniqueLineStrings.add(uniqueLineString);
					}
				}
				// asign the pvalue and qvalue to the peptide
				final String pvalue = split[indexByColumn.get(PVALUE)];
				quantPeptide.addScore(new ScoreEx(pvalue, PVALUE,
						"p-value at peptide level between " + columnsByExperiments.size() + " experiments", null));
				final String qvalue = split[indexByColumn.get(QVALUE)];
				quantPeptide.addScore(new ScoreEx(qvalue, QVALUE,
						"q-value at peptide level between " + columnsByExperiments.size() + " experiments", null));
				// create protein(s)
				String proteinAccs = split[indexByColumn.get(PROTEIN)];
				proteinAccs = removeQuotes(proteinAccs);
				final List<String> accs = new ArrayList<String>();
				if (proteinAccs.contains(",")) {
					final String[] split2 = proteinAccs.split(",");
					for (final String acc : split2) {
						accs.add(acc);
					}
				} else {
					accs.add(proteinAccs);
				}
				String proteinDescriptions = split[indexByColumn.get(PROTEIN_DESCRIPTION)];
				proteinDescriptions = removeQuotes(proteinDescriptions);
				final List<String> descriptions = new ArrayList<String>();
				if (proteinDescriptions.contains(",")) {
					final String[] split2 = proteinDescriptions.split(",");
					for (final String description : split2) {
						descriptions.add(description);
					}
				} else {
					descriptions.add(proteinDescriptions);
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
		for (final String header : split) {
			if (header.startsWith(EXP_HEADER_PREFIX)) {
				exp = Integer.valueOf(header.substring(EXP_HEADER_PREFIX.length()));
				columnsByExperiments.put(exp, new TObjectIntHashMap<String>());
				String conditionName = REPLICATE_PREFIX;
				if (header.length() > EXP_HEADER_PREFIX.length()) {
					conditionName += header.substring(EXP_HEADER_PREFIX.length());
				}
				conditionByExp.put(exp, new QuantCondition(conditionName));
			} else if (header.startsWith(NORM_INTENSITY)) {
				final int exp2 = Integer.valueOf(header.substring(NORM_INTENSITY.length()));
				normIntensityColumnPerExperiment.put(exp2, index);
			} else if (header.equals(PVALUE) || header.equals(QVALUE) || header.equals(PROTEIN)
					|| header.equals(PROTEIN_DESCRIPTION)) {
				indexByColumn.put(header, index);
			} else {
				if (exp > 0) {
					columnsByExperiments.get(exp).put(header, index);
				}
			}
			index++;
		}

	}

	private int getIndexByColumnAndExperiment(int experiment, String columnName) {
		return columnsByExperiments.get(experiment).get(columnName);
	}

}
