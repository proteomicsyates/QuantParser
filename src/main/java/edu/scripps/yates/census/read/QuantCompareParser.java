package edu.scripps.yates.census.read;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.QuantifiedPSM;
import edu.scripps.yates.census.read.model.QuantifiedPeptide;
import edu.scripps.yates.census.read.model.QuantifiedProtein;
import edu.scripps.yates.census.read.model.StaticQuantMaps;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
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

public class QuantCompareParser extends AbstractQuantParser {
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

	private final File file;
	private final TIntObjectMap<TObjectIntMap<String>> columnsByExperiments = new TIntObjectHashMap<TObjectIntMap<String>>();
	private final TObjectIntMap<String> indexByColumn = new TObjectIntHashMap<String>();

	private final TIntIntMap normIntensityColumnPerExperiment = new TIntIntHashMap();
	private boolean ignoreTaxonomies = true;
	private final TIntObjectMap<QuantCondition> conditionByExp = new TIntObjectHashMap<QuantCondition>();

	@Override
	public boolean isIgnoreTaxonomies() {
		return ignoreTaxonomies;
	}

	@Override
	public void setIgnoreTaxonomies(boolean ignoreTaxonomies) {
		this.ignoreTaxonomies = ignoreTaxonomies;
	}

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

		final List<String> lines = Files.readAllLines(this.file.toPath());
		for (final String line : lines) {
			final String[] split = line.split("\t");

			if (line.startsWith(H)) {
				if (split[1].equals(GROUP_SAMPLE)) {
					// do nothing yet. I have to figure out what is that
				}
			} else if (line.startsWith(SLINE)) {
				processColumns(split);
			} else if (split[0].equals(S)) {
				// this is a peptide line
				// take the sequence
				final String rawSequence = split[getIndexByColumnAndExperiment(1, SEQUENCE)];
				QuantifiedPeptideInterface peptide = null;
				for (int exp = 1; exp <= columnsByExperiments.size(); exp++) {
					// create a PSM per experiment
					final String scanNumberString = split[getIndexByColumnAndExperiment(exp, SCAN)];
					int scanNumber = -1;
					try {
						scanNumber = Integer.valueOf(scanNumberString);
					} catch (final NumberFormatException e) {

					}
					int chargeState = -1;
					final String chargeStateString = split[getIndexByColumnAndExperiment(exp, CSTATE)];
					try {
						chargeState = Integer.valueOf(chargeStateString);
					} catch (final NumberFormatException e) {

					}
					String rawFileName = split[getIndexByColumnAndExperiment(exp, FILENAME)];
					if ("NA".equals(rawFileName)) {
						rawFileName = "NA_" + exp;
					}
					final boolean singleton = false;// not used here
					QuantifiedPSMInterface quantPSM = new QuantifiedPSM(rawSequence, null, peptideToSpectraMap,
							scanNumber, chargeState, rawFileName, singleton);
					if (!localPsmMap.containsKey(quantPSM.getKey())) {
						localPsmMap.put(quantPSM.getKey(), quantPSM);
					}
					if (StaticQuantMaps.psmMap.containsKey(quantPSM.getKey())) {
						quantPSM = StaticQuantMaps.psmMap.getItem(quantPSM.getKey());
					}
					StaticQuantMaps.psmMap.addItem(quantPSM);
					if (peptide == null) {
						final String peptideKey = KeyUtils.getInstance().getSequenceKey(quantPSM, true);
						if (StaticQuantMaps.peptideMap.containsKey(peptideKey)) {
							peptide = StaticQuantMaps.peptideMap.getItem(peptideKey);
						} else {
							peptide = new QuantifiedPeptide(quantPSM, ignoreTaxonomies);
						}
						StaticQuantMaps.peptideMap.addItem(peptide);
						if (!localPeptideMap.containsKey(peptideKey)) {
							localPeptideMap.put(peptide.getKey(), peptide);
						}
					} else {
						peptide.addPSM(quantPSM, true);
					}
					// add the intensities to the PSM
					// INTENSITY
					final String intensityString = split[getIndexByColumnAndExperiment(exp, INTENSITY)];
					if (intensityString != null && !"NA".equals(intensityString)) {
						final double intensity = Double.valueOf(intensityString);
						final Amount intensityAmount = new AmountEx(intensity, AmountType.INTENSITY,
								conditionByExp.get(exp));
						quantPSM.addAmount(intensityAmount);
					}
					// NORM_INTENSITY
					final String normIntensityString = split[normIntensityColumnPerExperiment.get(exp)];
					if (normIntensityString != null && !"NA".equals(normIntensityString)) {
						final double normIntensity = Double.valueOf(normIntensityString);
						final Amount normIntensityAmount = new AmountEx(normIntensity, AmountType.NORMALIZED_INTENSITY,
								conditionByExp.get(exp));
						quantPSM.addAmount(normIntensityAmount);
					}
					// add the scores to the PSM
					// XCorr
					final String xcorrString = split[getIndexByColumnAndExperiment(exp, XCORR)];
					if (xcorrString != null && !"NA".equals(xcorrString)) {
						final Score score = new ScoreEx(xcorrString, XCORR, scoreType, null);
						quantPSM.addScore(score);
					}
					// XCorr
					final String dcnString = split[getIndexByColumnAndExperiment(exp, DCN)];
					if (dcnString != null && !"NA".equals(dcnString)) {
						final Score score = new ScoreEx(dcnString, DCN, scoreType, null);
						quantPSM.addScore(score);
					}
					// Retention time
					final String rtString = split[getIndexByColumnAndExperiment(exp, RETENTIONTIME)];
					if (rtString != null && !"NA".equals(rtString)) {
						try {
							final float rt = Float.valueOf(rtString);
							if (quantPSM instanceof AbstractPSM) {
								((AbstractPSM) quantPSM).setRtInMinutes(rt);
							}
						} catch (final NumberFormatException e) {
							// do nothing
						}
					}
				}
				// asign the pvalue and qvalue to the peptide
				final String pvalue = split[indexByColumn.get(PVALUE)];
				peptide.addScore(new ScoreEx(pvalue, PVALUE,
						"p-value at peptide level between " + columnsByExperiments.size() + " experiments", null));
				final String qvalue = split[indexByColumn.get(QVALUE)];
				peptide.addScore(new ScoreEx(qvalue, QVALUE,
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
						protein = new QuantifiedProtein(proteinACC, ignoreTaxonomies);

					}
					StaticQuantMaps.proteinMap.addItem(protein);

					if (!localProteinMap.containsKey(protein.getKey())) {
						localProteinMap.put(protein.getKey(), protein);
					}
					if (descriptions.size() > i) {
						protein.setDescription(descriptions.get(i));
					}
					protein.addPeptide(peptide, true);

				}
			}
		}
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
			if (header.startsWith("EXP_")) {
				exp = Integer.valueOf(header.substring("EXP_".length()));
				columnsByExperiments.put(exp, new TObjectIntHashMap<String>());
				conditionByExp.put(exp, new QuantCondition(header));
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
