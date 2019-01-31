package edu.scripps.yates.census.read.model;

import java.io.IOException;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantKeyUtils;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.proteomicsmodel.AbstractPSM;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.Ratio;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AggregationLevel;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class QuantifiedPSM extends AbstractPSM implements QuantifiedPSMInterface {
	private final Set<edu.scripps.yates.census.read.util.QuantificationLabel> labels = new THashSet<QuantificationLabel>();
	private final Set<String> rawFileNames = new THashSet<String>();
	private final Map<QuantificationLabel, QuantCondition> conditionsByLabels;
	private final Set<String> fileNames = new THashSet<String>();
	private boolean discarded;
	private boolean singleton;
	private String key;
	private Set<QuantRatio> quantRatios;

	public QuantifiedPSM(String sequence, Map<QuantCondition, QuantificationLabel> labelsByConditions,
			Map<String, Set<String>> peptideToSpectraMap, int scanNumber, int chargeState, String rawFileName,
			boolean singleton) throws IOException {
		setFullSequence(FastaParser.getSequenceInBetween(sequence));
		setSequence(FastaParser.cleanSequence(sequence));
		setScanNumber(String.valueOf(scanNumber));
		conditionsByLabels = new THashMap<QuantificationLabel, QuantCondition>();
		if (labelsByConditions != null) {
			for (final QuantCondition condition : labelsByConditions.keySet()) {
				final QuantificationLabel quantificationLabel = labelsByConditions.get(condition);
				conditionsByLabels.put(quantificationLabel, condition);
			}
		}
		setChargeState(chargeState);
		// remove the H of HEAVY
		if (rawFileName != null && rawFileName.startsWith("H")) {
			rawFileNames.add(rawFileName.substring(1));
		} else {
			rawFileNames.add(rawFileName);
		}
		this.singleton = singleton;
		final String peptideKey = QuantKeyUtils.getInstance().getSequenceKey(this, true);
		final String spectrumKey = QuantKeyUtils.getInstance().getSpectrumKey(this, true);
		addToMap(peptideKey, peptideToSpectraMap, spectrumKey);

	}

	@Override
	public String getKey() {
		if (key == null) {
			key = getIdentifier();
		}
		return key;
	}

	private void addToMap(String key, Map<String, Set<String>> map, String value) {
		if (map == null)
			return;
		if (map.containsKey(key)) {
			map.get(key).add(value);
		} else {
			final Set<String> set = new THashSet<String>();
			set.add(value);
			map.put(key, set);
		}

	}

	/**
	 * @return the fileName
	 */
	@Override
	public Set<String> getRawFileNames() {
		return rawFileNames;
	}

	@Override
	public boolean addProtein(Protein protein, boolean recursive) {
		super.addProtein(protein, recursive);
		// get the taxonomy
		final Set<String> taxonomies = protein.getTaxonomies();
		for (final String tax : taxonomies) {
			addTaxonomy(tax);
		}

		return true;
	}

	/**
	 * Gets the labels that this {@link QuantifiedPSM} has been labeled ONLY with
	 * some label.<br>
	 * So, may happen that contains any ratio and it is not labeled
	 *
	 * @return the labels
	 */
	public Set<QuantificationLabel> getLabels() {
		return labels;
	}

	@Override
	public String getIdentifier() {
		return QuantKeyUtils.getInstance().getSpectrumKey(this, true);
	}

	@Override
	public double getMeanRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {

		final Set<QuantRatio> ratioSet = getQuantRatios();
		final TDoubleArrayList ratios = new TDoubleArrayList();

		for (final QuantRatio isoRatio : ratioSet) {
			ratios.add(isoRatio.getLog2Ratio(quantConditionNumerator, quantConditionDenominator));
		}
		return Maths.mean(ratios);
	}

	@Override
	public double getSTDRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {

		final Set<QuantRatio> ratioSet = getQuantRatios();
		final TDoubleArrayList ratios = new TDoubleArrayList();
		for (final QuantRatio isoRatio : ratioSet) {
			ratios.add(isoRatio.getLog2Ratio(quantConditionNumerator, quantConditionDenominator));
		}
		return Maths.stddev(ratios);
	}

	@Override
	public boolean addRatio(Ratio ratio) {
		// dont add it if it is already one ratio with the same value and
		// description
		for (final Ratio ratio2 : getRatios()) {
			if (ratio2.getDescription().equals(ratio.getDescription())) {
				if (Double.compare(ratio2.getValue(), ratio.getValue()) == 0) {
					return false;
				}
			}
		}
		final boolean ret = addRatio(ratio);
		if (ret) {
			quantRatios = null;
		}
		return ret;
	}

	@Override
	public Double getMaxPeak() {

		double max = -Double.MAX_VALUE;
		for (final Amount amount : getAmounts()) {
			if (Double.compare(amount.getValue(), max) > 0) {
				max = amount.getValue();
			}
		}
		if (Double.compare(max, -Double.MAX_VALUE) != 0) {
			return max;
		}
		return null;
	}

	@Override
	public Set<QuantRatio> getNonInfinityRatios() {
		return QuantUtils.getNonInfinityRatios(getQuantRatios());
	}

	@Override
	public Set<String> getFileNames() {
		return fileNames;
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator) {
		return QuantUtils.getAverageRatio(QuantUtils.getNonInfinityRatios(getQuantRatios()), AggregationLevel.PSM);
	}

	@Override
	public boolean isDiscarded() {
		return discarded;
	}

	@Override
	public void setDiscarded(boolean discarded) {
		this.discarded = discarded;

	}

	/**
	 * @return the singleton
	 */
	@Override
	public boolean isSingleton() {
		return singleton;
	}

	@Override
	public void setSingleton(boolean singleton) {
		this.singleton = singleton;

	}

	@Override
	public String toString() {
		return getKey();
	}

	@Override
	public Set<QuantRatio> getQuantRatios() {
		if (quantRatios == null || quantRatios.isEmpty()) {
			quantRatios = new THashSet<QuantRatio>();

			final Set<Ratio> ratios = getRatios();
			for (final Ratio ratio : ratios) {
				if (ratio instanceof QuantRatio) {
					quantRatios.add((QuantRatio) ratio);
				}
			}
		}
		return quantRatios;
	}

	@Override
	public boolean addQuantRatio(QuantRatio ratio) {
		return getQuantRatios().add(ratio);
	}

	@Override
	public boolean isQuantified() {
		return true;
	}

	@Override
	public QuantifiedPeptideInterface getQuantifiedPeptide() {
		if (getPeptide() instanceof QuantifiedPeptideInterface) {
			return (QuantifiedPeptideInterface) getPeptide();
		}
		return null;
	}

	@Override
	public boolean setQuantifiedPeptide(QuantifiedPeptideInterface peptide, boolean recursively) {
		return setPeptide(peptide, recursively);
	}

	@Override
	public Set<QuantifiedProteinInterface> getQuantifiedProteins() {
		final Set<QuantifiedProteinInterface> ret = new THashSet<QuantifiedProteinInterface>();
		for (final Protein protein : getProteins()) {
			if (protein instanceof QuantifiedProteinInterface) {
				ret.add((QuantifiedProteinInterface) protein);
			}
		}
		return ret;
	}

	@Override
	public boolean addQuantifiedProtein(QuantifiedProteinInterface protein, boolean recursively) {
		return addProtein(protein, recursively);
	}

}
