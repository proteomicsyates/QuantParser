package edu.scripps.yates.census.read.model;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.analysis.util.KeyUtils;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.util.StringPosition;

public class QuantifiedPeptide extends AbstractContainsQuantifiedPSMs implements QuantifiedPeptideInterface {
	protected String sequenceKey;
	protected final Set<QuantifiedPSMInterface> psms = new HashSet<QuantifiedPSMInterface>();
	private final Set<Amount> amounts = new HashSet<Amount>();
	private boolean discarded;
	private List<StringPosition> ptms;
	private final String fullSequence;
	private final String sequence;

	/**
	 * Creates a {@link QuantifiedPeptide} object, adding the
	 * {@link QuantifiedPSMInterface} to its list of
	 * {@link QuantifiedPSMInterface}s
	 *
	 * @param quantPSM
	 * @param distinguishModifiedSequences
	 */
	public QuantifiedPeptide(QuantifiedPSMInterface quantPSM) {
		sequenceKey = KeyUtils.getSequenceKey(quantPSM, true);
		sequence = quantPSM.getSequence();
		fullSequence = quantPSM.getFullSequence();
		addQuantifiedPSM(quantPSM, true);

	}

	@Override
	public String getKey() {
		return sequenceKey;
	}

	@Override
	public void setKey(String key) {
		sequenceKey = key;
	}

	/**
	 * It assures that has the same sequence, taking into account the
	 * distinguishModifiedSequences of the instance
	 */
	@Override
	public boolean addQuantifiedPSM(QuantifiedPSMInterface quantPSM, boolean recursive) {
		if (sequenceKey.equals(KeyUtils.getSequenceKey(quantPSM, true))) {
			if (!psms.contains(quantPSM)) {
				psms.add(quantPSM);
				if (recursive) {
					quantPSM.setQuantifiedPeptide(this, false);
					final Set<QuantifiedProteinInterface> quantifiedProteins = quantPSM.getQuantifiedProteins();
					for (QuantifiedProteinInterface protein : quantifiedProteins) {
						protein.addPeptide(this, false);
					}
				}
			}
		}
		return false;
	}

	/**
	 * Create a Map of {@link QuantifiedPeptide} getting the peptides from the
	 * {@link QuantifiedPSMInterface}
	 *
	 * @param quantifiedPSMs
	 * @param distringuishModifiedPeptides
	 * @return
	 */
	public static Map<String, QuantifiedPeptideInterface> getQuantifiedPeptides(
			Collection<QuantifiedPSMInterface> quantifiedPSMs) {
		Map<String, QuantifiedPeptideInterface> peptideMap = new HashMap<String, QuantifiedPeptideInterface>();
		for (QuantifiedPSMInterface quantifiedPSM : quantifiedPSMs) {
			final QuantifiedPeptideInterface quantifiedPeptide = quantifiedPSM.getQuantifiedPeptide();
			if (!peptideMap.containsKey(quantifiedPeptide.getKey())) {
				peptideMap.put(quantifiedPeptide.getKey(), quantifiedPeptide);
			}
		}
		return peptideMap;
	}

	@Override
	public Set<QuantifiedPSMInterface> getQuantifiedPSMs() {
		return psms;
	}

	@Override
	public String getSequence() {
		return sequence;
	}

	/**
	 *
	 * @return NOte that the returned set is created everytime this method is
	 *         called, because proteins are taken from the psms of the peptide
	 */
	@Override
	public Set<QuantifiedProteinInterface> getQuantifiedProteins() {
		Set<QuantifiedProteinInterface> set = new HashSet<QuantifiedProteinInterface>();
		for (QuantifiedPSMInterface psm : psms) {
			for (QuantifiedProteinInterface quantProtein : psm.getQuantifiedProteins()) {
				set.add(quantProtein);
			}
		}
		return set;
	}

	@Override
	public Float getCalcMHplus() {
		if (!psms.isEmpty())
			return psms.iterator().next().getCalcMHplus();
		return null;
	}

	/**
	 * @return the fileNames
	 */
	@Override
	public Set<String> getRawFileNames() {
		Set<String> ret = new HashSet<String>();
		for (QuantifiedPSMInterface quantPSM : psms) {
			ret.addAll(quantPSM.getRawFileNames());
		}
		return ret;
	}

	@Override
	public String toString() {
		return getKey() + "(x" + getQuantifiedPSMs().size() + ")";
	}

	@Override
	public String getFullSequence() {
		return fullSequence;
	}

	@Override
	public Set<String> getTaxonomies() {
		Set<String> ret = new HashSet<String>();
		for (QuantifiedPSMInterface quantifiedPSM : psms) {
			ret.addAll(quantifiedPSM.getTaxonomies());
		}
		return ret;
	}

	@Override
	public Float getMHplus() {
		if (!psms.isEmpty()) {
			return psms.iterator().next().getMHplus();
		}
		return null;
	}

	@Override
	public Set<Amount> getAmounts() {
		if (amounts.isEmpty()) {
			for (QuantifiedPSMInterface psm : getQuantifiedPSMs()) {
				amounts.addAll(psm.getAmounts());
			}
		}
		return amounts;
	}

	@Override
	public void addAmount(Amount amount) {
		amounts.add(amount);

	}

	@Override
	public void addRatio(QuantRatio ratio) {
		ratios.add(ratio);
	}

	@Override
	public Set<QuantRatio> getNonInfinityRatios() {
		return QuantUtils.getNonInfinityRatios(getRatios());
	}

	@Override
	public Set<String> getFileNames() {
		Set<String> ret = new HashSet<String>();
		for (QuantifiedPSMInterface quantPSM : psms) {
			ret.addAll(quantPSM.getFileNames());
		}
		return ret;
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator) {

		return QuantUtils.getAverageRatio(QuantUtils.getNonInfinityRatios(getRatios()), AggregationLevel.PEPTIDE);
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator, String replicateName) {
		return QuantUtils.getAverageRatio(QuantUtils.getNonInfinityRatios(getRatios(replicateName)),
				AggregationLevel.PEPTIDE);
	}

	@Override
	public boolean isDiscarded() {

		return discarded;
	}

	@Override
	public void setDiscarded(boolean discarded) {
		this.discarded = discarded;
		final Set<QuantifiedPSMInterface> quantifiedPSMs = getQuantifiedPSMs();
		for (QuantifiedPSMInterface quantifiedPSMInterface : quantifiedPSMs) {
			quantifiedPSMInterface.setDiscarded(discarded);
		}
	}

	@Override
	public Set<QuantifiedProteinInterface> getNonDiscardedQuantifiedProteins() {
		Set<QuantifiedProteinInterface> ret = new HashSet<QuantifiedProteinInterface>();
		for (QuantifiedProteinInterface protein : getQuantifiedProteins()) {
			if (!protein.isDiscarded()) {
				ret.add(protein);
			}
		}
		return ret;

	}

	@Override
	public boolean containsPTMs() {

		return !getPtms().isEmpty();
	}

	/**
	 * @return the ptms
	 */
	@Override
	public List<StringPosition> getPtms() {
		if (ptms == null) {
			ptms = FastaParser.getInside(getFullSequence());
		}
		return ptms;
	}

	@Override
	public boolean addQuantifiedProtein(QuantifiedProteinInterface protein, boolean recursive) {
		for (QuantifiedPSMInterface psm : psms) {
			psm.addQuantifiedProtein(protein, false);
		}
		if (recursive) {
			protein.addPeptide(this, false);
		}
		return true;
	}

	@Override
	public boolean isQuantified() {
		for (QuantifiedPSMInterface psm : getQuantifiedPSMs()) {
			if (psm.isQuantified()) {
				return true;
			}
		}
		return false;
	}
}
