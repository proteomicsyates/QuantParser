package edu.scripps.yates.census.read.model;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantKeyUtils;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.proteomicsmodel.AbstractPeptide;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.Peptide;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.Ratio;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AggregationLevel;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.TCharHashSet;
import gnu.trove.set.hash.THashSet;

public class QuantifiedPeptide extends AbstractPeptide implements QuantifiedPeptideInterface {
	private boolean discarded;
	private Set<QuantRatio> quantRatios;

	/**
	 * Creates a {@link QuantifiedPeptide} object, adding the
	 * {@link QuantifiedPSMInterface} to its list of {@link QuantifiedPSMInterface}s
	 *
	 * @param quantPSM
	 * @param distinguishModifiedSequences
	 */
	public QuantifiedPeptide(QuantifiedPSMInterface quantPSM, boolean ignoreTaxonomies) {
		super(QuantKeyUtils.getInstance().getSequenceKey(quantPSM, true));
		super.setSequence(quantPSM.getSequence());
		setFullSequence(quantPSM.getFullSequence());
		setIgnoreTaxonomy(ignoreTaxonomies);
		addPSM(quantPSM, true);

	}

	/**
	 * It assures that has the same sequence, taking into account the
	 * distinguishModifiedSequences of the instance
	 */
	@Override
	public boolean addPSM(PSM quantPSM, boolean recursive) {
		if (getKey().equals(QuantKeyUtils.getInstance().getSequenceKey(quantPSM, true))) {
			if (!getPSMs().contains(quantPSM)) {
				final boolean ret = super.addPSM(quantPSM, recursive);
				return ret;
			}
			return false;
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
		final Map<String, QuantifiedPeptideInterface> peptideMap = new THashMap<String, QuantifiedPeptideInterface>();
		for (final QuantifiedPSMInterface quantifiedPSM : quantifiedPSMs) {
			final Peptide peptide = quantifiedPSM.getPeptide();
			final QuantifiedPeptideInterface quantifiedPeptide = (QuantifiedPeptideInterface) peptide;
			if (!peptideMap.containsKey(quantifiedPeptide.getKey())) {
				peptideMap.put(quantifiedPeptide.getKey(), quantifiedPeptide);
			}
		}
		return peptideMap;
	}

	/**
	 *
	 * @return NOte that the returned set is created everytime this method is
	 *         called, because proteins are taken from the psms of the peptide
	 */
	@Override
	public Set<Protein> getProteins() {
		final Set<Protein> set = new THashSet<Protein>();
		for (final PSM psm : getPSMs()) {
			for (final Protein quantProtein : psm.getProteins()) {
				set.add(quantProtein);
			}
		}
		return set;
	}

	/**
	 * @return the fileNames
	 */
	@Override
	public Set<String> getRawFileNames() {
		final Set<String> ret = new THashSet<String>();
		for (final PSM psm : getPSMs()) {
			final QuantifiedPSMInterface quantPSM = (QuantifiedPSMInterface) psm;
			ret.addAll(quantPSM.getRawFileNames());
		}
		return ret;
	}

	@Override
	public String toString() {
		return getKey() + "(x" + getPSMs().size() + ")";
	}

	@Override
	public Set<String> getTaxonomies() {
		final Set<String> ret = new THashSet<String>();
		if (!isIgnoreTaxonomy()) {
			for (final PSM psm : getPSMs()) {
				if (psm instanceof QuantifiedPSMInterface) {
					final QuantifiedPSMInterface quantPSM = (QuantifiedPSMInterface) psm;
					if (quantPSM.getTaxonomies() != null) {
						ret.addAll(quantPSM.getTaxonomies());
					}
				}
			}
		}
		return ret;
	}

	@Override
	public Set<Amount> getAmounts() {
		if (super.getAmounts().isEmpty()) {
			for (final PSM psm : getPSMs()) {
				if (psm.getAmounts() != null) {
					for (final Amount amount : psm.getAmounts()) {
						addAmount(amount);
					}
				}
			}
		}
		return super.getAmounts();
	}

	@Override
	public Set<QuantRatio> getNonInfinityRatios() {
		return QuantUtils.getNonInfinityRatios(getQuantRatios());
	}

	@Override
	public Set<String> getFileNames() {
		final Set<String> ret = new THashSet<String>();
		for (final PSM psm : getPSMs()) {
			final QuantifiedPSMInterface quantPSM = (QuantifiedPSMInterface) psm;
			ret.addAll(quantPSM.getFileNames());
		}
		return ret;
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator) {
		final List<QuantRatio> ratios = new ArrayList<QuantRatio>();
		final Set<PositionInPeptide> quantPositionsInPeptides = new THashSet<PositionInPeptide>();
		final TCharHashSet chars = new TCharHashSet();
		for (final PSM psm : getPSMs()) {
			if (psm instanceof QuantifiedPSMInterface) {
				final QuantRatio ratio = QuantUtils.getRepresentativeRatio((QuantifiedPSMInterface) psm);
				if (ratio != null) {
					if (ratio.getCondition1().getName().equals(quantConditionNumerator.getName())) {
						if (ratio.getCondition2().getName().equals(quantConditionDenominator.getName())) {
							ratios.add(ratio);
							if (ratio.getQuantifiedAA() != null) {
								chars.add(ratio.getQuantifiedAA());
							}
							if (ratio.getQuantifiedSitePositionInPeptide() != null) {
								quantPositionsInPeptides.addAll(ratio.getQuantifiedSitePositionInPeptide());
							}
						}
					}
				}
			}
		}

		final QuantRatio ratio = QuantUtils.getAverageRatio(QuantUtils.getNonInfinityRatios(ratios),
				AggregationLevel.PEPTIDE);
		for (final PositionInPeptide positionInPeptide : quantPositionsInPeptides) {
			ratio.addQuantifiedSitePositionInPeptide(positionInPeptide);
		}
		for (final char aa : chars.toArray()) {
			ratio.setQuantifiedAA(aa);
		}

		return ratio;
	}

	@Override
	public boolean isDiscarded() {

		return discarded;
	}

	@Override
	public void setDiscarded(boolean discarded) {
		this.discarded = discarded;
		final List<PSM> psms = getPSMs();
		for (final PSM psm : psms) {
			final QuantifiedPSMInterface quantifiedPSMInterface = (QuantifiedPSMInterface) psm;
			quantifiedPSMInterface.setDiscarded(discarded);
		}
	}

	@Override
	public Set<QuantifiedProteinInterface> getNonDiscardedQuantifiedProteins() {
		final Set<QuantifiedProteinInterface> ret = new THashSet<QuantifiedProteinInterface>();
		for (final Protein protein : getProteins()) {
			if (protein instanceof QuantifiedProteinInterface) {
				final QuantifiedProteinInterface quantProtein = (QuantifiedProteinInterface) protein;
				if (!quantProtein.isDiscarded()) {
					ret.add(quantProtein);
				}
			}
		}
		return ret;

	}

	@Override
	public boolean containsPTMs() {

		return !getPTMsInPeptide().isEmpty();
	}

	@Override
	public boolean isQuantified() {
		for (final PSM psm : getPSMs()) {
			if (psm instanceof QuantifiedPSMInterface) {
				final QuantifiedPSMInterface quantPSM = (QuantifiedPSMInterface) psm;
				if (quantPSM.isQuantified()) {
					return true;
				}
			}
		}
		return false;
	}

	@Override
	public boolean addQuantifiedPSM(QuantifiedPSMInterface quantifiedPSM, boolean b) {
		return addPSM(quantifiedPSM, b);
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

	@Override
	public Set<QuantifiedPSMInterface> getQuantifiedPSMs() {
		final Set<QuantifiedPSMInterface> ret = new THashSet<QuantifiedPSMInterface>();
		for (final PSM psm : getPSMs()) {
			if (psm instanceof QuantifiedPSMInterface) {
				ret.add((QuantifiedPSMInterface) psm);
			}
		}
		return ret;
	}

	@Override
	public Set<QuantRatio> getQuantRatios() {
		if (quantRatios == null || quantRatios.isEmpty()) {
			quantRatios = new THashSet<QuantRatio>();

			final Set<Ratio> ratios = getRatios();
			if (ratios.isEmpty()) {
				final Set<QuantifiedPSMInterface> quantifiedPSMs = getQuantifiedPSMs();
				for (final QuantifiedPSMInterface quantPSM : quantifiedPSMs) {
					final Set<QuantRatio> quantRatiosFromPSM = quantPSM.getQuantRatios();
					for (final QuantRatio ratio : quantRatiosFromPSM) {
						quantRatios.add(ratio);
					}
				}
			} else {
				for (final Ratio ratio : ratios) {
					if (ratio instanceof QuantRatio) {
						quantRatios.add((QuantRatio) ratio);
					}
				}
			}
		}
		return quantRatios;
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
		final boolean ret = super.addRatio(ratio);
		if (ret) {
			quantRatios = null;
		}
		return ret;
	}

	@Override
	public boolean addQuantRatio(QuantRatio ratio) {
		return getQuantRatios().add(ratio);
	}

	@Override
	public double getMeanRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {
		final TDoubleArrayList ratioValues = new TDoubleArrayList();
		for (final QuantifiedPSMInterface psm : getQuantifiedPSMs()) {
			final QuantRatio ratio = QuantUtils.getRepresentativeRatio(psm);
			ratioValues.add(ratio.getLog2Ratio(quantConditionNumerator, quantConditionDenominator));
		}
		return Maths.mean(ratioValues);
	}

	@Override
	public double getSTDRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {
		final TDoubleArrayList ratioValues = new TDoubleArrayList();

		for (final QuantifiedPSMInterface psm : getQuantifiedPSMs()) {
			final QuantRatio ratio = QuantUtils.getRepresentativeRatio(psm);
			ratioValues.add(ratio.getLog2Ratio(quantConditionNumerator, quantConditionDenominator));
		}

		return Maths.stddev(ratioValues);
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof QuantifiedPeptide) {
			return ((QuantifiedPeptide) obj).getKey().equals(getKey());
		}
		return super.equals(obj);
	}
}
