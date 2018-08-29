package edu.scripps.yates.census.read.model;

import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.builder.HashCodeBuilder;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.HasIsoRatios;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.maths.Maths;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class IsobaricQuantifiedPeptide extends QuantifiedPeptide implements QuantifiedPeptideInterface, HasIsoRatios {

	private final Map<String, IonCountRatio> countRatiosByConditionKey = new THashMap<String, IonCountRatio>();
	private Map<QuantCondition, Set<Ion>> ionsByConditions;

	/**
	 * Creates a {@link IsobaricQuantifiedPeptide} object, adding the
	 * {@link QuantifiedPSMInterface} to its list of
	 * {@link QuantifiedPSMInterface}s
	 *
	 * @param quantPSM
	 * @param distinguishModifiedSequences
	 */
	public IsobaricQuantifiedPeptide(IsobaricQuantifiedPSM quantPSM, boolean ignoreTaxonomies) {
		super(quantPSM, ignoreTaxonomies);
	}

	public Set<IsobaricQuantifiedPSM> getIsobaricQuantifiedPSMs() {
		final Set<IsobaricQuantifiedPSM> isoPsms = new THashSet<IsobaricQuantifiedPSM>();
		for (final QuantifiedPSMInterface quantifiedPSM : psms) {
			if (quantifiedPSM instanceof IsobaricQuantifiedPSM) {
				isoPsms.add((IsobaricQuantifiedPSM) quantifiedPSM);
			}
		}
		return isoPsms;
	}

	@Override
	public Set<IsoRatio> getIsoRatios() {
		final Set<IsoRatio> ret = new THashSet<IsoRatio>();

		final Set<IsobaricQuantifiedPSM> quantifiedPSMs = getIsobaricQuantifiedPSMs();
		for (final IsobaricQuantifiedPSM quantifiedPSM : quantifiedPSMs) {
			ret.addAll(quantifiedPSM.getIsoRatios());
		}

		return ret;
	}

	/**
	 * Returns true if the protein contains any {@link Ion} labeled with a
	 * certain {@link QuantificationLabel} not paired with any other label
	 *
	 * @param label
	 * @return
	 */
	@Override
	public boolean containsAnySingletonIon(QuantificationLabel label) {
		for (final IsobaricQuantifiedPSM quantifiedPSM : getIsobaricQuantifiedPSMs()) {
			if (quantifiedPSM.containsAnySingletonIon(label))
				return true;
		}
		return false;
	}

	/**
	 * Returns true if the protein contains any {@link Ion} labeled with a
	 * certain {@link QuantificationLabel} not matter if they are paired with
	 * any other label or not (getting ratios or not)
	 *
	 * @param label
	 * @return
	 */
	@Override
	public boolean containsAnyIon(QuantificationLabel label) {
		for (final IsobaricQuantifiedPSM quantifiedPeptide : getIsobaricQuantifiedPSMs()) {
			if (quantifiedPeptide.containsAnyIon(label))
				return true;
		}
		return false;
	}

	@Override
	public Map<QuantificationLabel, Set<Ion>> getIons(IonSerie ionSerie) {
		final Map<QuantificationLabel, Set<Ion>> ret = new THashMap<QuantificationLabel, Set<Ion>>();
		for (final IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			final Map<QuantificationLabel, Set<Ion>> ions = psm.getIons(ionSerie);
			mergeMaps(ret, ions);

		}
		return ret;
	}

	@Override
	public Set<Ion> getSingletonIonsByLabel(QuantificationLabel label) {
		final Set<Ion> ret = new THashSet<Ion>();
		for (final IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			ret.addAll(psm.getSingletonIonsByLabel(label));
		}
		return ret;
	}

	@Override
	public Set<Ion> getIonsByLabel(QuantificationLabel label) {
		final Set<Ion> ret = new THashSet<Ion>();
		for (final IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			ret.addAll(psm.getIonsByLabel(label));
		}
		return ret;
	}

	@Override
	public Map<QuantificationLabel, Set<Ion>> getSingletonIons(IonSerie ionSerie) {
		final Map<QuantificationLabel, Set<Ion>> ret = new THashMap<QuantificationLabel, Set<Ion>>();
		for (final IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			final Map<QuantificationLabel, Set<Ion>> singletonIons = psm.getSingletonIons(ionSerie);
			mergeMaps(ret, singletonIons);

		}
		return ret;
	}

	private void mergeMaps(Map<QuantificationLabel, Set<Ion>> receiverMap,
			Map<QuantificationLabel, Set<Ion>> donorMap) {
		for (final QuantificationLabel label : donorMap.keySet()) {
			if (receiverMap.containsKey(label)) {
				receiverMap.get(label).addAll(donorMap.get(label));
			} else {
				final Set<Ion> set = new THashSet<Ion>();
				set.addAll(donorMap.get(label));
				receiverMap.put(label, set);
			}
		}

	}

	@Override
	public Map<QuantCondition, Set<Ion>> getSingletonIonsByCondition() {
		final Map<QuantCondition, Set<Ion>> ret = new THashMap<QuantCondition, Set<Ion>>();
		for (final IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			final Map<QuantCondition, Set<Ion>> singletonIons = psm.getSingletonIonsByCondition();
			for (final QuantCondition condition : singletonIons.keySet()) {
				if (ret.containsKey(condition)) {
					ret.get(condition).addAll(singletonIons.get(condition));
				} else {
					final Set<Ion> set = new THashSet<Ion>();
					set.addAll(singletonIons.get(condition));
					ret.put(condition, set);
				}
			}
		}
		return ret;
	}

	@Override
	public Set<IsoRatio> getNonInfinityIsoRatios() {
		final Set<IsoRatio> ret = new THashSet<IsoRatio>();
		for (final IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			ret.addAll(psm.getNonInfinityIsoRatios());
		}
		return ret;
	}

	@Override
	public Double getMaxPeak() {
		double max = -Double.MAX_VALUE;
		for (final IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			if (max < psm.getMaxPeak()) {
				max = psm.getMaxPeak();
			}
		}
		return max;
	}

	@Override
	public Map<QuantificationLabel, Set<Ion>> getSingletonIonsByLabel() {
		final Map<QuantificationLabel, Set<Ion>> ret = new THashMap<QuantificationLabel, Set<Ion>>();
		for (final IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			final Map<QuantificationLabel, Set<Ion>> singletonIons = psm.getSingletonIonsByLabel();
			for (final QuantificationLabel label : singletonIons.keySet()) {
				if (ret.containsKey(label)) {
					ret.get(label).addAll(singletonIons.get(label));
				} else {
					final Set<Ion> set = new THashSet<Ion>();
					set.addAll(singletonIons.get(label));
					ret.put(label, set);
				}
			}
		}
		return ret;
	}

	@Override
	public Map<QuantificationLabel, Set<Ion>> getIonsByLabel() {
		final Map<QuantificationLabel, Set<Ion>> ret = new THashMap<QuantificationLabel, Set<Ion>>();
		for (final IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			final Map<QuantificationLabel, Set<Ion>> ions = psm.getIonsByLabel();
			for (final QuantificationLabel label : ions.keySet()) {
				if (ret.containsKey(label)) {
					ret.get(label).addAll(ions.get(label));
				} else {
					final Set<Ion> set = new THashSet<Ion>();
					set.addAll(ions.get(label));
					ret.put(label, set);
				}
			}
		}
		return ret;
	}

	@Override
	public double getMeanRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {
		final TDoubleArrayList ratioValues = new TDoubleArrayList();

		for (final IsoRatio ratio : getNonInfinityIsoRatios()) {
			ratioValues.add(ratio.getLog2Ratio(quantConditionNumerator, quantConditionDenominator));
		}

		return Maths.mean(ratioValues);
	}

	@Override
	public double getSTDRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {
		final TDoubleArrayList ratioValues = new TDoubleArrayList();

		for (final IsoRatio ratio : getNonInfinityIsoRatios()) {
			ratioValues.add(ratio.getLog2Ratio(quantConditionNumerator, quantConditionDenominator));
		}

		return Maths.stddev(ratioValues);
	}

	/**
	 *
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	@Override
	public IonCountRatio getIonCountRatio(QuantCondition cond1, QuantCondition cond2) {
		final String conditionKey = cond1.getName() + cond2.getName();
		if (countRatiosByConditionKey.containsKey(conditionKey)) {
			return countRatiosByConditionKey.get(conditionKey);
		} else {
			final IonCountRatio ratio = QuantUtils.getIonCountRatioForPeptide(this, cond1, cond2);

			countRatiosByConditionKey.put(conditionKey, ratio);
			return ratio;
		}
	}

	@Override
	public IonCountRatio getIonCountRatio(QuantCondition cond1, QuantCondition cond2, String replicateName) {

		return QuantUtils.getIonCountRatioForPeptide(this, cond1, cond2, replicateName);

	}

	@Override
	public Map<QuantCondition, Set<Ion>> getIonsByCondition() {
		if (ionsByConditions == null) {
			ionsByConditions = new THashMap<QuantCondition, Set<Ion>>();
			for (final IsobaricQuantifiedPSM quantPSM : getIsobaricQuantifiedPSMs()) {
				final Map<QuantCondition, Set<Ion>> ions = quantPSM.getIonsByCondition();
				for (final QuantCondition condition : ions.keySet()) {
					final Set<Ion> c = ions.get(condition);
					if (ionsByConditions.containsKey(condition)) {
						ionsByConditions.get(condition).addAll(c);
					} else {
						ionsByConditions.put(condition, c);
					}
				}
			}
		}
		return ionsByConditions;
	}

	@Override
	public Map<QuantCondition, Set<Ion>> getIonsByCondition(String replicateName) {
		if (replicateName == null) {
			if (ionsByConditions == null) {
				ionsByConditions = new THashMap<QuantCondition, Set<Ion>>();
				for (final IsobaricQuantifiedPSM quantPSM : getIsobaricQuantifiedPSMs()) {

					final Map<QuantCondition, Set<Ion>> ions = quantPSM.getIonsByCondition();
					for (final QuantCondition condition : ions.keySet()) {
						final Set<Ion> c = ions.get(condition);
						if (ionsByConditions.containsKey(condition)) {
							ionsByConditions.get(condition).addAll(c);
						} else {
							ionsByConditions.put(condition, c);
						}
					}

				}
			}
			return ionsByConditions;
		} else {
			final Map<QuantCondition, Set<Ion>> ionsByConditions = new THashMap<QuantCondition, Set<Ion>>();
			for (final IsobaricQuantifiedPSM quantPSM : getIsobaricQuantifiedPSMs()) {
				if (replicateName == null || quantPSM.getFileNames().contains(replicateName)) {
					final Map<QuantCondition, Set<Ion>> ions = quantPSM.getIonsByCondition();
					for (final QuantCondition condition : ions.keySet()) {
						final Set<Ion> c = ions.get(condition);
						if (ionsByConditions.containsKey(condition)) {
							ionsByConditions.get(condition).addAll(c);
						} else {
							ionsByConditions.put(condition, c);
						}
					}
				}
			}
			return ionsByConditions;
		}

	}

	@Override
	public Map<QuantCondition, Set<Ion>> getIonsByConditionForSites(String replicateName, char[] quantifiedAAs) {
		final Map<QuantCondition, Set<Ion>> ret = new THashMap<QuantCondition, Set<Ion>>();

		final Set<QuantifiedPSMInterface> quantifiedPSMs = getQuantifiedPSMs();
		for (final QuantifiedPSMInterface quantifiedPSM : quantifiedPSMs) {
			if (quantifiedPSM instanceof IsobaricQuantifiedPSM) {
				final IsobaricQuantifiedPSM isoPSM = (IsobaricQuantifiedPSM) quantifiedPSM;
				final Map<QuantCondition, Set<Ion>> ionsByConditionForSites = isoPSM
						.getIonsByConditionForSites(replicateName, quantifiedAAs);
				for (final QuantCondition cond : ionsByConditionForSites.keySet()) {
					if (ret.containsKey(cond)) {
						ret.get(cond).addAll(ionsByConditionForSites.get(cond));
					} else {
						final Set<Ion> ions = new THashSet<Ion>();
						ions.addAll(ionsByConditionForSites.get(cond));
						ret.put(cond, ions);
					}
				}
			}
		}

		return ret;
	}

	@Override
	public Map<QuantCondition, Set<Ion>> getIonsByConditionForSites(String replicateName, char[] quantifiedAAs,
			int positionInPeptide) {
		final Map<QuantCondition, Set<Ion>> ret = new THashMap<QuantCondition, Set<Ion>>();

		final Set<QuantifiedPSMInterface> quantifiedPSMs = getQuantifiedPSMs();
		for (final QuantifiedPSMInterface quantifiedPSM : quantifiedPSMs) {
			if (quantifiedPSM instanceof IsobaricQuantifiedPSM) {
				final IsobaricQuantifiedPSM isoPSM = (IsobaricQuantifiedPSM) quantifiedPSM;
				final Map<QuantCondition, Set<Ion>> ionsByConditionForSites = isoPSM
						.getIonsByConditionForSites(replicateName, quantifiedAAs, positionInPeptide);
				for (final QuantCondition cond : ionsByConditionForSites.keySet()) {
					if (ret.containsKey(cond)) {
						ret.get(cond).addAll(ionsByConditionForSites.get(cond));
					} else {
						final Set<Ion> ions = new THashSet<Ion>();
						ions.addAll(ionsByConditionForSites.get(cond));
						ret.put(cond, ions);
					}
				}
			}
		}

		return ret;
	}

	/**
	 * In case of OspbaricQuantifiedPeptide, the consensusRatio is the
	 * ionCountRatio
	 */
	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator) {
		return getIonCountRatio(quantConditionNumerator, quantConditionDenominator);
	}

	@Override
	public int hashCode() {
		return HashCodeBuilder.reflectionHashCode(getKey(), false);
	}
}
