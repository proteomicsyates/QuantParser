package edu.scripps.yates.census.read.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.HasIsoRatios;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AggregationLevel;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class IsobaricQuantifiedProteinGroup extends QuantifiedProteinGroup implements HasIsoRatios {
	private static final String SEPARATOR = " ## ";
	private Set<IsoRatio> ratios;
	private final Map<String, IonCountRatio> countRatiosByConditionKey = new THashMap<String, IonCountRatio>();
	private Map<QuantCondition, Set<Ion>> ionsByConditions;
	private StringBuilder accessionString;

	public IsobaricQuantifiedProteinGroup(ProteinGroup proteinGroup) {
		super(proteinGroup);
	}

	@Override
	public int size() {
		return getProteins().size();
	}

	public Set<IsobaricQuantifiedPSM> getIsobaricQuantifiedPSMs() {
		final Set<IsobaricQuantifiedPSM> ret = new THashSet<IsobaricQuantifiedPSM>();
		for (final IsobaricQuantifiedProtein quantifiedProtein : getIsobaricQuantifiedProteins()) {
			final Set<QuantifiedPSMInterface> quantifiedPSMs = quantifiedProtein.getQuantifiedPSMs();
			for (final QuantifiedPSMInterface quantifiedPSMInterface : quantifiedPSMs) {
				if (quantifiedPSMInterface instanceof IsobaricQuantifiedPSM) {
					ret.add((IsobaricQuantifiedPSM) quantifiedPSMInterface);
				}
			}
		}
		return ret;
	}

	/**
	 * NOTE THAT THIS RETURNED LIST IS NOT VALID FOR ADDING NEW PROTEINS TO THE
	 * GROUP
	 *
	 * @return the proteins
	 */
	public List<IsobaricQuantifiedProtein> getIsobaricQuantifiedProteins() {
		final List<IsobaricQuantifiedProtein> ret = new ArrayList<IsobaricQuantifiedProtein>();
		for (final QuantifiedProteinInterface protein : getProteins()) {
			if (protein instanceof IsobaricQuantifiedProtein) {
				ret.add((IsobaricQuantifiedProtein) protein);
			}
		}

		return ret;
	}

	@Override
	public Set<IsoRatio> getIsoRatios() {
		if (ratios == null || ratios.isEmpty()) {
			ratios = new THashSet<IsoRatio>();
			for (final IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
				ratios.addAll(psm.getIsoRatios());
			}
		}
		return ratios;
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
			final Set<Ion> ions1 = getIonsByCondition().get(cond1);
			int numIons1 = 0;
			if (ions1 != null) {
				numIons1 = ions1.size();
			}
			final Set<Ion> ions2 = getIonsByCondition().get(cond2);
			int numIons2 = 0;
			if (ions2 != null) {
				numIons2 = ions2.size();
			}
			final IonCountRatio ratio = new IonCountRatio(AggregationLevel.PROTEINGROUP);
			ratio.addIonCount(cond1, numIons1);
			ratio.addIonCount(cond2, numIons2);
			countRatiosByConditionKey.put(conditionKey, ratio);
			return ratio;
		}
	}

	@Override
	public IonCountRatio getIonCountRatio(QuantCondition cond1, QuantCondition cond2, String replicateName) {

		final Set<Ion> ions1 = getIonsByCondition(replicateName).get(cond1);
		int numIons1 = 0;
		if (ions1 != null) {
			numIons1 = ions1.size();
		}
		final Set<Ion> ions2 = getIonsByCondition(replicateName).get(cond2);
		int numIons2 = 0;
		if (ions2 != null) {
			numIons2 = ions2.size();
		}
		final IonCountRatio ratio = new IonCountRatio(AggregationLevel.PSM);
		ratio.addIonCount(cond1, numIons1);
		ratio.addIonCount(cond2, numIons2);
		return ratio;

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
		final Map<QuantCondition, Set<Ion>> ionsByConditions2 = new THashMap<QuantCondition, Set<Ion>>();
		for (final IsobaricQuantifiedPSM quantPSM : getIsobaricQuantifiedPSMs()) {
			if (quantPSM.getFileNames().contains(replicateName)) {
				final Map<QuantCondition, Set<Ion>> ions = quantPSM.getIonsByCondition();
				for (final QuantCondition condition : ions.keySet()) {
					final Set<Ion> c = ions.get(condition);
					if (ionsByConditions2.containsKey(condition)) {
						ionsByConditions2.get(condition).addAll(c);
					} else {
						ionsByConditions2.put(condition, c);
					}
				}
			}
		}
		return ionsByConditions2;
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

	@Override
	public int hashCode() {
		return Objects.hash(getAccessionString());
	}
}
