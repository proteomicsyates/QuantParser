package edu.scripps.yates.census.read.model;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.quant.xml.ProteinType;
import edu.scripps.yates.census.quant.xml.ProteinType.Peptide;
import edu.scripps.yates.census.read.model.interfaces.HasIsoRatios;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AggregationLevel;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class IsobaricQuantifiedProtein extends QuantifiedProtein implements HasIsoRatios {
	private static final Logger log = Logger.getLogger(IsobaricQuantifiedProtein.class);

	private final ProteinType protein;
	private boolean distinguishModifiedPeptides;
	private Map<QuantCondition, Set<Ion>> ionsByConditions;
	private final Map<String, IonCountRatio> countRatiosByConditionKey = new THashMap<String, IonCountRatio>();
	private Set<IsoRatio> isoRatios;

	public IsobaricQuantifiedProtein(ProteinType protein) throws IOException {
		super(protein.getLocus());
		this.protein = protein;
	}

	public IsobaricQuantifiedProtein(String proteinACC) throws IOException {
		super(proteinACC);
		protein = null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.quant.xml.RelexChro.Protein#getPeptide()
	 */

	public List<Peptide> getCensusChroPeptide() {
		if (protein != null)
			return protein.getPeptide();

		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.quant.xml.RelexChro.Protein#getSeqCt()
	 */

	public Integer getSeqCt() {
		if (protein != null)
			return protein.getSeqCt();
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.quant.xml.RelexChro.Protein#getSpecCt()
	 */

	@Override
	public Integer getSpectrumCount() {
		if (super.getSpectrumCount() == null && protein != null) {
			try {
				setSpectrumCount(Integer.valueOf(protein.getSpecCt()));
			} catch (final NumberFormatException e) {

			}
		}
		return super.getSpectrumCount();
	}

	@Override
	public Float getCoverage() {
		if (super.getCoverage() == null && protein != null) {
			try {
				setCoverage(Float.valueOf(protein.getSeqCov()));
			} catch (final NumberFormatException e) {

			}
		}
		return super.getCoverage();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.quant.xml.RelexChro.Protein#getLength()
	 */

	@Override
	public Integer getLength() {
		if (super.getLength() == null && protein != null) {
			try {
				setLength(Integer.valueOf(protein.getLength()));
			} catch (final NumberFormatException e) {

			}
		}
		return super.getLength();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.quant.xml.RelexChro.Protein#getMolwt()
	 */

	@Override
	public Float getMw() {
		if (super.getMw() == null && protein != null) {
			try {
				setMw(Float.valueOf(protein.getMolwt()));
			} catch (final NumberFormatException e) {

			}
		}
		return super.getMw();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.quant.xml.RelexChro.Protein#getPi()
	 */

	@Override
	public Float getPi() {

		if (super.getPi() == null && protein != null) {
			try {
				setPi(Float.valueOf(protein.getPi()));
			} catch (final NumberFormatException e) {

			}
		}
		return super.getPi();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.quant.xml.RelexChro.Protein#getVal()
	 */

	public String getVal() {
		if (protein != null)
			return protein.getVal();
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.quant.xml.RelexChro.Protein#getDesc()
	 */
	@Override
	public String getDescription() {
		if (super.getDescription() == null && protein != null) {
			setDescription(protein.getDesc());
		}
		return super.getDescription();
	}

	@Override
	public Set<String> getTaxonomies() {
		if (super.getTaxonomies() == null || super.getTaxonomies().isEmpty()) {
			String fastaHeader = null;
			final String accession = getAccession();
			if (protein != null) {
				fastaHeader = protein.getDesc();
			}
			final String organismNameFromFastaHeader = FastaParser.getOrganismNameFromFastaHeader(fastaHeader,
					accession);
			addTaxonomy(organismNameFromFastaHeader);
		}
		return super.getTaxonomies();
	}

	/**
	 * @return the distinguishModifiedPeptides
	 */
	public boolean isDistinguishModifiedPeptides() {
		return distinguishModifiedPeptides;
	}

	/**
	 * @param distinguishModifiedPeptides
	 *            the distinguishModifiedPeptides to set
	 */
	public void setDistinguishModifiedPeptides(boolean distinguishModifiedPeptides) {
		this.distinguishModifiedPeptides = distinguishModifiedPeptides;
	}

	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		sb.append(getAccession() + ": ");
		final List<QuantifiedPeptideInterface> list = new ArrayList<QuantifiedPeptideInterface>();
		final Set<QuantifiedPeptideInterface> quantifiedPeptides = getQuantifiedPeptides();
		for (final QuantifiedPeptideInterface peptide : quantifiedPeptides) {
			list.add(peptide);
		}

		Collections.sort(list, new Comparator<QuantifiedPeptideInterface>() {

			@Override
			public int compare(QuantifiedPeptideInterface o1, QuantifiedPeptideInterface o2) {
				return o1.getSequence().compareTo(o2.getSequence());
			}
		});
		final StringBuilder sb2 = new StringBuilder();
		for (final QuantifiedPeptideInterface quantifiedPeptide : list) {
			if (!"".equals(sb2.toString()))
				sb2.append(",");
			sb2.append(quantifiedPeptide.getSequence());
		}
		sb.append(sb2);
		return sb.toString();
	}

	public Set<IsobaricQuantifiedPSM> getIsobaricQuantifiedPSMs() {
		final Set<IsobaricQuantifiedPSM> ret = new THashSet<IsobaricQuantifiedPSM>();
		final Set<QuantifiedPSMInterface> quantifiedPSMs2 = getQuantifiedPSMs();
		for (final QuantifiedPSMInterface quantifiedPSMInterface : quantifiedPSMs2) {
			if (quantifiedPSMInterface instanceof IsobaricQuantifiedPSM) {
				ret.add((IsobaricQuantifiedPSM) quantifiedPSMInterface);
			}
		}
		return ret;
	}

	public Set<IsobaricQuantifiedPeptide> getIsobaricQuantifiedPeptides() {
		final Set<IsobaricQuantifiedPeptide> ret = new THashSet<IsobaricQuantifiedPeptide>();
		final Set<QuantifiedPeptideInterface> quantifiedPeptides2 = getQuantifiedPeptides();
		for (final QuantifiedPeptideInterface quantifiedPeptideInterface : quantifiedPeptides2) {
			if (quantifiedPeptideInterface instanceof IsobaricQuantifiedPeptide) {
				ret.add((IsobaricQuantifiedPeptide) quantifiedPeptideInterface);
			}
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
		Double max = Double.MIN_VALUE;
		for (final IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			if (max < psm.getMaxPeak()) {
				max = psm.getMaxPeak();
			}
		}
		if (!max.equals(Double.MIN_VALUE)) {
			return max;
		}
		return null;
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
	public Set<QuantifiedPeptideInterface> getNonDiscardedQuantifiedPeptides() {
		final Set<QuantifiedPeptideInterface> ret = new THashSet<QuantifiedPeptideInterface>();
		for (final QuantifiedPeptideInterface peptide : getQuantifiedPeptides()) {
			if (!peptide.isDiscarded()) {
				ret.add(peptide);
			}
		}
		return ret;

	}

	@Override
	public Set<IsoRatio> getIsoRatios() {
		if (isoRatios == null || isoRatios.isEmpty()) {
			isoRatios = new THashSet<IsoRatio>();
			for (final IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
				isoRatios.addAll(psm.getIsoRatios());
			}
		}
		return isoRatios;
	}

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
	public boolean isQuantified() {
		return true;
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

}
