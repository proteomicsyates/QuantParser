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
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.grouping.GroupablePSM;
import edu.scripps.yates.utilities.grouping.ProteinEvidence;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class IsobaricQuantifiedProtein extends AbstractContainsQuantifiedPSMs
		implements QuantifiedProteinInterface, HasIsoRatios {
	private static final Logger log = Logger.getLogger(IsobaricQuantifiedProtein.class);

	private final ProteinType protein;
	private final Set<QuantifiedPSMInterface> quantifiedPSMs = new THashSet<QuantifiedPSMInterface>();
	private boolean distinguishModifiedPeptides;
	private ProteinEvidence evidence;
	private String primaryAccession;
	private ProteinGroup proteinGroup;
	private String accessionType;
	private String description;
	private String taxonomy;
	private final Set<QuantRatio> ratios = new THashSet<QuantRatio>();
	private Map<QuantCondition, Set<Ion>> ionsByConditions;
	private final Set<Amount> amounts = new THashSet<Amount>();
	private final Map<String, IonCountRatio> countRatiosByConditionKey = new THashMap<String, IonCountRatio>();

	private boolean discarded;

	private Set<IsoRatio> isoRatios;

	private String key;

	public IsobaricQuantifiedProtein(ProteinType protein) throws IOException {
		this.protein = protein;
		final Pair<String, String> accPair = FastaParser.getACC(protein.getLocus());
		primaryAccession = accPair.getFirstelement();
		accessionType = accPair.getSecondElement();
	}

	public IsobaricQuantifiedProtein(String proteinACC) throws IOException {
		protein = null;
		final Pair<String, String> accPair = FastaParser.getACC(proteinACC);
		primaryAccession = accPair.getFirstelement();
		accessionType = accPair.getSecondElement();
	}

	@Override
	public String getKey() {
		if (key != null) {
			return key;
		}
		return getAccession();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.quant.xml.RelexChro.Protein#getPeptide()
	 */

	public List<Peptide> getPeptide() {
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

	public Integer getSpecCt() {
		if (protein != null)
			return protein.getSpecCt();
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.quant.xml.RelexChro.Protein#getSeqCov()
	 */

	public String getSeqCov() {
		if (protein != null)
			return protein.getSeqCov();
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.quant.xml.RelexChro.Protein#getLength()
	 */

	@Override
	public Integer getLength() {
		if (protein != null) {
			try {
				return Integer.valueOf(protein.getLength());
			} catch (NumberFormatException e) {

			}
		}
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.quant.xml.RelexChro.Protein#getMolwt()
	 */

	public Float getMolwt() {
		if (protein != null) {
			try {
				return Float.valueOf(protein.getMolwt());
			} catch (NumberFormatException e) {

			}
		}
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.quant.xml.RelexChro.Protein#getPi()
	 */

	public Float getPi() {
		if (protein != null) {
			try {
				return Float.valueOf(protein.getPi());
			} catch (NumberFormatException e) {

			}
		}
		return null;
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
		if (description == null) {
			if (protein != null) {
				description = protein.getDesc();
			}
		}
		return description;
	}

	/**
	 * @return the quantifiedPSMs
	 */
	@Override
	public Set<QuantifiedPSMInterface> getQuantifiedPSMs() {
		return quantifiedPSMs;
	}

	/**
	 * @return the quantifiedPeptides
	 */
	@Override
	public Set<QuantifiedPeptideInterface> getQuantifiedPeptides() {
		Set<QuantifiedPeptideInterface> ret = new THashSet<QuantifiedPeptideInterface>();
		for (QuantifiedPSMInterface psm : getQuantifiedPSMs()) {
			if (psm.getQuantifiedPeptide() != null) {
				ret.add(psm.getQuantifiedPeptide());
			}
		}
		return ret;
	}

	@Override
	public boolean addPSM(QuantifiedPSMInterface quantifiedPSM, boolean recursive) {
		if (quantifiedPSMs.contains(quantifiedPSM)) {
			return false;
		}
		quantifiedPSMs.add(quantifiedPSM);
		if (recursive) {
			quantifiedPSM.addQuantifiedProtein(this, false);
		}
		return true;
	}

	@Override
	public int getDBId() {
		return hashCode();
	}

	@Override
	public String getAccession() {
		return primaryAccession;
	}

	@Override
	public String getAccessionType() {
		return accessionType;
	}

	@Override
	public void setEvidence(ProteinEvidence evidence) {
		this.evidence = evidence;

	}

	@Override
	public void setProteinGroup(ProteinGroup proteinGroup) {
		this.proteinGroup = proteinGroup;

	}

	@Override
	public List<GroupablePSM> getGroupablePSMs() {
		List<GroupablePSM> list = new ArrayList<GroupablePSM>();
		list.addAll(getQuantifiedPSMs());
		return list;
	}

	@Override
	public ProteinGroup getProteinGroup() {
		return proteinGroup;
	}

	@Override
	public Set<String> getTaxonomies() {
		if (taxonomy == null) {
			String fastaHeader = null;
			final String accession = getAccession();
			if (protein != null)
				fastaHeader = protein.getDesc();

			final String organismNameFromFastaHeader = FastaParser.getOrganismNameFromFastaHeader(fastaHeader,
					accession);
			taxonomy = organismNameFromFastaHeader;
		}
		Set<String> set = new THashSet<String>();
		set.add(taxonomy);
		return set;
	}

	@Override
	public ProteinEvidence getEvidence() {
		return evidence;
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

	/**
	 * @return the fileNames
	 */
	@Override
	public Set<String> getRawFileNames() {
		Set<String> ret = new THashSet<String>();
		final Set<QuantifiedPSMInterface> quantifiedPSMs2 = getQuantifiedPSMs();
		for (QuantifiedPSMInterface quantifiedPSMInterface : quantifiedPSMs2) {
			ret.addAll(quantifiedPSMInterface.getRawFileNames());
		}
		return ret;
	}

	@Override
	public boolean addPeptide(QuantifiedPeptideInterface peptide, boolean recursive) {
		final Set<QuantifiedPSMInterface> quantifiedPSMs2 = peptide.getQuantifiedPSMs();
		if (quantifiedPSMs.containsAll(quantifiedPSMs2)) {
			return false;
		}
		quantifiedPSMs.addAll(quantifiedPSMs2);
		if (recursive) {
			peptide.addQuantifiedProtein(this, false);
		}
		return true;
	}

	/**
	 * @return the primaryAccession
	 */
	public String getPrimaryAccession() {
		return primaryAccession;
	}

	/**
	 * @param primaryAccession
	 *            the primaryAccession to set
	 */
	public void setPrimaryAccession(String primaryAccession) {
		this.primaryAccession = primaryAccession;
	}

	/**
	 * @param accessionType
	 *            the accessionType to set
	 */
	public void setAccessionType(String accessionType) {
		this.accessionType = accessionType;
	}

	@Override
	public void setAccession(String primaryAccession) {
		this.primaryAccession = primaryAccession;

	}

	@Override
	public void setDescription(String description) {
		this.description = description;

	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getAccession() + ": ");
		List<QuantifiedPeptideInterface> list = new ArrayList<QuantifiedPeptideInterface>();
		final Set<QuantifiedPeptideInterface> quantifiedPeptides = getQuantifiedPeptides();
		for (QuantifiedPeptideInterface peptide : quantifiedPeptides) {
			list.add(peptide);
		}

		Collections.sort(list, new Comparator<QuantifiedPeptideInterface>() {

			@Override
			public int compare(QuantifiedPeptideInterface o1, QuantifiedPeptideInterface o2) {
				return o1.getSequence().compareTo(o2.getSequence());
			}
		});
		StringBuilder sb2 = new StringBuilder();
		for (QuantifiedPeptideInterface quantifiedPeptide : list) {
			if (!"".equals(sb2.toString()))
				sb2.append(",");
			sb2.append(quantifiedPeptide.getSequence());
		}
		sb.append(sb2);
		return sb.toString();
	}

	@Override
	public void setTaxonomy(String taxonomy) {
		this.taxonomy = taxonomy;

	}

	public Set<IsobaricQuantifiedPSM> getIsobaricQuantifiedPSMs() {
		Set<IsobaricQuantifiedPSM> ret = new THashSet<IsobaricQuantifiedPSM>();
		final Set<QuantifiedPSMInterface> quantifiedPSMs2 = getQuantifiedPSMs();
		for (QuantifiedPSMInterface quantifiedPSMInterface : quantifiedPSMs2) {
			if (quantifiedPSMInterface instanceof IsobaricQuantifiedPSM) {
				ret.add((IsobaricQuantifiedPSM) quantifiedPSMInterface);
			}
		}
		return ret;
	}

	public Set<IsobaricQuantifiedPeptide> getIsobaricQuantifiedPeptides() {
		Set<IsobaricQuantifiedPeptide> ret = new THashSet<IsobaricQuantifiedPeptide>();
		final Set<QuantifiedPeptideInterface> quantifiedPeptides2 = getQuantifiedPeptides();
		for (QuantifiedPeptideInterface quantifiedPeptideInterface : quantifiedPeptides2) {
			if (quantifiedPeptideInterface instanceof IsobaricQuantifiedPeptide) {
				ret.add((IsobaricQuantifiedPeptide) quantifiedPeptideInterface);
			}
		}
		return ret;
	}

	@Override
	public Set<QuantRatio> getRatios() {
		if (ratios.isEmpty()) {
			for (IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
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
		for (IsobaricQuantifiedPSM quantifiedPSM : getIsobaricQuantifiedPSMs()) {
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
		for (IsobaricQuantifiedPSM quantifiedPeptide : getIsobaricQuantifiedPSMs()) {
			if (quantifiedPeptide.containsAnyIon(label))
				return true;
		}
		return false;
	}

	@Override
	public Map<QuantificationLabel, Set<Ion>> getIons(IonSerie ionSerie) {
		Map<QuantificationLabel, Set<Ion>> ret = new THashMap<QuantificationLabel, Set<Ion>>();
		for (IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			final Map<QuantificationLabel, Set<Ion>> ions = psm.getIons(ionSerie);
			mergeMaps(ret, ions);

		}
		return ret;
	}

	@Override
	public Set<Ion> getSingletonIonsByLabel(QuantificationLabel label) {
		Set<Ion> ret = new THashSet<Ion>();
		for (IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			ret.addAll(psm.getSingletonIonsByLabel(label));
		}
		return ret;
	}

	@Override
	public Set<Ion> getIonsByLabel(QuantificationLabel label) {
		Set<Ion> ret = new THashSet<Ion>();
		for (IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			ret.addAll(psm.getIonsByLabel(label));
		}
		return ret;
	}

	@Override
	public Map<QuantificationLabel, Set<Ion>> getSingletonIons(IonSerie ionSerie) {
		Map<QuantificationLabel, Set<Ion>> ret = new THashMap<QuantificationLabel, Set<Ion>>();
		for (IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			final Map<QuantificationLabel, Set<Ion>> singletonIons = psm.getSingletonIons(ionSerie);
			mergeMaps(ret, singletonIons);

		}
		return ret;
	}

	private void mergeMaps(Map<QuantificationLabel, Set<Ion>> receiverMap,
			Map<QuantificationLabel, Set<Ion>> donorMap) {
		for (QuantificationLabel label : donorMap.keySet()) {
			if (receiverMap.containsKey(label)) {
				receiverMap.get(label).addAll(donorMap.get(label));
			} else {
				Set<Ion> set = new THashSet<Ion>();
				set.addAll(donorMap.get(label));
				receiverMap.put(label, set);
			}
		}

	}

	@Override
	public Map<QuantCondition, Set<Ion>> getSingletonIonsByCondition() {
		Map<QuantCondition, Set<Ion>> ret = new THashMap<QuantCondition, Set<Ion>>();
		for (IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			final Map<QuantCondition, Set<Ion>> singletonIons = psm.getSingletonIonsByCondition();
			for (QuantCondition condition : singletonIons.keySet()) {
				if (ret.containsKey(condition)) {
					ret.get(condition).addAll(singletonIons.get(condition));
				} else {
					Set<Ion> set = new THashSet<Ion>();
					set.addAll(singletonIons.get(condition));
					ret.put(condition, set);
				}
			}
		}
		return ret;
	}

	@Override
	public Set<IsoRatio> getNonInfinityIsoRatios() {
		Set<IsoRatio> ret = new THashSet<IsoRatio>();
		for (IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			ret.addAll(psm.getNonInfinityIsoRatios());
		}
		return ret;
	}

	@Override
	public Double getMaxPeak() {
		Double max = Double.MIN_VALUE;
		for (IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
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
		Map<QuantificationLabel, Set<Ion>> ret = new THashMap<QuantificationLabel, Set<Ion>>();
		for (IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			final Map<QuantificationLabel, Set<Ion>> singletonIons = psm.getSingletonIonsByLabel();
			for (QuantificationLabel label : singletonIons.keySet()) {
				if (ret.containsKey(label)) {
					ret.get(label).addAll(singletonIons.get(label));
				} else {
					Set<Ion> set = new THashSet<Ion>();
					set.addAll(singletonIons.get(label));
					ret.put(label, set);
				}
			}
		}
		return ret;
	}

	@Override
	public Map<QuantificationLabel, Set<Ion>> getIonsByLabel() {
		Map<QuantificationLabel, Set<Ion>> ret = new THashMap<QuantificationLabel, Set<Ion>>();
		for (IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
			final Map<QuantificationLabel, Set<Ion>> ions = psm.getIonsByLabel();
			for (QuantificationLabel label : ions.keySet()) {
				if (ret.containsKey(label)) {
					ret.get(label).addAll(ions.get(label));
				} else {
					Set<Ion> set = new THashSet<Ion>();
					set.addAll(ions.get(label));
					ret.put(label, set);
				}
			}
		}
		return ret;
	}

	@Override
	public double getMeanRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {
		List<Double> ratioValues = new ArrayList<Double>();

		for (IsoRatio ratio : getNonInfinityIsoRatios()) {
			ratioValues.add(ratio.getLog2Ratio(quantConditionNumerator, quantConditionDenominator));
		}

		return Maths.mean(ratioValues.toArray(new Double[0]));
	}

	@Override
	public double getSTDRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {
		List<Double> ratioValues = new ArrayList<Double>();

		for (IsoRatio ratio : getNonInfinityIsoRatios()) {
			ratioValues.add(ratio.getLog2Ratio(quantConditionNumerator, quantConditionDenominator));
		}

		return Maths.stddev(ratioValues.toArray(new Double[0]));
	}

	@Override
	public Map<QuantCondition, Set<Ion>> getIonsByCondition() {
		if (ionsByConditions == null) {
			ionsByConditions = new THashMap<QuantCondition, Set<Ion>>();
			for (IsobaricQuantifiedPSM quantPSM : getIsobaricQuantifiedPSMs()) {
				final Map<QuantCondition, Set<Ion>> ions = quantPSM.getIonsByCondition();
				for (QuantCondition condition : ions.keySet()) {
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
	public Set<Amount> getAmounts() {
		return amounts;
	}

	@Override
	public void addAmount(Amount amount) {
		amounts.add(amount);
	}

	@Override
	public void addRatio(QuantRatio ratio) {
		getRatios();
		ratios.add(ratio);

	}

	@Override
	public Set<QuantRatio> getNonInfinityRatios() {
		return QuantUtils.getNonInfinityRatios(getRatios());
	}

	@Override
	public Set<String> getFileNames() {
		Set<String> ret = new THashSet<String>();
		for (QuantifiedPSMInterface quantPSM : getQuantifiedPSMs()) {
			ret.addAll(quantPSM.getRawFileNames());
		}
		return ret;
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition cond1, QuantCondition cond2) {
		return QuantUtils.getAverageRatio(QuantUtils.getNonInfinityRatios(getRatios()), AggregationLevel.PROTEIN);
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
	public Set<QuantifiedPeptideInterface> getNonDiscardedQuantifiedPeptides() {
		Set<QuantifiedPeptideInterface> ret = new THashSet<QuantifiedPeptideInterface>();
		for (QuantifiedPeptideInterface peptide : getQuantifiedPeptides()) {
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
			for (IsobaricQuantifiedPSM psm : getIsobaricQuantifiedPSMs()) {
				isoRatios.addAll(psm.getIsoRatios());
			}
		}
		return isoRatios;
	}

	@Override
	public IonCountRatio getIonCountRatio(QuantCondition cond1, QuantCondition cond2) {
		String conditionKey = cond1.getName() + cond2.getName();
		if (countRatiosByConditionKey.containsKey(conditionKey)) {
			return countRatiosByConditionKey.get(conditionKey);
		} else {
			Set<Ion> ions1 = getIonsByCondition().get(cond1);
			int numIons1 = 0;
			if (ions1 != null) {
				numIons1 = ions1.size();
			}
			Set<Ion> ions2 = getIonsByCondition().get(cond2);
			int numIons2 = 0;
			if (ions2 != null) {
				numIons2 = ions2.size();
			}
			IonCountRatio ratio = new IonCountRatio(AggregationLevel.PROTEINGROUP);
			ratio.addIonCount(cond1, numIons1);
			ratio.addIonCount(cond2, numIons2);
			countRatiosByConditionKey.put(conditionKey, ratio);
			return ratio;
		}
	}

	@Override
	public IonCountRatio getIonCountRatio(QuantCondition cond1, QuantCondition cond2, String replicateName) {

		Set<Ion> ions1 = getIonsByCondition(replicateName).get(cond1);
		int numIons1 = 0;
		if (ions1 != null) {
			numIons1 = ions1.size();
		}
		Set<Ion> ions2 = getIonsByCondition(replicateName).get(cond2);
		int numIons2 = 0;
		if (ions2 != null) {
			numIons2 = ions2.size();
		}
		IonCountRatio ratio = new IonCountRatio(AggregationLevel.PSM);
		ratio.addIonCount(cond1, numIons1);
		ratio.addIonCount(cond2, numIons2);
		return ratio;

	}

	@Override
	public Map<QuantCondition, Set<Ion>> getIonsByCondition(String replicateName) {
		Map<QuantCondition, Set<Ion>> ionsByConditions2 = new THashMap<QuantCondition, Set<Ion>>();
		for (IsobaricQuantifiedPSM quantPSM : getIsobaricQuantifiedPSMs()) {
			if (quantPSM.getFileNames().contains(replicateName)) {
				final Map<QuantCondition, Set<Ion>> ions = quantPSM.getIonsByCondition();
				for (QuantCondition condition : ions.keySet()) {
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
		Map<QuantCondition, Set<Ion>> ret = new THashMap<QuantCondition, Set<Ion>>();

		Set<QuantifiedPSMInterface> quantifiedPSMs = getQuantifiedPSMs();
		for (QuantifiedPSMInterface quantifiedPSM : quantifiedPSMs) {
			if (quantifiedPSM instanceof IsobaricQuantifiedPSM) {
				IsobaricQuantifiedPSM isoPSM = (IsobaricQuantifiedPSM) quantifiedPSM;
				Map<QuantCondition, Set<Ion>> ionsByConditionForSites = isoPSM.getIonsByConditionForSites(replicateName,
						quantifiedAAs);
				for (QuantCondition cond : ionsByConditionForSites.keySet()) {
					if (ret.containsKey(cond)) {
						ret.get(cond).addAll(ionsByConditionForSites.get(cond));
					} else {
						Set<Ion> ions = new THashSet<Ion>();
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
		Map<QuantCondition, Set<Ion>> ret = new THashMap<QuantCondition, Set<Ion>>();

		Set<QuantifiedPSMInterface> quantifiedPSMs = getQuantifiedPSMs();
		for (QuantifiedPSMInterface quantifiedPSM : quantifiedPSMs) {
			if (quantifiedPSM instanceof IsobaricQuantifiedPSM) {
				IsobaricQuantifiedPSM isoPSM = (IsobaricQuantifiedPSM) quantifiedPSM;
				Map<QuantCondition, Set<Ion>> ionsByConditionForSites = isoPSM.getIonsByConditionForSites(replicateName,
						quantifiedAAs, positionInPeptide);
				for (QuantCondition cond : ionsByConditionForSites.keySet()) {
					if (ret.containsKey(cond)) {
						ret.get(cond).addAll(ionsByConditionForSites.get(cond));
					} else {
						Set<Ion> ions = new THashSet<Ion>();
						ions.addAll(ionsByConditionForSites.get(cond));
						ret.put(cond, ions);
					}
				}
			}
		}

		return ret;
	}
}
