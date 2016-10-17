package edu.scripps.yates.census.read.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.grouping.GroupablePSM;
import edu.scripps.yates.utilities.grouping.ProteinEvidence;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.util.Pair;

public class QuantifiedProtein extends AbstractContainsQuantifiedPSMs implements QuantifiedProteinInterface {
	private static final Logger log = Logger.getLogger(QuantifiedProtein.class);

	private final Set<QuantifiedPSMInterface> quantifiedPSMs = new HashSet<QuantifiedPSMInterface>();
	private ProteinEvidence evidence;
	private String accession;
	private ProteinGroup proteinGroup;
	private String accessionType;
	private final Set<Amount> amounts = new HashSet<Amount>();

	protected String description;

	private String taxonomy;

	private final Set<String> fileNames = new HashSet<String>();

	private boolean discarded;

	private String key;

	public QuantifiedProtein(String proteinACC) {
		final Pair<String, String> accPair = FastaParser.getACC(proteinACC);
		accession = accPair.getFirstelement();
		accessionType = accPair.getSecondElement();
	}

	@Override
	public String getKey() {
		if (key != null) {
			return key;
		}
		return getAccession();
	}

	@Override
	public void setKey(String key) {
		this.key = key;
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
		Set<QuantifiedPeptideInterface> ret = new HashSet<QuantifiedPeptideInterface>();
		for (QuantifiedPSMInterface psm : quantifiedPSMs) {
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
		return accession;
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
			String fastaHeader = getDescription();
			final String accession = getAccession();
			taxonomy = FastaParser.getOrganismNameFromFastaHeader(fastaHeader, accession);
		}
		Set<String> set = new HashSet<String>();
		set.add(taxonomy);
		return set;

	}

	@Override
	public void setTaxonomy(String taxonomy) {
		this.taxonomy = taxonomy;
	}

	@Override
	public ProteinEvidence getEvidence() {
		return evidence;
	}

	/**
	 * @return the rawfileNames
	 */

	@Override
	public Set<String> getRawFileNames() {
		Set<String> ret = new HashSet<String>();
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
	 * @param accession
	 *            the primaryAccession to set
	 */
	@Override
	public void setAccession(String accession) {
		this.accession = accession;
	}

	/**
	 * @param accessionType
	 *            the accessionType to set
	 */
	public void setAccessionType(String accessionType) {
		this.accessionType = accessionType;
	}

	@Override
	public String getDescription() {
		return description;
	}

	/**
	 * @param description
	 *            the description to set
	 */
	@Override
	public void setDescription(String description) {
		this.description = description;
	}

	@Override
	public Integer getLength() {
		return null;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getKey() + ": ");
		List<QuantifiedPeptideInterface> list = new ArrayList<QuantifiedPeptideInterface>();
		list.addAll(getQuantifiedPeptides());
		Collections.sort(list, new Comparator<QuantifiedPeptideInterface>() {

			@Override
			public int compare(QuantifiedPeptideInterface o1, QuantifiedPeptideInterface o2) {
				return o1.getFullSequence().compareTo(o2.getFullSequence());
			}
		});
		StringBuilder sb2 = new StringBuilder();
		for (QuantifiedPeptideInterface quantifiedPeptide : list) {
			if (!"".equals(sb2.toString()))
				sb2.append(",");
			sb2.append(quantifiedPeptide.getFullSequence());
		}
		sb.append(sb2);
		return sb.toString();
	}

	@Override
	public void addRatio(QuantRatio ratio) {
		ratios.add(ratio);

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
	public Set<QuantRatio> getNonInfinityRatios() {
		return QuantUtils.getNonInfinityRatios(getRatios());
	}

	@Override
	public void addFileName(String fileName) {
		fileNames.add(fileName);

	}

	@Override
	public Set<String> getFileNames() {
		return fileNames;
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator) {
		return QuantUtils.getAverageRatio(QuantUtils.getNonInfinityRatios(getRatios()), AggregationLevel.PROTEIN);
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator, String replicateName) {
		return QuantUtils.getAverageRatio(QuantUtils.getNonInfinityRatios(getRatios(replicateName)),
				AggregationLevel.PROTEIN);
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
		Set<QuantifiedPeptideInterface> ret = new HashSet<QuantifiedPeptideInterface>();
		for (QuantifiedPeptideInterface peptide : getQuantifiedPeptides()) {
			if (!peptide.isDiscarded()) {
				ret.add(peptide);
			}
		}
		return ret;

	}

	@Override
	public boolean isQuantified() {
		return true;
	}
}