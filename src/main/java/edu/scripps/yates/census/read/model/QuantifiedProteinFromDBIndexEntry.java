package edu.scripps.yates.census.read.model;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.analysis.util.KeyUtils;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.census.read.util.QuantifiedPSMCollectionObserver;
import edu.scripps.yates.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.grouping.GroupablePeptide;
import edu.scripps.yates.utilities.grouping.ProteinEvidence;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.trove.ObservableTHashSet;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.set.hash.THashSet;

public class QuantifiedProteinFromDBIndexEntry extends AbstractContainsQuantifiedPSMs
		implements QuantifiedProteinInterface {
	private static final Logger log = Logger.getLogger(QuantifiedProteinFromDBIndexEntry.class);

	private final ObservableTHashSet<QuantifiedPSMInterface> quantifiedPSMs = new ObservableTHashSet<QuantifiedPSMInterface>();
	private boolean distinguishModifiedPeptides;
	private IndexedProtein indexedProtein;
	private ProteinEvidence evidence;
	private String primaryAccession;
	private ProteinGroup proteinGroup;
	private String accessionType;
	private final Set<Amount> amounts = new THashSet<Amount>();

	private String description;

	private String taxonomy;

	private boolean discarded;

	private String key;

	private final boolean ignoreTaxonomies;

	private final boolean ignoreACCFormat;

	private THashSet<QuantifiedPeptideInterface> quantifiedPeptides = new THashSet<QuantifiedPeptideInterface>();

	public QuantifiedProteinFromDBIndexEntry(IndexedProtein indexedProtein, boolean ignoreTaxonomies,
			boolean ignoreACCFormat) throws IOException {
		this.indexedProtein = indexedProtein;
		final Pair<String, String> accPair = FastaParser.getACC(indexedProtein.getFastaDefLine());
		primaryAccession = accPair.getFirstelement();
		accessionType = accPair.getSecondElement();
		this.ignoreTaxonomies = ignoreTaxonomies;
		this.ignoreACCFormat = ignoreACCFormat;
		quantifiedPSMs.addCollectionObserver(new QuantifiedPSMCollectionObserver(quantifiedPeptides));
	}

	@Override
	public String getKey() {
		if (key == null) {
			key = KeyUtils.getProteinKey(indexedProtein, ignoreACCFormat);
		}
		return key;
	}

	@Override
	public String getDescription() {
		if (description == null) {
			if (indexedProtein != null) {
				description = FastaParser.getDescription(indexedProtein.getFastaDefLine());
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
		if (quantifiedPeptides == null) {
			quantifiedPeptides = new THashSet<QuantifiedPeptideInterface>();
			for (final QuantifiedPSMInterface psm : quantifiedPSMs) {
				if (psm.getQuantifiedPeptide() != null) {
					quantifiedPeptides.add(psm.getQuantifiedPeptide());
				}
			}
		}
		return quantifiedPeptides;
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
	public List<GroupablePeptide> getGroupablePeptides() {
		final List<GroupablePeptide> list = new ArrayList<GroupablePeptide>();
		list.addAll(getQuantifiedPSMs());
		return list;
	}

	@Override
	public ProteinGroup getProteinGroup() {
		return proteinGroup;
	}

	@Override
	public Set<String> getTaxonomies() {
		if (taxonomy == null && !ignoreTaxonomies) {
			String fastaHeader = null;

			if (indexedProtein != null)
				fastaHeader = indexedProtein.getFastaDefLine();
			final String accession = getAccession();

			final String organismNameFromFastaHeader = FastaParser.getOrganismNameFromFastaHeader(fastaHeader,
					accession);
			taxonomy = organismNameFromFastaHeader;
		}
		final Set<String> set = new THashSet<String>();
		if (taxonomy != null) {
			set.add(taxonomy);
		}
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

	@Override
	public Set<String> getRawFileNames() {
		final Set<String> ret = new THashSet<String>();
		for (final QuantifiedPSMInterface quantifiedPSMInterface : getQuantifiedPSMs()) {
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

	/**
	 * @param indexedProtein
	 *            the indexedProtein to set
	 */
	public void setIndexedProtein(IndexedProtein indexedProtein) {
		this.indexedProtein = indexedProtein;
	}

	@Override
	public Integer getLength() {
		return null;
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
		final StringBuilder sb = new StringBuilder();
		sb.append(getAccession() + ": ");
		final List<QuantifiedPeptideInterface> list = new ArrayList<QuantifiedPeptideInterface>();
		list.addAll(getQuantifiedPeptides());
		Collections.sort(list, new Comparator<QuantifiedPeptideInterface>() {

			@Override
			public int compare(QuantifiedPeptideInterface o1, QuantifiedPeptideInterface o2) {
				return o1.getFullSequence().compareTo(o2.getFullSequence());
			}
		});
		final StringBuilder sb2 = new StringBuilder();
		for (final QuantifiedPeptideInterface quantifiedPeptide : list) {
			if (!"".equals(sb2.toString()))
				sb2.append(",");
			sb2.append(quantifiedPeptide.getFullSequence());
		}
		sb.append(sb2);
		return sb.toString();
	}

	@Override
	public void setTaxonomy(String taxonomy) {
		if (taxonomy != null)
			this.taxonomy = taxonomy;
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
		ratios.add(ratio);

	}

	@Override
	public Set<QuantRatio> getNonInfinityRatios() {
		return QuantUtils.getNonInfinityRatios(getRatios());
	}

	@Override
	public Set<String> getFileNames() {
		final Set<String> ret = new THashSet<String>();
		for (final QuantifiedPSMInterface quantPSM : getQuantifiedPSMs()) {
			ret.addAll(quantPSM.getRawFileNames());
		}
		return ret;
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator) {
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
		for (final QuantifiedPSMInterface quantifiedPSMInterface : quantifiedPSMs) {
			quantifiedPSMInterface.setDiscarded(discarded);
		}
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
	public boolean isQuantified() {
		return true;
	}
}
