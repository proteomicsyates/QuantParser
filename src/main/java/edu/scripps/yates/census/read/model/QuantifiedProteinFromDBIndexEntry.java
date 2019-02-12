package edu.scripps.yates.census.read.model;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantKeyUtils;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.grouping.GroupablePeptide;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.proteomicsmodel.AbstractProtein;
import edu.scripps.yates.utilities.proteomicsmodel.Accession;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.Ratio;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AggregationLevel;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.set.hash.THashSet;

public class QuantifiedProteinFromDBIndexEntry extends AbstractProtein implements QuantifiedProteinInterface {
	private boolean distinguishModifiedPeptides;
	private final IndexedProtein indexedProtein;
	private boolean discarded;

	private String key;

	private final boolean ignoreACCFormat;

	private THashSet<QuantifiedPeptideInterface> quantifiedPeptides = new THashSet<QuantifiedPeptideInterface>();

	private Set<QuantifiedPSMInterface> quantifiedPSMs;
	private Set<QuantRatio> quantRatios;

	public QuantifiedProteinFromDBIndexEntry(IndexedProtein indexedProtein, boolean ignoreTaxonomies,
			boolean ignoreACCFormat) throws IOException {
		this.indexedProtein = indexedProtein;
		final Accession primaryAcc = FastaParser.getACC(indexedProtein.getFastaDefLine());
		setPrimaryAccession(primaryAcc);
		setIgnoreTaxonomy(ignoreTaxonomies);
		this.ignoreACCFormat = ignoreACCFormat;
	}

	@Override
	public String getKey() {
		if (key == null) {
			key = QuantKeyUtils.getInstance().getProteinKey(indexedProtein, ignoreACCFormat);
		}
		return key;
	}

	@Override
	public String getDescription() {
		if (super.getDescription() == null) {
			if (indexedProtein != null) {
				setDescription(FastaParser.getDescription(indexedProtein.getFastaDefLine()));
			}
		}
		return super.getDescription();
	}

	/**
	 * @return the quantifiedPSMs
	 */
	@Override
	public Set<QuantifiedPSMInterface> getQuantifiedPSMs() {
		if (quantifiedPSMs == null) {
			quantifiedPSMs = new THashSet<QuantifiedPSMInterface>();
			for (final PSM psm : getPSMs()) {
				if (psm instanceof QuantifiedPSMInterface) {
					final QuantifiedPSMInterface quantPSM = (QuantifiedPSMInterface) psm;
					quantifiedPSMs.add(quantPSM);
				}
			}
		}
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
	public boolean addQuantifiedPSM(QuantifiedPSMInterface quantifiedPSM, boolean recursive) {
		return addPSM(quantifiedPSM, recursive);
	}

	@Override
	public int getUniqueID() {
		return hashCode();
	}

	@Override
	public List<GroupablePeptide> getGroupablePeptides() {
		final List<GroupablePeptide> list = new ArrayList<GroupablePeptide>();
		list.addAll(getQuantifiedPSMs());
		return list;
	}

	@Override
	public Set<String> getTaxonomies() {
		if ((super.getTaxonomies() == null || super.getTaxonomies().isEmpty()) && !isIgnoreTaxonomy()) {
			String fastaHeader = null;

			if (indexedProtein != null) {
				fastaHeader = indexedProtein.getFastaDefLine();
			}
			final String accession = getAccession();
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
	 * @param distinguishModifiedPeptides the distinguishModifiedPeptides to set
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
	public boolean addQuantifiedPeptide(QuantifiedPeptideInterface peptide, boolean recursive) {

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

	@Override
	public void setPrimaryAccession(String primaryAccession) {
		super.setPrimaryAccession(primaryAccession);
		// because the key depends on primaryAccession
		key = primaryAccession;
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
	public boolean addQuantRatio(QuantRatio ratio) {
		return getQuantRatios().add(ratio);
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
	public Set<QuantRatio> getNonInfinityRatios() {
		return QuantUtils.getNonInfinityRatios(getQuantRatios());
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
		return QuantUtils.getAverageRatio(QuantUtils.getNonInfinityRatios(getQuantRatios()), AggregationLevel.PROTEIN);
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

	@Override
	public Set<QuantRatio> getQuantRatios() {
		if (quantRatios == null) {
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
}
