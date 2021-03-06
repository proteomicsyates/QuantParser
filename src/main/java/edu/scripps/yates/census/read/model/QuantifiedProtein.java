package edu.scripps.yates.census.read.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
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
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.proteomicsmodel.AbstractProtein;
import edu.scripps.yates.utilities.proteomicsmodel.Accession;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.Ratio;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AggregationLevel;
import edu.scripps.yates.utilities.proteomicsmodel.factories.AccessionEx;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.set.hash.THashSet;

public class QuantifiedProtein extends AbstractProtein implements QuantifiedProteinInterface {
	private static final Logger log = Logger.getLogger(QuantifiedProtein.class);

	private boolean discarded;

	private Set<QuantifiedPSMInterface> quantifiedPSMs;
	private Set<QuantifiedPeptideInterface> quantifiedPeptides;
	private Set<QuantRatio> quantRatios;

	public QuantifiedProtein(String proteinAcc) {
		this(proteinAcc, false);
	}

	public QuantifiedProtein(String proteinACC, boolean ignoreTaxonomies) {
		this(proteinACC, FastaParser.getACC(proteinACC).getAccession(), ignoreTaxonomies);

	}

	public QuantifiedProtein(String proteinACC, String key, boolean ignoreTaxonomy, boolean takeAccAsItIs) {
		setKey(key);
		Accession acc = null;
		if (!takeAccAsItIs) {
			acc = FastaParser.getACC(proteinACC);
		} else {
			final Accession tmp = FastaParser.getACC(proteinACC);
			acc = new AccessionEx(proteinACC, tmp.getAccessionType());
		}
		super.setPrimaryAccession(acc);
		setIgnoreTaxonomy(ignoreTaxonomy);

	}

	public QuantifiedProtein(String proteinACC, String key, boolean ignoreTaxonomy) {
		this(proteinACC, key, ignoreTaxonomy, false);
	}

	/**
	 * @return the quantifiedPSMs
	 */
	@Override
	public Set<QuantifiedPSMInterface> getQuantifiedPSMs() {
		if (quantifiedPSMs == null || quantifiedPSMs.isEmpty()) {
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
		if (quantifiedPeptides == null || quantifiedPeptides.isEmpty()) {
			quantifiedPeptides = new THashSet<QuantifiedPeptideInterface>();
			for (final QuantifiedPSMInterface psm : getQuantifiedPSMs()) {
				if (psm.getQuantifiedPeptide() != null) {
					quantifiedPeptides.add(psm.getQuantifiedPeptide());
				}
			}
		}
		return quantifiedPeptides;
	}

	@Override
	public boolean addPSM(PSM psm, boolean recursive) {
		final boolean ret = super.addPSM(psm, recursive);
		if (ret) {
			quantifiedPSMs = null;
			quantifiedPeptides = null;
		}
		return ret;
	}

	/**
	 * @return the rawfileNames
	 */

	@Override
	public Set<String> getRawFileNames() {
		final Set<String> ret = new THashSet<String>();
		final Set<QuantifiedPSMInterface> quantifiedPSMs2 = getQuantifiedPSMs();
		for (final QuantifiedPSMInterface quantifiedPSMInterface : quantifiedPSMs2) {
			ret.addAll(quantifiedPSMInterface.getRawFileNames());
		}
		return ret;
	}

	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		sb.append(getKey() + ": ");
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
	public boolean addQuantifiedPeptide(QuantifiedPeptideInterface peptide, boolean recursively) {
		return addPeptide(peptide, recursively);
	}

	@Override
	public boolean addQuantifiedPSM(QuantifiedPSMInterface psm, boolean recursively) {
		return addPSM(psm, recursively);
	}

}
