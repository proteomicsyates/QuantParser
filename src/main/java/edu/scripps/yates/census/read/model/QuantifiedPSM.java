package edu.scripps.yates.census.read.model;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.annotations.util.UniprotEntryUtil;
import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.analysis.util.KeyUtils;
import edu.scripps.yates.census.read.model.interfaces.HasRatios;
import edu.scripps.yates.census.read.model.interfaces.PeptideSequenceInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.grouping.GroupablePeptide;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.PeptideRelation;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.sequence.PTMInPeptide;
import edu.scripps.yates.utilities.sequence.PTMInProtein;
import edu.scripps.yates.utilities.strings.StringUtils;
import edu.scripps.yates.utilities.util.StringPosition;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class QuantifiedPSM implements GroupablePeptide, PeptideSequenceInterface, HasRatios, QuantifiedPSMInterface {
	private static final Logger log = Logger.getLogger(QuantifiedPSM.class);
	private final Set<QuantifiedProteinInterface> quantifiedProteins = new THashSet<QuantifiedProteinInterface>();
	private final Set<edu.scripps.yates.census.read.util.QuantificationLabel> labels = new THashSet<QuantificationLabel>();
	private final Set<String> taxonomies = new THashSet<String>();

	private PeptideRelation relation;
	private final Set<QuantRatio> ratios = new THashSet<QuantRatio>();
	private final Set<String> rawFileNames = new THashSet<String>();
	private final String scan;
	private final String sequence;
	private final Map<QuantificationLabel, QuantCondition> conditionsByLabels;
	private QuantifiedPeptideInterface quantifiedPeptide;
	private final String fullSequence;
	private final Integer charge;
	private Float deltaCN;
	private Float xcorr;
	private final Set<Amount> amounts = new THashSet<Amount>();
	private final Set<String> fileNames = new THashSet<String>();
	private boolean discarded;
	private boolean singleton;
	private List<PTMInPeptide> ptms;
	private String key;

	public QuantifiedPSM(String sequence, Map<QuantCondition, QuantificationLabel> labelsByConditions,
			Map<String, Set<String>> peptideToSpectraMap, int scanNumber, int chargeState, String rawFileName,
			boolean singleton) throws IOException {
		fullSequence = FastaParser.getSequenceInBetween(sequence);
		this.sequence = FastaParser.cleanSequence(sequence);
		scan = String.valueOf(scanNumber);
		conditionsByLabels = new THashMap<QuantificationLabel, QuantCondition>();
		if (labelsByConditions != null) {
			for (final QuantCondition condition : labelsByConditions.keySet()) {
				final QuantificationLabel quantificationLabel = labelsByConditions.get(condition);
				conditionsByLabels.put(quantificationLabel, condition);
			}
		}
		charge = chargeState;
		// remove the H of HEAVY
		if (rawFileName != null && rawFileName.startsWith("H")) {
			rawFileNames.add(rawFileName.substring(1));
		} else {
			rawFileNames.add(rawFileName);
		}
		this.singleton = singleton;
		final String peptideKey = KeyUtils.getSequenceKey(this, true);
		final String spectrumKey = KeyUtils.getSpectrumKey(this, true);
		addToMap(peptideKey, peptideToSpectraMap, spectrumKey);

	}

	@Override
	public String getKey() {
		if (key == null) {
			key = getIdentifier();
		}
		return key;
	}

	private void addToMap(String key, Map<String, Set<String>> map, String value) {
		if (map == null)
			return;
		if (map.containsKey(key)) {
			map.get(key).add(value);
		} else {
			final Set<String> set = new THashSet<String>();
			set.add(value);
			map.put(key, set);
		}

	}

	/**
	 * @return the fileName
	 */
	@Override
	public Set<String> getRawFileNames() {
		return rawFileNames;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * edu.scripps.yates.census.quant.xml.RelexChro.Protein.Peptide#getScan()
	 */
	@Override
	public String getScan() {

		return scan;

	}

	@Override
	public boolean addQuantifiedProtein(QuantifiedProteinInterface quantifiedProtein, boolean recursive) {
		if (quantifiedProteins.contains(quantifiedProtein)) {
			return false;
		}
		quantifiedProteins.add(quantifiedProtein);
		if (recursive) {
			quantifiedProtein.addPeptide(getQuantifiedPeptide(), false);
			quantifiedProtein.addPSM(this, false);
		}
		// get the taxonomy
		final Set<String> taxonomies = quantifiedProtein.getTaxonomies();
		taxonomies.addAll(taxonomies);

		return true;
	}

	/**
	 * @return the quantifiedProteins
	 */
	@Override
	public Set<QuantifiedProteinInterface> getQuantifiedProteins() {
		return quantifiedProteins;
	}

	/**
	 * @return the ratios
	 */
	@Override
	public Set<QuantRatio> getRatios() {

		return ratios;
	}

	/**
	 * Gets the labels that this {@link QuantifiedPSM} has been labeled ONLY
	 * with some label.<br>
	 * So, may happen that contains any ratio and it is not labeled
	 *
	 * @return the labels
	 */
	public Set<QuantificationLabel> getLabels() {
		return labels;
	}

	/**
	 * @return the taxonomies
	 */
	@Override
	public Set<String> getTaxonomies() {
		return taxonomies;
	}

	@Override
	public String getSequence() {
		return sequence;
	}

	@Override
	public String getIdentifier() {
		return KeyUtils.getSpectrumKey(this, true);
	}

	@Override
	public void setRelation(PeptideRelation relation) {
		this.relation = relation;

	}

	@Override
	public PeptideRelation getRelation() {
		return relation;
	}

	@Override
	public List<GroupableProtein> getGroupableProteins() {
		final List<GroupableProtein> ret = new ArrayList<GroupableProtein>();
		ret.addAll(getQuantifiedProteins());
		return ret;
	}

	@Override
	public double getMeanRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {

		final Set<QuantRatio> ratioSet = getRatios();
		final TDoubleArrayList ratios = new TDoubleArrayList();

		for (final QuantRatio isoRatio : ratioSet) {
			ratios.add(isoRatio.getLog2Ratio(quantConditionNumerator, quantConditionDenominator));
		}
		return Maths.mean(ratios);
	}

	@Override
	public double getSTDRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {

		final Set<QuantRatio> ratioSet = getRatios();
		final TDoubleArrayList ratios = new TDoubleArrayList();
		for (final QuantRatio isoRatio : ratioSet) {
			ratios.add(isoRatio.getLog2Ratio(quantConditionNumerator, quantConditionDenominator));
		}
		return Maths.stddev(ratios);
	}

	@Override
	public void setQuantifiedPeptide(QuantifiedPeptideInterface quantifiedPeptide, boolean recursive) {
		this.quantifiedPeptide = quantifiedPeptide;
		if (quantifiedPeptide != null) {
			if (recursive) {
				quantifiedPeptide.addQuantifiedPSM(this, false);
				for (final QuantifiedProteinInterface protein : getQuantifiedProteins()) {
					quantifiedPeptide.addQuantifiedProtein(protein, false);
				}
			}
		}

	}

	/**
	 * @return the quantifiedPeptide
	 */
	@Override
	public QuantifiedPeptideInterface getQuantifiedPeptide() {
		return quantifiedPeptide;
	}

	@Override
	public String getFullSequence() {
		return fullSequence;
	}

	@Override
	public Integer getCharge() {
		return charge;
	}

	@Override
	public Float getCalcMHplus() {
		return null;
	}

	@Override
	public Float getMHplus() {
		return null;
	}

	@Override
	public void addRatio(QuantRatio ratio) {
		// dont add it if it is already one ratio with the same value and
		// description
		for (final QuantRatio ratio2 : ratios) {
			if (ratio2.getDescription().equals(ratio.getDescription())) {
				if (Double.compare(ratio2.getValue(), ratio.getValue()) == 0) {
					return;
				}
			}
		}
		ratios.add(ratio);
	}

	@Override
	public Float getDeltaCN() {
		return deltaCN;
	}

	@Override
	public Float getXcorr() {
		return xcorr;
	}

	@Override
	public Float getDeltaMass() {
		return null;
	}

	/**
	 * @param deltaCN
	 *            the deltaCN to set
	 */
	public void setDeltaCN(Float deltaCN) {
		this.deltaCN = deltaCN;
	}

	/**
	 * @param xcorr
	 *            the xcorr to set
	 */
	public void setXcorr(Float xcorr) {
		this.xcorr = xcorr;
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
	public Double getMaxPeak() {

		double max = -Double.MAX_VALUE;
		for (final Amount amount : getAmounts()) {
			if (Double.compare(amount.getValue(), max) > 0) {
				max = amount.getValue();
			}
		}
		if (Double.compare(max, -Double.MAX_VALUE) != 0) {
			return max;
		}
		return null;
	}

	@Override
	public Set<QuantRatio> getNonInfinityRatios() {
		return QuantUtils.getNonInfinityRatios(getRatios());
	}

	@Override
	public Set<String> getFileNames() {
		return fileNames;
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator) {
		return QuantUtils.getAverageRatio(QuantUtils.getNonInfinityRatios(getRatios()), AggregationLevel.PSM);
	}

	@Override
	public boolean isDiscarded() {

		return discarded;
	}

	@Override
	public void setDiscarded(boolean discarded) {
		this.discarded = discarded;

	}

	/**
	 * @return the singleton
	 */
	@Override
	public boolean isSingleton() {
		return singleton;
	}

	public void setSingleton(boolean singleton) {
		this.singleton = singleton;

	}

	@Override
	public boolean containsPTMs() {

		return !getPtms().isEmpty();
	}

	/**
	 * @return the ptms
	 */
	@Override
	public List<PTMInPeptide> getPtms() {
		if (ptms == null) {
			ptms = new ArrayList<PTMInPeptide>();
			final List<StringPosition> tmp = FastaParser.getInside(getFullSequence());
			for (final StringPosition stringPosition : tmp) {
				Double deltaMass = null;
				try {
					deltaMass = Double.valueOf(stringPosition.string);
				} catch (final NumberFormatException e) {

				}
				final PTMInPeptide ptmInPeptide = new PTMInPeptide(stringPosition.position,
						getSequence().charAt(stringPosition.position - 1), getSequence(), deltaMass);
				ptms.add(ptmInPeptide);
			}
		}
		return ptms;
	}

	@Override
	public String toString() {
		return getKey();
	}

	@Override
	public boolean isQuantified() {
		return true;
	}

	@Override
	public List<PTMInProtein> getPTMInProtein(UniprotProteinLocalRetriever uplr, Map<String, String> proteinSequences) {
		final List<PTMInProtein> ptmsInProtein = new ArrayList<PTMInProtein>();
		final List<PTMInPeptide> ptms = getPtms();
		if (ptms.isEmpty()) {
			return ptmsInProtein;
		}

		for (final PTMInPeptide ptmInPeptide : ptms) {
			final int positionInPeptide = ptmInPeptide.getPosition();

			final Set<QuantifiedProteinInterface> proteins = getQuantifiedProteins();
			for (final QuantifiedProteinInterface quantifiedProtein : proteins) {
				final String acc = quantifiedProtein.getAccession();
				String proteinSequence = null;
				// it is important that we look for any protein
				// CONTAINING the accession, so that we can search for
				// the isoforms and proteoforms
				if (proteinSequences != null && proteinSequences.containsKey(acc)) {
					proteinSequence = proteinSequences.get(acc);
				}
				if (proteinSequence == null && uplr != null) {
					final Map<String, Entry> annotatedProtein = uplr.getAnnotatedProtein(null, acc);
					final Entry entry = annotatedProtein.get(acc);
					if (entry != null) {
						proteinSequence = UniprotEntryUtil.getProteinSequence(entry);
					}
				}
				if (proteinSequence != null) {
					final TIntArrayList positionsInProteinSequence = StringUtils.allPositionsOf(proteinSequence,
							getSequence());
					for (final int positionInProteinSequence : positionsInProteinSequence.toArray()) {
						final int positionOfSiteInProtein = positionInProteinSequence + positionInPeptide - 1;
						final PTMInProtein ptmInProtein = new PTMInProtein(positionOfSiteInProtein,
								proteinSequence.charAt(positionOfSiteInProtein - 1), acc, ptmInPeptide.getDeltaMass());
						if (!ptmsInProtein.contains(ptmInProtein)) {
							ptmsInProtein.add(ptmInProtein);
						}
					}
				} else {
					throw new IllegalArgumentException("Protein sequence from protein " + acc
							+ " not found neither in the fasta file nor in Uniprot");
				}

			}

		}

		return ptmsInProtein;
	}
	// @Override
	// public int hashCode() {
	// return Objects.hash(getKey());
	// }

}
