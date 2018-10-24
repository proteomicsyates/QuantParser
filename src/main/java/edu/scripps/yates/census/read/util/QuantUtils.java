package edu.scripps.yates.census.read.util;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinRetriever;
import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.analysis.util.KeyUtils;
import edu.scripps.yates.census.read.CensusOutParser;
import edu.scripps.yates.census.read.model.CensusRatio;
import edu.scripps.yates.census.read.model.Ion;
import edu.scripps.yates.census.read.model.IonCountRatio;
import edu.scripps.yates.census.read.model.IsoRatio;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPSM;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPeptide;
import edu.scripps.yates.census.read.model.QuantifiedPSM;
import edu.scripps.yates.census.read.model.QuantifiedPeptide;
import edu.scripps.yates.census.read.model.RatioScore;
import edu.scripps.yates.census.read.model.interfaces.HasRatios;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.model.enums.AmountType;
import edu.scripps.yates.utilities.model.enums.CombinationType;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.sequence.PTMInPeptide;
import edu.scripps.yates.utilities.sequence.PTMInProtein;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import edu.scripps.yates.utilities.strings.StringUtils;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class QuantUtils {
	private static final Logger log = Logger.getLogger(QuantUtils.class);
	private static UniprotProteinRetriever uplr;

	public static void addToPeptideMap(QuantifiedPSMInterface quantifiedPSM, Map<String, QuantifiedPeptide> map,
			boolean ignoreTaxonomies) {
		final String sequenceKey = KeyUtils.getSequenceKey(quantifiedPSM, true);
		QuantifiedPeptide quantifiedPeptide = null;
		if (map.containsKey(sequenceKey)) {
			quantifiedPeptide = map.get(sequenceKey);
			quantifiedPeptide.addQuantifiedPSM(quantifiedPSM, true);
		} else {
			quantifiedPeptide = new QuantifiedPeptide(quantifiedPSM, ignoreTaxonomies);
			map.put(sequenceKey, quantifiedPeptide);
		}

		quantifiedPSM.setQuantifiedPeptide(quantifiedPeptide, true);
	}

	public static void addToIsobaricPeptideMap(IsobaricQuantifiedPSM quantifiedPSM,
			Map<String, IsobaricQuantifiedPeptide> map, boolean ignoreTaxonomies) {
		final String sequenceKey = KeyUtils.getSequenceKey(quantifiedPSM, true);
		IsobaricQuantifiedPeptide quantifiedPeptide = null;
		if (map.containsKey(sequenceKey)) {
			quantifiedPeptide = map.get(sequenceKey);
			quantifiedPeptide.addQuantifiedPSM(quantifiedPSM, true);
		} else {
			quantifiedPeptide = new IsobaricQuantifiedPeptide(quantifiedPSM, ignoreTaxonomies);
			map.put(sequenceKey, quantifiedPeptide);
		}

		quantifiedPSM.setQuantifiedPeptide(quantifiedPeptide, true);
	}

	/**
	 * Create a Map of {@link IsobaricQuantifiedPeptide} getting the peptides
	 * from the {@link QuantifiedPSMInterface}
	 *
	 * @param quantifiedPSMs
	 * @param distringuishModifiedPeptides
	 * @return
	 */
	public static Map<String, IsobaricQuantifiedPeptide> getIsobaricQuantifiedPeptides(
			Collection<IsobaricQuantifiedPSM> quantifiedPSMs) {
		final Map<String, IsobaricQuantifiedPeptide> peptideMap = new THashMap<String, IsobaricQuantifiedPeptide>();

		for (final IsobaricQuantifiedPSM quantifiedPSM : quantifiedPSMs) {
			final QuantifiedPeptideInterface quantifiedPeptide = quantifiedPSM.getQuantifiedPeptide();
			if (quantifiedPeptide instanceof IsobaricQuantifiedPeptide) {
				if (!peptideMap.containsKey(quantifiedPeptide.getKey())) {
					peptideMap.put(quantifiedPeptide.getKey(), (IsobaricQuantifiedPeptide) quantifiedPeptide);
				}
			}
		}

		return peptideMap;
	}

	public static Set<QuantRatio> getNonInfinityRatios(Collection<QuantRatio> ratios) {
		final Set<QuantRatio> set = new THashSet<QuantRatio>();
		for (final QuantRatio ratio : ratios) {
			final double log2Ratio = ratio.getLog2Ratio(ratio.getQuantCondition1(), ratio.getQuantCondition2());
			if (Double.isInfinite(log2Ratio) || Double.compare(log2Ratio, Double.MAX_VALUE) == 0
					|| Double.compare(log2Ratio, -Double.MAX_VALUE) == 0) {
				continue;
			}
			if (!Double.isInfinite(log2Ratio) && !Double.isNaN(log2Ratio)) {
				set.add(ratio);
			}
		}
		return set;
	}

	/**
	 * Creates a {@link QuantRatio} being the average of the provided ratios
	 *
	 * @param nonInfinityRatios
	 * @return
	 */
	public static QuantRatio getAverageRatio(Set<QuantRatio> nonInfinityRatios, AggregationLevel aggregationLevel) {
		QuantCondition cond1 = null;
		QuantCondition cond2 = null;
		final TDoubleArrayList ratioValues = new TDoubleArrayList();
		if (nonInfinityRatios.size() == 1) {
			return nonInfinityRatios.iterator().next();
		}
		for (final QuantRatio quantRatio : nonInfinityRatios) {
			if (cond1 == null) {
				cond1 = quantRatio.getQuantCondition1();
			} else {
				if (!cond1.equals(quantRatio.getQuantCondition1())) {
					continue;
				}
			}
			if (cond2 == null) {
				cond2 = quantRatio.getQuantCondition2();
			} else {
				if (!cond2.equals(quantRatio.getQuantCondition2())) {
					continue;
				}
			}
			final double log2Ratio = quantRatio.getLog2Ratio(cond1, cond2);
			if (!Double.isInfinite(log2Ratio) && !Double.isNaN(log2Ratio)) {
				ratioValues.add(log2Ratio);
			}
		}
		if (ratioValues.isEmpty()) {
			final CensusRatio ret = new CensusRatio(Double.NaN, false, cond1, cond2, aggregationLevel,
					"Peptide with no PSM ratios");

			return ret;
		}
		final double mean = Maths.mean(ratioValues);
		double stdev = Maths.stddev(ratioValues);
		if (ratioValues.size() == 1) {
			stdev = 0;
		}
		String ratioDescription = nonInfinityRatios.iterator().next().getDescription();
		if (ratioValues.size() > 1) {
			ratioDescription = "Average of " + ratioDescription;
		}
		final CensusRatio ret = new CensusRatio(mean, true, cond1, cond2, aggregationLevel, ratioDescription);
		ret.setCombinationType(CombinationType.AVERAGE);
		ret.setNumMeasurements(ratioValues.size());
		final RatioScore ratioScore = new RatioScore(String.valueOf(stdev), "STDEV",
				"Standard deviation of log2 ratios", "Standard deviation of multiple log2 ratios");
		ret.setRatioScore(ratioScore);
		return ret;
	}

	public static Double getMaxAmountValueByAmountType(Set<Amount> amounts, AmountType amountType) {
		double max = -Double.MAX_VALUE;
		if (amounts != null) {
			for (final Amount quantAmount : amounts) {
				if (quantAmount.getAmountType() == amountType) {
					if (quantAmount.getValue() > max) {
						max = quantAmount.getValue();
					}
				}
			}
		}
		if (max != -Double.MAX_VALUE) {
			return max;
		}
		return null;
	}

	public static QuantRatio getRepresentativeRatio(QuantifiedPSMInterface quantifiedPSM) {
		// RATIO for non singletons and AREA_RATIO for singletons
		// UPDATE ON 28th Nov 2017: Use just AREA_RATIO
		if (quantifiedPSM instanceof QuantifiedPSM) {
			// if (quantifiedPSM.isSingleton()) {
			final QuantRatio ret = getRatioByName(quantifiedPSM, CensusOutParser.AREA_RATIO);
			if (ret != null) {
				return ret;
			}
			return getRatioByName(quantifiedPSM, CensusOutParser.RATIO);

			// } else {
			// return getRatioByName(quantifiedPSM, CensusOutParser.RATIO);
			// }
		} else {
			if (quantifiedPSM instanceof IsobaricQuantifiedPSM) {
				throw new IllegalArgumentException("This shouldnt be an isobaric psm");
			} else {
				// QuantRatio ret = getRatioByName(quantifiedPSM,
				// CensusOutParser.RATIO);
				// if (ret == null) {
				// ret = getRatioByName(quantifiedPSM,
				// CensusOutParser.AREA_RATIO);
				// }
				QuantRatio ret = getRatioByName(quantifiedPSM, CensusOutParser.AREA_RATIO);
				if (ret == null) {
					ret = getRatioByName(quantifiedPSM, CensusOutParser.RATIO);
				}
				if (ret == null) {
					ret = getRatioByName(quantifiedPSM, CensusOutParser.NORM_RATIO);
				}
				return ret;
			}
		}
	}

	public static QuantRatio getRatioValidForIntegrationAnalysis(QuantifiedPSMInterface quantifiedPSM) {
		// RATIO for non singletons and AREA_RATIO for singletons
		// UPDATE ON 28th Nov 2017: Use just AREA_RATIO
		if (quantifiedPSM instanceof QuantifiedPSM) {
			// if (quantifiedPSM.isSingleton()) {
			return getRatioByName(quantifiedPSM, CensusOutParser.AREA_RATIO);
			// } else {
			// return getRatioByName(quantifiedPSM, CensusOutParser.RATIO);
			// }
		} else {
			if (quantifiedPSM instanceof IsobaricQuantifiedPSM) {
				throw new IllegalArgumentException(
						"This shouldnt be an isobaric psm. In that case, there is more than one ratio per psm for the integration analysis.");
			} else {
				// QuantRatio ret = getRatioByName(quantifiedPSM,
				// CensusOutParser.RATIO);
				// if (ret == null) {
				// ret = getRatioByName(quantifiedPSM,
				// CensusOutParser.AREA_RATIO);
				// }
				QuantRatio ret = getRatioByName(quantifiedPSM, CensusOutParser.AREA_RATIO);
				if (ret == null) {
					ret = getRatioByName(quantifiedPSM, CensusOutParser.RATIO);
				}
				if (ret == null) {
					ret = getRatioByName(quantifiedPSM, CensusOutParser.NORM_RATIO);
				}
				return ret;
			}
		}
	}

	public static QuantRatio getRatioByName(HasRatios quantifiedPSMorPeptideOrProtein, String ratioDescription) {
		final List<QuantRatio> ratiosByName = getRatiosByName(quantifiedPSMorPeptideOrProtein, ratioDescription);
		if (ratiosByName.isEmpty()) {
			return null;
		}
		return ratiosByName.get(0);
	}

	public static List<QuantRatio> getRatiosByName(HasRatios quantifiedPSMorPeptideOrProtein, String ratioDescription) {
		if (quantifiedPSMorPeptideOrProtein != null && quantifiedPSMorPeptideOrProtein.getRatios() != null) {
			return getRatiosByName(quantifiedPSMorPeptideOrProtein.getRatios(), ratioDescription);
		}
		return Collections.emptyList();
	}

	public static List<QuantRatio> getRatiosByName(Collection<QuantRatio> quantRatios, String ratioDescription) {
		final List<QuantRatio> ratios = new ArrayList<QuantRatio>();
		if (quantRatios != null) {
			for (final QuantRatio ratio : quantRatios) {
				if (ratio.getDescription().equals(ratioDescription)) {
					ratios.add(ratio);
				}
			}
		}
		return ratios;
	}

	public static List<QuantifiedPeptideInterface> getSortedPeptidesByFullSequence(
			Collection<QuantifiedPeptideInterface> peptides) {
		final List<QuantifiedPeptideInterface> ret = new ArrayList<QuantifiedPeptideInterface>();
		ret.addAll(peptides);
		Collections.sort(ret, new Comparator<QuantifiedPeptideInterface>() {

			@Override
			public int compare(QuantifiedPeptideInterface o1, QuantifiedPeptideInterface o2) {
				return o1.getFullSequence().compareTo(o2.getFullSequence());
			}
		});
		return ret;
	}

	/**
	 * Get a CVS list of peptide sequences after sorting them alphabetically by
	 * the sequence
	 *
	 * @param peptides
	 * @return
	 */
	public static String getPeptidesFullSequenceString(Collection<QuantifiedPeptideInterface> peptides) {

		final StringBuilder sb = new StringBuilder();
		for (final QuantifiedPeptideInterface peptide : QuantUtils.getSortedPeptidesByFullSequence(peptides)) {
			if (!"".equals(sb.toString()))
				sb.append("_");
			sb.append(peptide.getFullSequence());
		}
		return sb.toString();
	}

	public static List<QuantifiedProteinInterface> getSortedQuantifiedProteinsByAcc(
			Collection<QuantifiedProteinInterface> proteinsToSort) {
		final List<QuantifiedProteinInterface> list = new ArrayList<QuantifiedProteinInterface>();
		list.addAll(proteinsToSort);
		Collections.sort(list, new Comparator<QuantifiedProteinInterface>() {
			@Override
			public int compare(QuantifiedProteinInterface o1, QuantifiedProteinInterface o2) {
				return o1.getAccession().compareTo(o2.getAccession());
			}
		});
		return list;
	}

	public static int getIonCount(QuantifiedPSMInterface psm, QuantificationLabel label) {
		if (psm instanceof IsobaricQuantifiedPSM) {
			final IsobaricQuantifiedPSM isoPSM = (IsobaricQuantifiedPSM) psm;
			final Set<Ion> ionsByLabel = isoPSM.getIonsByLabel(label);
			if (ionsByLabel != null) {
				return ionsByLabel.size();
			}
		}
		return 0;
	}

	/**
	 * Removes the peptide from its proteins and all its psms from irs proteins
	 *
	 * @param quantifiedPeptide
	 */
	public static void discardPeptide(QuantifiedPeptideInterface quantifiedPeptide) {

		// get all psms of that peptide
		final Set<QuantifiedPSMInterface> quantifiedPSMs = quantifiedPeptide.getQuantifiedPSMs();
		for (final QuantifiedPSMInterface quantifiedPSM : quantifiedPSMs) {
			// remove this psms from its proteins
			final Iterator<QuantifiedProteinInterface> proteinsIterator = quantifiedPSM.getQuantifiedProteins()
					.iterator();
			while (proteinsIterator.hasNext()) {
				final QuantifiedProteinInterface quantifiedProtein = proteinsIterator.next();
				quantifiedProtein.getQuantifiedPSMs().remove(quantifiedPSM);
				proteinsIterator.remove();
				final int numPSMs = quantifiedProtein.getQuantifiedPSMs().size();
				if (quantifiedProtein.getAccession().equals("J3QTB2")) {
					log.info(quantifiedProtein + " " + numPSMs);
				}
			}
			// remove this psm from the parser
			// parser.getPSMMap().remove(quantifiedPSM.getPSMIdentifier());
		}
	}

	/**
	 * Removes the protein from its peptides and from its psms
	 *
	 * @param quantifiedProtein
	 */
	public static void discardProtein(QuantifiedProteinInterface quantifiedProtein) {

		// remove this protein from its peptides
		final Set<QuantifiedPeptideInterface> quantifiedPeptides = quantifiedProtein.getQuantifiedPeptides();
		for (final QuantifiedPeptideInterface quantifiedPeptide : quantifiedPeptides) {
			quantifiedPeptide.getQuantifiedProteins().remove(quantifiedProtein);
		}
		// remove the protein from its psms
		final Set<QuantifiedPSMInterface> quantifiedPSMs = quantifiedProtein.getQuantifiedPSMs();
		for (final QuantifiedPSMInterface quantifiedPSM : quantifiedPSMs) {
			quantifiedPSM.getQuantifiedProteins().remove(quantifiedProtein);
		}

	}

	/**
	 * Gets a consensus {@link IonCountRatio} from a set of
	 * {@link IsobaricQuantifiedPeptide} where the ions from each
	 * {@link QuantCondition} are pulled together and normalized by the number
	 * of {@link QuantifiedPSMInterface} per {@link QuantifiedPeptideInterface}
	 *
	 * @param isobaricQuantifiedPeptides
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	public static IonCountRatio getNormalizedIonCountRatioForPeptides(
			Set<IsobaricQuantifiedPeptide> isobaricQuantifiedPeptides, QuantCondition cond1, QuantCondition cond2) {
		return getNormalizedIonCountRatioForPeptides(isobaricQuantifiedPeptides, cond1, cond2, null);
	}

	/**
	 * Gets a consensus {@link IonCountRatio} from a set of
	 * {@link IsobaricQuantifiedPeptide} where the ions from each
	 * {@link QuantCondition} are pulled together and normalized by the number
	 * of {@link QuantifiedPSMInterface} per {@link QuantifiedPeptideInterface}
	 *
	 * @param isobaricQuantifiedPeptides
	 * @param cond1
	 * @param cond2
	 * @param replicateName
	 * @return
	 */
	public static IonCountRatio getNormalizedIonCountRatioForPeptides(
			Collection<IsobaricQuantifiedPeptide> isobaricQuantifiedPeptides, QuantCondition cond1,
			QuantCondition cond2, String replicateName) {

		final IonCountRatio ratio = new IonCountRatio(AggregationLevel.PEPTIDE);
		for (final IsobaricQuantifiedPeptide isoPeptide : isobaricQuantifiedPeptides) {
			final int numPSMs = isoPeptide.getQuantifiedPSMs().size();
			// get number of ions in one condition, and normalize by the
			// number of PSMs
			int peakCount1 = 0;
			final Map<QuantCondition, Set<Ion>> ionsByCondition = isoPeptide.getIonsByCondition(replicateName);
			if (ionsByCondition.containsKey(cond1)) {
				peakCount1 = ionsByCondition.get(cond1).size();
			}
			final double normalizedPeakCount1 = peakCount1 * 1.0 / numPSMs;
			int peakCount2 = 0;
			if (ionsByCondition.containsKey(cond2)) {
				peakCount2 = ionsByCondition.get(cond2).size();
			}
			final double normalizedPeakCount2 = peakCount2 * 1.0 / numPSMs;
			ratio.addIonCount(cond1, normalizedPeakCount1);
			ratio.addIonCount(cond2, normalizedPeakCount2);

		}
		ratio.setAsNormalizedIonCountRatio();
		return ratio;
	}

	/**
	 * Gets a consensus {@link IonCountRatio} from a set of
	 * {@link IsobaricQuantifiedPeptide} where the ions from each
	 * {@link QuantCondition} are pulled together and normalized by the number
	 * of {@link QuantifiedPSMInterface} per
	 * {@link QuantifiedPeptideInterface}.<br>
	 * The ions are used are only the ones that can distinguish the particular
	 * aminoacids that are quantified (aas)
	 *
	 * @param isobaricQuantifiedPeptides
	 * @param positionsInProteins
	 * @param cond1
	 * @param cond2
	 * @param replicateName
	 * @param quantifiedAAs
	 * @return
	 */
	public static IonCountRatio getNormalizedIonCountRatioForPeptidesForQuantifiedSites(
			Collection<Pair<IsobaricQuantifiedPeptide, List<PositionInPeptide>>> peptidesAndPositionInPeptides,
			QuantCondition cond1, QuantCondition cond2, String replicateName, char[] quantifiedAAs) {

		final IonCountRatio ratio = new IonCountRatio(AggregationLevel.PEPTIDE);
		for (final Pair<IsobaricQuantifiedPeptide, List<PositionInPeptide>> peptideAndPositionInPeptide : peptidesAndPositionInPeptides) {
			final IsobaricQuantifiedPeptide isoPeptide = peptideAndPositionInPeptide.getFirstelement();
			for (final PositionInPeptide positionInPeptideObj : peptideAndPositionInPeptide.getSecondElement()) {

				final int positionInPeptide = positionInPeptideObj.getPosition();
				final int numPSMs = isoPeptide.getQuantifiedPSMs().size();
				// get number of ions in one condition, and normalize by the
				// number of PSMs
				int peakCount1 = 0;
				final Map<QuantCondition, Set<Ion>> ionsByConditionForSites = isoPeptide
						.getIonsByConditionForSites(replicateName, quantifiedAAs, positionInPeptide);
				if (ionsByConditionForSites.containsKey(cond1)) {
					peakCount1 = ionsByConditionForSites.get(cond1).size();
				}
				final double normalizedPeakCount1 = peakCount1 * 1.0 / numPSMs;
				int peakCount2 = 0;
				if (ionsByConditionForSites.containsKey(cond2)) {
					peakCount2 = ionsByConditionForSites.get(cond2).size();
				}
				final double normalizedPeakCount2 = peakCount2 * 1.0 / numPSMs;
				ratio.addIonCount(cond1, normalizedPeakCount1);
				ratio.addIonCount(cond2, normalizedPeakCount2);
			}
		}
		ratio.setAsNormalizedIonCountRatio();
		return ratio;

	}

	/**
	 * Gets a consensus {@link IonCountRatio} from a set of
	 * {@link IsobaricQuantifiedPeptide} where the ions from each
	 * {@link QuantCondition} are pulled together and normalized by the number
	 * of {@link QuantifiedPSMInterface} per {@link QuantifiedPeptideInterface}
	 *
	 * @param isobaricQuantifiedPeptides
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	public static IonCountRatio getIonCountRatioForPeptide(IsobaricQuantifiedPeptide isoPeptide, QuantCondition cond1,
			QuantCondition cond2) {
		return getIonCountRatioForPeptide(isoPeptide, cond1, cond2, null);
	}

	/**
	 * Gets a consensus {@link IonCountRatio} from a set of
	 * {@link IsobaricQuantifiedPeptide} where the ions from each
	 * {@link QuantCondition} are pulled together and normalized by the number
	 * of {@link QuantifiedPSMInterface} per {@link QuantifiedPeptideInterface}
	 *
	 * @param isobaricQuantifiedPeptides
	 * @param cond1
	 * @param cond2
	 * @return
	 */
	public static IonCountRatio getIonCountRatioForPeptide(IsobaricQuantifiedPeptide isoPeptide, QuantCondition cond1,
			QuantCondition cond2, String replicateName) {

		final IonCountRatio ratio = new IonCountRatio(AggregationLevel.PEPTIDE);

		// get number of ions in one condition, and normalize by the
		// number of PSMs
		int peakCount1 = 0;
		if (isoPeptide.getIonsByCondition().containsKey(cond1)) {
			peakCount1 = isoPeptide.getIonsByCondition(replicateName).get(cond1).size();
		}

		int peakCount2 = 0;
		if (isoPeptide.getIonsByCondition().containsKey(cond2)) {
			peakCount2 = isoPeptide.getIonsByCondition(replicateName).get(cond2).size();
		}

		ratio.addIonCount(cond1, peakCount1);
		ratio.addIonCount(cond2, peakCount2);

		return ratio;
	}

	public static List<PTMInProtein> getPTMPositionsInProtein(String accession, QuantifiedPeptideInterface peptide,
			String uniprotVersion, File uniprotAnnotationsFolder) {
		final List<PTMInProtein> ptmPositionsInProtein = new ArrayList<PTMInProtein>();
		final UniprotProteinRetriever uplr2 = getUniprotProteinLocalRetriever(uniprotVersion, uniprotAnnotationsFolder);
		// get protein sequence from uniprot
		final Map<String, String> annotatedProteinSequence = uplr2.getAnnotatedProteinSequence(accession);
		if (annotatedProteinSequence.containsKey(accession)) {
			final String proteinSequence = annotatedProteinSequence.get(accession);
			final TIntArrayList positionsInProteinSequence = StringUtils.allPositionsOf(proteinSequence,
					peptide.getSequence());
			if (!positionsInProteinSequence.isEmpty()) {
				for (final int startingPosition : positionsInProteinSequence.toArray()) {
					final List<PTMInPeptide> ptms = peptide.getPtms();
					for (final PTMInPeptide ptmInPeptide : ptms) {
						final int positionInPeptide = ptmInPeptide.getPosition();
						final PTMInProtein ptmInProtein = new PTMInProtein(positionInPeptide + startingPosition - 1,
								peptide.getSequence().charAt(positionInPeptide - 1), accession,
								ptmInPeptide.getDeltaMass());
						ptmPositionsInProtein.add(ptmInProtein);
					}
				}

			}
		}

		return ptmPositionsInProtein;
	}

	private static UniprotProteinRetriever getUniprotProteinLocalRetriever(String uniprotVersion,
			File uniprotAnnotationsFolder) {
		if (uplr == null) {
			uplr = new UniprotProteinRetriever(uniprotVersion, uniprotAnnotationsFolder, true);
		}
		return uplr;
	}

	public static double parseCensusRatioValue(String stringValue) {
		if (stringValue == null) {
			return Double.NaN;
		}
		if (stringValue.toUpperCase().contains("INF")) {
			if (stringValue.toUpperCase().contains("-")) {
				return Double.NEGATIVE_INFINITY;
			} else {
				return Double.POSITIVE_INFINITY;
			}
		}
		if (stringValue.equals("1000")) {
			return Double.POSITIVE_INFINITY;
		}
		if (stringValue.equals("0.001")) {
			return Double.NEGATIVE_INFINITY;
		}

		try {
			return Double.valueOf(stringValue);
		} catch (final NumberFormatException e) {
			log.warn("Census ratio is not recognized: '" + stringValue + "'");
			return Double.NaN;
		}
	}

	public static Set<IsoRatio> getIsobaricRatiosForSiteFromPeptide(IsobaricQuantifiedPeptide peptide,
			int positionInPeptide) {
		final Set<IsoRatio> ret = new HashSet<IsoRatio>();

		final Set<IsoRatio> individualIsoRatios = peptide.getIsoRatios();
		for (final IsoRatio isoRatio : individualIsoRatios) {
			if (containsPosition(isoRatio.getQuantifiedSitePositionInPeptide(), positionInPeptide)) {
				ret.add(isoRatio);
			}
		}
		return ret;
	}

	public static boolean containsPosition(Collection<PositionInPeptide> positionInPeptides, int position) {
		if (positionInPeptides != null) {
			for (final PositionInPeptide positionInPeptide : positionInPeptides) {
				if (positionInPeptide.getPosition() == position) {
					return true;
				}
			}
		}
		return false;
	}

	public static String printPositionsInPeptideInOrder(Collection<PositionInPeptide> positionsInPeptides) {
		final List<PositionInPeptide> list = new ArrayList<PositionInPeptide>();
		list.addAll(positionsInPeptides);
		list.sort(new Comparator<PositionInPeptide>() {

			@Override
			public int compare(PositionInPeptide o1, PositionInPeptide o2) {
				final int ret = o1.getProteinACC().compareTo(o2.getProteinACC());
				if (ret == 0) {
					return Integer.compare(o1.getPosition(), o2.getPosition());
				} else {
					return ret;
				}
			}
		});
		final StringBuilder sb = new StringBuilder();
		for (final PositionInPeptide positionInPeptide : list) {
			if (!"".equals(sb.toString())) {
				sb.append("-");
			}
			sb.append(positionInPeptide.toString());
		}
		return sb.toString();
	}
}
