package edu.scripps.yates.census.read.util;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
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
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPSM;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedPeptide;
import edu.scripps.yates.census.read.model.QuantifiedPSM;
import edu.scripps.yates.census.read.model.QuantifiedPeptide;
import edu.scripps.yates.census.read.model.RatioScore;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.model.enums.AmountType;
import edu.scripps.yates.utilities.model.enums.CombinationType;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.strings.StringUtils;
import edu.scripps.yates.utilities.util.StringPosition;

public class QuantUtils {
	private static final Logger log = Logger.getLogger(QuantUtils.class);
	private static UniprotProteinRetriever uplr;

	public static void addToPeptideMap(QuantifiedPSMInterface quantifiedPSM, Map<String, QuantifiedPeptide> map) {
		final String sequenceKey = KeyUtils.getSequenceKey(quantifiedPSM, true);
		QuantifiedPeptide quantifiedPeptide = null;
		if (map.containsKey(sequenceKey)) {
			quantifiedPeptide = map.get(sequenceKey);
			quantifiedPeptide.addQuantifiedPSM(quantifiedPSM, true);
		} else {
			quantifiedPeptide = new QuantifiedPeptide(quantifiedPSM);
			map.put(sequenceKey, quantifiedPeptide);
		}

		quantifiedPSM.setQuantifiedPeptide(quantifiedPeptide, true);
	}

	public static void addToIsobaricPeptideMap(IsobaricQuantifiedPSM quantifiedPSM,
			Map<String, IsobaricQuantifiedPeptide> map) {
		final String sequenceKey = KeyUtils.getSequenceKey(quantifiedPSM, true);
		IsobaricQuantifiedPeptide quantifiedPeptide = null;
		if (map.containsKey(sequenceKey)) {
			quantifiedPeptide = map.get(sequenceKey);
			quantifiedPeptide.addPSM(quantifiedPSM);
		} else {
			quantifiedPeptide = new IsobaricQuantifiedPeptide(quantifiedPSM);
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
		Map<String, IsobaricQuantifiedPeptide> peptideMap = new HashMap<String, IsobaricQuantifiedPeptide>();

		for (IsobaricQuantifiedPSM quantifiedPSM : quantifiedPSMs) {
			final QuantifiedPeptideInterface quantifiedPeptide = quantifiedPSM.getQuantifiedPeptide();
			if (quantifiedPeptide instanceof IsobaricQuantifiedPeptide) {
				if (!peptideMap.containsKey(quantifiedPeptide.getKey())) {
					peptideMap.put(quantifiedPeptide.getKey(), (IsobaricQuantifiedPeptide) quantifiedPeptide);
				}
			}
		}

		return peptideMap;
	}

	public static Set<QuantRatio> getNonInfinityRatios(Set<QuantRatio> ratios) {
		Set<QuantRatio> set = new HashSet<QuantRatio>();
		for (QuantRatio ratio : ratios) {
			final double log2Ratio = ratio.getLog2Ratio(ratio.getQuantCondition1(), ratio.getQuantCondition2());
			if (Double.compare(log2Ratio, Double.MAX_VALUE) == 0 || Double.compare(log2Ratio, -Double.MAX_VALUE) == 0) {
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
		List<Double> ratioValues = new ArrayList<Double>();
		for (QuantRatio quantRatio : nonInfinityRatios) {
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
			final Double log2Ratio = quantRatio.getLog2Ratio(cond1, cond2);
			if (log2Ratio != null && !Double.isInfinite(log2Ratio) && !Double.isNaN(log2Ratio)) {
				ratioValues.add(log2Ratio);
			}
		}
		if (ratioValues.isEmpty()) {
			CensusRatio ret = new CensusRatio(Double.NaN, false, cond1, cond2, aggregationLevel,
					"Peptide with no PSM ratios");
			return ret;
		}
		final Double[] ratioValuesArray = ratioValues.toArray(new Double[0]);
		final double mean = Maths.mean(ratioValuesArray);
		final double stdev = Maths.stddev(ratioValuesArray);
		CensusRatio ret = new CensusRatio(mean, true, cond1, cond2, aggregationLevel, "Average of ratios");
		ret.setCombinationType(CombinationType.AVERAGE);
		RatioScore ratioScore = new RatioScore(String.valueOf(stdev), "STDEV", "Standard deviation of log2 ratios",
				"Standard deviation of multiple log2 ratios");
		ret.setRatioScore(ratioScore);
		return ret;
	}

	/**
	 * Gets a {@link Set} of {@link QuantRatio} as consensus ratios of the
	 * {@link Collection}of {@link QuantifiedPeptideInterface}, calling to
	 * getConsensusRatio() in each {@link QuantifiedPeptideInterface}
	 *
	 * @param peptides
	 * @param cond1
	 * @param cond2
	 * @param replicateName
	 *            if not null, get only the consensus {@link QuantRatio} for
	 *            that replicate
	 * @return
	 */
	public static Set<QuantRatio> getConsensusRatios(Collection<QuantifiedPeptideInterface> peptides,
			QuantCondition cond1, QuantCondition cond2, String replicateName) {
		Set<QuantRatio> ret = new HashSet<QuantRatio>();
		for (QuantifiedPeptideInterface peptide : peptides) {
			if (replicateName != null) {
				final QuantRatio consensusRatio = peptide.getConsensusRatio(cond1, cond2, replicateName);
				if (consensusRatio != null) {
					ret.add(consensusRatio);
				}
			} else {
				final QuantRatio consensusRatio = peptide.getConsensusRatio(cond1, cond2);
				if (consensusRatio != null) {
					ret.add(consensusRatio);
				}
			}
		}
		return ret;
	}

	public static Double getMaxAmountValueByAmountType(Set<Amount> amounts, AmountType amountType) {
		double max = -Double.MAX_VALUE;
		if (amounts != null) {
			for (Amount quantAmount : amounts) {
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

	/**
	 * For a {@link QuantifiedPSM}, depending on if it is singleton or not, we
	 * get a different ratio.<br>
	 * For singletons, we get the AREA_RATIO and for non singletons, we get the
	 * RATIO
	 *
	 * @param quantifiedPSM
	 * @return
	 */
	public static QuantRatio getRatioValidForAnalysis(QuantifiedPSM quantifiedPSM) {

		// RATIO for non singletons and AREA_RATIO for singletons
		if (quantifiedPSM.isSingleton()) {
			return getRatioByName(quantifiedPSM, CensusOutParser.AREA_RATIO);
		} else {
			return getRatioByName(quantifiedPSM, CensusOutParser.RATIO);
		}
	}

	private static QuantRatio getRatioByName(QuantifiedPSM quantifiedPSM, String ratioDescription) {
		if (quantifiedPSM != null && quantifiedPSM.getRatios() != null) {
			for (QuantRatio ratio : quantifiedPSM.getRatios()) {
				if (ratio.getDescription().equals(ratioDescription)) {
					return ratio;
				}
			}
		}
		return null;
	}

	public static List<QuantifiedPeptideInterface> getSortedPeptidesByFullSequence(
			Collection<QuantifiedPeptideInterface> peptides) {
		List<QuantifiedPeptideInterface> ret = new ArrayList<QuantifiedPeptideInterface>();
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

		StringBuilder sb = new StringBuilder();
		for (QuantifiedPeptideInterface peptide : QuantUtils.getSortedPeptidesByFullSequence(peptides)) {
			if (!"".equals(sb.toString()))
				sb.append("_");
			sb.append(peptide.getFullSequence());
		}
		return sb.toString();
	}

	public static List<QuantifiedProteinInterface> getSortedQuantifiedProteinsByAcc(
			Collection<QuantifiedProteinInterface> proteinsToSort) {
		List<QuantifiedProteinInterface> list = new ArrayList<QuantifiedProteinInterface>();
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
			IsobaricQuantifiedPSM isoPSM = (IsobaricQuantifiedPSM) psm;
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
		for (QuantifiedPSMInterface quantifiedPSM : quantifiedPSMs) {
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
		for (QuantifiedPeptideInterface quantifiedPeptide : quantifiedPeptides) {
			quantifiedPeptide.getQuantifiedProteins().remove(quantifiedProtein);
		}
		// remove the protein from its psms
		final Set<QuantifiedPSMInterface> quantifiedPSMs = quantifiedProtein.getQuantifiedPSMs();
		for (QuantifiedPSMInterface quantifiedPSM : quantifiedPSMs) {
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
			Set<IsobaricQuantifiedPeptide> isobaricQuantifiedPeptides, QuantCondition cond1, QuantCondition cond2,
			String replicateName) {

		IonCountRatio ratio = new IonCountRatio(AggregationLevel.PEPTIDE);
		for (IsobaricQuantifiedPeptide isoPeptide : isobaricQuantifiedPeptides) {
			final int numPSMs = isoPeptide.getQuantifiedPSMs().size();
			// get number of ions in one condition, and normalize by the
			// number of PSMs
			int peakCount1 = 0;
			if (isoPeptide.getIonsByCondition(replicateName).containsKey(cond1)) {
				peakCount1 = isoPeptide.getIonsByCondition(replicateName).get(cond1).size();
			}
			double normalizedPeakCount1 = peakCount1 * 1.0 / numPSMs;
			int peakCount2 = 0;
			if (isoPeptide.getIonsByCondition(replicateName).containsKey(cond2)) {
				peakCount2 = isoPeptide.getIonsByCondition(replicateName).get(cond2).size();
			}
			double normalizedPeakCount2 = peakCount2 * 1.0 / numPSMs;
			ratio.addIonCount(cond1, normalizedPeakCount1);
			ratio.addIonCount(cond2, normalizedPeakCount2);

		}
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

		IonCountRatio ratio = new IonCountRatio(AggregationLevel.PEPTIDE);

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

	public static List<StringPosition> getPTMPositionsInProtein(String accession, QuantifiedPeptideInterface peptide,
			String uniprotVersion, File uniprotAnnotationsFolder) {
		List<StringPosition> ptmPositionsInProtein = new ArrayList<StringPosition>();
		UniprotProteinRetriever uplr2 = getUniprotProteinLocalRetriever(uniprotVersion, uniprotAnnotationsFolder);
		// get protein sequence from uniprot
		final Map<String, String> annotatedProteinSequence = uplr2.getAnnotatedProteinSequence(accession);
		if (annotatedProteinSequence.containsKey(accession)) {
			final String proteinSequence = annotatedProteinSequence.get(accession);
			List<Integer> positionsInProteinSequence = StringUtils.allPositionsOf(proteinSequence,
					peptide.getSequence());
			if (!positionsInProteinSequence.isEmpty()) {
				for (Integer startingPosition : positionsInProteinSequence) {
					final List<StringPosition> ptms = peptide.getPtms();
					for (StringPosition ptmPositionInPeptide : ptms) {
						final int positionInPeptide = ptmPositionInPeptide.position;
						StringPosition ptmPositionInProtein = new StringPosition(ptmPositionInPeptide.string,
								positionInPeptide + startingPosition - 1);
						ptmPositionsInProtein.add(ptmPositionInProtein);
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
}