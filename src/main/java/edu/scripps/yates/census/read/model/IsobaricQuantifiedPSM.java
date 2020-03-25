package edu.scripps.yates.census.read.model;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.quant.xml.ProteinType.Peptide;
import edu.scripps.yates.census.quant.xml.ProteinType.Peptide.Frag;
import edu.scripps.yates.census.read.AbstractQuantParser;
import edu.scripps.yates.census.read.model.IonSerie.IonSerieType;
import edu.scripps.yates.census.read.model.interfaces.HasIsoRatios;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.IonExclusion;
import edu.scripps.yates.census.read.util.QuantKeyUtils;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.proteomicsmodel.AbstractPSM;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.Ratio;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AggregationLevel;
import edu.scripps.yates.utilities.proteomicsmodel.factories.AmountEx;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import edu.scripps.yates.utilities.strings.StringUtils;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.THashSet;

public class IsobaricQuantifiedPSM extends AbstractPSM implements QuantifiedPSMInterface, HasIsoRatios {
	private static final Logger log = Logger.getLogger(IsobaricQuantifiedPSM.class);
	private final Peptide censusChroPeptide;
	private final Set<edu.scripps.yates.census.read.util.QuantificationLabel> labels = new THashSet<QuantificationLabel>();
	private IonSerie serieYHeavy;
	private IonSerie serieYLight;
	private IonSerie serieBHeavy;
	private IonSerie serieBLight;
	private List<IsoRatio> ratiosSerieY;
	private List<IsoRatio> ratiosSerieB;
	private Set<QuantRatio> quantRatios;
	private final Collection<IonExclusion> ionExclusions;
	private static int scanNum = 0;
	private final Map<QuantCondition, QuantificationLabel> labelsByConditions;
	private final Map<QuantificationLabel, QuantCondition> conditionsByLabels;
	private final Map<String, IonCountRatio> countRatiosByConditionKey = new THashMap<String, IonCountRatio>();
	private final Set<String> fileNames = new THashSet<String>();
	private boolean discarded;
	private final Set<Character> quantifiedSites;

	/**
	 *
	 * @param peptide
	 * @param labelNumerator      for writing data file. Otherwise, it can be null
	 * @param labelDenominator    for writing data file. Otherwise it can be null
	 * @param spectrumToIonsMap
	 * @param peptideToSpectraMap
	 * @param ionKeys
	 * @param ionExclusions
	 * @param dataFileWriter
	 * @param chargeStateSensible
	 * @param string
	 * @throws IOException
	 */
	public IsobaricQuantifiedPSM(Peptide peptide, Map<QuantCondition, QuantificationLabel> labelsByConditions,
			Collection<IonExclusion> ionExclusions, Set<Character> quantifiedSites, boolean distinguishModifiedSequence,
			boolean chargeStateSensible) throws IOException {
		super(distinguishModifiedSequence, chargeStateSensible);
		censusChroPeptide = peptide;
		this.ionExclusions = ionExclusions;
		setScanNumber(peptide.getScan());
		setSequence(peptide.getSeq());
		this.labelsByConditions = labelsByConditions;
		this.quantifiedSites = quantifiedSites;
		conditionsByLabels = new THashMap<QuantificationLabel, QuantCondition>();
		for (final QuantCondition condition : labelsByConditions.keySet()) {
			final QuantificationLabel quantificationLabel = labelsByConditions.get(condition);
			conditionsByLabels.put(quantificationLabel, condition);
		}

		process();
		setKey(QuantKeyUtils.getInstance().getSpectrumKey(this, distinguishModifiedSequence, chargeStateSensible));
	}

	/**
	 * @return the rawfileName
	 */
	@Override
	public Set<String> getRawFileNames() {
		final Set<String> ret = new THashSet<String>();
		if (censusChroPeptide != null) {
			ret.add(censusChroPeptide.getFile());
		}
		return ret;
	}

	private void process() throws IOException {

		final Frag frag = censusChroPeptide.getFrag();
		final String yr = frag.getYr();
		final String ys = frag.getYs();

		// SERIE Y
		serieYHeavy = new IonSerie(QuantificationLabel.HEAVY, IonSerieType.Y, yr, ionExclusions);
		serieYLight = new IonSerie(QuantificationLabel.LIGHT, IonSerieType.Y, ys, ionExclusions);
		// check the ions and remove the ones that has the same intensities in
		// the two labels, which means that cannot be distinguished
		checkIons(serieYLight, serieYHeavy);
		ratiosSerieY = new ArrayList<IsoRatio>();
		if (quantifiedSites != null && !quantifiedSites.isEmpty()) {
			// ratios are referring to a particular site
			final TIntObjectHashMap<Set<IsoRatio>> siteSpecificRatios = getSiteSpecificRatiosSeriesY(quantifiedSites);
			for (final Set<IsoRatio> isoRatios : siteSpecificRatios.valueCollection()) {
				for (final IsoRatio isoRatio : isoRatios) {
					ratiosSerieY.add(isoRatio);
					addRatio(isoRatio);
				}
			}
		} else {
			// ratios are abundance of the peptide
			ratiosSerieY = getRatiosFromSeries(serieYLight, serieYHeavy);
			for (final IsoRatio isoRatio : ratiosSerieY) {
				addRatio(isoRatio);
			}
		}
		// SERIE B
		final String br = frag.getBr();
		final String bs = frag.getBs();
		serieBHeavy = new IonSerie(QuantificationLabel.HEAVY, IonSerieType.B, br, ionExclusions);
		serieBLight = new IonSerie(QuantificationLabel.LIGHT, IonSerieType.B, bs, ionExclusions);
		// check the ions and remove the ones that has the same intensities in
		// the two labels, which means that cannot be distinguished
		checkIons(serieBLight, serieBHeavy);
		ratiosSerieB = new ArrayList<IsoRatio>();
		if (quantifiedSites != null && !quantifiedSites.isEmpty()) {
			// ratios are referring to a particular site
			final TIntObjectHashMap<Set<IsoRatio>> siteSpecificRatios = getSiteSpecificRatiosSeriesB(quantifiedSites);
			for (final Set<IsoRatio> isoRatios : siteSpecificRatios.valueCollection()) {
				for (final IsoRatio isoRatio : isoRatios) {
					ratiosSerieB.add(isoRatio);
					addRatio(isoRatio);
				}
			}
		} else {
			// ratios are abundance of the peptide
			ratiosSerieB = getRatiosFromSeries(serieBLight, serieBHeavy);
			for (final IsoRatio isoRatio : ratiosSerieB) {
				addRatio(isoRatio);
			}
		}

		// create ion amounts
		createIonAmounts();
	}

	public void addSpectrumToIonsMaps(String spectrumKey, Map<String, Set<String>> spectrumToIonsMap,
			Set<String> ionKeys) {
		for (final IsoRatio ratio : ratiosSerieY) {
			addRatio(ratio);
			final String ionKey = QuantKeyUtils.getInstance().getIonKey(ratio, censusChroPeptide, true);
			if (ionKeys.contains(ionKey)) {
				continue;
			}
			AbstractQuantParser.addToMap(spectrumKey, spectrumToIonsMap, ionKey);
		}
		for (final IsoRatio ratio : ratiosSerieB) {
			addRatio(ratio);
			final String ionKey = QuantKeyUtils.getInstance().getIonKey(ratio, censusChroPeptide, true);
			if (ionKeys.contains(ionKey)) {
				continue;
			}
			AbstractQuantParser.addToMap(spectrumKey, spectrumToIonsMap, ionKey);
		}
	}

	private void createIonAmounts() {
		final Set<QuantCondition> conditions = getIonsByCondition().keySet();
		for (final QuantCondition condition : conditions) {
			final Set<Ion> ions = getIonsByCondition().get(condition);
			if (ions != null) {
				for (final Ion ion : ions) {
					final AmountEx amount = new AmountEx(ion.getIntensity(),
							edu.scripps.yates.utilities.proteomicsmodel.enums.AmountType.INTENSITY, condition);
					// singleton or not
					final boolean singleton = ion.isSingleton();
					amount.setSingleton(singleton);
					addAmount(amount);
				}
			}
		}

	}

	private void checkIons(IonSerie lightSerie, IonSerie heavySerie) {
		final TIntObjectHashMap<Ion> lightMap = lightSerie.getIonMap();
		final TIntObjectHashMap<Ion> heavyMap = heavySerie.getIonMap();
		int max = lightSerie.getMaxNumberIon();
		if (heavySerie.getMaxNumberIon() > max)
			max = heavySerie.getMaxNumberIon();

		for (int numIon = 1; numIon <= max; numIon++) {
			final Ion lightIon = lightMap.get(numIon);
			final Ion heavyIon = heavyMap.get(numIon);
			if (lightIon != null && heavyIon != null) {
				if (lightIon.getIntensity() == heavyIon.getIntensity()) {
					lightSerie.removeIon(numIon);
					heavySerie.removeIon(numIon);
				}
			}
		}
	}

	public String getChro() {

		return censusChroPeptide.getChro();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.quant.xml.RelexChro.Protein.Peptide#getFrag()
	 */

	public Frag getFrag() {

		return censusChroPeptide.getFrag();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.scripps.yates.census.quant.xml.RelexChro.Protein.Peptide#getUnique()
	 */

	public String getUnique() {

		return censusChroPeptide.getUnique();
	}

	@Override
	public Float getXCorr() {
		return censusChroPeptide.getXcorr();
	}

	@Override
	public Float getCalcMH() {
		return censusChroPeptide.getCalcMHplus();
	}

	@Override
	public Float getExperimentalMH() {
		return censusChroPeptide.getMHplus();
	}

	@Override
	public Float getTotalIntensity() {
		return censusChroPeptide.getTotalIntensity();
	}

	@Override
	public Integer getSpr() {
		return censusChroPeptide.getSpRank();
	}

	public Float getSpScore() {
		return censusChroPeptide.getSpScore();
	}

	public Integer getRedundancy() {
		return censusChroPeptide.getRedundancy();
	}

	@Override
	public Float getDeltaCn() {

		return censusChroPeptide.getDeltaCN();
	}

	@Override
	public Float getMassErrorPPM() {
		return censusChroPeptide.getDeltaMass();
	}

	@Override
	public Integer getChargeState() {

		return censusChroPeptide.getCharge();
	}

	public Integer getSpectrumCounts() {

		return censusChroPeptide.getSpC();
	}

	@Override
	public boolean addQuantifiedProtein(QuantifiedProteinInterface quantifiedProtein, boolean recursive) {
		final boolean ret = addProtein(quantifiedProtein, recursive);
		// get the taxonomy
		final Set<String> taxonomies = quantifiedProtein.getTaxonomies();
		if (taxonomies != null) {
			for (final String taxonomy : taxonomies) {
				addTaxonomy(taxonomy);
			}
		}
		return ret;
	}

	/**
	 * @return the quantifiedProteins
	 */
	@Override
	public Set<QuantifiedProteinInterface> getQuantifiedProteins() {
		final Set<QuantifiedProteinInterface> ret = new THashSet<QuantifiedProteinInterface>();
		for (final Protein protein : getProteins()) {
			if (protein instanceof QuantifiedProteinInterface) {
				ret.add((QuantifiedProteinInterface) protein);
			}
		}
		return ret;
	}

	/**
	 * Returns true if the peptide contains any {@link Ion} labeled with a certain
	 * {@link QuantificationLabel} and not composing any ratio.
	 *
	 * @param label
	 * @return
	 */
	@Override
	public boolean containsAnySingletonIon(QuantificationLabel label) {
		return containsAnySingletonIon(label, serieBHeavy) || containsAnySingletonIon(label, serieBLight)
				|| containsAnySingletonIon(label, serieYHeavy) || containsAnySingletonIon(label, serieYLight);
	}

	/**
	 * Returns true if the peptide contains any {@link Ion} labeled with a certain
	 * {@link QuantificationLabel}, not matter if they are composing any ratio or
	 * not.
	 *
	 * @param label
	 * @return
	 */
	@Override
	public boolean containsAnyIon(QuantificationLabel label) {
		return containsAnyIon(label, serieBHeavy) || containsAnyIon(label, serieBLight)
				|| containsAnyIon(label, serieYHeavy) || containsAnyIon(label, serieYLight);
	}

	/**
	 * Returns true if the peptide contains any {@link Ion} labeled with a certain
	 * {@link QuantificationLabel}, no matter if they are composing any ratio or
	 * not.
	 *
	 * @param label
	 * @return
	 */
	private boolean containsAnyIon(QuantificationLabel label, IonSerie serie) {
		if (serie != null) {
			return label.equals(serie.getNonNullLabel()) && !serie.getNonNullIons().isEmpty();
		}
		return false;
	}

	/**
	 * Returns true if the peptide contains any {@link Ion} labeled with a certain
	 * {@link QuantificationLabel} and not composing any ratio.
	 *
	 * @param label
	 * @return
	 */
	private boolean containsAnySingletonIon(QuantificationLabel label, IonSerie serie) {
		if (serie != null) {
			return label.equals(serie.getNonNullLabel()) && serie.isSingletonLabeled();
		}
		return false;
	}

	private List<IsoRatio> getRatiosFromSeries(IonSerie serier, IonSerie series) {

		final List<IsoRatio> ret = new ArrayList<IsoRatio>();
		if (series != null) {
			int max = serier.getMaxNumberIon();
			if (series.getMaxNumberIon() > max)
				max = series.getMaxNumberIon();
			for (int ionNumber = 1; ionNumber <= max; ionNumber++) {
				final Ion ionr = serier.getIon(ionNumber);
				final Ion ions = series.getIon(ionNumber);

				// ignore the ions in ionExclusion
				boolean excludeThisIon = false;
				if (ionExclusions != null) {
					for (final IonExclusion ionExlusion : ionExclusions) {
						if (ionExlusion.getIonNumber() == ionNumber) {
							if (ionExlusion.getIonSerieType() == serier.getIonSerieType()) {
								excludeThisIon = true;
							}
							if (ionExlusion.getIonSerieType() == series.getIonSerieType()) {
								excludeThisIon = true;
							}
						}
					}
				}
				if (excludeThisIon) {
					continue;
				}

				// if (ionr != null && ions != null) {
				// ignore if they are not null but they are the same. That is a
				// Not
				// Determined peak (ND)

				if (ionr != null && ions != null && ionr.getIntensity() == ions.getIntensity()) {
					continue;
				}
				// create the ratio if at least one of the is not null
				if (ionr == null && ions == null)
					continue;

				final QuantificationLabel label1 = series.getNonNullLabel();
				final QuantificationLabel label2 = serier.getNonNullLabel();
				final QuantCondition condition1 = conditionsByLabels.get(label1);
				final QuantCondition condition2 = conditionsByLabels.get(label2);
				ret.add(new IsoRatio(ions, label1, condition1, ionr, label2, condition2, ionNumber,
						serier.getIonSerieType(), AggregationLevel.PSM));

			}
		}
		return ret;
	}

	/**
	 * @return the ratios
	 */
	@Override
	public Set<IsoRatio> getIsoRatios() {
		final Set<IsoRatio> isoRatios = new THashSet<IsoRatio>();
		for (final QuantRatio ratio : getQuantRatios()) {
			if (ratio instanceof IsoRatio) {
				isoRatios.add((IsoRatio) ratio);
			}
		}
		return isoRatios;
	}

	/**
	 * @return the ratios
	 */
	@Override
	public Set<IsoRatio> getNonInfinityIsoRatios() {
		final Set<IsoRatio> ret = new THashSet<IsoRatio>();
		for (final IsoRatio isoRatio : getIsoRatios()) {
			final Double log2Ratio = isoRatio.getLog2Ratio(isoRatio.getLabel1(), isoRatio.getLabel2());
			if (Double.isNaN(log2Ratio) || Double.isInfinite(log2Ratio)
					|| Double.compare(log2Ratio, Double.MAX_VALUE) == 0
					|| Double.compare(log2Ratio, -Double.MAX_VALUE) == 0) {
				continue;
			}
			ret.add(isoRatio);
		}
		return ret;
	}

	/**
	 * @return the singletonIons
	 */
	@Override
	public Map<QuantificationLabel, Set<Ion>> getSingletonIonsByLabel() {
		final Map<QuantificationLabel, Set<Ion>> ret = new THashMap<QuantificationLabel, Set<Ion>>();
		final Map<QuantificationLabel, Set<Ion>> singletonIons = getSingletonIons(serieBHeavy);
		// if (!singletonIons.isEmpty())
		ret.putAll(singletonIons);
		final Map<QuantificationLabel, Set<Ion>> singletonIons2 = getSingletonIons(serieBLight);
		// if (!singletonIons2.isEmpty()) {
		addToMapByLabel(ret, singletonIons2);
		// }
		final Map<QuantificationLabel, Set<Ion>> singletonIons3 = getSingletonIons(serieYHeavy);
		// if (!singletonIons3.isEmpty()) {
		addToMapByLabel(ret, singletonIons3);
		// }
		final Map<QuantificationLabel, Set<Ion>> singletonIons4 = getSingletonIons(serieYLight);
		// if (!singletonIons4.isEmpty()) {
		addToMapByLabel(ret, singletonIons4);
		// }
		return ret;
	}

	/**
	 * @return the singletonIons
	 */
	@Override
	public Map<QuantificationLabel, Set<Ion>> getIonsByLabel() {
		final Map<QuantificationLabel, Set<Ion>> ret = new THashMap<QuantificationLabel, Set<Ion>>();
		final Map<QuantificationLabel, Set<Ion>> singletonIons = getIons(serieBHeavy);
		if (!singletonIons.isEmpty())
			addToMapByLabel(ret, singletonIons);
		final Map<QuantificationLabel, Set<Ion>> singletonIons2 = getIons(serieBLight);
		if (!singletonIons2.isEmpty()) {
			addToMapByLabel(ret, singletonIons2);
		}
		final Map<QuantificationLabel, Set<Ion>> singletonIons3 = getIons(serieYHeavy);
		if (!singletonIons3.isEmpty()) {
			addToMapByLabel(ret, singletonIons3);
		}
		final Map<QuantificationLabel, Set<Ion>> singletonIons4 = getIons(serieYLight);
		if (!singletonIons4.isEmpty()) {
			addToMapByLabel(ret, singletonIons4);
		}

		return ret;
	}

	private void addToMapByLabel(Map<QuantificationLabel, Set<Ion>> receiver,
			Map<QuantificationLabel, Set<Ion>> donor) {
		for (final QuantificationLabel label : donor.keySet()) {
			final Set<Ion> ions = donor.get(label);
			if (receiver.containsKey(label)) {
				receiver.get(label).addAll(ions);
			} else {
				final Set<Ion> set = new THashSet<Ion>();
				set.addAll(ions);
				receiver.put(label, set);
			}
		}

	}

	private void addToMapByCondition(Map<QuantCondition, Set<Ion>> receiver, Map<QuantificationLabel, Set<Ion>> donor) {
		for (final QuantificationLabel label : donor.keySet()) {
			final QuantCondition condition = conditionsByLabels.get(label);
			final Set<Ion> ions = donor.get(label);
			if (receiver.containsKey(condition)) {
				receiver.get(condition).addAll(ions);
			} else {
				final Set<Ion> set = new THashSet<Ion>();
				set.addAll(ions);
				receiver.put(condition, set);
			}
		}

	}

	/**
	 * @return the singletonIons
	 */
	@Override
	public Map<QuantificationLabel, Set<Ion>> getSingletonIons(IonSerie serie) {
		final Map<QuantificationLabel, Set<Ion>> ret = new THashMap<QuantificationLabel, Set<Ion>>();
		final Set<Ion> singletonIons = serie.getSingletonIons();
		// if (!singletonIons.isEmpty()) {
		ret.put(serie.getNonNullLabel(), singletonIons);
		// }
		return ret;
	}

	@Override
	public Map<QuantificationLabel, Set<Ion>> getIons(IonSerie serie) {
		final Map<QuantificationLabel, Set<Ion>> ret = new THashMap<QuantificationLabel, Set<Ion>>();
		if (serie != null) {
			final Set<Ion> ions = serie.getNonNullIons();
			if (!ions.isEmpty()) {
				ret.put(serie.getNonNullLabel(), ions);
			}
		}
		return ret;
	}

	@Override
	public Set<Ion> getIonsByLabel(QuantificationLabel label) {
		final Set<Ion> ret = new THashSet<Ion>();
		final Set<Ion> ions = getIonsByLabel().get(label);
		if (ions != null && !ions.isEmpty()) {
			ret.addAll(ions);
		}
		return ret;
	}

	@Override
	public Set<Ion> getSingletonIonsByLabel(QuantificationLabel label) {
		final Set<Ion> list = new THashSet<Ion>();
		final Set<Ion> singletonIons = getSingletonIonsByLabel().get(label);
		if (singletonIons != null && !singletonIons.isEmpty()) {
			list.addAll(singletonIons);
		}
		return list;
	}

	/**
	 * Gets the labels that this {@link IsobaricQuantifiedPSM} has been labeled ONLY
	 * with some label.<br>
	 * So, may happen that contains any ratio and it is not labeled
	 *
	 * @return the labels
	 */
	public Set<QuantificationLabel> getLabels() {
		return labels;
	}

	@Override
	public String getSequence() {
		if (super.getSequence() == null) {
			setFullSequence(censusChroPeptide.getSeq());
		}
		return super.getSequence();
	}

	@Override
	public String getIdentifier() {
		if (censusChroPeptide != null) {
			return QuantKeyUtils.getInstance().getSpectrumKey(censusChroPeptide, true);
		} else {
			return getRawFileNames().iterator().next() + "-" + getScanNumber();
		}
	}

	@Override
	public List<GroupableProtein> getGroupableProteins() {
		final List<GroupableProtein> ret = new ArrayList<GroupableProtein>();
		ret.addAll(getQuantifiedProteins());
		return ret;
	}

	@Override
	public Double getMaxPeak() {
		double max = -Double.MAX_VALUE;
		final Map<QuantificationLabel, Set<Ion>> ionsBylabel = getIonsByLabel();
		for (final Set<Ion> ions : ionsBylabel.values()) {
			for (final Ion ion : ions) {
				if (ion.getIntensity() > max)
					max = ion.getIntensity();
			}
		}
		return max;
	}

	public double getMaxPeakInRatio() {
		double max = -Long.MAX_VALUE;
		final Set<IsoRatio> ratios = getIsoRatios();
		for (final IsoRatio ratio : ratios) {
			if (ratio.getMaxIntensity() > max) {
				max = ratio.getMaxIntensity();
			}
		}
		return max;
	}

	@Override
	public double getMeanRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {

		final Set<IsoRatio> isobaricRatios = getNonInfinityIsoRatios();
		final TDoubleArrayList ratios = new TDoubleArrayList();

		for (final IsoRatio isoRatio : isobaricRatios) {
			ratios.add(isoRatio.getLog2Ratio(quantConditionNumerator, quantConditionDenominator));
		}
		return Maths.mean(ratios);
	}

	@Override
	public double getSTDRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {

		final Set<IsoRatio> isobaricRatios = getNonInfinityIsoRatios();
		final TDoubleArrayList ratios = new TDoubleArrayList();
		for (final IsoRatio isoRatio : isobaricRatios) {
			ratios.add(isoRatio.getLog2Ratio(quantConditionNumerator, quantConditionDenominator));
		}
		return Maths.stddev(ratios);
	}

	@Override
	public Map<QuantCondition, Set<Ion>> getSingletonIonsByCondition() {
		final Map<QuantCondition, Set<Ion>> ret = new THashMap<QuantCondition, Set<Ion>>();
		final Map<QuantificationLabel, Set<Ion>> singletonIons = getSingletonIons(serieBHeavy);
		// if (!singletonIons.isEmpty())
		addToMapByCondition(ret, singletonIons);

		final Map<QuantificationLabel, Set<Ion>> singletonIons2 = getSingletonIons(serieBLight);
		// if (!singletonIons2.isEmpty()) {
		addToMapByCondition(ret, singletonIons2);
		// }
		final Map<QuantificationLabel, Set<Ion>> singletonIons3 = getSingletonIons(serieYHeavy);
		// if (!singletonIons3.isEmpty()) {
		addToMapByCondition(ret, singletonIons3);

		// }
		final Map<QuantificationLabel, Set<Ion>> singletonIons4 = getSingletonIons(serieYLight);
		// if (!singletonIons4.isEmpty()) {
		addToMapByCondition(ret, singletonIons4);

		// }
		return ret;
	}

	@Override
	public Map<QuantCondition, Set<Ion>> getIonsByCondition() {
		final Map<QuantCondition, Set<Ion>> ret = new THashMap<QuantCondition, Set<Ion>>();
		final Map<QuantificationLabel, Set<Ion>> singletonIons = getIons(serieBHeavy);
		if (!singletonIons.isEmpty()) {
			addToMapByCondition(ret, singletonIons);
		}
		final Map<QuantificationLabel, Set<Ion>> singletonIons2 = getIons(serieBLight);
		if (!singletonIons2.isEmpty()) {
			addToMapByCondition(ret, singletonIons2);
		}
		final Map<QuantificationLabel, Set<Ion>> singletonIons3 = getIons(serieYHeavy);
		if (!singletonIons3.isEmpty()) {
			addToMapByCondition(ret, singletonIons3);
		}
		final Map<QuantificationLabel, Set<Ion>> singletonIons4 = getIons(serieYLight);
		if (!singletonIons4.isEmpty()) {
			addToMapByCondition(ret, singletonIons4);
		}
		return ret;
	}

	@Override
	public Map<QuantCondition, Set<Ion>> getIonsByCondition(String replicateName) {
		if (replicateName == null || getFileNames().contains(replicateName)) {
			return getIonsByCondition();
		}
		final Map<QuantCondition, Set<Ion>> ret = new THashMap<QuantCondition, Set<Ion>>();
		return ret;
	}

	@Override
	public boolean setQuantifiedPeptide(QuantifiedPeptideInterface quantifiedPeptide, boolean recursive) {
		return setPeptide(quantifiedPeptide, recursive);
	}

	/**
	 * @return the quantifiedPeptide
	 */
	@Override
	public QuantifiedPeptideInterface getQuantifiedPeptide() {
		final edu.scripps.yates.utilities.proteomicsmodel.Peptide peptide = getPeptide();
		if (peptide instanceof QuantifiedPeptideInterface) {
			return (QuantifiedPeptideInterface) peptide;
		}
		return null;
	}

	@Override
	public String getFullSequence() {
		if (super.getFullSequence() == null) {
			setFullSequence(censusChroPeptide.getSeq());
		}
		return super.getFullSequence();
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
			final IonCountRatio ratio = new IonCountRatio(AggregationLevel.PSM);
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

	public Peptide getCensusChroPeptide() {
		return censusChroPeptide;
	}

	@Override
	public Set<String> getFileNames() {
		return fileNames;
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator) {
		return QuantUtils.getAverageRatio(QuantUtils.getNonInfinityRatios(getQuantRatios()), AggregationLevel.PSM);
	}

	@Override
	public boolean isDiscarded() {
		return discarded;
	}

	@Override
	public void setDiscarded(boolean discarded) {
		this.discarded = discarded;

	}

	@Override
	public boolean isSingleton() {
		return getNonInfinityIsoRatios().isEmpty();
	}

	@Override
	public String toString() {
		return getKey();
	}

	@Override
	public boolean isQuantified() {
		return true;
	}

	public TIntObjectHashMap<Set<IsoRatio>> getSiteSpecificRatiosSeriesY(Set<Character> aas) {
		final String sequence = getSequence();
		final TIntObjectHashMap<Set<IsoRatio>> ret = new TIntObjectHashMap<Set<IsoRatio>>();
		// series Y, the positions are the opposite as the ion numbers
		final List<IsoRatio> ratiosFromSeries = getRatiosFromSeries(serieYLight, serieYHeavy);
		final TIntObjectHashMap<IsoRatio> isoRatiosByPosition = new TIntObjectHashMap<IsoRatio>();

		for (final IsoRatio isoRatio : ratiosFromSeries) {
			isoRatiosByPosition.put(isoRatio.getNumIon(), isoRatio);
		}
		for (int position = 1; position <= sequence.length(); position++) {
			if (isoRatiosByPosition.containsKey(position)) {
				final IsoRatio isoRatio = isoRatiosByPosition.get(position);
				final int sequencePosition = sequence.length() - position + 1;
				// go from there to the start of the sequence seeing how many aa
				// are from aas Set
				final List<Integer> tmpPositions = new ArrayList<Integer>();
				for (int j = sequencePosition; j <= sequence.length(); j++) {
					if (aas.contains(sequence.charAt(j - 1))) {
						tmpPositions.add(j);
					}
				}
				// only if num==1 we have a valid ratio
				if (tmpPositions.size() == 1) {
					final Integer quantifiedPosition = tmpPositions.get(0);
					final Character quantifiedAA = sequence.charAt(quantifiedPosition - 1);
					isoRatio.setQuantifiedAA(quantifiedAA);
					isoRatio.addQuantifiedSitePositionInPeptide(
							new PositionInPeptide(quantifiedPosition, quantifiedAA, sequence));
					if (ret.containsKey(quantifiedPosition)) {
						ret.get(quantifiedPosition).add(isoRatio);
					} else {
						final Set<IsoRatio> set = new THashSet<IsoRatio>();
						set.add(isoRatio);
						ret.put(quantifiedPosition, set);
					}
				} else if (!Double.isNaN(isoRatio.getValue())
						&& Double.isFinite(isoRatio.getLog2Ratio(serieBLight.getLabel(), serieBHeavy.getLabel()))) {
					// also, if the ratio is infinite, it doesn't matters
					// I will assign the ratio to all positions
					for (final Integer quantifiedPosition : tmpPositions) {
						final Character quantifiedAA = sequence.charAt(quantifiedPosition - 1);
						final IsoRatio clonedIsoRatio = cloneIsoRatio(isoRatio);
						clonedIsoRatio.setQuantifiedAA(quantifiedAA);
						clonedIsoRatio.addQuantifiedSitePositionInPeptide(
								new PositionInPeptide(quantifiedPosition, quantifiedAA, sequence));
						if (ret.containsKey(quantifiedPosition)) {
							ret.get(quantifiedPosition).add(clonedIsoRatio);
						} else {
							final Set<IsoRatio> set = new THashSet<IsoRatio>();
							set.add(clonedIsoRatio);
							ret.put(quantifiedPosition, set);
						}
					}
				}
			}
		}
		return ret;
	}

	public TIntObjectHashMap<Set<IsoRatio>> getSiteSpecificRatiosSeriesB(Set<Character> aas) {
		final String sequence = getSequence();
		final TIntObjectHashMap<Set<IsoRatio>> ret = new TIntObjectHashMap<Set<IsoRatio>>();
		// series B, the positions are the same as the ion numbers
		final List<IsoRatio> ratiosFromSeries = getRatiosFromSeries(serieBLight, serieBHeavy);
		final TIntObjectHashMap<IsoRatio> isoRatiosByPosition = new TIntObjectHashMap<IsoRatio>();

		for (final IsoRatio isoRatio : ratiosFromSeries) {
			isoRatiosByPosition.put(isoRatio.getNumIon(), isoRatio);
		}
		for (int position = 1; position <= sequence.length(); position++) {

			if (isoRatiosByPosition.containsKey(position)) {
				final IsoRatio isoRatio = isoRatiosByPosition.get(position);
				// go from there to the start of the sequence seeing how many aa
				// are from aas Set
				final List<Integer> tmpPositions = new ArrayList<Integer>();
				for (int j = position; j >= 1; j--) {
					if (aas.contains(sequence.charAt(j - 1))) {
						tmpPositions.add(j);
					}
				}
				// only if num==1 we have a valid ratio
				if (tmpPositions.size() == 1) {
					final Integer quantifiedPosition = tmpPositions.get(0);
					final Character quantifiedAA = sequence.charAt(quantifiedPosition - 1);
					isoRatio.setQuantifiedAA(quantifiedAA);
					isoRatio.addQuantifiedSitePositionInPeptide(
							new PositionInPeptide(quantifiedPosition, quantifiedAA, sequence));
					if (ret.containsKey(quantifiedPosition)) {
						ret.get(quantifiedPosition).add(isoRatio);
					} else {
						final Set<IsoRatio> set = new THashSet<IsoRatio>();
						set.add(isoRatio);
						ret.put(quantifiedPosition, set);
					}
				} else if (!Double.isNaN(isoRatio.getValue())
						&& Double.isInfinite(isoRatio.getLog2Ratio(serieBLight.getLabel(), serieBHeavy.getLabel()))) {
					// also, if the ratio is infinite, it doesn't matters
					// I will assign the ratio to all positions
					for (final Integer quantifiedPosition : tmpPositions) {
						final Character quantifiedAA = sequence.charAt(quantifiedPosition - 1);
						final IsoRatio clonedIsoRatio = cloneIsoRatio(isoRatio);
						clonedIsoRatio.setQuantifiedAA(quantifiedAA);
						clonedIsoRatio.addQuantifiedSitePositionInPeptide(
								new PositionInPeptide(quantifiedPosition, quantifiedAA, sequence));
						if (ret.containsKey(quantifiedPosition)) {
							ret.get(quantifiedPosition).add(clonedIsoRatio);
						} else {
							final Set<IsoRatio> set = new THashSet<IsoRatio>();
							set.add(clonedIsoRatio);
							ret.put(quantifiedPosition, set);
						}
					}
				}
			}
		}
		return ret;
	}

	private IsoRatio cloneIsoRatio(IsoRatio isoRatio) {
		final IsoRatio clonedRatio = new IsoRatio(isoRatio.getIon(isoRatio.getLabel1()), isoRatio.getLabel1(),
				isoRatio.getQuantCondition1(), isoRatio.getIon(isoRatio.getLabel2()), isoRatio.getLabel2(),
				isoRatio.getQuantCondition2(), isoRatio.getNumIon(), isoRatio.getIonSerieType(),
				isoRatio.getAggregationLevel());
		return clonedRatio;
	}

	public TIntObjectHashMap<Set<IsoRatio>> getSiteSpecificRatios(Set<Character> aas) {
		final TIntObjectHashMap<Set<IsoRatio>> ret = new TIntObjectHashMap<Set<IsoRatio>>();
		final TIntObjectHashMap<Set<IsoRatio>> siteSpecificRatiosSeriesB = getSiteSpecificRatiosSeriesB(aas);
		for (final int position : siteSpecificRatiosSeriesB.keys()) {
			final Set<IsoRatio> isoRatios = siteSpecificRatiosSeriesB.get(position);
			if (ret.containsKey(position)) {
				ret.get(position).addAll(isoRatios);
			} else {
				final Set<IsoRatio> set = new THashSet<IsoRatio>();
				set.addAll(isoRatios);
				ret.put(position, set);
			}

		}
		final TIntObjectHashMap<Set<IsoRatio>> siteSpecificRatiosSeriesY = getSiteSpecificRatiosSeriesY(aas);
		for (final int position : siteSpecificRatiosSeriesY.keys()) {
			final Set<IsoRatio> isoRatios = siteSpecificRatiosSeriesY.get(position);
			if (ret.containsKey(position)) {
				ret.get(position).addAll(isoRatios);
			} else {
				final Set<IsoRatio> set = new THashSet<IsoRatio>();
				set.addAll(isoRatios);
				ret.put(position, set);
			}

		}

		return ret;
	}

	private TIntObjectHashMap<IsoRatio> getIsoRatiosSortedByIonNumber(IonSerieType serie) {
		final TIntObjectHashMap<IsoRatio> ret = new TIntObjectHashMap<IsoRatio>();

		for (final QuantRatio ratio : getQuantRatios()) {
			if (ratio instanceof IsoRatio) {
				final IsoRatio isoRatio = (IsoRatio) ratio;
				if (isoRatio.getIonSerieType() == serie) {
					ret.put(isoRatio.getNumIon(), isoRatio);
				}
			}
		}

		return ret;

	}

	@Override
	public Map<QuantCondition, Set<Ion>> getIonsByConditionForSites(String replicateName, char[] quantifiedAAs) {
		final Map<QuantCondition, Set<Ion>> ret = new THashMap<QuantCondition, Set<Ion>>();
		final Map<QuantCondition, Set<Ion>> ionsByCondition2 = getIonsByCondition(replicateName);
		final TIntArrayList quantifiedPositions = StringUtils.getPositions(getSequence(), quantifiedAAs);
		if (ionsByCondition2 != null) {
			for (final QuantCondition condition : ionsByCondition2.keySet()) {
				final Set<Ion> ions = ionsByCondition2.get(condition);
				for (final Ion ion : ions) {
					final int ionNumber = ion.getIonNumber();
					TIntArrayList quantPositions = null;
					switch (ion.getIonSerieType()) {
					case B:
						quantPositions = quantifiedPositions;
						break;
					case Y:
						// we reverse positions. So for a peptide with lenth 9,
						// having a quantified AA at 3 means that now we are
						// going to have a quantified AA of
						quantPositions = reversePositions(quantifiedPositions, getSequence().length());
						break;
					default:
						break;
					}
					// we need that the ion number is >= to the
					// position of the aas is quantified and that there is
					// only one quantified aa with ion number >= to this ion
					if (getNumberOfNumbersEqualOrLessThanValue(quantPositions, ionNumber) == 1) {
						if (ret.containsKey(condition)) {
							ret.get(condition).add(ion);
						} else {
							final Set<Ion> set = new THashSet<Ion>();
							set.add(ion);
							ret.put(condition, set);
						}
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
		final Map<QuantCondition, Set<Ion>> ionsByCondition2 = getIonsByCondition(replicateName);
		final TIntArrayList quantifiedPositions = StringUtils.getPositions(getSequence(), quantifiedAAs);
		if (ionsByCondition2 != null) {
			for (final QuantCondition condition : ionsByCondition2.keySet()) {
				final Set<Ion> ions = ionsByCondition2.get(condition);
				for (final Ion ion : ions) {
					final int ionNumber = ion.getIonNumber();
					TIntArrayList quantPositions = null;
					int quantifiedPosition = -1;
					switch (ion.getIonSerieType()) {
					case B:
						quantPositions = quantifiedPositions;
						quantifiedPosition = positionInPeptide;
						break;
					case Y:
						// we reverse positions. So for a peptide with lenth 9,
						// having a quantified AA at 3 means that now we are
						// going to have a quantified AA of
						quantPositions = reversePositions(quantifiedPositions, getSequence().length());
						quantifiedPosition = reversePosition(positionInPeptide, getSequence().length());
						break;
					default:
						break;
					}
					// we need that the ion number is >= to the
					// position of the aas is quantified and that there is
					// only one quantified aa with ion number >= to this ion
					if (ionNumber >= quantifiedPosition) {
						// if we are here means that the ion has information
						// about the position of interest
						// let's see if if it only has information about one
						// site
						if (getNumberOfNumbersEqualOrLessThanValue(quantPositions, ionNumber) == 1) {
							if (ret.containsKey(condition)) {
								ret.get(condition).add(ion);
							} else {
								final Set<Ion> set = new THashSet<Ion>();
								set.add(ion);
								ret.put(condition, set);
							}
						}
					}

				}
			}
		}

		return ret;
	}

	private TIntArrayList reversePositions(TIntArrayList quantifiedPositions, int length) {
		final TIntArrayList ret = new TIntArrayList();
		for (final int quantifiedPosition : quantifiedPositions.toArray()) {
			ret.add(length - quantifiedPosition + 1);
		}
		return ret;
	}

	private int reversePosition(int quantifiedPosition, int length) {
		return length - quantifiedPosition + 1;
	}

	private int getNumberOfNumbersEqualOrLessThanValue(TIntArrayList numberlist, int number) {
		int ret = 0;
		for (final int integer : numberlist.toArray()) {
			if (integer <= number) {
				ret++;
			}
		}
		return ret;
	}

	@Override
	public void setSingleton(boolean singleton) {
		// TODO Auto-generated method stub

	}

}
