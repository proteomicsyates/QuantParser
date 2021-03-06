package edu.scripps.yates.census.read.model;

import java.util.Map;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.CensusChroParser;
import edu.scripps.yates.census.read.model.IonSerie.IonSerieType;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AggregationLevel;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.THashMap;

public class IsoRatio extends CensusRatio {

	private final int numIon;
	private final IonSerieType ionSerieType;
	private final Map<QuantificationLabel, Ion> ionsByLabel = new THashMap<QuantificationLabel, Ion>();
	private final QuantificationLabel quantificationLabel1;
	private final QuantificationLabel quantificationLabel2;
	private final Map<QuantificationLabel, Double> massesByLabel = new THashMap<QuantificationLabel, Double>();
	private final Map<QuantCondition, QuantificationLabel> labelsByConditions;
	private final Map<QuantificationLabel, QuantCondition> conditionsByLabels;
	private RatioScore ratioScore;

	/**
	 *
	 * @param ion1
	 *            can be null
	 * @param quantificationLabel1
	 * @param ion2
	 *            can be null
	 * @param quantificationLabel2
	 * @param numIon
	 * @param ionSerieType
	 */
	public IsoRatio(Ion ion1, QuantificationLabel quantificationLabel1, QuantCondition condition1, Ion ion2,
			QuantificationLabel quantificationLabel2, QuantCondition condition2, int numIon, IonSerieType ionSerieType,
			AggregationLevel aggregationLevel) {
		super(null, null, null, condition1, condition2, quantificationLabel1, quantificationLabel2, aggregationLevel,
				CensusChroParser.ISOBARIC_INTENSITY_RATIO);

		if (ion1 == null && ion2 == null)
			throw new IllegalArgumentException("Ions and Ionr cannot be null at the same time");

		if (ion1 != null) {
			ionsByLabel.put(quantificationLabel1, ion1);
			ion1.setRatio(this);
			massesByLabel.put(quantificationLabel1, ion1.getMass());
			if (ion2 == null) {
				ion1.setSingleton(true);
			}
		}
		this.quantificationLabel1 = quantificationLabel1;
		if (ion2 != null) {
			ionsByLabel.put(quantificationLabel2, ion2);
			ion2.setRatio(this);
			massesByLabel.put(quantificationLabel2, ion2.getMass());
			if (ion1 == null) {
				ion2.setSingleton(true);
			}
		}
		this.quantificationLabel2 = quantificationLabel2;
		this.numIon = numIon;
		this.ionSerieType = ionSerieType;

		labelsByConditions = new THashMap<QuantCondition, QuantificationLabel>();
		labelsByConditions.put(condition1, quantificationLabel1);
		labelsByConditions.put(condition2, quantificationLabel2);
		conditionsByLabels = new THashMap<QuantificationLabel, QuantCondition>();
		conditionsByLabels.put(quantificationLabel1, condition1);
		conditionsByLabels.put(quantificationLabel2, condition2);
	}

	/**
	 * Gets the log2 of the ratio
	 *
	 * @return the ratio
	 */
	@Override
	public double getLog2Ratio(QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) {
		if (ionsByLabel.containsKey(labelNumerator) && ionsByLabel.containsKey(labelDenominator)) {
			final Ion ion1 = ionsByLabel.get(labelNumerator);
			final Ion ion2 = ionsByLabel.get(labelDenominator);
			final double value = Double.valueOf(ion1.getIntensity()) / Double.valueOf(ion2.getIntensity());
			return Math.log(value) / Math.log(2);
		} else if (ionsByLabel.containsKey(labelNumerator) && !ionsByLabel.containsKey(labelDenominator)) {
			return Double.POSITIVE_INFINITY;
		} else if (!ionsByLabel.containsKey(labelNumerator) && ionsByLabel.containsKey(labelDenominator)) {
			// ion1==null && ion2!=null
			return Double.NEGATIVE_INFINITY;
		} else {
			throw new IllegalArgumentException(
					"No ratio value found for labels: " + labelNumerator + "/" + labelDenominator);
		}
	}

	/**
	 * @return the maxPeak
	 */
	public double getMaxIntensity() {
		double max = Double.MIN_VALUE;
		if (ionsByLabel.containsKey(quantificationLabel1) && ionsByLabel.get(quantificationLabel1).getIntensity() > max)
			max = ionsByLabel.get(quantificationLabel1).getIntensity();
		if (ionsByLabel.containsKey(quantificationLabel2) && ionsByLabel.get(quantificationLabel2).getIntensity() > max)
			max = ionsByLabel.get(quantificationLabel2).getIntensity();
		return max;
	}

	/**
	 * @return the average of the intensities of the peaks of the ratio
	 */
	public double getAverageIntensityPeak() {
		final TDoubleArrayList values = new TDoubleArrayList();
		if (ionsByLabel.containsKey(quantificationLabel1)) {
			values.add(Double.valueOf(ionsByLabel.get(quantificationLabel1).getIntensity()));
		}
		if (ionsByLabel.containsKey(quantificationLabel2)) {
			values.add(Double.valueOf(ionsByLabel.get(quantificationLabel2).getIntensity()));
		}
		if (!values.isEmpty())
			return Maths.mean(values);
		return Double.NaN;
	}

	/**
	 * @return the ionSerieType
	 */
	public IonSerieType getIonSerieType() {
		return ionSerieType;
	}

	/**
	 * @return the numIon
	 */
	public int getNumIon() {
		return numIon;
	}

	/**
	 * @return the intensity1
	 */
	public double getIntensity(QuantificationLabel label) {
		if (ionsByLabel.containsKey(label))
			return ionsByLabel.get(label).getIntensity();
		return Double.NaN;
	}

	/**
	 * gets the ratio. No logarithmic transformation
	 *
	 * @return the norLogRatio
	 */
	@Override
	public double getNonLogRatio(QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) {
		if (ionsByLabel.containsKey(labelNumerator) && ionsByLabel.containsKey(labelDenominator)) {
			final Ion ion1 = ionsByLabel.get(labelNumerator);
			final Ion ion2 = ionsByLabel.get(labelDenominator);
			final double value = Double.valueOf(ion1.getIntensity()) / Double.valueOf(ion2.getIntensity());
			return value;

		} else if (ionsByLabel.containsKey(labelNumerator) && !ionsByLabel.containsKey(labelDenominator)) {
			return Double.POSITIVE_INFINITY; // n/0
		} else if (!ionsByLabel.containsKey(labelNumerator) && ionsByLabel.containsKey(labelDenominator)) {
			// ion1==null && ion2!=null
			return 0.0; // 0/n
		} else {
			throw new IllegalArgumentException(
					"No ratio value found for labels: " + labelNumerator + "/" + labelDenominator);
		}
	}

	@Override
	public QuantificationLabel getLabel1() {
		return quantificationLabel1;
	}

	@Override
	public QuantificationLabel getLabel2() {
		return quantificationLabel2;

	}

	public static void main(String[] args) {
		double num = Double.NEGATIVE_INFINITY;
		System.out.println(num);
		num = Double.POSITIVE_INFINITY;
		System.out.println(num);
		num = Double.NaN;
		System.out.println(num);

		System.out.println(1.0 / 0);

		final TDoubleArrayList ratioValues = new TDoubleArrayList();
		ratioValues.add(Double.POSITIVE_INFINITY);
		ratioValues.add(4.0);
		ratioValues.add(6.0);
		final double mean = Maths.mean(ratioValues);
		System.out.println("mean=" + mean);

		System.out.println("OK");
	}

	public Ion getIon(QuantificationLabel label) {
		return ionsByLabel.get(label);
	}

	/**
	 * Gets the mass of the peak labeled with a certain
	 * {@link QuantificationLabel}
	 *
	 * @param label
	 * @return
	 */
	public Double getMass(QuantificationLabel label) {
		if (massesByLabel.containsKey(label)) {
			return massesByLabel.get(label);
		}
		for (final QuantificationLabel label2 : massesByLabel.keySet()) {
			if (label2.isLight() == label.isLight()) {
				return massesByLabel.get(label);
			}
		}
		return null;
	}

	@Override
	public double getLog2Ratio(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {
		return getLog2Ratio(labelsByConditions.get(quantConditionNumerator),
				labelsByConditions.get(quantConditionDenominator));
	}

	public Double getMass(QuantCondition condition) {
		return massesByLabel.get(labelsByConditions.get(condition));
	}

	public double getIntensity(QuantCondition conditionDenominator) {
		if (ionsByLabel.containsKey(labelsByConditions.get(conditionDenominator)))
			return ionsByLabel.get(labelsByConditions.get(conditionDenominator)).getIntensity();
		return Double.NaN;
	}

}
