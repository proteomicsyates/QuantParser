package edu.scripps.yates.census.read.model;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;

public class IonCountRatio extends CensusRatio {

	public IonCountRatio(int numIons1, int numIons2, QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator, AggregationLevel aggregationLevel) {
		this(numIons1, numIons2, quantConditionNumerator, quantConditionDenominator, null, null, aggregationLevel);
	}

	public IonCountRatio(int numIons1, int numIons2, QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator, AggregationLevel aggregationLevel) {
		super(getNonLogValue(numIons1, numIons2), null, false, quantConditionNumerator, quantConditionDenominator,
				labelNumerator, labelDenominator, aggregationLevel, "Ion Count Ratio");
	}

	public static double getNonLogValue(int numIons1, int numIons2) {
		if (numIons1 == 0 && numIons2 == 0) {
			return Double.NaN;
		} else {
			if (numIons2 == 0) {
				return Double.POSITIVE_INFINITY;
			} else {
				return Double.valueOf(numIons1) / Double.valueOf(numIons2);
			}
		}

	}
}
