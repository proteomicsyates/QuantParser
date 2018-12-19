package edu.scripps.yates.census.read.model.interfaces;

import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;

public interface HasQuantRatios extends edu.scripps.yates.utilities.proteomicsmodel.HasRatios {

	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator);

	public Set<QuantRatio> getNonInfinityRatios();

	public Set<QuantRatio> getQuantRatios();

	public boolean addQuantRatio(QuantRatio quantRatio);

	public double getMeanRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator);

	public double getSTDRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator);

}
