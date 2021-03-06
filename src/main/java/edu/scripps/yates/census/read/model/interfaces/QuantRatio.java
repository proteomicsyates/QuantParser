package edu.scripps.yates.census.read.model.interfaces;

import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;

public interface QuantRatio extends edu.scripps.yates.utilities.proteomicsmodel.Ratio {

	public double getLog2Ratio(QuantificationLabel labelNumerator, QuantificationLabel labelDenominator);

	public double getNonLogRatio(QuantificationLabel labelNumerator, QuantificationLabel labelDenominator);

	public QuantificationLabel getLabel1();

	public QuantificationLabel getLabel2();

	public double getLog2Ratio(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator);

	public double getNonLogRatio(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator);

	public double getLog2Ratio(String condition1Name, String condition2Name);

	public double getNonLogRatio(String condition1Name, String condition2Name);

	public QuantCondition getQuantCondition1();

	public QuantCondition getQuantCondition2();

	public Character getQuantifiedAA();

	public void setQuantifiedAA(Character c);

	public void addQuantifiedSitePositionInPeptide(PositionInPeptide quantifiedSitePositionInPeptide);

	public Set<PositionInPeptide> getQuantifiedSitePositionInPeptide();

	/**
	 * The number of measurements that have been used to calculate this ratio.
	 * In case of being an average ratio over 3 ratios, this function will
	 * return a 3.
	 * 
	 * @return
	 */
	public int getNumMeasurements();

	public void setNumMeasurements(int numMeasurements);
}
