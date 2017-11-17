package edu.scripps.yates.census.read.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.HasRatios;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.utilities.maths.Maths;
import gnu.trove.set.hash.THashSet;

public abstract class AbstractContainsQuantifiedPSMs implements HasRatios {

	protected final Set<QuantRatio> ratios = new THashSet<QuantRatio>();

	public abstract Set<QuantifiedPSMInterface> getQuantifiedPSMs();

	@Override
	public Set<QuantRatio> getRatios() {
		if (ratios == null || ratios.isEmpty()) {
			for (QuantifiedPSMInterface psm : getQuantifiedPSMs()) {
				ratios.addAll(psm.getRatios());
			}
		}
		return ratios;
	}

	@Override
	public double getMeanRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {
		List<Double> ratioValues = new ArrayList<Double>();
		for (QuantifiedPSMInterface psm : getQuantifiedPSMs()) {
			QuantRatio ratio = QuantUtils.getRepresentativeRatio(psm);
			ratioValues.add(ratio.getLog2Ratio(quantConditionNumerator, quantConditionDenominator));
		}
		return Maths.mean(ratioValues.toArray(new Double[0]));
	}

	@Override
	public double getSTDRatios(QuantCondition quantConditionNumerator, QuantCondition quantConditionDenominator) {
		List<Double> ratioValues = new ArrayList<Double>();

		for (QuantifiedPSMInterface psm : getQuantifiedPSMs()) {
			QuantRatio ratio = QuantUtils.getRepresentativeRatio(psm);
			ratioValues.add(ratio.getLog2Ratio(quantConditionNumerator, quantConditionDenominator));
		}

		return Maths.stddev(ratioValues.toArray(new Double[0]));
	}

}
