package edu.scripps.yates.census.read.model;

import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.HasRatios;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.utilities.maths.Maths;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.set.hash.THashSet;

public abstract class AbstractContainsQuantifiedPSMs implements HasRatios {

	protected final Set<QuantRatio> ratios = new THashSet<QuantRatio>();

	public abstract Set<QuantifiedPSMInterface> getQuantifiedPSMs();

	@Override
	public Set<QuantRatio> getRatios() {
		if (ratios == null || ratios.isEmpty()) {
			for (final QuantifiedPSMInterface psm : getQuantifiedPSMs()) {
				ratios.addAll(psm.getRatios());
			}
		}
		return ratios;
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

}
