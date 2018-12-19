package edu.scripps.yates.census.read.model;

import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.HasQuantRatios;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.proteomicsmodel.Ratio;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.set.hash.THashSet;

public abstract class AbstractContainsQuantifiedPSMs implements HasQuantRatios {

	protected final Set<QuantRatio> quantRatios = new THashSet<QuantRatio>();
	private Set<Ratio> ratios;

	public abstract Set<QuantifiedPSMInterface> getQuantifiedPSMs();

	@Override
	public Set<QuantRatio> getQuantRatios() {
		if (quantRatios == null || quantRatios.isEmpty()) {
			for (final QuantifiedPSMInterface psm : getQuantifiedPSMs()) {
				quantRatios.addAll(psm.getQuantRatios());
			}
		}
		return quantRatios;
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

	@Override
	public Set<Ratio> getRatios() {
		if (ratios == null || ratios.isEmpty()) {
			ratios = new THashSet<Ratio>();
			for (final QuantifiedPSMInterface psm : getQuantifiedPSMs()) {
				ratios.addAll(psm.getRatios());
			}
		}
		return ratios;
	}

	@Override
	public boolean addRatio(Ratio ratio) {
		return getRatios().add(ratio);
	}

	@Override
	public Set<QuantRatio> getNonInfinityRatios() {
		return QuantUtils.getNonInfinityRatios(getQuantRatios());
	}

	@Override
	public boolean addQuantRatio(QuantRatio ratio) {
		return addRatio(ratio);
	}
}
