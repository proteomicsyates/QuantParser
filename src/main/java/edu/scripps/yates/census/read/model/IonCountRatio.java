package edu.scripps.yates.census.read.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.model.enums.CombinationType;
import edu.scripps.yates.utilities.proteomicsmodel.Condition;
import edu.scripps.yates.utilities.proteomicsmodel.Score;
import gnu.trove.map.hash.THashMap;

public class IonCountRatio implements QuantRatio {
	public static final IonCountRatio NAN_RATIO = new IonCountRatio(null);
	private final Map<QuantCondition, Double> ionCountMap = new THashMap<QuantCondition, Double>();
	private Score score;
	private CombinationType combinationType;
	private final AggregationLevel aggregationLevel;
	private String description = "Rc";
	private Condition condition2;
	private Integer quantifiedSitePositionInPeptide;
	private Character quantifiedAA;

	public IonCountRatio(AggregationLevel aggregationLevel) {
		this.aggregationLevel = aggregationLevel;
	}

	public void addIonCount(QuantCondition cond, double peakCount) {
		if (ionCountMap.containsKey(cond)) {
			ionCountMap.put(cond, ionCountMap.get(cond) + peakCount);
		} else {
			ionCountMap.put(cond, peakCount);
		}
	}

	public void setAsNormalizedIonCountRatio() {
		setDescription("NormRc");
	}

	private void setDescription(String description) {
		this.description = description;
	}

	@Override
	public double getNonLogRatio(QuantCondition cond1, QuantCondition cond2) {
		if (ionCountMap.containsKey(cond1) && ionCountMap.containsKey(cond2)) {
			return ionCountMap.get(cond1) / ionCountMap.get(cond2);
		}
		return Double.NaN;
	}

	@Override
	public double getLog2Ratio(QuantCondition cond1, QuantCondition cond2) {
		final Double countRatio = getNonLogRatio(cond1, cond2);
		if (countRatio.isNaN()) {
			return countRatio;
		}
		return Math.log(countRatio) / Math.log(2.0);
	}

	public Double getIonCount(QuantCondition cond) {
		if (ionCountMap.containsKey(cond)) {
			return ionCountMap.get(cond);
		}
		return 0.0;
	}

	@Override
	public double getValue() {
		throw new IllegalArgumentException("Not supported in a consensus isobaric count ratio");
	}

	@Override
	public Condition getCondition1() {
		final List<QuantCondition> conditionList = getSortedConditionlist();
		if (conditionList.size() > 0)
			return conditionList.get(0);
		return null;
	}

	private List<QuantCondition> getSortedConditionlist() {
		final Set<QuantCondition> keySet = ionCountMap.keySet();
		final List<QuantCondition> conditionList = new ArrayList<QuantCondition>();
		conditionList.addAll(keySet);
		Collections.sort(conditionList, new Comparator<QuantCondition>() {

			@Override
			public int compare(QuantCondition o1, QuantCondition o2) {
				return o1.getName().compareTo(o2.getName());
			}
		});
		return conditionList;
	}

	@Override
	public Condition getCondition2() {
		final List<QuantCondition> conditionList = getSortedConditionlist();
		if (conditionList.size() > 1)
			return conditionList.get(1);
		return null;
	}

	@Override
	public String getDescription() {

		return description;
	}

	@Override
	public Score getAssociatedConfidenceScore() {
		return score;
	}

	public void setAssociatedConfidenceScore(Score score) {
		this.score = score;
	}

	@Override
	public CombinationType getCombinationType() {
		return combinationType;
	}

	public void setCombinationType(CombinationType combinationType) {
		this.combinationType = combinationType;
	}

	@Override
	public AggregationLevel getAggregationLevel() {
		return aggregationLevel;
	}

	@Override
	public double getLog2Ratio(QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) {
		throw new IllegalArgumentException(
				"This is not supported. This is a consensus ratio that may come from ions counts with different labels");
	}

	@Override
	public double getNonLogRatio(QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) {
		throw new IllegalArgumentException(
				"This is not supported. This is a consensus ratio that may come from ions counts with different labels");

	}

	@Override
	public QuantificationLabel getLabel1() {
		throw new IllegalArgumentException(
				"This is not supported. This is a consensus ratio that may come from ions counts with different labels");

	}

	@Override
	public QuantificationLabel getLabel2() {
		throw new IllegalArgumentException(
				"This is not supported. This is a consensus ratio that may come from ions counts with different labels");

	}

	@Override
	public double getLog2Ratio(String condition1Name, String condition2Name) {
		final QuantCondition condition1 = getConditionByName(condition1Name);
		final QuantCondition condition2 = getConditionByName(condition1Name);
		if (condition1 != null && condition2 != null) {
			return getLog2Ratio(condition1, condition2);
		}
		return Double.NaN;
	}

	private QuantCondition getConditionByName(String conditionName) {
		final Set<QuantCondition> conditionSet = ionCountMap.keySet();
		for (final QuantCondition quantCondition : conditionSet) {
			if (quantCondition.getName().equalsIgnoreCase(conditionName)) {
				return quantCondition;
			}
		}
		return null;
	}

	@Override
	public double getNonLogRatio(String condition1Name, String condition2Name) {
		final QuantCondition condition1 = getConditionByName(condition1Name);
		final QuantCondition condition2 = getConditionByName(condition1Name);
		if (condition1 != null && condition2 != null) {
			return getNonLogRatio(condition1, condition2);
		}
		return Double.NaN;
	}

	@Override
	public QuantCondition getQuantCondition1() {
		return (QuantCondition) getCondition1();
	}

	@Override
	public QuantCondition getQuantCondition2() {
		return (QuantCondition) getCondition2();
	}

	@Override
	public Integer getQuantifiedSitePositionInPeptide() {
		return quantifiedSitePositionInPeptide;
	}

	@Override
	public void setQuantifiedSitePositionInPeptide(int quantifiedSitePositionInPeptide) {
		this.quantifiedSitePositionInPeptide = quantifiedSitePositionInPeptide;
	}

	@Override
	public Character getQuantifiedAA() {
		return quantifiedAA;
	}

	@Override
	public void setQuantifiedAA(Character quantifiedAA) {
		this.quantifiedAA = quantifiedAA;
	}
}
