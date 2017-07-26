package edu.scripps.yates.census.read.model;

import java.util.HashMap;
import java.util.Map;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.util.QuantificationLabel;

public class RatioDescriptor {

	private final QuantificationLabel label1;
	private final QuantificationLabel label2;
	private final QuantCondition condition1;
	private final QuantCondition condition2;

	public RatioDescriptor(QuantificationLabel label1, QuantificationLabel label2, QuantCondition condition1,
			QuantCondition condition2) {
		super();
		this.label1 = label1;
		this.label2 = label2;
		this.condition1 = condition1;
		this.condition2 = condition2;
	}

	public QuantificationLabel getLabel1() {
		return label1;
	}

	public QuantificationLabel getLabel2() {
		return label2;
	}

	public QuantCondition getQuantCondition1() {
		return condition1;
	}

	public QuantCondition getQuantCondition2() {
		return condition2;
	}

	/**
	 * This returns _L_H for a light over heavy ratio, _L_M for a medium over
	 * heavy ratio, and so on
	 *
	 * @return
	 */
	public String getRatioSuffix() {
		String letter1 = getLetterForLabel(label1);
		String letter2 = getLetterForLabel(label2);

		return "_" + letter1 + "_" + letter2;
	}

	/**
	 * This returns a map of a suffix for each quantcondition, so [_L,
	 * lightcondition] for a light condition
	 * 
	 * @return
	 */
	public Map<String, QuantCondition> getConditionsByIndividualRatioSuffixes() {
		Map<String, QuantCondition> ret = new HashMap<String, QuantCondition>();
		ret.put("_" + getLetterForLabel(label1), condition1);
		ret.put("_" + getLetterForLabel(label2), condition2);
		return ret;
	}

	private String getLetterForLabel(QuantificationLabel label) {
		if (label == QuantificationLabel.LIGHT) {
			return "L";
		} else if (label == QuantificationLabel.MEDIUM) {
			return "M";
		} else if (label == QuantificationLabel.HEAVY) {
			return "H";
		}
		return label.name().substring(0, 1);
	}

}
