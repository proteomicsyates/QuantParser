package edu.scripps.yates.census.read.model;

import org.apache.commons.lang3.builder.HashCodeBuilder;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.proteomicsmodel.Condition;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AmountType;
import edu.scripps.yates.utilities.proteomicsmodel.enums.CombinationType;

public class QuantAmount implements Amount {
	private double value;
	private AmountType amountType;

	private CombinationType combinationType;

	private Boolean singleton = false;

	private Boolean manualSpc = false;
	private QuantCondition quantCondition;
	private QuantificationLabel label;
	private Condition condition;
	private int hashCode = -1;

	public QuantAmount(double value, AmountType amountType, QuantCondition condition, QuantificationLabel label) {

		this.value = value;
		this.amountType = amountType;
		quantCondition = condition;
		this.setLabel(label);
	}

	/**
	 * @return the value
	 */

	@Override
	public double getValue() {
		return value;
	}

	/**
	 * @param value the value to set
	 */
	public void setValue(double value) {
		this.value = value;
	}

	/**
	 * @return the amountType
	 */

	@Override
	public AmountType getAmountType() {
		return amountType;
	}

	/**
	 * @param amountType the amountType to set
	 */
	public void setAmountType(AmountType amountType) {
		this.amountType = amountType;
	}

	/**
	 * @return the combinationType
	 */

	@Override
	public CombinationType getCombinationType() {
		return combinationType;
	}

	/**
	 * @param combinationType the combinationType to set
	 */
	public void setCombinationType(CombinationType combinationType) {
		this.combinationType = combinationType;
	}

	/**
	 * @return the singleton
	 */

	@Override
	public Boolean isSingleton() {
		return singleton;
	}

	/**
	 * @param singleton the singleton to set
	 */
	public void setSingleton(boolean singleton) {
		this.singleton = singleton;
	}

	/**
	 * @return the manualSpc
	 */

	@Override
	public Boolean isManualSpc() {
		return manualSpc;
	}

	/**
	 * @param manualSpc the manualSpc to set
	 */
	public void setManualSpc(boolean manualSpc) {
		this.manualSpc = manualSpc;
	}

	/**
	 * @param condition the condition to set
	 */
	public void setQuantCondition(QuantCondition condition) {
		quantCondition = condition;
	}

	@Override
	public Condition getCondition() {
		if (condition != null) {
			return condition;
		}
		return quantCondition;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof QuantAmount) {
			final boolean equals = toString().equals(obj.toString());
			return equals;
		}
		return super.equals(obj);
	}

	@Override
	public int hashCode() {
		if (hashCode == -1) {
			hashCode = HashCodeBuilder.reflectionHashCode(toString(), false);
		}
		return hashCode;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "QuantAmount [value=" + value + ", label=" + label + ", amountType=" + amountType + ", combinationType="
				+ combinationType + ", singleton=" + singleton + ", manualSpc=" + manualSpc + ", quantCondition="
				+ quantCondition + "]";
	}

	/**
	 * @param condition the condition to set
	 */
	@Override
	public void setCondition(Condition condition) {
		this.condition = condition;
	}

	public QuantificationLabel getLabel() {
		return label;
	}

	public void setLabel(QuantificationLabel label) {
		this.label = label;
	}

}
