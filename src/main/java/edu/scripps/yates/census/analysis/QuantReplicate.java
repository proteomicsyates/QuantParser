package edu.scripps.yates.census.analysis;

import java.io.FileNotFoundException;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.census.read.model.IonSerie.IonSerieType;
import edu.scripps.yates.census.read.model.interfaces.IsobaricQuantParser;
import edu.scripps.yates.census.read.model.interfaces.QuantParser;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import gnu.trove.map.hash.THashMap;

public class QuantReplicate {
	private final Map<QuantificationLabel, QuantCondition> conditionsByLabels = new THashMap<QuantificationLabel, QuantCondition>();
	private final String name;
	private final QuantParser parser;

	public QuantReplicate(String name, QuantParser quantParser,
			Map<QuantificationLabel, QuantCondition> conditionsByLabels) throws FileNotFoundException {
		this.name = name;
		parser = quantParser;

		this.conditionsByLabels.putAll(conditionsByLabels);
	}

	public void addIonExclusion(IonSerieType serieType, int ionNumber) {
		if (parser instanceof IsobaricQuantParser) {
			((IsobaricQuantParser) parser).addIonExclusion(serieType, ionNumber);
		} else {
			throw new IllegalArgumentException("addIonExclusion is only valid for IsobaricQuantParser");
		}
	}

	/**
	 * @return the labelsByConditions
	 */
	public Map<QuantCondition, Set<QuantificationLabel>> getLabelsByConditions() {
		return QuantUtils.getLabelsByConditions(conditionsByLabels);
	}

	/**
	 * @return the conditionsByLabels
	 */
	public Map<QuantificationLabel, QuantCondition> getConditionsByLabels() {
		return conditionsByLabels;
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @return the parser
	 */
	public QuantParser getParser() {
		return parser;
	}

}
