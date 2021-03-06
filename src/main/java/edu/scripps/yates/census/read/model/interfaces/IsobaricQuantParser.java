package edu.scripps.yates.census.read.model.interfaces;

import java.util.Collection;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.census.read.QuantParserException;
import edu.scripps.yates.census.read.model.IonSerie.IonSerieType;
import edu.scripps.yates.census.read.util.IonExclusion;

public interface IsobaricQuantParser extends QuantParser {

	Map<String, Set<String>> getSpectrumToIonsMap() throws QuantParserException;

	void addIonExclusion(IonSerieType serieType, int ionNumber);

	void addIonExclusions(Collection<IonExclusion> ionExclusions);

}
