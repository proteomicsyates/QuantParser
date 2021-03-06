package edu.scripps.yates.census.read;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.IonSerie.IonSerieType;
import edu.scripps.yates.census.read.model.interfaces.IsobaricQuantParser;
import edu.scripps.yates.census.read.util.IonExclusion;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public abstract class AbstractIsobaricQuantParser extends AbstractQuantParser implements IsobaricQuantParser {
	// key=spectrumKeys, values=ionKeys
	protected final Map<String, Set<String>> spectrumToIonsMap = new THashMap<String, Set<String>>();
	protected final Set<String> ionKeys = new THashSet<String>();

	// ion exclusions
	protected final List<IonExclusion> ionExclusions = new ArrayList<IonExclusion>();

	public AbstractIsobaricQuantParser() {
		super();
	}

	public AbstractIsobaricQuantParser(List<RemoteSSHFileReference> remoteSSHServers,
			List<Map<QuantificationLabel, QuantCondition>> conditionsByLabels, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) {
		super(remoteSSHServers, conditionsByLabels, labelNumerator, labelDenominator);
	}

	public AbstractIsobaricQuantParser(Map<QuantificationLabel, QuantCondition> conditionsByLabels,
			Collection<RemoteSSHFileReference> remoteSSHServers, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) {
		super(conditionsByLabels, remoteSSHServers, labelNumerator, labelDenominator);
	}

	public AbstractIsobaricQuantParser(RemoteSSHFileReference remoteSSHServer,
			Map<QuantificationLabel, QuantCondition> conditionsByLabels, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) throws FileNotFoundException {
		super(remoteSSHServer, conditionsByLabels, labelNumerator, labelDenominator);
	}

	public AbstractIsobaricQuantParser(File xmlFile, Map<QuantificationLabel, QuantCondition> conditionsByLabels,
			QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) throws FileNotFoundException {
		super(xmlFile, conditionsByLabels, labelNumerator, labelDenominator);
	}

	public AbstractIsobaricQuantParser(File[] xmlFiles, Map<QuantificationLabel, QuantCondition> conditionsByLabels,
			QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) throws FileNotFoundException {
		super(xmlFiles, conditionsByLabels, labelNumerator, labelDenominator);
	}

	public AbstractIsobaricQuantParser(File[] xmlFiles, Map<QuantificationLabel, QuantCondition>[] conditionsByLabels,
			QuantificationLabel[] labelNumerator, QuantificationLabel[] labelDenominator) throws FileNotFoundException {
		super(xmlFiles, conditionsByLabels, labelNumerator, labelDenominator);
	}

	public AbstractIsobaricQuantParser(Collection<File> xmlFiles,
			Map<QuantificationLabel, QuantCondition> conditionsByLabels, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) throws FileNotFoundException {
		super(xmlFiles, conditionsByLabels, labelNumerator, labelDenominator);
	}

	public AbstractIsobaricQuantParser(RemoteSSHFileReference remoteServer, QuantificationLabel label1,
			QuantCondition cond1, QuantificationLabel label2, QuantCondition cond2) {
		super(remoteServer, label1, cond1, label2, cond2);
	}

	public AbstractIsobaricQuantParser(File inputFile, QuantificationLabel label1, QuantCondition cond1,
			QuantificationLabel label2, QuantCondition cond2) throws FileNotFoundException {
		super(inputFile, label1, cond1, label2, cond2);
	}

	/**
	 * Gets the set of spectrum to ions map
	 *
	 * @return the spectrumMap
	 * @throws IOException
	 */
	@Override
	public Map<String, Set<String>> getSpectrumToIonsMap() throws QuantParserException {
		if (!processed)
			process();
		return spectrumToIonsMap;
	}

	/**
	 * @param excludeY1Ions the excludeY1Ions to set
	 */
	@Override
	public void addIonExclusion(IonSerieType serieType, int ionNumber) {
		ionExclusions.add(new IonExclusion(serieType, ionNumber));
	}

	/**
	 * @param excludeY1Ions the excludeY1Ions to set
	 */
	public void addIonExclusion(IonExclusion ionExclusion) {
		ionExclusions.add(ionExclusion);
	}

	/**
	 * @param excludeY1Ions the excludeY1Ions to set
	 */
	@Override
	public void addIonExclusions(Collection<IonExclusion> ionExclusions) {
		this.ionExclusions.addAll(ionExclusions);
	}

	/**
	 * @return the ionExclusions
	 */
	public List<IonExclusion> getIonExclusions() {
		return ionExclusions;
	}

	@Override
	protected abstract void process() throws QuantParserException;

}
