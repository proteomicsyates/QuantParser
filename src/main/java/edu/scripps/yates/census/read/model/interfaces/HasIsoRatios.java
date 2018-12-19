package edu.scripps.yates.census.read.model.interfaces;

import java.util.Map;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.Ion;
import edu.scripps.yates.census.read.model.IonCountRatio;
import edu.scripps.yates.census.read.model.IonSerie;
import edu.scripps.yates.census.read.model.IsoRatio;
import edu.scripps.yates.census.read.util.QuantificationLabel;

public interface HasIsoRatios extends HasQuantRatios {

	public Set<IsoRatio> getIsoRatios();

	public Set<IsoRatio> getNonInfinityIsoRatios();

	public Map<QuantificationLabel, Set<Ion>> getSingletonIonsByLabel();

	public Map<QuantificationLabel, Set<Ion>> getIonsByLabel();

	public Set<Ion> getIonsByLabel(QuantificationLabel label);

	public Map<QuantCondition, Set<Ion>> getSingletonIonsByCondition();

	public Map<QuantificationLabel, Set<Ion>> getSingletonIons(IonSerie ionSerie);

	public Map<QuantCondition, Set<Ion>> getIonsByCondition();

	public Map<QuantCondition, Set<Ion>> getIonsByCondition(String replicateName);

	public Map<QuantCondition, Set<Ion>> getIonsByConditionForSites(String replicateName, char[] quantifiedAAs);

	public Map<QuantCondition, Set<Ion>> getIonsByConditionForSites(String replicateName, char[] quantifiedAAs,
			int positionInPeptide);

	public Map<QuantificationLabel, Set<Ion>> getIons(IonSerie ionSerie);

	public boolean containsAnySingletonIon(QuantificationLabel label);

	public boolean containsAnyIon(QuantificationLabel label);

	public Set<Ion> getSingletonIonsByLabel(QuantificationLabel label);

	public IonCountRatio getIonCountRatio(QuantCondition condition1, QuantCondition condition2);

	public IonCountRatio getIonCountRatio(QuantCondition cond1, QuantCondition cond2, String replicateName);

	Double getMaxPeak();

}
