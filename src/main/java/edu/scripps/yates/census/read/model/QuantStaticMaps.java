package edu.scripps.yates.census.read.model;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.model.interfaces.StaticItemStorage;

public class QuantStaticMaps {
	/**
	 * The map that stores the {@link QuantifiedPSMInterface} by the spectrum
	 * key
	 */
	public static StaticItemStorage<QuantifiedPSMInterface> psmMap = new StaticItemStorage<QuantifiedPSMInterface>();

	// public static Map<String, QuantifiedPSMInterface> psmMap = new
	// HashMap<String, QuantifiedPSMInterface>();

	/**
	 * The map that stores the {@link QuantifiedProteinInterface} by the protein
	 * key
	 */
	public static StaticItemStorage<QuantifiedProteinInterface> proteinMap = new StaticItemStorage<QuantifiedProteinInterface>();

	/**
	 * The map that stores the {@link QuantifiedPeptide} by the peptide key
	 */
	public static StaticItemStorage<QuantifiedPeptideInterface> peptideMap = new StaticItemStorage<QuantifiedPeptideInterface>();

	public static void clearInfo() {
		proteinMap.clear();
		psmMap.clear();
		peptideMap.clear();
	}

}
