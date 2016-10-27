package edu.scripps.yates.census.read.model;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.model.interfaces.StaticItemStorage;

public class StaticQuantMaps {
	private final static Logger log = Logger.getLogger(StaticQuantMaps.class);
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
		log.info("Clearing static quantitative maps");
		if (!proteinMap.isEmpty()) {
			log.info("Clearing quant protein map");
		}
		proteinMap.clear();
		if (!psmMap.isEmpty()) {
			log.info("Clearing quant PSM map");
		}
		psmMap.clear();
		if (!peptideMap.isEmpty()) {
			log.info("Clearing quant peptide map");
		}
		peptideMap.clear();
	}

}
