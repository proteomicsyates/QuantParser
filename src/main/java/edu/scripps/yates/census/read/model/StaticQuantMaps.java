package edu.scripps.yates.census.read.model;

import java.util.Set;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock.ReadLock;
import java.util.concurrent.locks.ReentrantReadWriteLock.WriteLock;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.model.interfaces.StaticItemStorage;
import gnu.trove.set.hash.THashSet;

public class StaticQuantMaps {
	private final static Logger log = Logger.getLogger(StaticQuantMaps.class);

	private static Set<String> rawFileNames = new THashSet<String>();
	private final static ReentrantReadWriteLock rawFileNamesLock = new ReentrantReadWriteLock();
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
			log.info("Clearing quant protein map of " + proteinMap.size());
		}
		proteinMap.clear();

		if (!psmMap.isEmpty()) {
			log.info("Clearing quant PSM map of " + psmMap.size());
		}
		psmMap.clear();

		if (!peptideMap.isEmpty()) {
			log.info("Clearing quant peptide map " + peptideMap.size());
		}
		peptideMap.clear();

		WriteLock writeLock = rawFileNamesLock.writeLock();
		try {
			writeLock.lock();
			rawFileNames.clear();
		} finally {
			writeLock.unlock();
		}
	}

	public static boolean rawFileNamesContains(String rawFileName) {
		ReadLock readLock = rawFileNamesLock.readLock();
		try {
			readLock.lock();
			return rawFileNames.contains(rawFileName);
		} finally {
			readLock.unlock();
		}
	}

	public static void addRawFileName(String rawFileName) {
		WriteLock writeLock = rawFileNamesLock.writeLock();
		try {
			writeLock.lock();
			rawFileNames.add(rawFileName);
		} finally {
			writeLock.unlock();
		}
	}

}
