package edu.scripps.yates.census.read.util;

import java.util.HashMap;

import org.apache.log4j.Logger;

/**
 * Class wrapper of a map of protein sequences, that will be needed to check
 * whether
 * 
 * @author salvador
 *
 */
public class ProteinSequences extends HashMap<String, String> {
	private final static Logger log = Logger.getLogger(ProteinSequences.class);
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	@Override
	public String put(String key, String value) {
		// if (key.contains("contaminant")) {
		// log.info(key + "\t" + value);
		// }
		return super.put(key, value);
	}

}
