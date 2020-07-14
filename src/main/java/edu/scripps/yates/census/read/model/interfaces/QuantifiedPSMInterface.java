package edu.scripps.yates.census.read.model.interfaces;

import edu.scripps.yates.utilities.proteomicsmodel.HasKey;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;

/**
 * A quantified psm that is quantified in a particular MS Run (one fileName)
 *
 * @author Salva
 *
 */
public interface QuantifiedPSMInterface
		extends PSM, HasKey, QuantifiedItem, HasQuantRatios, HasQuantifiedPeptide, HasQuantifiedProteins, Discardable {

	/**
	 * Gets the peak twith more intensity
	 *
	 * @return
	 */
	public Double getMaxPeak();

	public boolean isSingleton();

	public void setSingleton(boolean singleton);

}
