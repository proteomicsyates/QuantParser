package edu.scripps.yates.census.read.model.interfaces;

import java.util.Set;

import edu.scripps.yates.utilities.proteomicsmodel.HasKey;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;

/**
 * A quantified protein that can be quantified in different msruns and
 * experiments.
 *
 * @author Salva
 *
 */
public interface QuantifiedProteinInterface
		extends Protein, HasQuantRatios, HasKey, QuantifiedItem, HasQuantifiedPeptides, HasQuantifiedPSMs, Discardable {

	public Set<QuantifiedPeptideInterface> getNonDiscardedQuantifiedPeptides();

}
