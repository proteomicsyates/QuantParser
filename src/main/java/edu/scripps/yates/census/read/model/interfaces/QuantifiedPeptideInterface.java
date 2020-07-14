package edu.scripps.yates.census.read.model.interfaces;

import java.util.Set;

import edu.scripps.yates.utilities.proteomicsmodel.HasKey;
import edu.scripps.yates.utilities.proteomicsmodel.Peptide;

/**
 * A quantified peptide that can be coming from different msruns (filenames),
 * and can contain psms from different msruns
 *
 * @author Salva
 *
 */
public interface QuantifiedPeptideInterface
		extends Peptide, HasQuantifiedProteins, HasQuantifiedPSMs, HasQuantRatios, QuantifiedItem, HasKey, Discardable {

	public Set<QuantifiedProteinInterface> getNonDiscardedQuantifiedProteins();

}
