package edu.scripps.yates.census.read.model.interfaces;

import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.utilities.proteomicsmodel.HasAmounts;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import edu.scripps.yates.utilities.sequence.PositionInProtein;

/**
 * A quantified peptide that can be coming from different msruns (filenames),
 * and can contain psms from different msruns
 *
 * @author Salva
 *
 */
public interface QuantifiedPeptideInterface extends PeptideSequenceInterface, HasRatios, HasAmounts, QuantifiedItem {

	public Set<QuantifiedProteinInterface> getQuantifiedProteins();

	public boolean addQuantifiedProtein(QuantifiedProteinInterface protein, boolean recursive);

	public Set<QuantifiedProteinInterface> getNonDiscardedQuantifiedProteins();

	public Set<QuantifiedPSMInterface> getQuantifiedPSMs();

	public boolean addQuantifiedPSM(QuantifiedPSMInterface psm, boolean recursive);

	/**
	 * Get a list of {@link PositionInProtein} for each quantified site in the
	 * peptide sequence (represented as a {@link PositionInPeptide}).<br>
	 * Examples:<br>
	 * "PEPTIDE#4 {PROTEIN1#234, PROTEIN2#123}
	 * 
	 * @param quantifiedAAs
	 * @param uplr
	 *            used in order to get the protein sequence
	 * @return
	 */
	public Map<PositionInPeptide, List<PositionInProtein>> getProteinKeysByPeptideKeysForQuantifiedAAs(
			char[] quantifiedAAs, UniprotProteinLocalRetriever uplr);

	public Map<Character, List<PositionInProtein>> getPositionInProteinForSites(char[] quantifiedAAs,
			UniprotProteinLocalRetriever uplr);

}
