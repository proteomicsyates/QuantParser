package edu.scripps.yates.census.read.model.interfaces;

import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.utilities.proteomicsmodel.HasAmounts;
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
	 * Get a list of {@link PositionInProtein}.<br>
	 * Examples:<br>
	 * {"PROTEIN1_K234, PROTEIN2_K123"}
	 * 
	 * @param quantifiedAAs
	 * @param uplr
	 *            used in order to get the protein sequence
	 * @return
	 */
	public List<PositionInProtein> getKeysForQuantifiedAAs(char[] quantifiedAAs, UniprotProteinLocalRetriever uplr);

	public Map<Character, List<PositionInProtein>> getPositionInProteinForSites(char[] quantifiedAAs,
			UniprotProteinLocalRetriever uplr);

}
