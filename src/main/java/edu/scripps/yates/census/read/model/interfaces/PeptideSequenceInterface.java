package edu.scripps.yates.census.read.model.interfaces;

import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.utilities.sequence.PTMInPeptide;
import edu.scripps.yates.utilities.sequence.PTMInProtein;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import edu.scripps.yates.utilities.sequence.PositionInProtein;

public interface PeptideSequenceInterface extends HasKey {
	public String getSequence();

	public String getFullSequence();

	public Float getCalcMHplus();

	public Set<String> getTaxonomies();

	public Float getMHplus();

	public boolean containsPTMs();

	public List<PTMInPeptide> getPtms();

	/**
	 * Get a list of {@link PositionInProtein} for each ptm site in the peptide
	 * sequence (represented as a {@link PositionInPeptide}).<br>
	 * Examples:<br>
	 * "PEPTIDE#4 {PROTEIN1#234#238, PROTEIN2#123#127}
	 * 
	 * @param uplr
	 *            used in order to get the protein sequence
	 * @param proteinSequences
	 *            map of protein sequences
	 * @param dBIndex
	 * @return
	 */
	public List<PTMInProtein> getPTMInProtein(UniprotProteinLocalRetriever uplr, Map<String, String> proteinSequences);
}
