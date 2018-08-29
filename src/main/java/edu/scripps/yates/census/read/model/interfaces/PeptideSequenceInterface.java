package edu.scripps.yates.census.read.model.interfaces;

import java.util.List;
import java.util.Set;

import edu.scripps.yates.utilities.sequence.PTMInPeptide;

public interface PeptideSequenceInterface extends HasKey {
	public String getSequence();

	public String getFullSequence();

	public Float getCalcMHplus();

	public Set<String> getTaxonomies();

	public Float getMHplus();

	public boolean containsPTMs();

	public List<PTMInPeptide> getPtms();
}
