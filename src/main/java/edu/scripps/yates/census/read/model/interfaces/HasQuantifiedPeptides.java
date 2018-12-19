package edu.scripps.yates.census.read.model.interfaces;

import java.util.Set;

public interface HasQuantifiedPeptides {
	public Set<QuantifiedPeptideInterface> getQuantifiedPeptides();

	public boolean addQuantifiedPeptide(QuantifiedPeptideInterface peptide, boolean recursively);
}
