package edu.scripps.yates.census.read.model.interfaces;

public interface HasQuantifiedPeptide {
	public QuantifiedPeptideInterface getQuantifiedPeptide();

	public boolean setQuantifiedPeptide(QuantifiedPeptideInterface peptide, boolean recursively);
}
