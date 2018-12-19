package edu.scripps.yates.census.read.model.interfaces;

import java.util.Set;

public interface HasQuantifiedProteins {
	public Set<QuantifiedProteinInterface> getQuantifiedProteins();

	public boolean addQuantifiedProtein(QuantifiedProteinInterface protein, boolean recursively);
}
