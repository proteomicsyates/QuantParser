package edu.scripps.yates.census.read.model.interfaces;

import java.util.Set;

public interface HasQuantifiedPSMs {
	public Set<QuantifiedPSMInterface> getQuantifiedPSMs();

	public boolean addQuantifiedPSM(QuantifiedPSMInterface psm, boolean recursively);
}
