package edu.scripps.yates.census.read.util;

import java.util.Collection;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.utilities.trove.CollectionObserver;

public class QuantifiedPSMCollectionObserver extends CollectionObserver<QuantifiedPSMInterface> {
	private final Collection<QuantifiedPeptideInterface> quantifiedPeptides;

	public QuantifiedPSMCollectionObserver(Collection<QuantifiedPeptideInterface> quantifiedPeptides) {
		this.quantifiedPeptides = quantifiedPeptides;
	}

	@Override
	public boolean add(QuantifiedPSMInterface obj) {
		return quantifiedPeptides.add(obj.getQuantifiedPeptide());
	}

	@Override
	public void clear() {
		quantifiedPeptides.clear();
	}

	@Override
	public boolean remove(Object obj) {
		if (obj instanceof QuantifiedPSMInterface) {
			return quantifiedPeptides.remove(((QuantifiedPSMInterface) obj).getQuantifiedPeptide());
		}
		return false;
	}

	@Override
	public boolean addAll(Collection<? extends QuantifiedPSMInterface> collection) {
		boolean ret = false;
		for (final QuantifiedPSMInterface quantifiedPSMInterface : collection) {
			final boolean b = quantifiedPeptides.add(quantifiedPSMInterface.getQuantifiedPeptide());
			if (b) {
				ret = true;
			}
		}
		return ret;
	}

	@Override
	public boolean removeAll(Collection<?> collection) {
		boolean ret = true;
		for (final Object obj : collection) {
			if (obj instanceof QuantifiedPSMInterface) {
				final boolean b = quantifiedPeptides.remove(((QuantifiedPSMInterface) obj).getQuantifiedPeptide());
				if (!b) {
					ret = false;
				}
			}
		}
		return ret;
	}

}
