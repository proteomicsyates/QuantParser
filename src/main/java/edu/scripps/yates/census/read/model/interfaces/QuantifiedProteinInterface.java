package edu.scripps.yates.census.read.model.interfaces;

import java.util.Set;

import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.proteomicsmodel.HasAmounts;

/**
 * A quantified protein that can be quantified in different msruns and
 * experiments.
 *
 * @author Salva
 *
 */
public interface QuantifiedProteinInterface extends GroupableProtein, HasRatios, HasAmounts, HasKey, QuantifiedItem {

	public String getDescription();

	public Set<String> getTaxonomies();

	public Set<QuantifiedPeptideInterface> getQuantifiedPeptides();

	public Set<QuantifiedPeptideInterface> getNonDiscardedQuantifiedPeptides();

	public Set<QuantifiedPSMInterface> getQuantifiedPSMs();

	public boolean addPeptide(QuantifiedPeptideInterface peptide, boolean recursive);

	public String getAccessionType();

	public boolean addPSM(QuantifiedPSMInterface quantifiedPSM, boolean recursive);

	public void setTaxonomy(String taxonomy);

	public Integer getLength();

	public void setAccession(String primaryAccession);

	public void setDescription(String description);

}
