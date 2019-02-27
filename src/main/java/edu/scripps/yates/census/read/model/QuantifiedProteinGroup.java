package edu.scripps.yates.census.read.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AggregationLevel;
import gnu.trove.set.hash.THashSet;

public class QuantifiedProteinGroup extends AbstractContainsQuantifiedPSMs {
	private static final String SEPARATOR = " ## ";
	protected final Set<QuantifiedProteinInterface> proteins = new THashSet<QuantifiedProteinInterface>();
	private StringBuilder accessionString;

	public QuantifiedProteinGroup(ProteinGroup proteinGroup) {
		for (final GroupableProtein groupableProtein : proteinGroup) {
			if (groupableProtein instanceof IsobaricQuantifiedProtein) {
				proteins.add((IsobaricQuantifiedProtein) groupableProtein);
			}
		}

	}

	public int size() {
		return getProteins().size();
	}

	@Override
	public Set<QuantifiedPSMInterface> getQuantifiedPSMs() {
		final Set<QuantifiedPSMInterface> ret = new THashSet<QuantifiedPSMInterface>();
		for (final QuantifiedProteinInterface quantifiedProtein : getProteins()) {
			ret.addAll(quantifiedProtein.getQuantifiedPSMs());
		}

		return ret;
	}

	public String getAccessionString() {
		if (accessionString == null) {
			accessionString = new StringBuilder();
			for (final QuantifiedProteinInterface quantifiedProtein : getProteins()) {
				if (!"".equals(accessionString.toString()))
					accessionString.append(SEPARATOR);
				accessionString.append(quantifiedProtein.getAccession());
			}
		}
		return accessionString.toString();
	}

	public String getDescriptionString() {
		final StringBuilder sb = new StringBuilder();
		for (final QuantifiedProteinInterface quantifiedProtein : getProteins()) {
			if (!"".equals(sb.toString()))
				sb.append(SEPARATOR);
			sb.append(quantifiedProtein.getDescription());
		}
		return sb.toString();
	}

	public List<String> getTaxonomies() {
		final List<String> ret = new ArrayList<String>();
		for (final QuantifiedProteinInterface quantifiedProtein : getProteins()) {
			if (quantifiedProtein.getTaxonomies() != null) {
				ret.addAll(quantifiedProtein.getTaxonomies());
			}
		}
		return ret;
	}

	/**
	 * NOTE THAT THIS RETURNED LIST IS NOT VALID FOR ADDING NEW PROTEINS TO THE
	 * GROUP
	 *
	 * @return the proteins
	 */
	public List<QuantifiedProteinInterface> getProteins() {
		final List<QuantifiedProteinInterface> ret = new ArrayList<QuantifiedProteinInterface>();
		ret.addAll(proteins);
		Collections.sort(ret, new Comparator<QuantifiedProteinInterface>() {

			@Override
			public int compare(QuantifiedProteinInterface o1, QuantifiedProteinInterface o2) {
				final String accession1 = o1.getAccession();
				final String accession2 = o2.getAccession();
				return accession1.compareTo(accession2);
			}
		});
		return ret;
	}

	public String getGeneNameString() {
		final StringBuilder sb = new StringBuilder();
		for (final QuantifiedProteinInterface quantifiedProtein : getProteins()) {
			if (!"".equals(sb.toString()))
				sb.append(SEPARATOR);
			String geneFromFastaHeader = FastaParser.getGeneFromFastaHeader(quantifiedProtein.getAccession());
			if (geneFromFastaHeader == null) {
				geneFromFastaHeader = FastaParser.getGeneFromFastaHeader(quantifiedProtein.getDescription());
			}
			sb.append(geneFromFastaHeader);
		}
		return sb.toString();
	}

	public Set<String> getFileNames() {
		final Set<String> ret = new THashSet<String>();
		for (final QuantifiedProteinInterface quantprotein : proteins) {
			ret.addAll(quantprotein.getFileNames());
		}
		return ret;
	}

	@Override
	public QuantRatio getConsensusRatio(QuantCondition quantConditionNumerator,
			QuantCondition quantConditionDenominator) {
		return QuantUtils.getAverageRatio(QuantUtils.getNonInfinityRatios(getQuantRatios()),
				AggregationLevel.PROTEINGROUP);
	}

}
