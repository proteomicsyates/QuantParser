package edu.scripps.yates.census.analysis.util;

import java.util.Collections;
import java.util.Comparator;

import edu.scripps.yates.census.quant.xml.ProteinType;
import edu.scripps.yates.census.quant.xml.ProteinType.Peptide;
import edu.scripps.yates.census.read.model.IsoRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.ProteinGroup;

public class KeyUtils {
	public static String getSequenceChargeKey(QuantifiedPSMInterface psm, boolean chargeStateSensible) {
		final String seq = FastaParser.cleanSequence(psm.getFullSequence());
		if (chargeStateSensible) {
			return seq + "-" + psm.getCharge();
		} else {
			return seq;
		}
	}

	// public static String getProteinKey(DTASelectProtein protein) {
	// final String locus = protein.getLocus();
	// return FastaParser.getACC(locus).getFirstelement();
	// }

	// public static String getProteinKey(String locus) {
	// return FastaParser.getACC(locus).getFirstelement();
	// }
	public static String getProteinKey(ProteinType protein, boolean ignoreACCFormat) {
		if (ignoreACCFormat) {
			if (protein.getLocus().contains(" ")) {
				return protein.getLocus().substring(0, protein.getLocus().indexOf(" "));
			} else {
				return protein.getLocus();
			}
		}
		return FastaParser.getACC(protein.getLocus()).getFirstelement();
	}

	public static String getProteinKey(IndexedProtein indexedProtein, boolean ignoreACCFormat) {
		String fastaDefLine = indexedProtein.getFastaDefLine();

		if (ignoreACCFormat) {
			if (fastaDefLine.startsWith(">")) {
				fastaDefLine = fastaDefLine.substring(1);
			}
			if (fastaDefLine.contains(" ")) {
				return fastaDefLine.substring(0, fastaDefLine.indexOf(" "));
			} else {
				return fastaDefLine;
			}
		}
		return FastaParser.getACC(fastaDefLine).getFirstelement();
	}

	/**
	 *
	 * @param psm
	 * @param chargeSensible
	 *            if true, then, the charge will be considered for
	 *            differentiating peptides with different charge states. If
	 *            false, peptides with different charge states will have the
	 *            same key
	 * @return
	 */
	public static String getSpectrumKey(QuantifiedPSMInterface psm, boolean chargeSensible) {

		final StringBuilder sb = new StringBuilder();
		if (psm.getRawFileNames() != null && !psm.getRawFileNames().isEmpty())
			sb.append(psm.getRawFileNames().iterator().next());
		if (!"".equals(sb.toString())) {
			sb.append("-");
		}
		if (psm.getScan() != null) {
			sb.append(psm.getScan());
		}
		if (!"".equals(sb.toString())) {
			sb.append("-");
		}
		if (psm.getFullSequence() != null) {
			sb.append(psm.getFullSequence());
		}

		if (chargeSensible) {
			if (!"".equals(sb.toString())) {
				sb.append("-");
			}
			sb.append(psm.getCharge());
		}
		return sb.toString();
	}

	/**
	 *
	 * @param psm
	 * @param chargeSensible
	 *            if true, then, the charge will be considered for
	 *            differentiating peptides with different charge states. If
	 *            false, peptides with different charge states will have the
	 *            same key
	 * @return
	 */
	public static String getSpectrumKey(Peptide peptide, boolean chargeSensible) {

		final String string = peptide.getFile() + "-" + peptide.getScan() + "-" + peptide.getSeq();
		if (chargeSensible) {
			return string + "-" + peptide.getCharge();
		} else {
			return string;
		}
	}

	// /**
	// * returns ionSerieTypeName+numIon+spectrumKey
	// *
	// * @param ratio
	// * @param peptide
	// * @param chargeSensible
	// * if true, then, the charge will be considered for
	// * differentiating peptides with different charge states. If
	// * false, peptides with different charge states will have the
	// * same key
	// * @return
	// */
	// public static String getIonKey(IsoRatio ratio, QuantifiedPSMInterface
	// quantPSM, boolean chargeSensible) {
	// String ionSerieTypeName = "";
	// if (ratio.getIonSerieType() != null) {
	// ionSerieTypeName = ratio.getIonSerieType().name();
	// }
	// String numIon = "";
	// if (ratio.getNumIon() > 0)
	// numIon = String.valueOf(ratio.getNumIon());
	// String ret = ionSerieTypeName + numIon;
	// if (!"".equals(ret)) {
	// ret += "-";
	// }
	// ret += getSpectrumKey(quantPSM, chargeSensible);
	// return ret;
	//
	// }

	/**
	 * returns ionSerieTypeName+numIon+spectrumKey
	 *
	 * @param ratio
	 * @param peptide
	 * @param chargeSensible
	 *            if true, then, the charge will be considered for
	 *            differentiating peptides with different charge states. If
	 *            false, peptides with different charge states will have the
	 *            same key
	 * @return
	 */
	public static String getIonKey(IsoRatio ratio, Peptide peptide, boolean chargeSensible) {

		return ratio.getIonSerieType().name() + ratio.getNumIon() + "-" + getSpectrumKey(peptide, chargeSensible);

	}

	public static String getGroupKey(ProteinGroup group) {
		Collections.sort(group, new Comparator<GroupableProtein>() {

			@Override
			public int compare(GroupableProtein o1, GroupableProtein o2) {
				return o1.getAccession().compareTo(o2.getAccession());
			}
		});
		String key = "";
		for (final GroupableProtein protein : group) {
			if (!"".equals(key))
				key += ",";
			key += protein.getAccession();
		}
		return key;
	}

	public static String getSequenceKey(QuantifiedPSMInterface quantPSM, boolean distinguishModifiedSequence) {
		if (distinguishModifiedSequence) {
			return quantPSM.getFullSequence();
		} else {
			return FastaParser.cleanSequence(quantPSM.getSequence());
		}
	}
}
