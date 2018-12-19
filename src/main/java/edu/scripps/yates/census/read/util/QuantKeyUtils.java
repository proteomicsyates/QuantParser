package edu.scripps.yates.census.read.util;

import edu.scripps.yates.census.quant.xml.ProteinType;
import edu.scripps.yates.census.read.model.IsoRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.Peptide;
import edu.scripps.yates.utilities.proteomicsmodel.utils.KeyUtils;

public class QuantKeyUtils extends KeyUtils {
	private static QuantKeyUtils instance;

	private QuantKeyUtils() {
		super();
	}

	public static QuantKeyUtils getInstance() {
		if (instance == null) {
			instance = new QuantKeyUtils();
		}
		return instance;
	}

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
	public String getIonKey(IsoRatio ratio, ProteinType.Peptide peptide, boolean chargeSensible) {

		return ratio.getIonSerieType().name() + ratio.getNumIon() + "-" + getSpectrumKey(peptide, chargeSensible);

	}

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
	public String getIonKey(IsoRatio ratio, Peptide peptide) {

		return ratio.getIonSerieType().name() + ratio.getNumIon() + "-" + peptide.getFullSequence();

	}

	public String getProteinKey(ProteinType protein, boolean ignoreACCFormat) {
		if (ignoreACCFormat) {
			if (protein.getLocus().contains(" ")) {
				return protein.getLocus().substring(0, protein.getLocus().indexOf(" "));
			} else {
				return protein.getLocus();
			}
		}
		return FastaParser.getACC(protein.getLocus()).getAccession();
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
	@Override
	public String getSpectrumKey(PSM psm, boolean chargeSensible) {

		final StringBuilder sb = new StringBuilder();
		if (psm instanceof QuantifiedPSMInterface) {
			if (((QuantifiedPSMInterface) psm).getRawFileNames() != null
					&& !((QuantifiedPSMInterface) psm).getRawFileNames().isEmpty())
				sb.append(((QuantifiedPSMInterface) psm).getRawFileNames().iterator().next());
			if (!"".equals(sb.toString())) {
				sb.append("-");
			}
		}
		if (psm.getScanNumber() != null) {
			sb.append(psm.getScanNumber());
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
			sb.append(psm.getChargeState());
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
	public String getSpectrumKey(ProteinType.Peptide peptide, boolean chargeSensible) {

		final String string = peptide.getFile() + "-" + peptide.getScan() + "-" + peptide.getSeq();
		if (chargeSensible) {
			return string + "-" + peptide.getCharge();
		} else {
			return string;
		}
	}
}
