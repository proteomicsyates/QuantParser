package edu.scripps.yates.census.read;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.CensusRatio;
import edu.scripps.yates.census.read.model.QuantifiedPSM;
import edu.scripps.yates.census.read.model.QuantifiedPeptide;
import edu.scripps.yates.census.read.model.QuantifiedProtein;
import edu.scripps.yates.census.read.model.QuantifiedProteinFromDBIndexEntry;
import edu.scripps.yates.census.read.model.RatioScore;
import edu.scripps.yates.census.read.model.StaticQuantMaps;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.dbindex.util.PeptideNotFoundInDBIndexException;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AggregationLevel;
import edu.scripps.yates.utilities.proteomicsmodel.utils.KeyUtils;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import edu.scripps.yates.utilities.strings.StringUtils;
import gnu.trove.list.array.TIntArrayList;

public class SeparatedValuesParser extends AbstractQuantParser {
	private final static Logger log = Logger.getLogger(SeparatedValuesParser.class);
	private final String separator;
	private static int scanNumberTMP = 0;
	private final static int PSM_ID_COL = 0;
	private final static int SEQ_COL = 1;
	private final static int RATIO_COL = 2;
	private final static int RATIO_WEIGHT_COL = 3;
	private final static int PROTEIN_ACC_COL = 4;
	public static final String RATIO_WEIGHT = "Ratio initial weigh";
	public static final String RATIO = "RATIO";

	public SeparatedValuesParser(String separator) {
		super();
		this.separator = separator;
	}

	public SeparatedValuesParser(List<RemoteSSHFileReference> remoteSSHServers, String separator,
			List<Map<QuantificationLabel, QuantCondition>> conditionsByLabels, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator, boolean ignoreTaxonomies) {
		super(remoteSSHServers, conditionsByLabels, labelNumerator, labelDenominator);
		this.separator = separator;
		setIgnoreTaxonomies(ignoreTaxonomies);
	}

	public SeparatedValuesParser(String separator, Map<QuantificationLabel, QuantCondition> conditionsByLabels,
			Collection<RemoteSSHFileReference> remoteSSHServers, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator, boolean ignoreTaxonomies) {
		super(conditionsByLabels, remoteSSHServers, labelNumerator, labelDenominator);
		this.separator = separator;
		setIgnoreTaxonomies(ignoreTaxonomies);
	}

	public SeparatedValuesParser(RemoteSSHFileReference remoteSSHServer, String separator,
			Map<QuantificationLabel, QuantCondition> conditionsByLabels, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator, boolean ignoreTaxonomies) throws FileNotFoundException {
		super(remoteSSHServer, conditionsByLabels, labelNumerator, labelDenominator);
		this.separator = separator;
		setIgnoreTaxonomies(ignoreTaxonomies);
	}

	public SeparatedValuesParser(File xmlFile, String separator,
			Map<QuantificationLabel, QuantCondition> conditionsByLabels, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator, boolean ignoreTaxonomies) throws FileNotFoundException {
		super(xmlFile, conditionsByLabels, labelNumerator, labelDenominator);
		this.separator = separator;
		setIgnoreTaxonomies(ignoreTaxonomies);
	}

	public SeparatedValuesParser(File[] xmlFiles, String separator,
			Map<QuantificationLabel, QuantCondition> conditionsByLabels, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator, boolean ignoreTaxonomies) throws FileNotFoundException {
		super(xmlFiles, conditionsByLabels, labelNumerator, labelDenominator);
		this.separator = separator;
		setIgnoreTaxonomies(ignoreTaxonomies);
	}

	public SeparatedValuesParser(File[] xmlFiles, String separator,
			Map<QuantificationLabel, QuantCondition>[] conditionsByLabels, QuantificationLabel[] labelNumerator,
			QuantificationLabel[] labelDenominator, boolean ignoreTaxonomies) throws FileNotFoundException {
		super(xmlFiles, conditionsByLabels, labelNumerator, labelDenominator);
		this.separator = separator;
		setIgnoreTaxonomies(ignoreTaxonomies);
	}

	public SeparatedValuesParser(Collection<File> xmlFiles, String separator,
			Map<QuantificationLabel, QuantCondition> conditionsByLabels, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator, boolean ignoreTaxonomies) throws FileNotFoundException {
		super(xmlFiles, conditionsByLabels, labelNumerator, labelDenominator);
		this.separator = separator;
		setIgnoreTaxonomies(ignoreTaxonomies);
	}

	public SeparatedValuesParser(RemoteSSHFileReference remoteServer, String separator, QuantificationLabel label1,
			QuantCondition cond1, QuantificationLabel label2, QuantCondition cond2, boolean ignoreTaxonomies) {
		super(remoteServer, label1, cond1, label2, cond2);
		this.separator = separator;
		setIgnoreTaxonomies(ignoreTaxonomies);
	}

	public SeparatedValuesParser(File inputFile, String separator, QuantificationLabel label1, QuantCondition cond1,
			QuantificationLabel label2, QuantCondition cond2, boolean ignoreTaxonomies) throws FileNotFoundException {
		super(inputFile, label1, cond1, label2, cond2);
		this.separator = separator;
		setIgnoreTaxonomies(ignoreTaxonomies);
	}

	@Override
	protected void process() throws QuantParserException {
		processed = false;
		log.info("Processing file...");

		int numDecoy = 0;
		boolean someValidFile = false;
		for (final RemoteSSHFileReference remoteFileRetriever : remoteFileRetrievers) {
			final Map<QuantificationLabel, QuantCondition> conditionsByLabels = conditionsByLabelsByFile
					.get(remoteFileRetriever);

			final QuantificationLabel labelNumerator = ratioDescriptorsByFile.get(remoteFileRetriever).get(0)
					.getLabel1();
			final QuantificationLabel labelDenominator = ratioDescriptorsByFile.get(remoteFileRetriever).get(0)
					.getLabel2();

			final String experimentKey = FilenameUtils
					.getBaseName(remoteFileRetriever.getOutputFile().getAbsolutePath());
			log.info(experimentKey);
			// get all the Quantified PSMs first
			// Set<QuantifiedPSMInterface> psms = new
			// HashSet<QuantifiedPSMInterface>();
			final File remoteFile = remoteFileRetriever.getRemoteFile();
			if (remoteFile == null || !remoteFile.exists())
				continue;
			log.info("Reading " + remoteFile.getAbsolutePath());
			someValidFile = true;
			BufferedReader br = null;
			try {
				br = new BufferedReader(new FileReader(remoteFile));

				String line;

				int numLine = 0;
				while ((line = br.readLine()) != null) {
					numLine++;
					if ("".equals(line.trim())) {
						continue;
					}
					if (numLine > 1 && line.contains(separator)) {
						String psmID = null;
						String seq = null;
						Double ratio = null;
						Double ratioWeigth = null;
						String proteinAcc = null;
						final String[] split = line.split(separator);
						if (split.length > PSM_ID_COL) {
							psmID = split[PSM_ID_COL].trim();
						}
						if (split.length > SEQ_COL) {
							seq = split[SEQ_COL].trim();
						}
						if (split.length > RATIO_COL && !"".equals(split[RATIO_COL])) {
							try {
								ratio = Double.valueOf(split[RATIO_COL]);
							} catch (final NumberFormatException e) {
								e.printStackTrace();
								throw new IllegalArgumentException(
										"Error in line:" + numLine + ", col:" + RATIO_COL + 1 + "\t" + e.getMessage());
							}
						}
						if (split.length > RATIO_WEIGHT_COL && !"".equals(split[RATIO_WEIGHT_COL])) {
							try {
								ratioWeigth = Double.valueOf(split[RATIO_WEIGHT_COL]);
							} catch (final NumberFormatException e) {
								e.printStackTrace();
								throw new IllegalArgumentException("Error in line:" + numLine + ", col:"
										+ RATIO_WEIGHT_COL + 1 + "\t" + e.getMessage());
							}
						}
						if (split.length > PROTEIN_ACC_COL) {
							proteinAcc = split[PROTEIN_ACC_COL].trim();
						}
						// apply the pattern if available
						if (decoyPattern != null) {
							final Matcher matcher = decoyPattern.matcher(proteinAcc);
							if (matcher.find()) {
								log.debug("Discarding decoy: " + proteinAcc);
								numDecoy++;
								continue;
							}
						}
						processPSMLine(psmID, seq, ratio, ratioWeigth, proteinAcc, conditionsByLabels, labelNumerator,
								labelDenominator, experimentKey, remoteFileRetriever);
					}

				}
				br.close();
			} catch (final PeptideNotFoundInDBIndexException e) {
				if (!super.ignoreNotFoundPeptidesInDB) {
					throw e;
				}
			} catch (final IOException e) {
				e.printStackTrace();
				log.error(e.getMessage());
				throw new QuantParserException(e);
			} catch (final DBIndexStoreException e) {

				e.printStackTrace();
				log.error(e.getMessage());
				throw new QuantParserException(e);
			} catch (final Exception e) {
				e.printStackTrace();
				log.error(e.getMessage());
				throw e;
			} finally {
				if (br != null) {
					try {
						br.close();
					} catch (final IOException e) {
						throw new QuantParserException(e);
					}
				}
			}

			log.info("(" + experimentKey + ") " + localPsmMap.size() + " PSMs from this parser. "
					+ StaticQuantMaps.psmMap.size() + " PSMs in the system");
			log.info("(" + experimentKey + ") " + localProteinMap.size() + " Proteins from this parser. "
					+ StaticQuantMaps.proteinMap.size() + " Proteins in the system");
			log.info("(" + experimentKey + ") " + localPeptideMap.size() + " Peptides from this parser. "
					+ StaticQuantMaps.peptideMap.size() + " Peptides in the system");
			if (decoyPattern != null) {
				log.info(numDecoy + " decoy Proteins were discarded  in " + experimentKey);
			}

		}
		if (!someValidFile)
			throw new IllegalArgumentException("some error occurred while reading the files");

		processed = true;

	}

	private void processPSMLine(String psmId, String sequence, Double nonLogRatioValue, Double ratioWeigth,
			String proteinACC, Map<QuantificationLabel, QuantCondition> conditionsByLabels,
			QuantificationLabel labelNumerator, QuantificationLabel labelDenominator, String experimentKey,
			RemoteSSHFileReference remoteFileRetriever) throws IOException, DBIndexStoreException {

		// new psm

		// dont look into the QuantifiedPSM.map because each
		// line is always a new PSM
		final String inputFileName = FilenameUtils.getName(remoteFileRetriever.getOutputFile().getAbsolutePath());
		String rawFileName = FastaParser.getFileNameFromPSMIdentifier(psmId);
		if (rawFileName == null) {
			rawFileName = inputFileName;
		}
		StaticQuantMaps.addRawFileName(rawFileName);
		// scan number
		String scanNumber = "0";
		try {
			scanNumber = String.valueOf(Integer.valueOf(FastaParser.getScanFromPSMIdentifier(psmId)));
		} catch (final Exception e) {
			// get next available number
			scanNumber = String.valueOf(++scanNumberTMP);
		}

		// charge state
		int chargeState = 0;
		try {
			chargeState = Double.valueOf(FastaParser.getChargeStateFromPSMIdentifier(psmId)).intValue();
		} catch (final Exception e) {
		}
		QuantifiedPSMInterface quantifiedPSM = null;
		// if (!isGetPTMInProteinMap()) {
		// quantifiedPSM = new QuantifiedPSM(sequence, conditionsByLabels,
		// peptideToSpectraMap, scanNumber,
		// chargeState, rawFileName, false);
		// } else {
		quantifiedPSM = new QuantifiedPSM(sequence, peptideToSpectraMap, scanNumber, chargeState, rawFileName, false,
				isDistinguishModifiedSequences(), isChargeSensible());
		// }
		quantifiedPSM.getFileNames().add(inputFileName);
		final String psmKey = KeyUtils.getInstance().getSpectrumKey(quantifiedPSM, isDistinguishModifiedSequences(),
				isChargeSensible());
		// in case of TMT, the psm may have been created before
		if (StaticQuantMaps.psmMap.containsKey(psmKey)) {
			quantifiedPSM = StaticQuantMaps.psmMap.getItem(psmKey);
		}
		StaticQuantMaps.psmMap.addItem(quantifiedPSM);

		// psms.add(quantifiedPSM);
		// add to map
		if (!localPsmMap.containsKey(quantifiedPSM.getKey())) {
			localPsmMap.put(quantifiedPSM.getKey(), quantifiedPSM);
		}
		// PSM regular ratio
		// add PSM ratios from census out
		if (nonLogRatioValue != null) {
			try {
				final CensusRatio ratio = new CensusRatio(nonLogRatioValue, false, conditionsByLabels, labelNumerator,
						labelDenominator, AggregationLevel.PSM, RATIO);
				// set singleton
				if (nonLogRatioValue == 0 || Double.compare(Double.POSITIVE_INFINITY, nonLogRatioValue) == 0) {
					if (quantifiedPSM instanceof QuantifiedPSM) {
						((QuantifiedPSM) quantifiedPSM).setSingleton(true);
					}
				}
				if (ratioWeigth != null) {
					final RatioScore ratioScore = new RatioScore(String.valueOf(ratioWeigth), RATIO_WEIGHT,
							"PSM-level quantification confidence metric", RATIO_WEIGHT);
					ratio.setRatioScore(ratioScore);
				}

				// add ratio to PSM
				quantifiedPSM.addRatio(ratio);
				if (!getQuantifiedAAs().isEmpty()) {
					for (final Character c : getQuantifiedAAs()) {
						if (quantifiedPSM.getSequence().contains(String.valueOf(c))) {
							ratio.setQuantifiedAA(c);
						}
					}
					// check for ambiguity on the quantified site
					int numSites = 0;
					PositionInPeptide quantifiedSitePositionInPeptide = null;
					for (final Character c : getQuantifiedAAs()) {
						final TIntArrayList allPositionsOf = StringUtils.allPositionsOf(quantifiedPSM.getSequence(), c);
						numSites = +allPositionsOf.size();
						if (allPositionsOf.size() == 1) {
							quantifiedSitePositionInPeptide = new PositionInPeptide(allPositionsOf.get(0), c,
									quantifiedPSM.getSequence());
							ratio.setQuantifiedAA(c);
						}
					}
					// if no ambiguities
					if (numSites == 1) {
						ratio.addQuantifiedSitePositionInPeptide(quantifiedSitePositionInPeptide);
					}
				}
			} catch (final NumberFormatException e) {
				// skip this
			}
		}

		// create the peptide
		QuantifiedPeptideInterface quantifiedPeptide = null;
		final String peptideKey = KeyUtils.getInstance().getSequenceChargeKey(quantifiedPSM,
				isDistinguishModifiedSequences(), isChargeSensible());
		if (StaticQuantMaps.peptideMap.containsKey(peptideKey)) {
			quantifiedPeptide = StaticQuantMaps.peptideMap.getItem(peptideKey);
		} else {
			quantifiedPeptide = new QuantifiedPeptide(quantifiedPSM, isIgnoreTaxonomies(),
					isDistinguishModifiedSequences(), isChargeSensible());
		}
		StaticQuantMaps.peptideMap.addItem(quantifiedPeptide);
		quantifiedPSM.setQuantifiedPeptide(quantifiedPeptide, true);
		// add peptide to map
		if (!localPeptideMap.containsKey(peptideKey)) {
			localPeptideMap.put(peptideKey, quantifiedPeptide);
		}

		if (dbIndex != null) {
			final String cleanSeq = quantifiedPSM.getSequence();
			final Set<IndexedProtein> indexedProteins = dbIndex.getProteins(cleanSeq);
			if (indexedProteins.isEmpty()) {
				if (!super.ignoreNotFoundPeptidesInDB) {
					throw new PeptideNotFoundInDBIndexException("The peptide " + cleanSeq
							+ " is not found in Fasta DB.\nReview the default indexing parameters such as the number of allowed misscleavages.");
				}
				log.warn("The peptide " + cleanSeq + " is not found in Fasta DB with current digestion parameters.");
				// continue;
			}
			// create a new Quantified Protein for each
			// indexedProtein
			for (final IndexedProtein indexedProtein : indexedProteins) {
				final String proteinKey = KeyUtils.getInstance().getProteinKey(indexedProtein, isIgnoreACCFormat());
				QuantifiedProteinInterface quantifiedProtein = null;
				if (StaticQuantMaps.proteinMap.containsKey(proteinKey)) {
					quantifiedProtein = StaticQuantMaps.proteinMap.getItem(proteinKey);
				} else {
					quantifiedProtein = new QuantifiedProteinFromDBIndexEntry(indexedProtein, isIgnoreTaxonomies(),
							isIgnoreACCFormat());
				}
				StaticQuantMaps.proteinMap.addItem(quantifiedProtein);

				// add psm to the proteins
				quantifiedProtein.addPSM(quantifiedPSM, true);
				// add protein to the psm
				quantifiedPSM.addQuantifiedProtein(quantifiedProtein, true);
				// add peptide to the protein
				quantifiedProtein.addPeptide(quantifiedPeptide, true);
				// add to the map (if it was already there
				// is not a problem, it will be only once)
				addToMap(proteinKey, proteinToPeptidesMap, KeyUtils.getInstance().getSequenceChargeKey(quantifiedPSM,
						isDistinguishModifiedSequences(), isChargeSensible()));
				// add protein to protein map
				localProteinMap.put(proteinKey, quantifiedProtein);
				// add to protein-experiment map
				addToMap(experimentKey, experimentToProteinsMap, proteinKey);

			}
		}
		if (proteinACC != null) {
			final String proteinKey = proteinACC;
			QuantifiedProteinInterface quantifiedProtein = null;
			if (StaticQuantMaps.proteinMap.containsKey(proteinKey)) {
				quantifiedProtein = StaticQuantMaps.proteinMap.getItem(proteinKey);
			} else {
				quantifiedProtein = new QuantifiedProtein(proteinKey, isIgnoreTaxonomies());
			}
			StaticQuantMaps.proteinMap.addItem(quantifiedProtein);

			// add psm to the proteins
			quantifiedProtein.addPSM(quantifiedPSM, true);
			// add protein to the psm
			quantifiedPSM.addQuantifiedProtein(quantifiedProtein, true);
			// add peptide to the protein
			quantifiedProtein.addPeptide(quantifiedPeptide, true);
			// add to the map (if it was already there
			// is not a problem, it will be only once)
			addToMap(proteinKey, proteinToPeptidesMap, KeyUtils.getInstance().getSequenceChargeKey(quantifiedPSM,
					isDistinguishModifiedSequences(), isChargeSensible()));
			// add protein to protein map
			localProteinMap.put(proteinKey, quantifiedProtein);
			// add to protein-experiment map
			addToMap(experimentKey, experimentToProteinsMap, proteinKey);

		}
		if (proteinACC == null && dbIndex == null) {
			throw new IllegalArgumentException("Protein missing for peptide  " + quantifiedPeptide.getKey() + " ("
					+ psmId + "). Either provide a protein column or a Fasta file");
		}

	}

}
