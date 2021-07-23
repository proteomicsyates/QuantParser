package edu.scripps.yates.census.read;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.QuantAmount;
import edu.scripps.yates.census.read.model.QuantifiedPSM;
import edu.scripps.yates.census.read.model.QuantifiedPeptide;
import edu.scripps.yates.census.read.model.QuantifiedProtein;
import edu.scripps.yates.census.read.model.QuantifiedProteinFromDBIndexEntry;
import edu.scripps.yates.census.read.model.StaticQuantMaps;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantKeyUtils;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.dbindex.util.PeptideNotFoundInDBIndexException;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.proteomicsmodel.PTM;
import edu.scripps.yates.utilities.proteomicsmodel.PTMPosition;
import edu.scripps.yates.utilities.proteomicsmodel.Score;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AmountType;
import edu.scripps.yates.utilities.proteomicsmodel.factories.GeneEx;
import edu.scripps.yates.utilities.proteomicsmodel.factories.PTMEx;
import edu.scripps.yates.utilities.proteomicsmodel.factories.PTMSiteEx;
import edu.scripps.yates.utilities.proteomicsmodel.factories.ScoreEx;
import edu.scripps.yates.utilities.proteomicsmodel.utils.KeyUtils;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import edu.scripps.yates.utilities.util.StringPosition;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;

public class MaxLFQuantParser extends AbstractQuantParser {
	private static final Logger log = Logger.getLogger(MaxLFQuantParser.class);

	private static final String H = "H";
	private static final String P = "P";
	private static final String S = "S";

	private static final String SEQUENCE = "Sequence";
	private static final String RAW_FILE = "Raw file";
	// synonym from previous versions
	public static final String FILE_name = "Filename";
	private static final String DESCRIPTION = "DESCRIPTION";

	private static final String MODIFIED_SEQUENCE = "Modified sequence";
	private static final String PIF = "PIF"; // Parent Ion Fraction (Interference)
	private static final String MSMS_IDS = "MS/MS IDs";
	private static final String SCORE = "Score";
	private static final String DELTA_SCORE = "Delta score";
	private static final String INTENSITY = "Intensity";
	private static final String CHARGE = "Charge";
	private static final String PRECURSOR_SCAN = "Precursor full scan number";
	private static final String ID = "id";
	private static final String LEADING_PROTEINS = "Leading proteins";

	private static final String LEADING_RAZOR_PROTEIN = "Leading razor protein";

	private static final String GENES = "Gene names";

	public MaxLFQuantParser() {
		super();
	}

	public MaxLFQuantParser(File maxQuantFolder, QuantificationLabel label, QuantCondition cond)
			throws FileNotFoundException {
		super(maxQuantFolder, label, cond, null, null);
	}

	/**
	 * Returns a list of headers that start by "Intensity_"
	 * 
	 * @param headersList
	 * @return
	 */
	private List<String> getIntensityHeaders(List<String> headersList) {
		final List<String> ret = new ArrayList<String>();
		for (final String header : headersList) {
			if (header.startsWith("Intensity")) {
				ret.add(header);
			}
		}
		return ret;
	}

	/**
	 *
	 * 
	 * @throws QuantParserException
	 */
	@Override
	protected void process() throws QuantParserException {
		processed = false;
		log.info("Processing quant file...");

		try {
			final int numDecoy = 0;
			boolean someValidFile = false;
			for (final RemoteSSHFileReference remoteFileRetriever : remoteFileRetrievers) {

				final Map<QuantCondition, Set<QuantificationLabel>> labelsByConditions = QuantUtils
						.getLabelsByConditions(conditionsByLabelsByFile.get(remoteFileRetriever));
				// it should only contain 1 condition and 1 label:
				final QuantCondition condition = labelsByConditions.keySet().iterator().next();
				final QuantificationLabel label = labelsByConditions.get(condition).iterator().next();

				final String experimentKey = FilenameUtils
						.getBaseName(remoteFileRetriever.getOutputFile().getAbsolutePath());
				log.info(experimentKey);

				log.info("Reading " + remoteFileRetriever.getRemoteFileName() + " from "
						+ remoteFileRetriever.getRemotePath());
				someValidFile = true;

				final File maxQuantFolder = remoteFileRetriever.getOutputFile();
				// get ms/ms file
				final File msmsFile = getFileFromMaxQuantFolder(maxQuantFolder, "msms.txt");
				final TIntObjectMap<QuantifiedPSMInterface> psmsById = readPSMsFromMSMSFile(msmsFile);
				// get mzTab.mzTab file to read the ptm mass shifts
				final File mzTabFile = getFileFromMaxQuantFolder(maxQuantFolder, "mzTab.mzTab");
				final Map<String, Double> ptmNameToMassShiftMap = readPTMNameToMassShiftFromMzTab(mzTabFile);

				// get evidence.txt file
				final File[] files = maxQuantFolder.listFiles(new FilenameFilter() {

					@Override
					public boolean accept(File dir, String name) {
						if (name.equalsIgnoreCase("evidence.txt")) {
							return true;
						}
						return false;
					}
				});
				if (files == null || files.length == 0) {
					throw new IllegalArgumentException("Evidence.txt file not found in MaxQuant folder '"
							+ maxQuantFolder.getAbsolutePath() + "'");
				}
				final File evidenceFile = files[0];

				BufferedReader br = null;
				try {
					br = new BufferedReader(
							new InputStreamReader(new BufferedInputStream(new FileInputStream(evidenceFile))));

					String line;
					final List<String> headersList = new ArrayList<String>();

					int numLine = 1;
					while ((line = br.readLine()) != null) {
						try {
							if (numLine == 1) {
								final String[] split = line.split("\t");
								for (int i = 0; i < split.length; i++) {
									headersList.add(split[i]);
								}
							} else {
								processEvidenceLine(line, headersList, condition, label, experimentKey,
										remoteFileRetriever, psmsById, ptmNameToMassShiftMap);
							}
						} finally {
							numLine++;
						}
					}

					br.close();
					numLine++;
				} catch (final PeptideNotFoundInDBIndexException e) {
					e.printStackTrace();
					if (!super.ignoreNotFoundPeptidesInDB) {
						throw e;
					}
				} catch (final Exception e) {
					e.printStackTrace();
					throw e;
				} finally {
					if (br != null) {
						br.close();
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
		} catch (final IOException e) {
			e.printStackTrace();
		} finally {
			// if (processed) {
			// // to create the peptides at the end
			// peptideMap.clear();
			// peptideMap.putAll(
			// QuantifiedPeptide.getQuantifiedPeptides(getPSMMap().values(),
			// distinguishModifiedPeptides));
			// }
		}
	}

	private TIntObjectMap<QuantifiedPSMInterface> readPSMsFromMSMSFile(File msmsFile) throws IOException {
		final TIntObjectMap<QuantifiedPSMInterface> ret = new TIntObjectHashMap<QuantifiedPSMInterface>();

		BufferedReader br = null;
		try {
			final FileReader fr = new FileReader(msmsFile);
			br = new BufferedReader(fr);
			String line = br.readLine();
			int numLine = 1;
			final List<String> headersList = new ArrayList<String>();
			while (line != null) {
				String sequence = null;
				String scanNumber = null;
				int chargeState = -1;
				String rawFileName = null;
				int id = -1;
				try {
					if (numLine == 1) {
						final String[] split = line.split("\t");
						for (int i = 0; i < split.length; i++) {
							headersList.add(split[i]);
						}
					} else {
						final Map<String, String> mapValues = getMapFromSLine(headersList, line);

						if (mapValues.containsKey(SEQUENCE)) {
							sequence = mapValues.get(SEQUENCE).trim();
						} else {
							throw new IllegalArgumentException("Column " + SEQUENCE + " not found in msms.txt file");
						}
						if (mapValues.containsKey(RAW_FILE)) {
							rawFileName = mapValues.get(RAW_FILE).trim();
						} else {
							throw new IllegalArgumentException("Column " + RAW_FILE + " not found in msms.txt file");
						}
						if (mapValues.containsKey(CHARGE)) {
							chargeState = Integer.valueOf(mapValues.get(CHARGE).trim());
						} else {
							throw new IllegalArgumentException("Column " + CHARGE + " not found in msms.txt file");
						}
						if (mapValues.containsKey(PRECURSOR_SCAN)) {
							scanNumber = mapValues.get(PRECURSOR_SCAN).trim();
						} else {
							throw new IllegalArgumentException(
									"Column " + PRECURSOR_SCAN + " not found in msms.txt file");
						}
						if (mapValues.containsKey(ID)) {
							id = Integer.valueOf(mapValues.get(ID).trim());
						} else {
							throw new IllegalArgumentException("Column " + ID + " not found in msms.txt file");
						}
						final QuantifiedPSMInterface psm = new QuantifiedPSM(sequence, peptideToSpectraMap, scanNumber,
								chargeState, rawFileName, false, isDistinguishModifiedSequences(), isChargeSensible());
						ret.put(id, psm);
					}
				} finally {
					line = br.readLine();
					numLine++;
				}
			}
		} finally {
			br.close();
			log.info(ret.size() + " PSMs read from file " + msmsFile.getAbsolutePath());
		}
		return ret;
	}

	/**
	 * Reads mzTab File and grabs the mass shifts associated to the name of the ptms
	 * 
	 * @param mzTabFile
	 * @return
	 * @throws IOException
	 */
	private Map<String, Double> readPTMNameToMassShiftFromMzTab(File mzTabFile) throws IOException {
		BufferedReader br = null;
		final Map<String, Double> ret = new THashMap<String, Double>();
		try {
			final FileReader fr = new FileReader(mzTabFile);
			br = new BufferedReader(fr);
			String line = br.readLine();

			while (line != null) {
				try {
					if (line.startsWith("COM")) {
						// line like this:
						// COM [, CHEMMOD:203.079372533, GlcNAc (N),]
						final String regexp = ".*CHEMMOD:(-?\\d+\\.\\d+)\\s*,\\s*(.+)\\s*,.*\\].*";
						final Pattern pattern = Pattern.compile(regexp);
						final Matcher m = pattern.matcher(line);
						if (m.find()) {
							final String massShiftString = m.group(1).trim();
							final double massShift = Double.valueOf(massShiftString);
							final String ptmName = m.group(2).trim();
							ret.put(ptmName, massShift);
						}
					}
				} finally {
					line = br.readLine();
				}
			}
		} finally {
			if (br != null) {
				br.close();
			}
		}
		return ret;
	}

	private File getFileFromMaxQuantFolder(File maxQuantFolder, String fileName) {
		final File[] files = maxQuantFolder.listFiles(new FilenameFilter() {

			@Override
			public boolean accept(File dir, String name) {
				if (name.equalsIgnoreCase(fileName)) {
					return true;
				}
				return false;
			}
		});
		if (files == null || files.length == 0) {
			throw new IllegalArgumentException(
					"File " + fileName + " not found in maxQuant results folder " + maxQuantFolder.getAbsolutePath());
		}
		return files[0];
	}

	private void processEvidenceLine(String line, List<String> sLineHeaderList,

			QuantCondition condition, QuantificationLabel label, String experimentKey,
			RemoteSSHFileReference remoteFileRetriever, TIntObjectMap<QuantifiedPSMInterface> psmsById,
			Map<String, Double> ptmNameToMassMap) throws IOException, WrongTMTLabels {

		// new psm
		try {

			final Map<String, String> mapValues = getMapFromSLine(sLineHeaderList, line);

			String rawSequence = mapValues.get(MODIFIED_SEQUENCE);
			// remove _ as prefix or suffix
			if (rawSequence.startsWith("_")) {
				rawSequence = rawSequence.substring(1);
			}
			if (rawSequence.endsWith("_")) {
				rawSequence = rawSequence.substring(0, rawSequence.length() - 1);
			}
			final String sequence = FastaParser.cleanSequence(rawSequence);
			final List<StringPosition> modifications = FastaParser.getInside(rawSequence);
			final List<PTM> ptms = new ArrayList<PTM>();
			if (modifications != null) {
				for (final StringPosition stringPosition : modifications) {
					final int positionInPeptide = stringPosition.position;
					final String modificationString = stringPosition.string;
					final Double massShift = ptmNameToMassMap.get(modificationString);
					final PTMEx ptm = new PTMEx(modificationString, massShift);
					if (positionInPeptide == 0) {
						ptm.addPtmSite(new PTMSiteEx(null, positionInPeptide, PTMPosition.NTERM));
						ptms.add(ptm);
					} else if (positionInPeptide == sequence.length() + 1) {
						ptm.addPtmSite(new PTMSiteEx(null, positionInPeptide, PTMPosition.CTERM));
						ptms.add(ptm);
					} else {
						final String aa = String.valueOf(sequence.charAt(positionInPeptide - 1));
						ptm.addPtmSite(new PTMSiteEx(aa, positionInPeptide, PTMPosition.NONE));
						ptms.add(ptm);
					}
				}

				// dont look into the QuantifiedPSM.map because each
				// line is always a new PSM
				final String inputFileName = FilenameUtils
						.getName(remoteFileRetriever.getOutputFile().getAbsolutePath());
				String rawFileName = null;
				if (mapValues.containsKey(RAW_FILE)) {
					rawFileName = mapValues.get(RAW_FILE);
				} else {
					rawFileName = inputFileName;
				}
				StaticQuantMaps.addRawFileName(rawFileName);

				// see the psm references
				final List<QuantifiedPSMInterface> psmsOfPeptide = new ArrayList<QuantifiedPSMInterface>();
				final TIntList psmIDs = new TIntArrayList();
				if (mapValues.containsKey(MSMS_IDS)) {
					final String ids = mapValues.get(MSMS_IDS);
					if (ids.contains(";")) {
						final String[] split = ids.split(";");
						for (final String string : split) {
							psmIDs.add(Integer.valueOf(string));
						}
					} else {
						psmIDs.add(Integer.valueOf(ids));
					}
				} else {
					throw new IllegalArgumentException("Column " + MSMS_IDS + " is not found in file");
				}
				for (final int psmId : psmIDs.toArray()) {
					final QuantifiedPSMInterface psm = psmsById.get(psmId);
					psmsOfPeptide.add(psm);
				}

				QuantifiedPeptideInterface peptide = new QuantifiedPeptide(psmsOfPeptide.get(0), isIgnoreTaxonomies(),
						isDistinguishModifiedSequences(), isChargeSensible());
				// add ptms
				for (final PTM ptm : ptms) {
					peptide.addPTM(ptm);
				}
				// add psms to peptide and ptms to psms
				if (psmsOfPeptide.size() > 1) {
					for (int i = 1; i < psmsOfPeptide.size(); i++) {
						final QuantifiedPSMInterface psm = psmsOfPeptide.get(i);
						peptide.addPSM(psm, true);
						for (final PTM ptm : ptms) {
							psm.addPTM(ptm);
						}
					}
				}
				// PIF
				if (mapValues.containsKey(PIF)) {
					try {
						final Float pif = Float.valueOf(mapValues.get(PIF));
						final Score score = new ScoreEx(String.valueOf(pif), PIF, PIF, PIF);
						peptide.addScore(score);
					} catch (final NumberFormatException e) {

					}
				}
				// Score
				if (mapValues.containsKey(SCORE)) {
					try {
						final Float value = Float.valueOf(mapValues.get(SCORE));
						final Score score = new ScoreEx(String.valueOf(value), "Andromeda Score", null, null);
						peptide.addScore(score);
					} catch (final NumberFormatException e) {

					}
				}
				// Delta score
				if (mapValues.containsKey(DELTA_SCORE)) {
					try {
						final Float value = Float.valueOf(mapValues.get(DELTA_SCORE));
						final Score score = new ScoreEx(String.valueOf(value), DELTA_SCORE, null, null);
						peptide.addScore(score);
					} catch (final NumberFormatException e) {

					}
				}

				// localization score
				final String localizationScore = null; // TODO

				final String peptideKey = KeyUtils.getInstance().getSequenceChargeKey(peptide,
						isDistinguishModifiedSequences(), isChargeSensible());
				// in case of TMT, the psm may have been created before
				if (StaticQuantMaps.peptideMap.containsKey(peptideKey)) {
					peptide = StaticQuantMaps.peptideMap.getItem(peptideKey);
				}
				StaticQuantMaps.peptideMap.addItem(peptide);

				// psms.add(quantifiedPSM);
				// add to map
				if (!localPeptideMap.containsKey(peptide.getKey())) {
					localPeptideMap.put(peptide.getKey(), peptide);
				}

				// Peptide intensity: amounts

				// SAM_INT
				// light peptide peak area from reconstructed
				// chromatogram
				if (mapValues.containsKey(INTENSITY)) {
					try {
						final double value = Double.valueOf(mapValues.get(INTENSITY));

						final QuantAmount amount = new QuantAmount(value, AmountType.XIC, condition, null);

						// add amount to peptide
						peptide.addAmount(amount);
					} catch (final NumberFormatException e) {
						// skip this
					}
				}
				final List<QuantifiedProteinInterface> quantifiedProteins = new ArrayList<QuantifiedProteinInterface>();

				if (dbIndex != null) {
					final String cleanSeq = peptide.getSequence();
					final Set<IndexedProtein> indexedProteins = dbIndex.getProteins(cleanSeq);
					if (indexedProteins.isEmpty()) {
						if (!ignoreNotFoundPeptidesInDB) {
							throw new PeptideNotFoundInDBIndexException("The peptide " + cleanSeq
									+ " is not found in Fasta DB.\nReview the default indexing parameters such as the number of allowed misscleavages.");
						}
						// log.warn("The peptide " + cleanSeq +
						// " is not found in Fasta DB.");
						// continue;
					}
					// create a new Quantified Protein for each
					// indexedProtein
					for (final IndexedProtein indexedProtein : indexedProteins) {
						final String proteinKey = QuantKeyUtils.getInstance().getProteinKey(indexedProtein,
								isIgnoreACCFormat());
						QuantifiedProteinInterface newQuantifiedProtein = null;
						if (StaticQuantMaps.proteinMap.containsKey(proteinKey)) {
							newQuantifiedProtein = StaticQuantMaps.proteinMap.getItem(proteinKey);

						} else {
							newQuantifiedProtein = new QuantifiedProteinFromDBIndexEntry(indexedProtein,
									isIgnoreTaxonomies(), isIgnoreACCFormat());

						}
						quantifiedProteins.add(newQuantifiedProtein);
						StaticQuantMaps.proteinMap.addItem(newQuantifiedProtein);
						// add protein to protein map
						if (newQuantifiedProtein.getTaxonomies() != null) {
							taxonomies.addAll(newQuantifiedProtein.getTaxonomies());
						}
						final QuantifiedProteinInterface tmp = localProteinMap.put(proteinKey, newQuantifiedProtein);
						// add to protein-experiment map
						addToMap(experimentKey, experimentToProteinsMap, proteinKey);
						// add psm to the protein
						for (final QuantifiedPSMInterface quantifiedPSM : psmsOfPeptide) {
							newQuantifiedProtein.addPSM(quantifiedPSM, true);
						}
						// add peptide to the protein
						newQuantifiedProtein.addPeptide(peptide, true);
						// add protein to the psms
						for (final QuantifiedPSMInterface quantifiedPSM : psmsOfPeptide) {
							quantifiedPSM.addQuantifiedProtein(newQuantifiedProtein, true);
						}
						// add to the map (if it was already
						// there is not a problem, it will be
						// only once)
						addToMap(proteinKey, proteinToPeptidesMap, KeyUtils.getInstance().getSequenceChargeKey(peptide,
								isDistinguishModifiedSequences(), isChargeSensible()));

					}
				} else {
					// proteins column
					if (mapValues.containsKey(LEADING_RAZOR_PROTEIN)) {
						final String leadingProteinsAccessions = mapValues.get(LEADING_RAZOR_PROTEIN);
						final List<String> accs = new ArrayList<String>();
						if (leadingProteinsAccessions.contains(";")) {
							final String[] split = leadingProteinsAccessions.split(";");
							for (final String leadingProteinAcc : split) {
								accs.add(leadingProteinAcc);
							}
						} else {
							accs.add(leadingProteinsAccessions);
						}
						for (final String acc : accs) {
							final QuantifiedProteinInterface protein = new QuantifiedProtein(acc);
							quantifiedProteins.add(protein);

						}
						if (mapValues.containsKey(GENES)) {
							final String geneNames = mapValues.get(GENES);
							if (!"".equals(geneNames)) {
								final List<String> genes = new ArrayList<String>();
								if (geneNames.contains(";")) {
									final String[] split = geneNames.split(";");
									for (final String gene : split) {
										genes.add(gene);
									}
								} else {
									genes.add(geneNames);
								}
								for (final String gene : genes) {
									for (final QuantifiedProteinInterface protein : quantifiedProteins) {
										protein.addGene(new GeneEx(gene));
									}
								}
							}
						}
					} else {
						throw new IllegalArgumentException("Column " + LEADING_RAZOR_PROTEIN + " is not found");
					}

				}
				// use the already created quantified
				// protein

				for (final QuantifiedProteinInterface quantifiedProtein : quantifiedProteins) {

					for (final QuantifiedPSMInterface quantifiedPSM : psmsOfPeptide) {
						// add psm to the proteins
						quantifiedProtein.addPSM(quantifiedPSM, true);
						// add protein to the psm
						quantifiedPSM.addQuantifiedProtein(quantifiedProtein, true);
					}
					// add peptide to the protein
					quantifiedProtein.addPeptide(peptide, true);
					// add to the map (if it was already there
					// is not a problem, it will be only once)
					final String proteinKey = quantifiedProtein.getKey();
					addToMap(proteinKey, proteinToPeptidesMap, KeyUtils.getInstance().getSequenceChargeKey(peptide,
							isDistinguishModifiedSequences(), isChargeSensible()));
					// add protein to protein map
					localProteinMap.put(proteinKey, quantifiedProtein);
					// add to protein-experiment map
					addToMap(experimentKey, experimentToProteinsMap, proteinKey);
				}
			}
		} catch (final IllegalArgumentException e) {
			e.printStackTrace();
			log.warn(e);
			log.info("Error reading line '" + line + "' from file. Skipping it...");

		} catch (final NullPointerException e) {
			e.printStackTrace();
			log.warn(e);
			log.warn("Error reading line '" + line + "' from file. Skipping it...");
		} catch (final DBIndexStoreException e) {

			e.printStackTrace();
			log.warn(e);
			log.warn("Error reading line '" + line + "' from file. Skipping it...");
		}

	}

	private Map<String, String> getMapFromSLine(List<String> headerList, String line) {
		final Map<String, String> map = new THashMap<String, String>();
		final String[] split = line.split("\t");
		for (int i = 0; i < split.length; i++) {
			final String value = split[i];
			final String header = headerList.get(i);
			map.put(header, value);
		}
		return map;
	}

	private String[] removeElements(String[] split, String elementToRemove, int upToThisIndex) {
		final List<String> list = new ArrayList<String>();
		int index = -1;
		for (final String splitElement : split) {
			index++;
			if (splitElement.equals(elementToRemove) && index < upToThisIndex) {
				continue;
			}
			list.add(splitElement);

		}
		return list.toArray(new String[0]);
	}

	@Override
	public boolean canRead() {
		try {

			for (final RemoteSSHFileReference remoteFileRetriever : this.remoteFileRetrievers) {
				final File file = remoteFileRetriever.getOutputFile();

				List<String> lines = null;

				// check whether it is an excel file
				if (FileUtils.isExcelFile(file)) {
					lines = FileUtils.readLinesFromXLSX(file, "\t", 0);
				} else {
					lines = FileUtils.readFirstLines(file, 1);
				}

				final String line = lines.get(0);

				final String[] split = line.split("\t");

				if (!line.startsWith(H) || !split[1].startsWith("Census version")) {
					return false;
				}

			}

		} catch (final Exception e) {
			return false;
		}
		return true;
	}
}
