package edu.scripps.yates.census.read;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.RatioDescriptor;
import edu.scripps.yates.census.read.model.StaticQuantMaps;
import edu.scripps.yates.census.read.model.interfaces.QuantParser;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.dbindex.DBIndexInterface;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.ipi.IPI2UniprotACCMap;
import edu.scripps.yates.utilities.model.enums.AccessionType;
import edu.scripps.yates.utilities.model.factories.AccessionEx;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import edu.scripps.yates.utilities.proteomicsmodel.Accession;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import edu.scripps.yates.utilities.util.Pair;

public abstract class AbstractQuantParser implements QuantParser {
	private static final Logger log = Logger.getLogger(AbstractQuantParser.class);
	protected DBIndexInterface dbIndex;

	protected final List<RemoteSSHFileReference> remoteFileRetrievers = new ArrayList<RemoteSSHFileReference>();

	protected Pattern decoyPattern;
	public static Set<String> peptidesMissingInDB = new HashSet<String>();
	protected boolean ignoreNotFoundPeptidesInDB;

	// MAPS
	// key=experimentkey, values=proteinKeys
	protected final HashMap<String, Set<String>> experimentToProteinsMap = new HashMap<String, Set<String>>();
	// key=proteinkey, values=peptidekeys
	protected final HashMap<String, Set<String>> proteinToPeptidesMap = new HashMap<String, Set<String>>();
	// key=peptideKey, values=spectrumKeys
	protected final HashMap<String, Set<String>> peptideToSpectraMap = new HashMap<String, Set<String>>();

	// the key is the protein key
	protected final Map<String, QuantifiedProteinInterface> localProteinMap = new HashMap<String, QuantifiedProteinInterface>();
	// the key is the spectrum key
	protected final Map<String, QuantifiedPSMInterface> localPsmMap = new HashMap<String, QuantifiedPSMInterface>();
	// the key is the peptide key (the peptide sequence, distinguising between
	// modified or not, depending on 'distinguishModifiedPeptides' variable
	protected final Map<String, QuantifiedPeptideInterface> localPeptideMap = new HashMap<String, QuantifiedPeptideInterface>();

	protected final Set<String> taxonomies = new HashSet<String>();
	protected boolean processed = false;

	protected final Map<RemoteSSHFileReference, Map<QuantCondition, QuantificationLabel>> labelsByConditionsByFile = new HashMap<RemoteSSHFileReference, Map<QuantCondition, QuantificationLabel>>();

	// protected final Map<RemoteSSHFileReference, QuantificationLabel>
	// numeratorLabelByFile = new HashMap<RemoteSSHFileReference,
	// QuantificationLabel>();
	// protected final Map<RemoteSSHFileReference, QuantificationLabel>
	// denominatorLabelByFile = new HashMap<RemoteSSHFileReference,
	// QuantificationLabel>();
	protected final Map<RemoteSSHFileReference, List<RatioDescriptor>> ratioDescriptorsByFile = new HashMap<RemoteSSHFileReference, List<RatioDescriptor>>();
	private UniprotProteinLocalRetriever uplr;
	private String uniprotVersion;
	protected boolean clearStaticMapsBeforeReading = true;
	private boolean retrieveFastaIsoforms;

	public AbstractQuantParser() {

	}

	public AbstractQuantParser(List<RemoteSSHFileReference> remoteSSHServers,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditions, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) {
		int index = 0;
		for (RemoteSSHFileReference remoteSSHServer : remoteSSHServers) {
			addFile(remoteSSHServer, labelsByConditions.get(index), labelNumerator, labelDenominator);
			index++;
		}
	}

	public AbstractQuantParser(Map<QuantCondition, QuantificationLabel> labelsByConditions,
			Collection<RemoteSSHFileReference> remoteSSHServers, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) {

		for (RemoteSSHFileReference remoteSSHServer : remoteSSHServers) {
			addFile(remoteSSHServer, labelsByConditions, labelNumerator, labelDenominator);
		}
	}

	public AbstractQuantParser(RemoteSSHFileReference remoteSSHServer,
			Map<QuantCondition, QuantificationLabel> labelsByConditions, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) throws FileNotFoundException {
		addFile(remoteSSHServer, labelsByConditions, labelNumerator, labelDenominator);
	}

	public AbstractQuantParser(File inputFile, Map<QuantCondition, QuantificationLabel> labelsByConditions,
			QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) throws FileNotFoundException {
		addFile(inputFile, labelsByConditions, labelNumerator, labelDenominator);
	}

	public AbstractQuantParser(File inputFile, Map<QuantCondition, QuantificationLabel> labelsByConditions,
			QuantificationLabel labelL, QuantificationLabel labelM, QuantificationLabel labelH)
			throws FileNotFoundException {
		addFile(inputFile, labelsByConditions, labelL, labelM, labelH);
	}

	public AbstractQuantParser(File[] inputFiles, Map<QuantCondition, QuantificationLabel> labelsByConditions,
			QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) throws FileNotFoundException {
		this(Arrays.asList(inputFiles), labelsByConditions, labelNumerator, labelDenominator);
	}

	public AbstractQuantParser(File[] inputFiles, Map<QuantCondition, QuantificationLabel>[] labelsByConditions,
			QuantificationLabel[] labelNumerator, QuantificationLabel[] labelDenominator) throws FileNotFoundException {
		for (int i = 0; i < inputFiles.length; i++) {
			addFile(inputFiles[i], labelsByConditions[i], labelNumerator[i], labelDenominator[i]);
		}
	}

	public AbstractQuantParser(Collection<File> inputFiles, Map<QuantCondition, QuantificationLabel> labelsByConditions,
			QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) throws FileNotFoundException {
		for (File xmlFile : inputFiles) {
			addFile(xmlFile, labelsByConditions, labelNumerator, labelDenominator);
		}
	}

	public AbstractQuantParser(RemoteSSHFileReference remoteServer, QuantificationLabel label1, QuantCondition cond1,
			QuantificationLabel label2, QuantCondition cond2) {
		Map<QuantCondition, QuantificationLabel> map = new HashMap<QuantCondition, QuantificationLabel>();
		map.put(cond1, label1);
		map.put(cond2, label2);
		addFile(remoteServer, map, label1, label2);
	}

	public AbstractQuantParser(RemoteSSHFileReference remoteServer, QuantificationLabel label1, QuantCondition cond1,
			QuantificationLabel label2, QuantCondition cond2, QuantificationLabel label3, QuantCondition cond3) {
		Map<QuantCondition, QuantificationLabel> map = new HashMap<QuantCondition, QuantificationLabel>();
		map.put(cond1, label1);
		map.put(cond2, label2);
		map.put(cond3, label3);
		addFile(remoteServer, map, label1, label2);
	}

	public AbstractQuantParser(File inputFile, QuantificationLabel label1, QuantCondition cond1,
			QuantificationLabel label2, QuantCondition cond2) throws FileNotFoundException {
		Map<QuantCondition, QuantificationLabel> map = new HashMap<QuantCondition, QuantificationLabel>();
		map.put(cond1, label1);
		map.put(cond2, label2);
		addFile(inputFile, map, label1, label2);
	}

	public AbstractQuantParser(File inputFile, QuantificationLabel label1, QuantCondition cond1,
			QuantificationLabel label2, QuantCondition cond2, QuantificationLabel label3, QuantCondition cond3)
			throws FileNotFoundException {
		Map<QuantCondition, QuantificationLabel> map = new HashMap<QuantCondition, QuantificationLabel>();
		map.put(cond1, label1);
		map.put(cond2, label2);
		map.put(cond3, label3);
		addFile(inputFile, map, label1, label2, label3);
	}

	@Override
	public void addFile(File xmlFile, Map<QuantCondition, QuantificationLabel> labelsByConditions,
			QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) throws FileNotFoundException {
		if (!xmlFile.exists()) {
			throw new FileNotFoundException(xmlFile.getAbsolutePath() + " is not found in the file system");
		}
		final RemoteSSHFileReference remoteFileReference = new RemoteSSHFileReference(xmlFile);
		addFile(remoteFileReference, labelsByConditions, labelNumerator, labelDenominator);
	}

	@Override
	public void addFile(File xmlFile, Map<QuantCondition, QuantificationLabel> labelsByConditions,
			QuantificationLabel labelL, QuantificationLabel labelM, QuantificationLabel labelH)
			throws FileNotFoundException {
		if (!xmlFile.exists()) {
			throw new FileNotFoundException(xmlFile.getAbsolutePath() + " is not found in the file system");
		}
		final RemoteSSHFileReference remoteFileReference = new RemoteSSHFileReference(xmlFile);
		addFile(remoteFileReference, labelsByConditions, labelL, labelM, labelH);
	}

	@Override
	public void addFile(RemoteSSHFileReference remoteFileReference,
			Map<QuantCondition, QuantificationLabel> labelsByConditions, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) {
		labelsByConditionsByFile.put(remoteFileReference, labelsByConditions);
		QuantCondition condition1 = null;
		QuantCondition condition2 = null;
		for (QuantCondition condition : labelsByConditions.keySet()) {
			final QuantificationLabel quantificationLabel = labelsByConditions.get(condition);
			if (quantificationLabel == labelNumerator) {
				condition1 = condition;
			} else if (quantificationLabel == labelDenominator) {
				condition2 = condition;
			}
		}

		RatioDescriptor ratioDescriptor = new RatioDescriptor(labelNumerator, labelDenominator, condition1, condition2);
		if (ratioDescriptorsByFile.containsKey(remoteFileReference)) {
			ratioDescriptorsByFile.get(remoteFileReference).add(ratioDescriptor);
		} else {
			List<RatioDescriptor> list = new ArrayList<RatioDescriptor>();
			list.add(ratioDescriptor);
			ratioDescriptorsByFile.put(remoteFileReference, list);
		}

		// numeratorLabelByFile.put(remoteFileReference, labelNumerator);
		// denominatorLabelByFile.put(remoteFileReference, labelDenominator);
		remoteFileRetrievers.add(remoteFileReference);
		// clearStaticInfo();
		checkParameters();
	}

	@Override
	public void addFile(RemoteSSHFileReference remoteFileReference,
			Map<QuantCondition, QuantificationLabel> labelsByConditions, QuantificationLabel labelL,
			QuantificationLabel labelM, QuantificationLabel labelH) {
		labelsByConditionsByFile.put(remoteFileReference, labelsByConditions);
		QuantCondition conditionL = null;
		QuantCondition conditionM = null;
		QuantCondition conditionH = null;
		for (QuantCondition condition : labelsByConditions.keySet()) {
			final QuantificationLabel quantificationLabel = labelsByConditions.get(condition);
			if (quantificationLabel == labelL) {
				conditionL = condition;
			} else if (quantificationLabel == labelM) {
				conditionM = condition;
			} else if (quantificationLabel == labelH) {
				conditionH = condition;
			}
		}
		// L/M
		RatioDescriptor ratioDescriptorLM = new RatioDescriptor(labelL, labelM, conditionL, conditionM);
		addRatioDescriptor(remoteFileReference, ratioDescriptorLM);
		// L/H
		RatioDescriptor ratioDescriptorLH = new RatioDescriptor(labelL, labelH, conditionL, conditionH);
		addRatioDescriptor(remoteFileReference, ratioDescriptorLH);
		// M/H
		RatioDescriptor ratioDescriptorMH = new RatioDescriptor(labelM, labelH, conditionM, conditionH);
		addRatioDescriptor(remoteFileReference, ratioDescriptorMH);

		// numeratorLabelByFile.put(remoteFileReference, labelNumerator);
		// denominatorLabelByFile.put(remoteFileReference, labelDenominator);
		remoteFileRetrievers.add(remoteFileReference);
		// clearStaticInfo();
		checkParameters();
	}

	private void addRatioDescriptor(RemoteSSHFileReference remoteFileReference, RatioDescriptor ratioDescriptor) {
		if (ratioDescriptorsByFile.containsKey(remoteFileReference)) {
			ratioDescriptorsByFile.get(remoteFileReference).add(ratioDescriptor);
		} else {
			List<RatioDescriptor> list = new ArrayList<RatioDescriptor>();
			list.add(ratioDescriptor);
			ratioDescriptorsByFile.put(remoteFileReference, list);
		}
	}

	/**
	 * It clears the static maps by protein keys and psm keys in
	 * {@link QuantifiedProteinInterface} and {@link QuantifiedPSMInterface}
	 * classes.<br>
	 * This method should be called at the beggining of an analysis in order to
	 * create just one {@link QuantifiedProteinInterface} for all the replicates
	 * and experiments for a given accession.
	 */
	private void clearStaticInfo() {
	}

	protected void checkParameters() {
		if (remoteFileRetrievers == null)
			throw new IllegalArgumentException("Input stream is null");

	}

	@Override
	public void setDecoyPattern(String patternString) throws PatternSyntaxException {
		if (patternString != null) {
			decoyPattern = Pattern.compile(patternString);
		} else {
			decoyPattern = null;
		}
	}

	/**
	 * @return the remoteFileRetrievers
	 */
	@Override
	public List<RemoteSSHFileReference> getRemoteFileRetrievers() {
		return remoteFileRetrievers;
	}

	/**
	 * @return the proteinMap @
	 */
	@Override
	public HashMap<String, Set<String>> getProteinToPeptidesMap() {
		if (!processed) {
			startProcess();
		}
		return proteinToPeptidesMap;
	}

	/**
	 * @return the peptideMap @
	 */
	@Override
	public HashMap<String, Set<String>> getPeptideToSpectraMap() {
		if (!processed) {
			startProcess();
		}
		return peptideToSpectraMap;
	}

	/**
	 * Gets the Protein by Protein key.
	 *
	 * @return the proteinMap @
	 */
	@Override
	public final Map<String, QuantifiedProteinInterface> getProteinMap() {
		if (!processed) {
			startProcess();
		}
		return localProteinMap;
	}

	@Override
	public final Map<String, QuantifiedPSMInterface> getPSMMap() {
		if (!processed) {
			startProcess();
		}
		return localPsmMap;
	}

	/**
	 * @param dbIndex
	 *            the dbIndex to set
	 */
	@Override
	public void setDbIndex(DBIndexInterface dbIndex) {
		this.dbIndex = dbIndex;
	}

	@Override
	public Set<String> getTaxonomies() {
		if (!processed) {
			startProcess();
		}
		return taxonomies;
	}

	/**
	 * @return the peptideMap @
	 */
	@Override
	public Map<String, QuantifiedPeptideInterface> getPeptideMap() {
		if (!processed) {
			startProcess();
		}
		return localPeptideMap;
	}

	private void startProcess() {
		if (clearStaticMapsBeforeReading) {
			// clear information in static maps
			StaticQuantMaps.clearInfo();
		}
		// first process
		process();
		// remove psms assigned to decoy proteins that were discarded
		removeDecoyPSMs();
		// second expand protein map
		mapIPI2Uniprot();
		// third merge proteins with secondary accessions
		mergeProteinsWithSecondaryAccessionsInParser();
	}

	protected abstract void process();

	public static void addToMap(String key, HashMap<String, Set<String>> map, String value) {
		if (map.containsKey(key)) {
			map.get(key).add(value);
		} else {
			Set<String> set = new HashSet<String>();
			set.add(value);
			map.put(key, set);
		}

	}

	/**
	 * @return the ignoreNotFoundPeptidesInDB
	 */
	@Override
	public boolean isIgnoreNotFoundPeptidesInDB() {
		return ignoreNotFoundPeptidesInDB;
	}

	/**
	 * @param ignoreNotFoundPeptidesInDB
	 *            the ignoreNotFoundPeptidesInDB to set
	 */
	@Override
	public void setIgnoreNotFoundPeptidesInDB(boolean ignoreNotFoundPeptidesInDB) {
		this.ignoreNotFoundPeptidesInDB = ignoreNotFoundPeptidesInDB;
	}

	/**
	 * Gets a set of Uniprot Accessions from the protein set in the parser. If
	 * the accessions are not Uniprot formatted, they are not retrieved here.
	 * <br>
	 * This function is used for getting annotations in uniprot for the proteins
	 * that are actually from uniprot.
	 *
	 * @return @
	 */
	@Override
	public Set<String> getUniprotAccSet() {
		Set<String> ret = new HashSet<String>();
		final Set<String> keySet = getProteinMap().keySet();
		for (String acc : keySet) {
			if (FastaParser.getUniProtACC(acc) != null) {
				ret.add(acc);
			}
		}
		return ret;
	}

	/**
	 * To be called after process().<br>
	 * If proteins have IPI accessions, look for the mapping from IPI 2 Uniprot.
	 * It adds new entries to the map, but it doesn't create any new
	 * {@link QuantifiedProteinInterface}
	 */
	private void mapIPI2Uniprot() {
		if (!localProteinMap.isEmpty()) {
			int originalNumberOfEntries = localProteinMap.size();
			Map<String, QuantifiedProteinInterface> newMap = new HashMap<String, QuantifiedProteinInterface>();
			for (String accession : localProteinMap.keySet()) {

				final Pair<String, String> acc = FastaParser.getACC(accession);
				if (acc.getSecondElement().equals("IPI")) {
					final QuantifiedProteinInterface quantProtein = localProteinMap.get(accession);
					Accession primaryAccession = new AccessionEx(accession, AccessionType.IPI);
					Pair<Accession, Set<Accession>> pair = IPI2UniprotACCMap.getInstance()
							.getPrimaryAndSecondaryAccessionsFromIPI(primaryAccession);
					if (pair.getFirstelement() != null) {
						primaryAccession = pair.getFirstelement();
						if (!newMap.containsKey(primaryAccession)) {
							newMap.put(primaryAccession.getAccession(), quantProtein);
						}
					}
					final Set<Accession> secondaryAccs = pair.getSecondElement();
					if (secondaryAccs != null) {
						for (Accession secondaryAcc : secondaryAccs) {
							if (!newMap.containsKey(secondaryAcc.getAccession())) {
								newMap.put(secondaryAcc.getAccession(), quantProtein);
							}
						}

					}
				}
			}
			for (String acc : newMap.keySet()) {
				if (!localProteinMap.containsKey(acc)) {
					localProteinMap.put(acc, newMap.get(acc));
				}
			}
			if (originalNumberOfEntries != localProteinMap.size()) {
				log.info("Protein Map expanded from " + originalNumberOfEntries + " to " + localProteinMap.size());
			}
		}
	}

	public void enableProteinMergingBySecondaryAccessions(UniprotProteinLocalRetriever uplr, String uniprotVersion) {
		this.uplr = uplr;
		this.uniprotVersion = uniprotVersion;
	}

	private void mergeProteinsWithSecondaryAccessionsInParser() {
		if (uplr == null) {
			return;
		}
		Set<String> accessions = new HashSet<String>();
		accessions.addAll(getProteinMap().keySet());
		String latestVersion = "latestVersion";
		if (uniprotVersion != null) {
			latestVersion = "version " + uniprotVersion;
		}
		// split into chunks of 500 accessions in order to show progress
		int chunckSize = 500000;
		List<Set<String>> listOfSets = new ArrayList<Set<String>>();
		Set<String> set = new HashSet<String>();
		for (String accession : accessions) {
			set.add(accession);
			if (set.size() == chunckSize) {
				listOfSets.add(set);
				set = new HashSet<String>();
			}
		}
		listOfSets.add(set);

		int numObsoletes = 0;
		log.info("Merging proteins that have secondary accessions according to Uniprot " + latestVersion + "...");
		ProgressCounter counter = new ProgressCounter(accessions.size(), ProgressPrintingType.PERCENTAGE_STEPS, 1);
		int initialSize = accessions.size();
		for (Set<String> accessionSet : listOfSets) {
			Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(uniprotVersion, accessionSet,
					retrieveFastaIsoforms);
			for (String accession : accessionSet) {
				counter.increment();
				String percentage = counter.printIfNecessary();
				if (!"".equals(percentage)) {
					log.info(percentage);
				}
				QuantifiedProteinInterface quantifiedProtein = StaticQuantMaps.proteinMap.getItem(accession);
				Entry entry = annotatedProteins.get(accession);
				if (entry != null && entry.getAccession() != null && !entry.getAccession().isEmpty()) {
					String primaryAccession = entry.getAccession().get(0);
					if (!accession.equals(primaryAccession) && !accession.contains(primaryAccession)) {
						log.info("Replacing Uniprot accession " + quantifiedProtein.getAccession() + " by "
								+ primaryAccession);
						quantifiedProtein.setAccession(primaryAccession);
						if (StaticQuantMaps.proteinMap.containsKey(primaryAccession)) {
							// there was already a protein with that
							// primaryAccession
							QuantifiedProteinInterface quantifiedProtein2 = StaticQuantMaps.proteinMap
									.getItem(primaryAccession);
							// merge quantifiedPRotein and quantifiedPRotein2
							mergeProteins(quantifiedProtein, quantifiedProtein2);

						} else {
							numObsoletes++;
						}
						// remove old/secondary accession
						getProteinMap().remove(accession);
						StaticQuantMaps.proteinMap.remove(accession);
						getProteinMap().put(primaryAccession, quantifiedProtein);

						StaticQuantMaps.proteinMap.addItem(quantifiedProtein);
					}
				} else {
					// // remove the protein because is obsolete
					// log.info(quantifiedProtein.getAccession());
					// parser.getProteinMap().remove(accession);
				}
			}
		}
		int finalSize = getProteinMap().size();
		if (initialSize != finalSize) {
			log.info(initialSize - finalSize
					+ " proteins with secondary accessions were merged with the corresponding protein with primary accession");
		}
		log.info("Obsolete accessions from " + numObsoletes + " proteins were changed to primary ones");
	}

	private static void mergeProteins(QuantifiedProteinInterface proteinReceiver,
			QuantifiedProteinInterface proteinDonor) {
		// PSMS
		for (QuantifiedPSMInterface psm : proteinDonor.getQuantifiedPSMs()) {
			proteinReceiver.addPSM(psm, true);
			psm.getQuantifiedProteins().remove(proteinDonor);
			psm.addQuantifiedProtein(proteinReceiver, true);
		}
		// Peptides
		for (QuantifiedPeptideInterface peptide : proteinDonor.getQuantifiedPeptides()) {
			proteinReceiver.addPeptide(peptide, true);

			peptide.getQuantifiedProteins().remove(proteinDonor);
		}
	}

	private void removeDecoyPSMs() {
		if (decoyPattern != null) {
			// in case of decoyPattern is enabled, we may have some PSMs
			// assigned to
			// those decoy proteins that have not been saved,
			// so we need to discard them
			// We iterate over the psms, and we will remove the ones with no
			// proteins
			Set<String> keysToDelete = new HashSet<String>();
			for (String key : localPsmMap.keySet()) {
				if (localPsmMap.get(key).getQuantifiedProteins().isEmpty()) {
					keysToDelete.add(key);
				}
			}
			log.info("Removing " + keysToDelete.size() + " PSMs not assigned to proteins");
			for (String key : keysToDelete) {
				final QuantifiedPSMInterface psm = localPsmMap.get(key);
				if (!psm.getQuantifiedProteins().isEmpty()) {
					throw new IllegalArgumentException("This should not happen");
				}

				// remove psmTableByPsmID
				localPsmMap.remove(key);
			}

			log.info(keysToDelete.size() + " PSMs discarded");
		}
	}

	public void setRetrieveFastaIsoforms(boolean retrieveFastaIsoforms) {
		this.retrieveFastaIsoforms = retrieveFastaIsoforms;
	}
}
