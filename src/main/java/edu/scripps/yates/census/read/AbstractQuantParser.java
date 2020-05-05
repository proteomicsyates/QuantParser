package edu.scripps.yates.census.read;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.RatioDescriptor;
import edu.scripps.yates.census.read.model.StaticQuantMaps;
import edu.scripps.yates.census.read.model.interfaces.QuantParser;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.ProteinSequences;
import edu.scripps.yates.census.read.util.QuantUtils;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.annotations.UniprotProteinLocalRetrieverInterface;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexInterface;
import edu.scripps.yates.utilities.ipi.IPI2UniprotACCMap;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import edu.scripps.yates.utilities.proteomicsmodel.Accession;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AccessionType;
import edu.scripps.yates.utilities.proteomicsmodel.factories.AccessionEx;
import edu.scripps.yates.utilities.proteomicsmodel.utils.KeyUtils;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import edu.scripps.yates.utilities.sequence.PTMInProtein;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public abstract class AbstractQuantParser implements QuantParser {
	private static final Logger log = Logger.getLogger(AbstractQuantParser.class);
	protected DBIndexInterface dbIndex;

	protected final List<RemoteSSHFileReference> remoteFileRetrievers = new ArrayList<RemoteSSHFileReference>();

	protected Pattern decoyPattern;
	public static Set<String> peptidesMissingInDB = new THashSet<String>();
	protected boolean ignoreNotFoundPeptidesInDB;
	private final Set<Character> quantifiedAAs = new THashSet<Character>();
	// MAPS
	// key=experimentkey, values=proteinKeys
	protected final Map<String, Set<String>> experimentToProteinsMap = new THashMap<String, Set<String>>();
	// key=proteinkey, values=peptidekeys
	protected final Map<String, Set<String>> proteinToPeptidesMap = new THashMap<String, Set<String>>();
	// key=peptideKey, values=spectrumKeys
	protected final Map<String, Set<String>> peptideToSpectraMap = new THashMap<String, Set<String>>();
	// key=siteInProtein(more than one can be), values=spectrumKeys
	protected final Map<String, Set<String>> ptmToSpectraMap = new THashMap<String, Set<String>>();
	// the key is the protein key
	protected final Map<String, QuantifiedProteinInterface> localProteinMap = new THashMap<String, QuantifiedProteinInterface>();
	// the key is the spectrum key
	protected final Map<String, QuantifiedPSMInterface> localPsmMap = new THashMap<String, QuantifiedPSMInterface>();
	// the key is the peptide key (the peptide sequence, distinguising between
	// modified or not, depending on 'distinguishModifiedPeptides' variable
	protected final Map<String, QuantifiedPeptideInterface> localPeptideMap = new THashMap<String, QuantifiedPeptideInterface>();

	protected final Set<String> taxonomies = new THashSet<String>();
	protected boolean processed = false;

	protected final Map<RemoteSSHFileReference, Map<QuantCondition, QuantificationLabel>> labelsByConditionsByFile = new THashMap<RemoteSSHFileReference, Map<QuantCondition, QuantificationLabel>>();
	protected final Map<RemoteSSHFileReference, Map<QuantificationLabel, QuantCondition>> conditionsByLabelsByFile = new THashMap<RemoteSSHFileReference, Map<QuantificationLabel, QuantCondition>>();

	// protected final Map<RemoteSSHFileReference, QuantificationLabel>
	// numeratorLabelByFile = new THashMap<RemoteSSHFileReference,
	// QuantificationLabel>();
	// protected final Map<RemoteSSHFileReference, QuantificationLabel>
	// denominatorLabelByFile = new THashMap<RemoteSSHFileReference,
	// QuantificationLabel>();
	protected final Map<RemoteSSHFileReference, List<RatioDescriptor>> ratioDescriptorsByFile = new THashMap<RemoteSSHFileReference, List<RatioDescriptor>>();
	protected UniprotProteinLocalRetrieverInterface uplr;
	protected String uniprotVersion;
	protected boolean clearStaticMapsBeforeReading = true;
	private boolean retrieveFastaIsoforms;
	private boolean ignoreTaxonomies;
	private boolean ignoreACCFormat;
	protected boolean getPTMInProteinMap = false;
	protected ProteinSequences proteinSequences;

	// PEPTIDE KEY GROUPING SETTINGS
	private boolean distinguishModifiedSequences = true;
	private boolean chargeSensible = true;
	private final Set<String> ratiosToCapture = new THashSet<String>();

	public boolean isDistinguishModifiedSequences() {
		return distinguishModifiedSequences;
	}

	public void setDistinguishModifiedSequences(boolean distinguishModifiedSequences) {
		this.distinguishModifiedSequences = distinguishModifiedSequences;
	}

	public boolean isChargeSensible() {
		return chargeSensible;
	}

	public void setChargeSensible(boolean chargeSensible) {
		this.chargeSensible = chargeSensible;
	}

	public AbstractQuantParser() {

	}

	public AbstractQuantParser(List<RemoteSSHFileReference> remoteSSHServers,
			List<Map<QuantCondition, QuantificationLabel>> labelsByConditions, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) {
		int index = 0;
		for (final RemoteSSHFileReference remoteSSHServer : remoteSSHServers) {
			addFile(remoteSSHServer, labelsByConditions.get(index), labelNumerator, labelDenominator);
			index++;
		}
	}

	public AbstractQuantParser(Map<QuantCondition, QuantificationLabel> labelsByConditions,
			Collection<RemoteSSHFileReference> remoteSSHServers, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator) {

		for (final RemoteSSHFileReference remoteSSHServer : remoteSSHServers) {
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

	public AbstractQuantParser(File inputFile, Map<QuantificationLabel, QuantCondition> conditionsByLabels)
			throws FileNotFoundException {

		if (!inputFile.exists()) {
			throw new FileNotFoundException(inputFile.getAbsolutePath() + " is not found in the file system");
		}
		final RemoteSSHFileReference remoteFileReference = new RemoteSSHFileReference(inputFile);
		conditionsByLabelsByFile.put(remoteFileReference, conditionsByLabels);

		remoteFileRetrievers.add(remoteFileReference);
		// clearStaticInfo();
		checkParameters();
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
		for (final File xmlFile : inputFiles) {
			addFile(xmlFile, labelsByConditions, labelNumerator, labelDenominator);
		}
	}

	public AbstractQuantParser(RemoteSSHFileReference remoteServer, QuantificationLabel label1, QuantCondition cond1,
			QuantificationLabel label2, QuantCondition cond2) {
		final Map<QuantCondition, QuantificationLabel> map = new THashMap<QuantCondition, QuantificationLabel>();
		map.put(cond1, label1);
		map.put(cond2, label2);
		addFile(remoteServer, map, label1, label2);
	}

	public AbstractQuantParser(RemoteSSHFileReference remoteServer, QuantificationLabel label1, QuantCondition cond1,
			QuantificationLabel label2, QuantCondition cond2, QuantificationLabel label3, QuantCondition cond3) {
		final Map<QuantCondition, QuantificationLabel> map = new THashMap<QuantCondition, QuantificationLabel>();
		map.put(cond1, label1);
		map.put(cond2, label2);
		map.put(cond3, label3);
		addFile(remoteServer, map, label1, label2);
	}

	public AbstractQuantParser(File inputFile, QuantificationLabel label1, QuantCondition cond1,
			QuantificationLabel label2, QuantCondition cond2) throws FileNotFoundException {
		final Map<QuantCondition, QuantificationLabel> map = new THashMap<QuantCondition, QuantificationLabel>();
		map.put(cond1, label1);
		map.put(cond2, label2);
		addFile(inputFile, map, label1, label2);
	}

	public AbstractQuantParser(File inputFile, QuantificationLabel label1, QuantCondition cond1,
			QuantificationLabel label2, QuantCondition cond2, QuantificationLabel label3, QuantCondition cond3)
			throws FileNotFoundException {
		final Map<QuantCondition, QuantificationLabel> map = new THashMap<QuantCondition, QuantificationLabel>();
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
		if (labelsByConditions != null) {
			for (final QuantCondition condition : labelsByConditions.keySet()) {
				final QuantificationLabel quantificationLabel = labelsByConditions.get(condition);
				if (quantificationLabel == labelNumerator) {
					condition1 = condition;
				} else if (quantificationLabel == labelDenominator) {
					condition2 = condition;
				}
			}
		}
		if (condition1 != null && condition2 != null) {
			final RatioDescriptor ratioDescriptor = new RatioDescriptor(labelNumerator, labelDenominator, condition1,
					condition2);
			addRatioDescriptor(remoteFileReference, ratioDescriptor);

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
		for (final QuantCondition condition : labelsByConditions.keySet()) {
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
		final RatioDescriptor ratioDescriptorLM = new RatioDescriptor(labelL, labelM, conditionL, conditionM);
		addRatioDescriptor(remoteFileReference, ratioDescriptorLM);
		// L/H
		final RatioDescriptor ratioDescriptorLH = new RatioDescriptor(labelL, labelH, conditionL, conditionH);
		addRatioDescriptor(remoteFileReference, ratioDescriptorLH);
		// M/H
		final RatioDescriptor ratioDescriptorMH = new RatioDescriptor(labelM, labelH, conditionM, conditionH);
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
			final List<RatioDescriptor> list = new ArrayList<RatioDescriptor>();
			list.add(ratioDescriptor);
			ratioDescriptorsByFile.put(remoteFileReference, list);
		}
	}

	public boolean isGetPTMInProteinMap() {
		return getPTMInProteinMap;
	}

	public void setGetPTMInProteinMap(boolean getPTMInProteinMap) {
		this.getPTMInProteinMap = getPTMInProteinMap;
	}

	/**
	 * It clears the static maps by protein keys and psm keys in
	 * {@link QuantifiedProteinInterface} and {@link QuantifiedPSMInterface}
	 * classes.<br>
	 * This method should be called at the beggining of an analysis in order to
	 * create just one {@link QuantifiedProteinInterface} for all the replicates and
	 * experiments for a given accession.
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
	 * @throws IOException
	 */
	@Override
	public Map<String, Set<String>> getProteinToPeptidesMap() throws IOException {
		if (!processed) {
			startProcess();
		}
		return proteinToPeptidesMap;
	}

	/**
	 * @return the peptideMap @
	 * @throws IOException
	 */
	@Override
	public Map<String, Set<String>> getPeptideToSpectraMap() throws IOException {
		if (!processed) {
			startProcess();
		}
		return peptideToSpectraMap;
	}

	/**
	 * @return the peptideMap @
	 * @throws IOException
	 */
	@Override
	public Map<String, Set<String>> getPTMToSpectraMap() throws IOException {
		if (!processed) {
			startProcess();
		}
		return ptmToSpectraMap;
	}

	/**
	 * Gets the Protein by Protein key.
	 *
	 * @return the proteinMap @
	 * @throws IOException
	 */
	@Override
	public final Map<String, QuantifiedProteinInterface> getProteinMap() throws IOException {
		if (!processed) {
			startProcess();
		}
		return localProteinMap;
	}

	@Override
	public final Map<String, QuantifiedPSMInterface> getPSMMap() throws IOException {
		if (!processed) {
			startProcess();
		}
		return localPsmMap;
	}

	/**
	 * @param dbIndexInterface the dbIndex to set
	 */
	@Override
	public void setDbIndex(DBIndexInterface dbIndexInterface) {
		dbIndex = dbIndexInterface;
	}

	@Override
	public Set<String> getTaxonomies() throws IOException {
		if (!processed) {
			startProcess();
		}
		return taxonomies;
	}

	/**
	 * @return the peptideMap @
	 * @throws IOException
	 */
	@Override
	public Map<String, QuantifiedPeptideInterface> getPeptideMap() throws IOException {
		if (!processed) {
			startProcess();
		}
		return localPeptideMap;
	}

	private void startProcess() throws IOException {
		if (clearStaticMapsBeforeReading) {
			// clear information in static maps
			StaticQuantMaps.clearInfo();
		}
		// first process
		process();
		// set processed to true
		processed = true;
		// remove psms assigned to decoy proteins that were discarded
		removeDecoyPSMs();
		// second expand protein map
		mapIPI2Uniprot();
		// third merge proteins with secondary accessions
		mergeProteinsWithSecondaryAccessionsInParser();
		// get ptmsInProteins, which implies to get uniprot annotations
		createPTMsInProteins();
	}

	private void createPTMsInProteins() throws IOException {

		if (getPTMInProteinMap) {
			final Set<String> accs = getProteinMap().keySet();
			uplr.getAnnotatedProteins(uniprotVersion, accs);

			for (final QuantifiedPSMInterface psm : getPSMMap().values()) {
				final String sequence = psm.getSequence();
				String ptmKey = sequence; // by default if no ptms
				if (!psm.getPTMsInPeptide().isEmpty()) {
					final List<PTMInProtein> ptmsInProtein = psm.getPTMsInProtein(uplr, proteinSequences);
					ptmKey = getPTMKeyFromPTMsInProtein(ptmsInProtein);

				}
				final String spectrumKey = KeyUtils.getInstance().getSpectrumKey(psm, isDistinguishModifiedSequences(),
						isChargeSensible());
				addToMap(ptmKey, ptmToSpectraMap, spectrumKey);
			}
		}

	}

	private String getPTMKeyFromPTMsInProtein(List<PTMInProtein> ptmsInProtein) {
		final StringBuilder sb = new StringBuilder();
		Collections.sort(ptmsInProtein, new Comparator<PTMInProtein>() {

			@Override
			public int compare(PTMInProtein o1, PTMInProtein o2) {
				final int ret = o1.getProteinACC().compareTo(o2.getProteinACC());
				if (ret != 0) {
					return ret;
				}
				return Integer.compare(o1.getPosition(), o2.getPosition());
			}
		});
		for (final PTMInProtein ptmInProtein : ptmsInProtein) {
			if (!"".equals(sb.toString())) {
				sb.append(QuantUtils.KEY_SEPARATOR);
			}
			sb.append(ptmInProtein.toString());
		}
		return sb.toString();
	}

	protected abstract void process() throws IOException;

	public static void addToMap(String key, Map<String, Set<String>> map, String value) {
		if (map.containsKey(key)) {
			map.get(key).add(value);
		} else {
			final Set<String> set = new THashSet<String>();
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
	 * @param ignoreNotFoundPeptidesInDB the ignoreNotFoundPeptidesInDB to set
	 */
	@Override
	public void setIgnoreNotFoundPeptidesInDB(boolean ignoreNotFoundPeptidesInDB) {
		this.ignoreNotFoundPeptidesInDB = ignoreNotFoundPeptidesInDB;
	}

	/**
	 * Gets a set of Uniprot Accessions from the protein set in the parser. If the
	 * accessions are not Uniprot formatted, they are not retrieved here. <br>
	 * This function is used for getting annotations in uniprot for the proteins
	 * that are actually from uniprot.
	 *
	 * @return @
	 * @throws IOException
	 */
	@Override
	public Set<String> getUniprotAccSet() throws IOException {
		final Set<String> ret = new THashSet<String>();
		final Set<String> keySet = getProteinMap().keySet();
		for (final String acc : keySet) {
			if (FastaParser.getUniProtACC(acc) != null) {
				ret.add(acc);
			}
		}
		return ret;
	}

	/**
	 * To be called after process().<br>
	 * If proteins have IPI accessions, look for the mapping from IPI 2 Uniprot. It
	 * adds new entries to the map, but it doesn't create any new
	 * {@link QuantifiedProteinInterface}
	 */
	private void mapIPI2Uniprot() {
		if (!localProteinMap.isEmpty()) {
			final int originalNumberOfEntries = localProteinMap.size();
			final Map<String, QuantifiedProteinInterface> newMap = new THashMap<String, QuantifiedProteinInterface>();
			for (final String accession : localProteinMap.keySet()) {

				final Accession acc = FastaParser.getACC(accession);
				if (acc.getAccessionType() == AccessionType.IPI) {
					final QuantifiedProteinInterface quantProtein = localProteinMap.get(accession);
					Accession primaryAccession = new AccessionEx(accession, AccessionType.IPI);
					final Pair<Accession, Set<Accession>> pair = IPI2UniprotACCMap.getInstance()
							.getPrimaryAndSecondaryAccessionsFromIPI(primaryAccession);
					if (pair.getFirstelement() != null) {
						primaryAccession = pair.getFirstelement();
						if (!newMap.containsKey(primaryAccession)) {
							newMap.put(primaryAccession.getAccession(), quantProtein);
						}
					}
					final Set<Accession> secondaryAccs = pair.getSecondElement();
					if (secondaryAccs != null) {
						for (final Accession secondaryAcc : secondaryAccs) {
							if (!newMap.containsKey(secondaryAcc.getAccession())) {
								newMap.put(secondaryAcc.getAccession(), quantProtein);
							}
						}

					}
				}
			}
			for (final String acc : newMap.keySet()) {
				if (!localProteinMap.containsKey(acc)) {
					localProteinMap.put(acc, newMap.get(acc));
				}
			}
			if (originalNumberOfEntries != localProteinMap.size()) {
				log.info("Protein Map expanded from " + originalNumberOfEntries + " to " + localProteinMap.size());
			}
		}
	}

	public void enableProteinMergingBySecondaryAccessions(UniprotProteinLocalRetrieverInterface uplr,
			String uniprotVersion) {
		this.uplr = uplr;
		this.uniprotVersion = uniprotVersion;
	}

	private void mergeProteinsWithSecondaryAccessionsInParser() throws IOException {
		if (uplr == null) {
			return;
		}
		final Set<String> accessions = new THashSet<String>();
		accessions.addAll(getProteinMap().keySet());
		String latestVersion = "latestVersion";
		if (uniprotVersion != null) {
			latestVersion = "version " + uniprotVersion;
		}
		// split into chunks of 500 accessions in order to show progress
		final int chunckSize = 500000;
		final List<Set<String>> listOfSets = new ArrayList<Set<String>>();
		Set<String> set = new THashSet<String>();
		for (final String accession : accessions) {
			set.add(accession);
			if (set.size() == chunckSize) {
				listOfSets.add(set);
				set = new THashSet<String>();
			}
		}
		listOfSets.add(set);

		int numObsoletes = 0;
		log.info("Merging proteins that have secondary accessions according to Uniprot " + latestVersion + "...");
		final ProgressCounter counter = new ProgressCounter(accessions.size(), ProgressPrintingType.PERCENTAGE_STEPS,
				0);
		final int initialSize = accessions.size();
		for (final Set<String> accessionSet : listOfSets) {
			final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(uniprotVersion, accessionSet,
					retrieveFastaIsoforms, false);
			for (final String accession : accessionSet) {
				counter.increment();
				final String percentage = counter.printIfNecessary();
				if (!"".equals(percentage)) {
					log.debug(percentage);
				}
				final QuantifiedProteinInterface quantifiedProtein = StaticQuantMaps.proteinMap.getItem(accession);
				final Entry entry = annotatedProteins.get(accession);
				if (entry != null && entry.getAccession() != null && !entry.getAccession().isEmpty()) {
					final String primaryAccession = entry.getAccession().get(0);
					if (!accession.equals(primaryAccession) && !accession.contains(primaryAccession)) {
						log.info("Replacing accession " + quantifiedProtein.getAccession() + " by primary accession "
								+ primaryAccession);
						quantifiedProtein.setPrimaryAccession(primaryAccession);
						if (StaticQuantMaps.proteinMap.containsKey(primaryAccession)) {
							// there was already a protein with that
							// primaryAccession
							final QuantifiedProteinInterface quantifiedProtein2 = StaticQuantMaps.proteinMap
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
		final int finalSize = getProteinMap().size();
		if (initialSize != finalSize) {
			log.info(initialSize - finalSize
					+ " proteins with secondary accessions were merged with the corresponding protein with primary accession");
		}
		log.info("Obsolete accessions from " + numObsoletes + " proteins were changed to primary ones");
	}

	private static void mergeProteins(QuantifiedProteinInterface proteinReceiver,
			QuantifiedProteinInterface proteinDonor) {
		// PSMS
		for (final QuantifiedPSMInterface psm : proteinDonor.getQuantifiedPSMs()) {
			proteinReceiver.addPSM(psm, true);
			psm.getQuantifiedProteins().remove(proteinDonor);
			psm.addQuantifiedProtein(proteinReceiver, true);
		}
		// Peptides
		for (final QuantifiedPeptideInterface peptide : proteinDonor.getQuantifiedPeptides()) {
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
			final Set<String> keysToDelete = new THashSet<String>();
			for (final String key : localPsmMap.keySet()) {
				if (localPsmMap.get(key).getQuantifiedProteins().isEmpty()) {
					keysToDelete.add(key);
				}
			}
			if (!keysToDelete.isEmpty()) {
				log.info("Removing " + keysToDelete.size() + " PSMs not assigned to proteins");
			}
			for (final String key : keysToDelete) {
				final QuantifiedPSMInterface psm = localPsmMap.get(key);
				if (!psm.getQuantifiedProteins().isEmpty()) {
					throw new IllegalArgumentException("This should not happen");
				}

				// remove psmTableByPsmID
				localPsmMap.remove(key);
			}
			if (!keysToDelete.isEmpty()) {
				log.info(keysToDelete.size() + " PSMs discarded");
			}
		}
	}

	public void setRetrieveFastaIsoforms(boolean retrieveFastaIsoforms) {
		this.retrieveFastaIsoforms = retrieveFastaIsoforms;
	}

	@Override
	public void addQuantifiedAA(char aa) {
		quantifiedAAs.add(aa);
	}

	@Override
	public Set<Character> getQuantifiedAAs() {
		return quantifiedAAs;
	}

	@Override
	public void setIgnoreTaxonomies(boolean ignoreTaxonomies) {
		this.ignoreTaxonomies = ignoreTaxonomies;
	}

	@Override
	public boolean isIgnoreTaxonomies() {
		return ignoreTaxonomies;
	}

	public boolean isIgnoreACCFormat() {
		return ignoreACCFormat;
	}

	public void setIgnoreACCFormat(boolean ignoreACCFormat) {
		this.ignoreACCFormat = ignoreACCFormat;
	}

	public UniprotProteinLocalRetrieverInterface getUplr() {
		return uplr;
	}

	public void setUplr(UniprotProteinLocalRetrieverInterface uplr) {
		this.uplr = uplr;
	}

	public ProteinSequences getProteinSequences() {
		return proteinSequences;
	}

	public void setProteinSequences(ProteinSequences proteinSequences) {
		this.proteinSequences = proteinSequences;
	}

	public String getUniprotVersion() {
		return uniprotVersion;
	}

	public void setUniprotVersion(String uniprotVersion) {
		this.uniprotVersion = uniprotVersion;
	}

	/**
	 * Adds ratio names to be captured if present in the input file
	 * 
	 * @param ratioName
	 */
	public void addRatioNameToCapture(String ratioName) {
		this.ratiosToCapture.add(ratioName.toLowerCase());
	}

	/**
	 * Adds ratio names to be captured if present in the input file
	 * 
	 * @param ratioNames Take them from the parser static properties, such as
	 *                   CensusOutParser.AREA_RATIO
	 */
	public void addRatioNamesToCapture(String... ratioNames) {
		for (final String ratioName : ratioNames) {
			this.ratiosToCapture.add(ratioName.toLowerCase());
		}

	}

	public Set<String> getRatiosNameToCapture() {
		return ratiosToCapture;
	}

	public boolean isCapturingRatioName(String ratioName) {
		return ratiosToCapture.contains(ratioName.toLowerCase());
	}
}
