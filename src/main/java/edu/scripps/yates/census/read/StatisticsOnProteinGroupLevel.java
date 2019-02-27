package edu.scripps.yates.census.read;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.read.model.IsobaricQuantifiedPSM;
import edu.scripps.yates.census.read.model.IsobaricQuantifiedProtein;
import edu.scripps.yates.census.read.model.interfaces.QuantRatio;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.CensusChroUtil;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.grouping.GroupablePeptide;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.PAnalyzer;
import edu.scripps.yates.utilities.grouping.PanalyzerStats;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.proteomicsmodel.utils.KeyUtils;
import edu.scripps.yates.utilities.taxonomy.UniprotOrganism;
import edu.scripps.yates.utilities.taxonomy.UniprotSpeciesCodeMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class StatisticsOnProteinGroupLevel {
	private static final Logger log = Logger.getLogger(StatisticsOnProteinGroupLevel.class);
	private final static String SEPARATOR = "\t\t";
	private final CensusChroParser census;
	private final Map<String, ProteinGroup> groupsInBothSpecies = new THashMap<String, ProteinGroup>();
	private final Map<String, Map<String, ProteinGroup>> groupsInOneSpecie = new THashMap<String, Map<String, ProteinGroup>>();
	private final Set<String> allTaxonomies = new THashSet<String>();
	private final Map<String, ProteinGroup> groupsContainingRatios = new THashMap<String, ProteinGroup>();
	private final Map<String, ProteinGroup> groupsNotContainingRatios = new THashMap<String, ProteinGroup>();
	private final DecimalFormat df = new DecimalFormat("#.#");
	private static final String YES = "1";
	private static final String NO = "0";
	private static final String PG = "PGR";
	private static final String PT = "PRT";
	private static final String PEP = "PEP";
	private static final String PSM = "PSM";
	private static final String LS = "LS";
	// private DBIndexInterface dbIndex;
	// private final File dbIndexParamFile;
	private List<ProteinGroup> proteinGroups;
	private PanalyzerStats groupingStats;
	private final boolean includePSMInformation;
	private ArrayList<String> taxonomyList;
	private final QuantificationLabel labelNumerator;
	private final QuantificationLabel labelDenominator;

	public StatisticsOnProteinGroupLevel(CensusChroParser census, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator, boolean includePSMInformation) throws IOException {
		this.census = census;
		this.includePSMInformation = includePSMInformation;
		this.labelDenominator = labelDenominator;
		this.labelNumerator = labelNumerator;
		processByGroups();
	}

	private void process() throws IOException {
		// grouping using panalyzer before doing nothing in order to build the
		// groups inside of the proteins
		log.info("Building groups...");
		getProteinGroups();
		log.info("Getting statistics from data...");

		log.info("Iterating over " + census.getPSMMap().size() + " quantified psms");
		int count = 0;
		for (final QuantifiedPSMInterface quantifiedPSM : census.getPSMMap().values()) {
			if (count++ % 500 == 0) {
				log.info(df.format(Double.valueOf(count) * 100 / census.getPSMMap().size()) + " % of PSMs...");
			}
			final Set<String> taxonomies = getTaxonomies(quantifiedPSM);
			allTaxonomies.addAll(taxonomies);
			if (!quantifiedPSM.getRatios().isEmpty()) {
				addGroupsFromPSM(quantifiedPSM, groupsContainingRatios);
			} else {
				addGroupsFromPSM(quantifiedPSM, groupsNotContainingRatios);
			}

			if (taxonomies.size() > 1) {
				// both
				addGroupsFromPSM(quantifiedPSM, groupsInBothSpecies);
			} else if (taxonomies.size() == 1) {
				final String taxonomy = taxonomies.iterator().next();
				// just in taxonomy
				if (groupsInOneSpecie.containsKey(taxonomy)) {
					addGroupsFromPSM(quantifiedPSM, groupsInOneSpecie.get(taxonomy));
				} else {
					final Map<String, ProteinGroup> map = new THashMap<String, ProteinGroup>();
					addGroupsFromPSM(quantifiedPSM, map);
					groupsInOneSpecie.put(taxonomy, map);
				}
			} else {
				log.warn("This cannot happen");
			}

		}
		final int numProteinsWithNoRatiosBefore = groupsNotContainingRatios.size();
		log.info("Consolidating " + numProteinsWithNoRatiosBefore + " Protein groups with no ratios...");
		count = 0;
		final Iterator<String> iterator = groupsNotContainingRatios.keySet().iterator();
		final int total = groupsNotContainingRatios.size();
		final Set<String> keysForDeletion = new THashSet<String>();
		while (iterator.hasNext()) {
			if (count++ % 100 == 0) {
				log.info(df.format(Double.valueOf(count) * 100 / total) + " % of Protein groups with no ratios...");
			}
			final String proteinKey = iterator.next();
			if (groupsContainingRatios.containsKey(proteinKey)) {
				keysForDeletion.add(proteinKey);
			}
		}
		for (final String proteinKey : keysForDeletion) {
			groupsNotContainingRatios.remove(proteinKey);
		}

		log.info("Removed " + (numProteinsWithNoRatiosBefore - groupsNotContainingRatios.size())
				+ " Protein groups from the list of proteins without ratios");

		log.info(census.getProteinMap().size() + " representative proteins quantified");
		log.info("Protein groups in both species: " + groupsInBothSpecies.size());
		for (final String taxonomy : groupsInOneSpecie.keySet()) {
			log.info("Protein groups just in " + taxonomy + ": " + groupsInOneSpecie.get(taxonomy).size());
		}
		log.info(groupsNotContainingRatios.size() + " Protein groups with no ratios");
		log.info(groupsContainingRatios.size() + " Protein groups containing ratios");
		log.info("Statistics done.");

	}

	private void processByGroups() throws IOException {
		// grouping using panalyzer before doing nothing in order to build the
		// groups inside of the proteins
		log.info("Building groups...");
		final List<ProteinGroup> groups = getProteinGroups();
		log.info("Getting statistics from data...");

		log.info("Iterating over " + groups.size() + " protein groups");
		int count = 0;
		for (final ProteinGroup proteinGroup : groups) {
			if (count++ % 500 == 0) {
				log.info(df.format(Double.valueOf(count) * 100 / groups.size()) + " % of Groups...");
			}
			final Set<String> taxonomies = getTaxonomiesFromGroup(proteinGroup);
			allTaxonomies.addAll(taxonomies);
			if (!CensusChroUtil.getRatiosFromProteinGroup(proteinGroup).isEmpty()) {
				addGroup(proteinGroup, groupsContainingRatios);
			} else {
				addGroup(proteinGroup, groupsNotContainingRatios);
			}

			if (taxonomies.size() > 1) {
				// both
				addGroup(proteinGroup, groupsInBothSpecies);
			} else if (taxonomies.size() == 1) {
				final String taxonomy = taxonomies.iterator().next();
				// just in taxonomy
				if (groupsInOneSpecie.containsKey(taxonomy)) {
					addGroup(proteinGroup, groupsInOneSpecie.get(taxonomy));
				} else {
					final Map<String, ProteinGroup> map = new THashMap<String, ProteinGroup>();
					addGroup(proteinGroup, map);
					groupsInOneSpecie.put(taxonomy, map);
				}
			} else {
				log.warn("This cannot happen");
			}

		}
		final int numProteinsWithNoRatiosBefore = groupsNotContainingRatios.size();
		log.info("Consolidating " + numProteinsWithNoRatiosBefore + " Protein groups with no ratios...");
		count = 0;
		final Iterator<String> iterator = groupsNotContainingRatios.keySet().iterator();
		final int total = groupsNotContainingRatios.size();
		final Set<String> keysForDeletion = new THashSet<String>();
		while (iterator.hasNext()) {
			if (count++ % 100 == 0) {
				log.info(df.format(Double.valueOf(count) * 100 / total) + " % of Protein groups with no ratios...");
			}
			final String proteinKey = iterator.next();
			if (groupsContainingRatios.containsKey(proteinKey)) {
				keysForDeletion.add(proteinKey);
			}
		}
		for (final String proteinKey : keysForDeletion) {
			groupsNotContainingRatios.remove(proteinKey);
		}

		log.info("Removed " + (numProteinsWithNoRatiosBefore - groupsNotContainingRatios.size())
				+ " Protein groups from the list of proteins without ratios");

		log.info("Protein groups in both species: " + groupsInBothSpecies.size());
		for (final String taxonomy : groupsInOneSpecie.keySet()) {
			log.info("Protein groups just in " + taxonomy + ": " + groupsInOneSpecie.get(taxonomy).size());
		}
		log.info(groupsNotContainingRatios.size() + " Protein groups with no ratios");
		log.info(groupsContainingRatios.size() + " Protein groups containing ratios");
		log.info("Statistics done.");

	}

	private void addGroup(ProteinGroup proteinGroup, Map<String, ProteinGroup> map) {
		map.put(KeyUtils.getInstance().getGroupKey(proteinGroup), proteinGroup);

	}

	private Set<String> getTaxonomiesFromGroup(ProteinGroup proteinGroup) {
		final Set<String> ret = new THashSet<String>();
		for (final GroupableProtein protein : proteinGroup) {
			if (protein instanceof IsobaricQuantifiedProtein) {
				final Set<String> taxonomies = ((IsobaricQuantifiedProtein) protein).getTaxonomies();
				if (taxonomies != null) {
					for (final String taxonomy : taxonomies) {
						final UniprotOrganism uniprotOrganism = UniprotSpeciesCodeMap.getInstance().get(taxonomy);
						if (uniprotOrganism != null) {
							ret.add(uniprotOrganism.getScientificName());
						} else {
							ret.add(taxonomy);
						}
					}
				}
			}
		}
		return ret;
	}

	private Set<String> getTaxonomies(QuantifiedPSMInterface quantifiedPSM) {

		final Set<String> taxonomies = quantifiedPSM.getTaxonomies();
		final Set<String> ret = new THashSet<String>();
		if (taxonomies != null) {
			for (final String taxonomy : taxonomies) {
				final UniprotOrganism uniprotOrganism = UniprotSpeciesCodeMap.getInstance().get(taxonomy);
				if (uniprotOrganism != null) {
					ret.add(uniprotOrganism.getCode());
				} else {
					ret.add(taxonomy);
				}
			}
		}
		return ret;

	}

	/**
	 * This function will add all the accessions of the proteins that the quantified
	 * psm belongs by using two different approaches depending on the constructor
	 * used for this class
	 *
	 * @param psm
	 * @param map
	 */
	private void addGroupsFromPSM(QuantifiedPSMInterface psm, Map<String, ProteinGroup> map) {
		final Set<QuantifiedProteinInterface> quantifiedProteins = psm.getQuantifiedProteins();

		for (final QuantifiedProteinInterface quantifiedProtein : quantifiedProteins) {
			final ProteinGroup group = quantifiedProtein.getProteinGroup();
			map.put(KeyUtils.getInstance().getGroupKey(group), group);
		}

	}

	/**
	 * @return the proteinsInBothSpecies
	 */
	public Map<String, ProteinGroup> getGroupsInBothSpecies() {
		return groupsInBothSpecies;
	}

	/**
	 * @return the proteinsInOneSpecie
	 */
	public Map<String, Map<String, ProteinGroup>> getGroupsInOneSpecie() {
		return groupsInOneSpecie;
	}

	/**
	 * Gets the taxonomies
	 *
	 * @return
	 */
	public Set<String> getTaxonomies() {
		return allTaxonomies;
	}

	/**
	 * Prints a table with the number of protein identified in each tag (light or
	 * heavy) for each species.
	 *
	 * @return
	 */
	public String printProteinGroupQuantificationTableBySpecies() {

		final StringBuilder sb = new StringBuilder();
		final StringBuilder allspecies = new StringBuilder();
		final StringBuilder allLabels = new StringBuilder();

		final Set<String> species = getTaxonomies();
		final List<String> sortedSpecies = new ArrayList<String>();
		sortedSpecies.addAll(species);
		Collections.sort(sortedSpecies);
		sb.append(SEPARATOR);

		for (final String specie : sortedSpecies) {
			sb.append(specie + SEPARATOR);
			if (!"".equals(allspecies.toString()))
				allspecies.append("+");
			allspecies.append(specie);
		}
		sb.append(allspecies.toString() + "\n");
		final QuantificationLabel[] labels = QuantificationLabel.values();
		for (int i = 0; i < labels.length; i++) {
			final QuantificationLabel label = labels[i];
			if (!"".equals(allLabels.toString()))
				allLabels.append("+");
			allLabels.append(label);
			sb.append(label + SEPARATOR);
			for (final String specie : sortedSpecies) {
				final Collection<ProteinGroup> proteinGroups = groupsInOneSpecie.get(specie).values();
				int proteinsInSpecieWithLabel = 0;
				if (proteinGroups != null) {
					for (final ProteinGroup proteinGroup : proteinGroups) {
						if (proteinGroup != null && CensusChroUtil.containsAnyIon(proteinGroup, label)
								&& !CensusChroUtil.containsAnyIonWithAnyOtherLabelThan(proteinGroup, label)) {
							proteinsInSpecieWithLabel++;
						}
					}
				}
				sb.append(proteinsInSpecieWithLabel + SEPARATOR);
			}
			// proteins in both species
			int proteinsInBothSpecieWithLabel = 0;
			for (final ProteinGroup proteinGroup : groupsInBothSpecies.values()) {
				if (CensusChroUtil.containsAnyIon(proteinGroup, label)
						&& !CensusChroUtil.containsAnyIonWithAnyOtherLabelThan(proteinGroup, label)) {
					proteinsInBothSpecieWithLabel++;
				}
			}
			sb.append(proteinsInBothSpecieWithLabel + "\n");
		}
		sb.append(allLabels.toString() + SEPARATOR);
		// both labels and one species
		for (final String specie : sortedSpecies) {
			if (groupsInOneSpecie.containsKey(specie)) {
				final Collection<ProteinGroup> proteinGroups = groupsInOneSpecie.get(specie).values();
				int numProteinsWithAllLabels = 0;
				for (final ProteinGroup proteinGroup : proteinGroups) {
					boolean found = true;
					for (final QuantificationLabel label : QuantificationLabel.values()) {
						if (!CensusChroUtil.containsAnyIon(proteinGroup, label)) {
							found = false;
							break;
						}
					}
					if (found)
						numProteinsWithAllLabels++;
				}

				sb.append(numProteinsWithAllLabels + SEPARATOR);
			}
		}
		// both labels and both species
		int numProteinsWithAllLabels = 0;
		for (final ProteinGroup proteinGroup : groupsInBothSpecies.values()) {
			boolean found = true;
			for (final QuantificationLabel label : QuantificationLabel.values()) {
				if (!CensusChroUtil.containsAnyIon(proteinGroup, label)) {
					found = false;
					break;
				}
			}
			if (found) {
				numProteinsWithAllLabels++;
			}
		}
		sb.append(numProteinsWithAllLabels + "\n");
		return sb.toString();
	}

	public String printDetailedProteinQuantificationTableBySpecies() throws IOException {
		final StringBuilder sb = new StringBuilder();
		final StringBuilder allspecies = new StringBuilder();
		final StringBuilder allLabels = new StringBuilder();

		final Set<String> species = getTaxonomies();
		taxonomyList = new ArrayList<String>();
		taxonomyList.addAll(species);
		Collections.sort(taxonomyList);

		for (final String specie : taxonomyList) {
			if (!"".equals(allspecies.toString()))
				allspecies.append("+");
			allspecies.append(specie);
		}

		final QuantificationLabel[] labels = QuantificationLabel.values();
		for (int i = 0; i < labels.length; i++) {
			final QuantificationLabel label = labels[i];
			if (!"".equals(allLabels.toString()))
				allLabels.append("+");
			allLabels.append(label);

			for (final String specie : taxonomyList) {
				sb.append(LS + SEPARATOR + label + SEPARATOR + specie + ": \n");
				final Collection<ProteinGroup> proteinGroups = groupsInOneSpecie.get(specie).values();
				int numProteinsInSpecieWithLabel = 0;
				if (proteinGroups != null) {
					for (final ProteinGroup proteinGroup : proteinGroups) {
						if (proteinGroup != null && CensusChroUtil.containsAnyIon(proteinGroup, label)
								&& !CensusChroUtil.containsAnyIonWithAnyOtherLabelThan(proteinGroup, label)) {
							numProteinsInSpecieWithLabel++;
							sb.append(PG + SEPARATOR + numProteinsInSpecieWithLabel + " " + proteinGroup);
							if (includePSMInformation) {
								printPSMInfo(sb, proteinGroup);
							}
							sb.append("\n");
						}
					}
				}
				sb.append("Count=" + numProteinsInSpecieWithLabel + "\n\n");
				sb.append("*********\n");
			}
			// proteins in both species
			int proteinsInBothSpecieWithLabel = 0;
			sb.append(LS + SEPARATOR + label + SEPARATOR + allspecies + ": \n");
			for (final ProteinGroup proteinGroup : groupsInBothSpecies.values()) {
				// boolean found = false;
				// for (GroupableProtein groupableProtein : proteinGroup) {
				// if (groupableProtein.getPrimaryAccession().getDescription()
				// .contains("GJ21323")) {
				// found = true;
				// }
				// }
				// if (found) {
				// for (GroupableProtein groupableProtein : proteinGroup) {
				// System.out.println(groupableProtein
				// .getPrimaryAccession().getAccession()
				// + " "
				// + groupableProtein.getPrimaryAccession()
				// .getDescription()
				// + " "
				// + groupableProtein.getGroupablePSMs().size());
				// for (GroupablePSM psm : groupableProtein
				// .getGroupablePSMs()) {
				// System.out.println(psm.getPSMIdentifier() + " "
				// + psm.getSequence());
				// }
				// }
				// }
				if (CensusChroUtil.containsAnyIon(proteinGroup, label)
						&& !CensusChroUtil.containsAnyIonWithAnyOtherLabelThan(proteinGroup, label)) {
					proteinsInBothSpecieWithLabel++;
					sb.append(proteinsInBothSpecieWithLabel + " " + proteinGroup);
					if (includePSMInformation) {
						printPSMInfo(sb, proteinGroup);
					}
					sb.append("\n");
				}
			}
			sb.append("Count=" + proteinsInBothSpecieWithLabel + "\n\n");
			sb.append("*********\n");
		}

		// both labels and one species
		for (final String specie : taxonomyList) {
			sb.append(LS + SEPARATOR + allLabels + SEPARATOR + specie + ": \n");
			if (groupsInOneSpecie.containsKey(specie)) {
				final Collection<ProteinGroup> proteinGroups = groupsInOneSpecie.get(specie).values();
				int numProteinsWithAllLabels = 0;
				for (final ProteinGroup proteinGroup : proteinGroups) {
					boolean found = true;
					for (final QuantificationLabel label : QuantificationLabel.values()) {
						if (!CensusChroUtil.containsAnyIon(proteinGroup, label)) {
							found = false;
							break;
						}
					}
					if (found) {
						numProteinsWithAllLabels++;
						sb.append(PG + SEPARATOR + numProteinsWithAllLabels + " " + proteinGroup + "\n");
						if (includePSMInformation)
							printPSMInfo(sb, proteinGroup);
					}
				}

				sb.append("Count=" + numProteinsWithAllLabels + "\n" + "\n");
				sb.append("*********\n");
			}
		}
		// both labels and both species
		int numProteinsWithAllLabels = 0;
		sb.append(LS + SEPARATOR + allLabels + SEPARATOR + allspecies + ": \n");
		for (final ProteinGroup proteinGroup : groupsInBothSpecies.values()) {
			boolean found = true;
			for (final QuantificationLabel label : QuantificationLabel.values()) {
				if (!CensusChroUtil.containsAnyIon(proteinGroup, label)) {
					found = false;
					break;
				}
			}
			if (found) {
				numProteinsWithAllLabels++;
				sb.append(PG + SEPARATOR + numProteinsWithAllLabels + " " + proteinGroup + "\n");
				if (includePSMInformation)
					printPSMInfo(sb, proteinGroup);
			}
		}
		sb.append("Count=" + numProteinsWithAllLabels + "\n" + "\n");
		sb.append("*********\n");
		return sb.toString();
	}

	private void printPSMInfo(StringBuilder psmTableString, ProteinGroup proteinGroup) throws IOException {
		final List<GroupablePeptide> psMs = proteinGroup.getPSMs();
		sortedPsmsBySequence(psMs);
		final boolean printPeptideTotal = true;
		int count = 0;
		String currentPeptideSequence = null;
		int peptideNumRatios = 0;
		int psmNumber = 0;
		final Map<QuantificationLabel, Integer> peptideIonsLabeled = new THashMap<QuantificationLabel, Integer>();
		final Map<String, Boolean> peptideIdentifiedBySpecie = new THashMap<String, Boolean>();
		final List<Double> peptideRatios = new ArrayList<Double>();
		for (final GroupablePeptide groupablePSM : psMs) {
			final QuantifiedPSMInterface quantifiedPSM = (QuantifiedPSMInterface) groupablePSM;
			final String seq = quantifiedPSM.getFullSequence();

			if (currentPeptideSequence != null && !seq.equals(currentPeptideSequence) && printPeptideTotal) {
				printPeptideTotal(psmTableString, currentPeptideSequence, peptideIonsLabeled, peptideRatios, psmNumber,
						peptideIdentifiedBySpecie, peptideNumRatios);
				psmNumber = 0;
			}
			currentPeptideSequence = seq;
			final String psmID = KeyUtils.getInstance().getSpectrumKey(quantifiedPSM, true);
			if (count++ % 500 == 0) {
				log.info(df.format(Double.valueOf(count) * 100 / census.getPSMMap().size()) + " % of PSMs...");
			}

			// PEP SEQ
			psmTableString.append(PSM + SEPARATOR + seq + SEPARATOR);
			// PSM ID
			psmTableString.append(psmID + SEPARATOR);
			// taxonomies
			for (final String taxonomy : taxonomyList) {
				final Set<String> taxonomies = quantifiedPSM.getTaxonomies();
				if (taxonomies != null && taxonomies.contains(taxonomy)) {
					psmTableString.append(YES + SEPARATOR);
					peptideIdentifiedBySpecie.put(taxonomy, true);
				} else {
					psmTableString.append(NO + SEPARATOR);
				}
			}
			// labeled
			for (final QuantificationLabel label : QuantificationLabel.values()) {
				if (quantifiedPSM instanceof IsobaricQuantifiedPSM) {
					final IsobaricQuantifiedPSM isoPsm = (IsobaricQuantifiedPSM) quantifiedPSM;
					final int numLabeledIons = isoPsm.getSingletonIonsByLabel(label).size();
					psmTableString.append(numLabeledIons + SEPARATOR);
					if (peptideIonsLabeled.containsKey(label)) {
						peptideIonsLabeled.put(label, peptideIonsLabeled.get(label) + numLabeledIons);
					} else {
						peptideIonsLabeled.put(label, numLabeledIons);
					}
				}
			}

			Double avgRatio = null;
			if (!quantifiedPSM.getRatios().isEmpty()) {
				double sum = 0.0;
				int valid = 0;
				for (final QuantRatio ratio : quantifiedPSM.getQuantRatios()) {
					if (!Maths.isMaxOrMinValue(ratio.getLog2Ratio(labelNumerator, labelDenominator))
							&& !Double.isNaN(ratio.getLog2Ratio(labelNumerator, labelDenominator))) {
						sum += ratio.getLog2Ratio(labelNumerator, labelDenominator);
						valid++;
					}
				}
				avgRatio = sum / valid;
				peptideRatios.add(avgRatio);
			}
			if (avgRatio == null) {
				psmTableString.append("-" + SEPARATOR);
			} else {
				psmTableString.append(avgRatio + SEPARATOR);
			}
			psmTableString.append(quantifiedPSM.getRatios().size() + SEPARATOR);
			peptideNumRatios += quantifiedPSM.getRatios().size();
			psmNumber++;
			psmTableString.append("\n");

		}
		printPeptideTotal(psmTableString, currentPeptideSequence, peptideIonsLabeled, peptideRatios, psmNumber,
				peptideIdentifiedBySpecie, peptideNumRatios);
	}

	private void printPeptideTotal(StringBuilder psmTableString, String currentPeptideSequence,
			Map<QuantificationLabel, Integer> peptideIonsLabeled, List<Double> peptideRatios, int psmNumber,
			Map<String, Boolean> peptideIdentifiedBySpecie, int peptideNumRatios) {
		// print PEPTIDE data
		psmTableString.append(PEP + SEPARATOR + currentPeptideSequence + SEPARATOR);
		psmTableString.append(psmNumber + SEPARATOR);
		for (final String taxonomy : taxonomyList) {
			if (peptideIdentifiedBySpecie.containsKey(taxonomy) && peptideIdentifiedBySpecie.get(taxonomy)) {
				psmTableString.append(YES + SEPARATOR);
			} else {
				psmTableString.append(NO + SEPARATOR);
			}
		}
		for (final QuantificationLabel label : QuantificationLabel.values()) {
			if (peptideIonsLabeled.containsKey(label)) {
				psmTableString.append(peptideIonsLabeled.get(label) + SEPARATOR);
			} else {
				psmTableString.append("0" + SEPARATOR);
			}
		}
		if (!peptideRatios.isEmpty()) {
			double sumOfAvgRatios = 0.0;
			for (final Double avratio : peptideRatios) {
				sumOfAvgRatios += avratio;
			}
			final double avgOfAvgRatios = sumOfAvgRatios / peptideRatios.size();
			psmTableString.append(avgOfAvgRatios + SEPARATOR);
		} else {
			psmTableString.append("-" + SEPARATOR);
		}
		psmTableString.append(peptideNumRatios + SEPARATOR);
		psmTableString.append("\n\n");

		// reset all peptide data
		peptideIdentifiedBySpecie.clear();
		peptideIonsLabeled.clear();
		peptideRatios.clear();
		peptideNumRatios = 0;
		psmNumber = 0;

	}

	private void sortedPsmsBySequence(List<GroupablePeptide> psms) {

		Collections.sort(psms, new Comparator<GroupablePeptide>() {

			@Override
			public int compare(GroupablePeptide o1, GroupablePeptide o2) {
				final String seq1 = o1.getSequence();
				final String seq2 = o2.getSequence();
				return seq1.compareTo(seq2);
			}

		});

	}

	public List<ProteinGroup> getProteinGroups() throws IOException {
		if (proteinGroups == null) {
			final PAnalyzer pa = new PAnalyzer(false);
			final List<GroupableProtein> set = new ArrayList<GroupableProtein>();
			set.addAll(census.getProteinMap().values());
			proteinGroups = pa.run(set);
			groupingStats = pa.getStats();
		}
		return proteinGroups;
	}

	/**
	 * @return the groupingStats
	 * @throws IOException
	 */
	public PanalyzerStats getGroupingStats() throws IOException {
		if (groupingStats == null) {
			getProteinGroups();
		}
		return groupingStats;
	}

	/**
	 * @return the allTaxonomies
	 */
	public Set<String> getAllTaxonomies() {
		return allTaxonomies;
	}

	/**
	 * @return the proteinsContainingRatios
	 */
	public Map<String, ProteinGroup> getGroupsContainingRatios() {
		return groupsContainingRatios;
	}

	/**
	 * @return the proteinsNotContainingRatios
	 */
	public Map<String, ProteinGroup> getGroupsNotContainingRatios() {
		return groupsNotContainingRatios;
	}

	public String printSummaryReport() throws IOException {
		final StringBuilder sb = new StringBuilder();
		sb.append(printProteinGroupQuantificationTableBySpecies() + "\n");
		sb.append(getGroupsContainingRatios().size()
				+ " protein groups containing at least one ratio in one peptide\n\n");
		sb.append(
				getGroupsNotContainingRatios().size() + " protein groups not containing ratios in their peptides\n\n");
		sb.append("Protein grouping statistics:\n");
		sb.append(getGroupingStats().toString());
		return sb.toString();
	}

	public String printFullReport() throws IOException {
		final StringBuilder sb = new StringBuilder();
		sb.append(printSummaryReport() + "\n");

		sb.append(printDetailedProteinQuantificationTableBySpecies() + "\n");

		return sb.toString();
	}
}
