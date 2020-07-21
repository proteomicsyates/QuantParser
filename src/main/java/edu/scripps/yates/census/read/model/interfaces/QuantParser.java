package edu.scripps.yates.census.read.model.interfaces;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.PatternSyntaxException;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.QuantParserException;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexInterface;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;

public interface QuantParser {

	void addFile(File xmlFile, Map<QuantificationLabel, QuantCondition> conditionsByLabels,
			QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) throws FileNotFoundException;

	void addFile(File xmlFile, Map<QuantificationLabel, QuantCondition> conditionsByLabels, QuantificationLabel labelL,
			QuantificationLabel labelM, QuantificationLabel labelH) throws FileNotFoundException;

	void addFile(RemoteSSHFileReference remoteFileReference,
			Map<QuantificationLabel, QuantCondition> conditionsByLabels, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator);

	void addFile(RemoteSSHFileReference remoteFileReference,
			Map<QuantificationLabel, QuantCondition> conditionsByLabels, QuantificationLabel labelL,
			QuantificationLabel labelM, QuantificationLabel labelH);

	Map<String, QuantifiedPeptideInterface> getPeptideMap() throws QuantParserException;

	Set<String> getTaxonomies() throws QuantParserException;

	void setDbIndex(DBIndexInterface dbIndex);

	Map<String, QuantifiedPSMInterface> getPSMMap() throws QuantParserException;

	Map<String, QuantifiedProteinInterface> getProteinMap() throws QuantParserException;

	Map<String, Set<String>> getPeptideToSpectraMap() throws QuantParserException;

	Map<String, Set<String>> getProteinToPeptidesMap() throws QuantParserException;

	List<RemoteSSHFileReference> getRemoteFileRetrievers() throws QuantParserException;

	void setDecoyPattern(String patternString) throws PatternSyntaxException;

	Set<String> getUniprotAccSet() throws QuantParserException;

	void setIgnoreNotFoundPeptidesInDB(boolean ignoreNotFoundPeptidesInDB);

	boolean isIgnoreNotFoundPeptidesInDB();

	void addQuantifiedAA(char aa);

	Set<Character> getQuantifiedAAs();

	void setIgnoreTaxonomies(boolean ignoreTaxonomies);

	boolean isIgnoreTaxonomies();

	Map<String, Set<String>> getPTMToSpectraMap() throws QuantParserException;

	/**
	 * Returns the IonCount, that is, the number of PSMs that are present with the
	 * peptide (full sequence) + charge in this parser
	 * 
	 * @param psm
	 * @return
	 */
	int getReCalculatedIonCount(QuantifiedPSMInterface psm);

	/**
	 * Regardless of how the peptides are formed (using distinguish by charge state
	 * or PTMs), this will return the ion count of the ion that represent this
	 * peptide, that is, the peptide's full sequence + charge
	 * 
	 * @param peptide
	 * @return
	 */
	int getReCalculatedIonCount(QuantifiedPeptideInterface peptide);
}
