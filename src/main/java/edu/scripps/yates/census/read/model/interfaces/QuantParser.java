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

	void addFile(File xmlFile, Map<QuantCondition, QuantificationLabel> labelsByConditions,
			QuantificationLabel labelNumerator, QuantificationLabel labelDenominator) throws FileNotFoundException;

	void addFile(File xmlFile, Map<QuantCondition, QuantificationLabel> labelsByConditions, QuantificationLabel labelL,
			QuantificationLabel labelM, QuantificationLabel labelH) throws FileNotFoundException;

	void addFile(RemoteSSHFileReference remoteFileReference,
			Map<QuantCondition, QuantificationLabel> labelsByConditions, QuantificationLabel labelNumerator,
			QuantificationLabel labelDenominator);

	void addFile(RemoteSSHFileReference remoteFileReference,
			Map<QuantCondition, QuantificationLabel> labelsByConditions, QuantificationLabel labelL,
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
}
