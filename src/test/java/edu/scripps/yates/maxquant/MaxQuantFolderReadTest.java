package edu.scripps.yates.maxquant;

import static org.junit.Assert.fail;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Map;

import org.junit.Test;

import edu.scripps.yates.census.read.MaxLFQuantParser;
import edu.scripps.yates.census.read.QuantParserException;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;

public class MaxQuantFolderReadTest {

	@Test
	public void readTest() {
		final File folder = new File("C:\\Users\\salvador\\Desktop\\Saby\\HIV project\\MaxQuant data\\MQ_Saby_txt");
		MaxLFQuantParser parser;
		try {
			parser = new MaxLFQuantParser(folder, null, null);

			final Map<String, QuantifiedPeptideInterface> peptideMap = parser.getPeptideMap();
			for (final String peptideKey : peptideMap.keySet()) {
				final QuantifiedPeptideInterface peptide = peptideMap.get(peptideKey);
				System.out.println(
						peptideKey + "\t" + peptide.getFullSequence() + "\t" + peptide.getQuantifiedPSMs().size() + "\t"
								+ peptide.getAmounts().iterator().next().getValue());
			}
		} catch (final FileNotFoundException | QuantParserException e) {

			e.printStackTrace();
			fail();
		}
	}
}
