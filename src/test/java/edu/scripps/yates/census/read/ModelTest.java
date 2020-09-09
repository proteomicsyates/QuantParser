package edu.scripps.yates.census.read;

import static org.junit.Assert.fail;

import java.util.List;

import org.junit.Test;

import edu.scripps.yates.census.read.model.QuantifiedPSM;
import edu.scripps.yates.census.read.model.QuantifiedProtein;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.util.StringPosition;

public class ModelTest {

	@Test
	public void usingFakePeptides() {
		try {
			FastaParser.setTolerantToFakePeptides(true);
			final String accs = "P05387,P08670,P08670";
			final String cleanSequence = FastaParser.cleanSequence(accs);
			System.out.println(cleanSequence);
			final List<StringPosition> inside = FastaParser.getInside(accs);
			for (final StringPosition stringPosition : inside) {
				System.out.println(stringPosition);
			}
			final QuantifiedPSMInterface psm = new QuantifiedPSM(accs, null, null, 0, "run", false, true, true);
			System.out.println(psm.getSequence());
			System.out.println(psm.getFullSequence());

			final QuantifiedProteinInterface protein = new QuantifiedProtein("1");
			protein.addPSM(psm, true);

			System.out.println(protein.getPeptides().size());
		} catch (final Exception e) {
			e.printStackTrace();
			fail();
		}
	}

}
