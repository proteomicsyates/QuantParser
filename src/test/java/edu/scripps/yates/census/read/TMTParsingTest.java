package edu.scripps.yates.census.read;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.springframework.core.io.ClassPathResource;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.QuantAmount;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AmountType;
import gnu.trove.map.hash.THashMap;
import junit.framework.Assert;

public class TMTParsingTest {
	@Test
	public void tmtParsingTest1() {
		try {
			final File tmtCensusOut = new ClassPathResource("bioTMT_b_census-out.txt").getFile();

			final Map<QuantificationLabel, QuantCondition> conditionsByLabels = new THashMap<QuantificationLabel, QuantCondition>();
			for (final QuantificationLabel label : QuantificationLabel.values()) {
				switch (label) {
				case TMT_10PLEX_126:
					conditionsByLabels.put(label, getCondition("AD1"));
					break;
				case TMT_10PLEX_127N:
					conditionsByLabels.put(label, getCondition("AD2"));
					break;
				case TMT_10PLEX_127C:
					conditionsByLabels.put(label, getCondition("AD3"));
					break;
				case TMT_10PLEX_128N:
					conditionsByLabels.put(label, getCondition("AD4"));
					break;
				case TMT_10PLEX_128C:
					conditionsByLabels.put(label, getCondition("AD5"));
					break;
				case TMT_10PLEX_129N:
					conditionsByLabels.put(label, getCondition("C1"));
					break;
				case TMT_10PLEX_129C:
					conditionsByLabels.put(label, getCondition("C2"));
					break;
				case TMT_10PLEX_130N:
					conditionsByLabels.put(label, getCondition("C3"));
					break;
				case TMT_10PLEX_130C:
					conditionsByLabels.put(label, getCondition("C4"));
					break;
				case TMT_10PLEX_131:
					conditionsByLabels.put(label, getCondition("C5"));
					break;

				default:
					break;
				}
			}

			int numAssertionsCorrect = 0;
			final CensusOutParser parser = new CensusOutParser(tmtCensusOut, conditionsByLabels, null, null);
			final Map<String, QuantifiedPSMInterface> psmMap = parser.getPSMMap();
			final QuantifiedPSMInterface psm = psmMap.get("bioTMT_S2_F07-17290-VDFSLAGALNAGFK(+339.162)ETR-3");
			System.out.println(psm.getKey());
			final Set<Amount> amounts = psm.getAmounts();
			for (final Amount amount : amounts) {
				final QuantAmount qamount = (QuantAmount) amount;
				if (amount.getAmountType() == AmountType.NORMALIZED_INTENSITY) {
					if (qamount.getLabel() == QuantificationLabel.TMT_10PLEX_126) {
						Assert.assertEquals(2175, amount.getValue(), 0);
						numAssertionsCorrect++;
					}
					if (qamount.getLabel() == QuantificationLabel.TMT_10PLEX_128C) {
						Assert.assertEquals(8811, amount.getValue(), 0);
						numAssertionsCorrect++;
					}
				}
			}
			final QuantifiedPSMInterface psm2 = psmMap
					.get("bioTMT_S2_F08-16111-GDVVPK(+36.075)DVNAAIATIK(+339.162)TK(+36.075)R-4");
			System.out.println(psm2.getKey());
			final Set<Amount> amounts2 = psm2.getAmounts();
			for (final Amount amount : amounts2) {
				final QuantAmount qamount = (QuantAmount) amount;
				if (amount.getAmountType() == AmountType.NORMALIZED_INTENSITY) {
					if (qamount.getLabel() == QuantificationLabel.TMT_10PLEX_128N) {
						Assert.assertEquals(109227, amount.getValue(), 0);
						numAssertionsCorrect++;
					}
					if (qamount.getLabel() == QuantificationLabel.TMT_10PLEX_131) {
						Assert.assertEquals(41781, amount.getValue(), 0);
						numAssertionsCorrect++;
					}
				}
			}
			final QuantifiedProteinInterface protein = parser.getProteinMap().get("Q71U36");
			System.out.println(protein.getKey());
			final Set<Amount> amounts3 = protein.getAmounts();
			for (final Amount amount : amounts3) {
				final QuantAmount qamount = (QuantAmount) amount;
				if (amount.getAmountType() == AmountType.NORMALIZED_INTENSITY) {
					if (qamount.getLabel() == QuantificationLabel.TMT_10PLEX_128N) {
						Assert.assertEquals(960166, amount.getValue(), 0);
						numAssertionsCorrect++;
					}
					if (qamount.getLabel() == QuantificationLabel.TMT_10PLEX_131) {
						Assert.assertEquals(410299, amount.getValue(), 0);
						numAssertionsCorrect++;
					}
				}
			}
			Assert.assertEquals(6, numAssertionsCorrect);
		} catch (final IOException | QuantParserException e) {
			e.printStackTrace();
			Assert.fail();
		}
	}

	private QuantCondition getCondition(String condName) {
		return new QuantCondition(condName);
	}
}
