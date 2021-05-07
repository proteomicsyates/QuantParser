package edu.scripps.yates.census.read.util;

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

public enum QuantificationLabel {

	LIGHT(true, false, false), //
	HEAVY(false, false, true), //
	MEDIUM(false, true, false), //
	//
	TMT_4PLEX_127(false, false, false), //
	TMT_4PLEX_128(false, false, false), //
	TMT_4PLEX_129(false, false, false), //
	TMT_4PLEX_130(false, false, false), //

	//
	TMT_6PLEX_126(true, false, false), //
	TMT_6PLEX_127(false, false, false), //
	TMT_6PLEX_128(false, false, false), //
	TMT_6PLEX_129(false, false, false), //
	TMT_6PLEX_130(false, false, false), //
	TMT_6PLEX_131(false, false, true), //
	//
	TMT_9PLEX_126(true, false, false), //
	TMT_9PLEX_127N(false, false, false), //
	TMT_9PLEX_127C(false, false, false), //
	TMT_9PLEX_128N(false, false, false), //
	TMT_9PLEX_128C(false, false, false), //
	TMT_9PLEX_129N(false, false, false), //
	TMT_9PLEX_129C(false, false, false), //
	TMT_9PLEX_130N(false, false, false), //
	TMT_9PLEX_130C(false, false, true), //
	//
	TMT_10PLEX_126(true, false, false), //
	TMT_10PLEX_127N(false, false, false), //
	TMT_10PLEX_127C(false, false, false), //
	TMT_10PLEX_128N(false, false, false), //
	TMT_10PLEX_128C(false, false, false), //
	TMT_10PLEX_129N(false, false, false), //
	TMT_10PLEX_129C(false, false, false), //
	TMT_10PLEX_130N(false, false, false), //
	TMT_10PLEX_130C(false, false, false), //
	TMT_10PLEX_131(false, false, true), //
	//
	TMT_11PLEX_126(true, false, false), //
	TMT_11PLEX_127N(false, false, false), //
	TMT_11PLEX_127C(false, false, false), //
	TMT_11PLEX_128N(false, false, false), //
	TMT_11PLEX_128C(false, false, false), //
	TMT_11PLEX_129N(false, false, false), //
	TMT_11PLEX_129C(false, false, false), //
	TMT_11PLEX_130N(false, false, false), //
	TMT_11PLEX_130C(false, false, false), //
	TMT_11PLEX_131(false, false, true), //
	TMT_11PLEX_131C(false, false, true), //
	//
	TMT_16PLEX_126(true, false, false), //
	TMT_16PLEX_127N(false, false, false), //
	TMT_16PLEX_127C(false, false, false), //
	TMT_16PLEX_128N(false, false, false), //
	TMT_16PLEX_128C(false, false, false), //
	TMT_16PLEX_129N(false, false, false), //
	TMT_16PLEX_129C(false, false, false), //
	TMT_16PLEX_130N(false, false, false), //
	TMT_16PLEX_130C(false, false, false), //
	TMT_16PLEX_131N(false, false, false), //
	TMT_16PLEX_131C(false, false, false), //
	TMT_16PLEX_132N(false, false, false), //
	TMT_16PLEX_132C(false, false, false), //
	TMT_16PLEX_133N(false, false, false), //
	TMT_16PLEX_133C(false, false, false), //
	TMT_16PLEX_134(false, false, true), //
	//
	N14(true, false, false), //
	N15(false, false, true), //

	TMT_CHANNEL(false, false, false);

	private final static Logger log = Logger.getLogger(QuantificationLabel.class);
	private final boolean isLight;
	private final boolean isMedium;
	private final boolean isHeavy;
	private String customNameForTMT;

	private QuantificationLabel(boolean isLight, boolean isMedium, boolean isHeavy) {
		this(isLight, isMedium, isHeavy, null);
	}

	private QuantificationLabel(boolean isLight, boolean isMedium, boolean isHeavy, String customName) {
		this.isLight = isLight;
		this.isMedium = isMedium;
		this.isHeavy = isHeavy;
		this.setCustomNameForTMT(customName);
	}

	public static QuantificationLabel getByName(String labelName) {
		if (labelName == null)
			return null;
		if (labelName.equalsIgnoreCase("l"))
			return LIGHT;
		if (labelName.equalsIgnoreCase("h"))
			return HEAVY;
		if (labelName.equalsIgnoreCase("m"))
			return MEDIUM;
		final QuantificationLabel[] values = values();
		for (final QuantificationLabel quantificationLabel : values) {
			if (quantificationLabel.name().equalsIgnoreCase(labelName))
				return quantificationLabel;

		}

		log.error("Quantification label '" + labelName + "' not recognized");
		return null;
	}

	/**
	 * @return the isLight
	 */
	public boolean isLight() {
		return isLight;
	}

	/**
	 * @return the isLight
	 */
	public boolean isMedium() {

		return isMedium;
	}

	/**
	 * @return the isHeavy
	 */
	public boolean isHeavy() {
		return isHeavy;
	}

	public static String getValuesString() {
		final StringBuilder sb = new StringBuilder();
		for (final QuantificationLabel label : values()) {
			if (!"".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(label.name());
		}
		return sb.toString();
	}

	public static boolean isN15(QuantificationLabel label) {
		if (label != null) {
			if (label == QuantificationLabel.N14 || label == QuantificationLabel.N15) {
				return true;
			}
		}
		return false;
	}

	public static boolean isTMT4PLEX(QuantificationLabel label) {
		return getTMT4PlexLabels().contains(label);
	}

	public static boolean isTMT6PLEX(QuantificationLabel label) {
		return getTMT6PlexLabels().contains(label);
	}

	public static boolean isTMT9PLEX(QuantificationLabel label) {
		return getTMT9PlexLabels().contains(label);
	}

	public static boolean isTMT10PLEX(QuantificationLabel label) {
		return getTMT10PlexLabels().contains(label);
	}

	public static boolean isTMT11PLEX(QuantificationLabel label) {
		return getTMT11PlexLabels().contains(label);
	}

	public static boolean isTMT16PLEX(QuantificationLabel label) {
		return getTMT16PlexLabels().contains(label);
	}

	public static List<QuantificationLabel> getTMT4PlexLabels() {
		final QuantificationLabel[] array = { QuantificationLabel.TMT_4PLEX_127, QuantificationLabel.TMT_4PLEX_128,
				QuantificationLabel.TMT_4PLEX_129, QuantificationLabel.TMT_4PLEX_130 };
		final List<QuantificationLabel> ret = new ArrayList<QuantificationLabel>();
		for (final QuantificationLabel quantificationLabel : array) {
			ret.add(quantificationLabel);
		}
		return ret;
	}

	public static List<QuantificationLabel> getTMT6PlexLabels() {
		final QuantificationLabel[] array = { QuantificationLabel.TMT_6PLEX_126, QuantificationLabel.TMT_6PLEX_127,
				QuantificationLabel.TMT_6PLEX_128, QuantificationLabel.TMT_6PLEX_129, QuantificationLabel.TMT_6PLEX_130,
				QuantificationLabel.TMT_6PLEX_131 };
		final List<QuantificationLabel> ret = new ArrayList<QuantificationLabel>();
		for (final QuantificationLabel quantificationLabel : array) {
			ret.add(quantificationLabel);
		}
		return ret;
	}

	public static List<QuantificationLabel> getTMT9PlexLabels() {
		final QuantificationLabel[] array = { QuantificationLabel.TMT_9PLEX_126, QuantificationLabel.TMT_9PLEX_127N,
				QuantificationLabel.TMT_9PLEX_127C, QuantificationLabel.TMT_9PLEX_128N,
				QuantificationLabel.TMT_9PLEX_128C, QuantificationLabel.TMT_9PLEX_129N,
				QuantificationLabel.TMT_9PLEX_129C, QuantificationLabel.TMT_9PLEX_130N,
				QuantificationLabel.TMT_9PLEX_130C };
		final List<QuantificationLabel> ret = new ArrayList<QuantificationLabel>();
		for (final QuantificationLabel quantificationLabel : array) {
			ret.add(quantificationLabel);
		}
		return ret;
	}

	public static List<QuantificationLabel> getTMT10PlexLabels() {
		final QuantificationLabel[] array = { QuantificationLabel.TMT_10PLEX_126, QuantificationLabel.TMT_10PLEX_127N,
				QuantificationLabel.TMT_10PLEX_127C, QuantificationLabel.TMT_10PLEX_128N,
				QuantificationLabel.TMT_10PLEX_128C, QuantificationLabel.TMT_10PLEX_129N,
				QuantificationLabel.TMT_10PLEX_129C, QuantificationLabel.TMT_10PLEX_130N,
				QuantificationLabel.TMT_10PLEX_130C, QuantificationLabel.TMT_10PLEX_131 };
		final List<QuantificationLabel> ret = new ArrayList<QuantificationLabel>();
		for (final QuantificationLabel quantificationLabel : array) {
			ret.add(quantificationLabel);
		}
		return ret;
	}

	public static List<QuantificationLabel> getTMT11PlexLabels() {
		final QuantificationLabel[] array = { QuantificationLabel.TMT_11PLEX_126, QuantificationLabel.TMT_11PLEX_127N,
				QuantificationLabel.TMT_11PLEX_127C, QuantificationLabel.TMT_11PLEX_128N,
				QuantificationLabel.TMT_11PLEX_128C, QuantificationLabel.TMT_11PLEX_129N,
				QuantificationLabel.TMT_11PLEX_129C, QuantificationLabel.TMT_11PLEX_130N,
				QuantificationLabel.TMT_11PLEX_130C, QuantificationLabel.TMT_11PLEX_131,
				QuantificationLabel.TMT_11PLEX_131C };
		final List<QuantificationLabel> ret = new ArrayList<QuantificationLabel>();
		for (final QuantificationLabel quantificationLabel : array) {
			ret.add(quantificationLabel);
		}
		return ret;
	}

	public static List<QuantificationLabel> getTMT16PlexLabels() {
		final QuantificationLabel[] array = { QuantificationLabel.TMT_16PLEX_126, QuantificationLabel.TMT_16PLEX_127N,
				QuantificationLabel.TMT_16PLEX_127C, QuantificationLabel.TMT_16PLEX_128N,
				QuantificationLabel.TMT_16PLEX_128C, QuantificationLabel.TMT_16PLEX_129N,
				QuantificationLabel.TMT_16PLEX_129C, QuantificationLabel.TMT_16PLEX_130N,
				QuantificationLabel.TMT_16PLEX_130C, QuantificationLabel.TMT_16PLEX_131N,
				QuantificationLabel.TMT_16PLEX_131C, QuantificationLabel.TMT_16PLEX_132N,
				QuantificationLabel.TMT_16PLEX_132C, QuantificationLabel.TMT_16PLEX_133N,
				QuantificationLabel.TMT_16PLEX_133C, QuantificationLabel.TMT_16PLEX_134 };
		final List<QuantificationLabel> ret = new ArrayList<QuantificationLabel>();
		for (final QuantificationLabel quantificationLabel : array) {
			ret.add(quantificationLabel);
		}
		return ret;
	}

	/**
	 * 
	 * @param plex
	 * @param tmtHeaders the actual headers we read from data file
	 * @return
	 */
	public static List<QuantificationLabel> getTMTPlexLabels(int plex) {
		if (plex == 4) {
			return getTMT4PlexLabels();
		} else if (plex == 6) {
			return getTMT6PlexLabels();
		} else if (plex == 9) {
			return getTMT9PlexLabels();
		} else if (plex == 10) {
			return getTMT10PlexLabels();
		} else if (plex == 11) {
			return getTMT11PlexLabels();
		} else if (plex == 16) {
			return getTMT16PlexLabels();
		}

		throw new IllegalArgumentException(
				"TMT plex " + plex + " is not yet supported. Contact salvador@scripps.edu to fix that!");
	}

	public static int getTMTPlex(QuantificationLabel label) {
		if (getTMT4PlexLabels().contains(label)) {
			return 4;
		} else if (getTMT6PlexLabels().contains(label)) {
			return 6;
		} else if (getTMT9PlexLabels().contains(label)) {
			return 9;
		} else if (getTMT10PlexLabels().contains(label)) {
			return 10;
		} else if (getTMT11PlexLabels().contains(label)) {
			return 11;
		} else if (getTMT16PlexLabels().contains(label)) {
			return 16;
		}
		throw new IllegalArgumentException("Label " + label + " is not detected as a TMT label");
	}

	public String getCustomNameForTMT() {
		return customNameForTMT;
	}

	public void setCustomNameForTMT(String customNameForTMT) {
		this.customNameForTMT = customNameForTMT;
	}
}
