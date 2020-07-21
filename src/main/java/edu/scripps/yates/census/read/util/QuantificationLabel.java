package edu.scripps.yates.census.read.util;

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

public enum QuantificationLabel {

	LIGHT(true, false, false), //
	HEAVY(false, false, true), //
	MEDIUM(false, true, false), //
	//
	TMT_6PLEX_126(true, false, false), //
	TMT_6PLEX_127(false, false, false), //
	TMT_6PLEX_128(false, false, false), //
	TMT_6PLEX_129(false, false, false), //
	TMT_6PLEX_130(false, false, false), //
	TMT_6PLEX_131(false, false, true), //
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
	N14(true, false, false), //
	N15(false, false, true);

	private final static Logger log = Logger.getLogger(QuantificationLabel.class);
	private final boolean isLight;
	private final boolean isMedium;
	private final boolean isHeavy;

	private QuantificationLabel(boolean isLight, boolean isMedium, boolean isHeavy) {
		this.isLight = isLight;
		this.isMedium = isMedium;
		this.isHeavy = isHeavy;
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

	public static boolean isTMT6PLEX(QuantificationLabel label) {
		return getTMT6PlexLabels().contains(label);
	}

	public static boolean isTMT10PLEX(QuantificationLabel label) {
		return getTMT10PlexLabels().contains(label);
	}

	public static boolean isTMT11PLEX(QuantificationLabel label) {
		return getTMT11PlexLabels().contains(label);
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

	public static List<QuantificationLabel> getTMT10PlexLabels() {
		final QuantificationLabel[] array = { QuantificationLabel.TMT_10PLEX_126,
				QuantificationLabel.TMT_10PLEX_127N, QuantificationLabel.TMT_10PLEX_127C,
				QuantificationLabel.TMT_10PLEX_128N, QuantificationLabel.TMT_10PLEX_128C,
				QuantificationLabel.TMT_10PLEX_129N, QuantificationLabel.TMT_10PLEX_129C,
				QuantificationLabel.TMT_10PLEX_130N, QuantificationLabel.TMT_10PLEX_130C,
				QuantificationLabel.TMT_10PLEX_131 };
		final List<QuantificationLabel> ret = new ArrayList<QuantificationLabel>();
		for (final QuantificationLabel quantificationLabel : array) {
			ret.add(quantificationLabel);
		}
		return ret;
	}

	public static List<QuantificationLabel> getTMT11PlexLabels() {
		final QuantificationLabel[] array = { QuantificationLabel.TMT_11PLEX_126,
				QuantificationLabel.TMT_11PLEX_127N, QuantificationLabel.TMT_11PLEX_127C,
				QuantificationLabel.TMT_11PLEX_128N, QuantificationLabel.TMT_11PLEX_128C,
				QuantificationLabel.TMT_11PLEX_129N, QuantificationLabel.TMT_11PLEX_129C,
				QuantificationLabel.TMT_11PLEX_130N, QuantificationLabel.TMT_11PLEX_130C,
				QuantificationLabel.TMT_11PLEX_131, QuantificationLabel.TMT_11PLEX_131C };
		final List<QuantificationLabel> ret = new ArrayList<QuantificationLabel>();
		for (final QuantificationLabel quantificationLabel : array) {
			ret.add(quantificationLabel);
		}
		return ret;
	}
}
