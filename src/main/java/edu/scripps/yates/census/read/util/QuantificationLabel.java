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
	TMT_10PLEX_126_127726(true, false, false), //
	TMT_10PLEX_127_124761(false, false, false), //
	TMT_10PLEX_127_131081(false, false, false), //
	TMT_10PLEX_128_128116(false, false, false), //
	TMT_10PLEX_128_134436(false, false, false), //
	TMT_10PLEX_129_131471(false, false, false), //
	TMT_10PLEX_129_13779(false, false, false), //
	TMT_10PLEX_130_134825(false, false, false), //
	TMT_10PLEX_130_141145(false, false, false), //
	TMT_10PLEX_131_13818(false, false, true), //
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
		final QuantificationLabel[] array = { QuantificationLabel.TMT_10PLEX_126_127726,
				QuantificationLabel.TMT_10PLEX_127_124761, QuantificationLabel.TMT_10PLEX_127_131081,
				QuantificationLabel.TMT_10PLEX_128_128116, QuantificationLabel.TMT_10PLEX_128_134436,
				QuantificationLabel.TMT_10PLEX_129_131471, QuantificationLabel.TMT_10PLEX_129_13779,
				QuantificationLabel.TMT_10PLEX_130_134825, QuantificationLabel.TMT_10PLEX_130_141145,
				QuantificationLabel.TMT_10PLEX_131_13818 };
		final List<QuantificationLabel> ret = new ArrayList<QuantificationLabel>();
		for (final QuantificationLabel quantificationLabel : array) {
			ret.add(quantificationLabel);
		}
		return ret;
	}
}
