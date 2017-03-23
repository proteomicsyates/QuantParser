package edu.scripps.yates.census.read.util;

import org.apache.log4j.Logger;

public enum QuantificationLabel {

	LIGHT(true, false, false), HEAVY(false, false, true), MEDIUM(false, true, false), //
	TMT_6PLEX_126(true, false, false), TMT_6PLEX_127(false, false, false), TMT_6PLEX_128(false, false,
			false), TMT_6PLEX_129(false, false,
					false), TMT_6PLEX_130(false, false, false), TMT_6PLEX_131(false, false, true), //
	N14(true, false, false), N15(false, false, true);
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
		for (QuantificationLabel quantificationLabel : values) {
			if (quantificationLabel.name().equalsIgnoreCase(labelName))
				return quantificationLabel;

		}
		if (labelName.contains("126")) {
			return QuantificationLabel.TMT_6PLEX_126;
		}
		if (labelName.contains("127")) {
			return QuantificationLabel.TMT_6PLEX_127;
		}
		if (labelName.contains("128")) {
			return QuantificationLabel.TMT_6PLEX_128;
		}
		if (labelName.contains("129")) {
			return QuantificationLabel.TMT_6PLEX_129;
		}
		if (labelName.contains("130")) {
			return QuantificationLabel.TMT_6PLEX_130;
		}
		if (labelName.contains("131")) {
			return QuantificationLabel.TMT_6PLEX_131;
		}
		log.warn("Quantification label '" + labelName + "' not recognized");
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
		StringBuilder sb = new StringBuilder();
		for (QuantificationLabel label : values()) {
			if (!"".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(label.name());
		}
		return sb.toString();
	}
}
