package edu.scripps.yates.census.analysis;

import java.io.File;

public class QuantParameters {
	public static final long DEFAULT_TIMEOUT = 1000 * 60 * 1;// 1 min
	public static final int DEFAULT_MAX_ITERATIONS = 300;

	private boolean performCalibration = true;
	private Double outlierRemovalFDR = 0.01;
	private File sanxotScriptsFolder;
	private long timeout;
	private int maxIterations;
	private boolean usePython = true; // by default
	private String ratioName;

	public QuantParameters(boolean calibration, Double outlierRemovalFDR, File sanxotScriptsFolder, long timeout,
			int maxIterations, boolean usePython) {
		performCalibration = calibration;
		this.outlierRemovalFDR = outlierRemovalFDR;
		this.sanxotScriptsFolder = sanxotScriptsFolder;
		setTimeout(timeout);
		setMaxIterations(maxIterations);
		setUsePython(usePython);
	}

	public QuantParameters() {

	}

	/**
	 * @return the performCalibration
	 */
	public boolean isPerformCalibration() {
		return performCalibration;
	}

	/**
	 * @return the fdr
	 */
	public Double getOutlierRemovalFDR() {
		return outlierRemovalFDR;
	}

	/**
	 * @param performCalibration
	 *            the performCalibration to set
	 */
	public void setPerformCalibration(boolean performCalibration) {
		this.performCalibration = performCalibration;
	}

	/**
	 * @param fdr
	 *            the fdr to set
	 */
	public void setOutlierRemovalFDR(Double fdr) {
		outlierRemovalFDR = fdr;
	}

	/**
	 * @return the sanxotScriptsFolder
	 */
	public File getSanxotScriptsFolder() {
		return sanxotScriptsFolder;
	}

	/**
	 * @param sanxotScriptsFolder
	 *            the sanxotScriptsFolder to set
	 */
	public void setSanxotScriptsFolder(File sanxotScriptsFolder) {
		this.sanxotScriptsFolder = sanxotScriptsFolder;
	}

	/**
	 * @return the timeout
	 */
	public long getTimeout() {
		return timeout;
	}

	/**
	 * @param timeout
	 *            the timeout to set
	 */
	public void setTimeout(long timeout) {
		this.timeout = timeout;
	}

	/**
	 * @return the maxIterations
	 */
	public int getMaxIterations() {
		return maxIterations;
	}

	/**
	 * @param maxIterations
	 *            the maxIterations to set
	 */
	public void setMaxIterations(int maxIterations) {
		this.maxIterations = maxIterations;
	}

	/**
	 * @return the usePython
	 */
	public boolean isUsePython() {
		return usePython;
	}

	/**
	 * @param usePython
	 *            the usePython to set
	 */
	public void setUsePython(boolean usePython) {
		this.usePython = usePython;
	}

	public String getRatioName() {
		return ratioName;
	}

	public void setRatioName(String ratioName) {
		this.ratioName = ratioName;
	}
}
