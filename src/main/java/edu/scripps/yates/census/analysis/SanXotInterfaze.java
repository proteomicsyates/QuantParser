package edu.scripps.yates.census.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;

import javax.swing.SwingWorker;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.census.analysis.wrappers.IntegrationResultWrapper;
import edu.scripps.yates.census.analysis.wrappers.KalibrateResultWrapper;
import edu.scripps.yates.census.analysis.wrappers.OutlierRemovalResultWrapper;
import edu.scripps.yates.census.analysis.wrappers.SanXotAnalysisResult;
import edu.scripps.yates.census.analysis.wrappers.SanxotQuantResult;
import edu.scripps.yates.census.read.util.FileSplitter;
import edu.scripps.yates.utilities.exec.CommandLineRunner;
import edu.scripps.yates.utilities.exec.ProcessExecutor;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class SanXotInterfaze extends SwingWorker<Object, Void> {
	private static final Logger log = Logger.getLogger(SanXotInterfaze.class);
	private static final String SANXOT_EXE = "sanxot.exe";
	private static final String KLIBRATE_EXE = "klibrate.exe";
	private static final String SANXOT_SIEVE_EXE = "sanxotsieve.exe";
	private static final String SANXOT_PY = "sanxot.py";
	private static final String KLIBRATE_PY = "klibrate.py";
	private static final String SANXOT_SIEVE_PY = "sanxotsieve.py";

	// fire change event
	public static final String CALIBRATING = "calibrating";
	public static final String CALIBRATING_DONE = "calibrating_done";
	public static final String CALIBRATING_ERROR = "calibrating_error";
	public static final String INTEGRATING = "integrating";
	public static final String INTEGRATING_DONE = "integrating_done";
	private static final String OUTLIER_REMOVAL = "outlier removal";
	private static final String OUTLIER_REMOVAL_DONE = "outlier removal done";
	public static final String STARTING_COMMAND = "starting_command";
	public static final String END_COMMAND = "end_command";

	public static final String END_ANALYSIS = "end analysis";
	private static final String PYTHON = "python";

	private final FileMappingResults fileMappingResults;

	private FileWriter logFileWriter;
	private final SanXotAnalysisResult result;
	private boolean keepExperimentsSeparated;
	private final QuantParameters quantParameters;

	public SanXotInterfaze(FileMappingResults fileMappingResults, QuantParameters quantParameters) {
		this.quantParameters = quantParameters;
		this.fileMappingResults = fileMappingResults;
		// TODO customize the FDR depending on the level
		result = new SanXotAnalysisResult(fileMappingResults);

	}

	/**
	 * @return the keepExperimentsSeparated
	 */
	public boolean isKeepExperimentsSeparated() {
		return keepExperimentsSeparated;
	}

	/**
	 * @param keepExperimentsSeparated
	 *            the keepExperimentsSeparated to set
	 */
	public void setKeepExperimentsSeparated(boolean keepExperimentsSeparated) {
		this.keepExperimentsSeparated = keepExperimentsSeparated;
	}

	public void setSanxotLocation(File folder) {
		quantParameters.setSanxotScriptsFolder(folder);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see javax.swing.SwingWorker#done()
	 */
	@Override
	protected void done() {
		log.info("Sanxot done");
		if (isCancelled())
			log.info("Sanxot cancelled");
		if (logFileWriter != null) {
			try {
				logFileWriter.close();
			} catch (final IOException e) {
				e.printStackTrace();
			}
		}

		firePropertyChange(END_ANALYSIS, null, result);

		super.done();
	}

	@Override
	protected Object doInBackground() throws Exception {
		analyze();
		return null;
	}

	public void analyze() throws IOException, InterruptedException, ExecutionException {
		logFileWriter = new FileWriter(new File(System.getProperty("user.dir") + File.separator + "sanxot_cmds.log"));

		log.info("Starting SanXot interfaze");
		log.info("Timeout set at " + quantParameters.getTimeout() / 1000 + " sg.");
		int lowLevel = 0;
		int upperLevel = 0;
		File dataFile = fileMappingResults.getDataFile();
		File infoFile = null;
		try {
			final Map<String, List<String>> experimentAndReplicateNames = fileMappingResults
					.getExperimentAndReplicateNames();
			Pair<Integer, File> lowLevelPair = fileMappingResults.getFirstLevel();
			lowLevel = lowLevelPair.getFirstelement();
			Pair<Integer, File> upperLevelPair = fileMappingResults.getNextAvailableLevel(lowLevel);
			upperLevel = upperLevelPair.getFirstelement();
			File relatFile = lowLevelPair.getSecondElement();

			for (final String experimentName : experimentAndReplicateNames.keySet()) {
				log.info("Experiment: " + experimentName);
				// get the names of the replicates and experiments in order to
				// split the relation files
				final List<String> dataSetNames = getDataSetNames(experimentAndReplicateNames, experimentName);

				// split relatFile in many files as replicates
				// Map<String, File> relatFiles = FileSplitter
				// .splitFiles(fileMappingResults.getFirstLevel().getSecondElement(),
				// dataSetNames);
				relatFile = fileMappingResults.getFirstLevel().getSecondElement();

				// split dataFile in many files as replicates
				Map<String, File> dataFiles = null;
				if (dataSetNames.size() == 1 && dataSetNames.get(0).equals("")) {
					dataFiles = new THashMap<String, File>();
					dataFiles.put("", fileMappingResults.getDataFile());
				} else {
					dataFiles = FileSplitter.splitFiles(fileMappingResults.getDataFile(), dataSetNames);
				}

				final Map<String, File> calibratedDataFiles = new THashMap<String, File>();
				for (final String datasetName : dataSetNames) {
					// relatFile = relatFiles.get(datasetName);
					dataFile = dataFiles.get(datasetName);
					// get lowLevels again
					lowLevel = fileMappingResults.getFirstLevel().getFirstelement();
					upperLevel = fileMappingResults.getNextAvailableLevel(lowLevel).getFirstelement();
					// CALIBRATION
					if (quantParameters.isPerformCalibration()) {
						final KalibrateResultWrapper calibrationResult = calibrate(lowLevel, upperLevel, relatFile,
								dataFile, "_" + datasetName, quantParameters.getTimeout());
						result.setKalibrationResult(calibrationResult);
						if (calibrationResult.getCalibratedDataFile() != null) {
							dataFile = calibrationResult.getCalibratedDataFile();
						}
						calibratedDataFiles.put(datasetName, dataFile);
					} else {
						calibratedDataFiles.put(datasetName, dataFile);
					}
				}

				final Map<String, File> lastDataFiles = new THashMap<String, File>();
				// loop
				boolean dataMergingNeeded = false;
				for (final String replicateName : dataSetNames) {
					try {
						log.info("Replicate: " + replicateName);
						// relatFile = relatFiles.get(replicateName);
						dataFile = calibratedDataFiles.get(replicateName);
						// get lowLevels again
						lowLevel = fileMappingResults.getFirstLevel().getFirstelement();
						upperLevel = fileMappingResults.getNextAvailableLevel(lowLevel).getFirstelement();
						dataMergingNeeded = false;
						while (!dataMergingNeeded && upperLevel <= fileMappingResults.getMaxLevel()) {
							dataMergingNeeded = fileMappingResults.isDataMergingNeeded(lowLevel);

							relatFile = fileMappingResults.getFileLevel(lowLevel);
							// split relatFile in many files as datasets
							// relatFiles = FileSplitter.splitFiles(relatFile,
							// dataSetNames);
							// relatFile = relatFiles.get(replicateName);
							lowLevelPair = fileMappingResults.getFilePairLevel(lowLevel);
							upperLevelPair = fileMappingResults.getNextAvailableLevel(lowLevelPair.getFirstelement());
							upperLevel = upperLevelPair.getFirstelement();
							log.info(lowLevel + " to " + upperLevel + " on " + replicateName + " "
									+ fileMappingResults.isDataMergingNeeded(lowLevel));
							IntegrationResultWrapper integrationResult = null;
							if (dataMergingNeeded) {
								log.info("Skipping this step. It is not necessary");
								// integrationResult = integrate(lowLevel,
								// upperLevel, null, dataFile, null,
								// datasetName,
								// null);
								//
								// result.addIntegrationResult(integrationResult);
							} else {
								// first integration, estimating the variance by
								// fitting
								// algorithm

								integrationResult = integrate(lowLevel, upperLevel, relatFile, dataFile, null,
										replicateName, null, true);

								result.addIntegrationResult(integrationResult);
								result.addReplicateExperimentIntegrationResult(integrationResult, experimentName,
										replicateName);
								if (dataSetNames.size() == 1) {
									// add this integration result as
									// experiment level
									result.addExperimentIntegrationResult(integrationResult, experimentName);
								}
								// remove outliers
								if (quantParameters.getOutlierRemovalFDR() != null) {
									// not perform in the last interation
									if (fileMappingResults.getNextAvailableLevel(upperLevel) == null) {
									} else {
										infoFile = integrationResult.getInfoFile();
										final OutlierRemovalResultWrapper removeOutliers = removeOutliers(
												lowLevelPair.getFirstelement(), upperLevelPair.getFirstelement(),
												lowLevelPair.getSecondElement(), dataFile, infoFile, quantParameters);
										relatFile = removeOutliers.getRelatFile();
										// second integration, using the
										// variance
										// calculated in previous one (use of
										// infoFile
										// for forcing the
										// variance to be that one)

										integrationResult = integrate(lowLevel, upperLevel, relatFile, dataFile,
												infoFile, replicateName, null, false);
										result.addIntegrationResult(integrationResult);
										integrationResult.setOutlierRemovalResult(removeOutliers);
										result.addReplicateExperimentIntegrationResult(integrationResult,
												experimentName, replicateName);
										if (dataSetNames.size() == 1) {
											// add this integration result as
											// experiment level
											result.addExperimentIntegrationResult(integrationResult, experimentName);
										}
									}
								}
								// get datafile from higherlevel data file
								dataFile = integrationResult.getHigherLevelDataFile();
								lastDataFiles.put(replicateName, dataFile);

								// next level
								lowLevel = upperLevel;
							}
						}
					} catch (final NextLevelException e) {
						log.debug(e);
						// do nothing
					}
				}
				// data merging
				if (dataMergingNeeded) {
					lowLevelPair = fileMappingResults.getFilePairLevel(lowLevel);
					upperLevelPair = fileMappingResults.getFilePairLevel(upperLevel);
					relatFile = lowLevelPair.getSecondElement();
					if (lastDataFiles.size() > 1) {
						// merge the latest data files into only one
						// Concatenate higherlevel result files in a single one
						final File mergedDataFile = new File(fileMappingResults.getWorkingFolder().getAbsolutePath()
								+ File.separator + experimentName + "_" + lowLevel + "_" + upperLevel + "_merged.tsv");
						FileUtils.mergeFiles(lastDataFiles.values(), mergedDataFile, true);

						// // up one level
						// lowLevelPair = upperLevelPair;
						// lowLevel = lowLevelPair.getFirstelement();
						// upperLevelPair = fileMappingResults
						// .getNextAvailableLevel(upperLevel);
						// upperLevel = upperLevelPair.getFirstelement();
						// relatFile = lowLevelPair.getSecondElement();
						final IntegrationResultWrapper integrationResult = integrate(lowLevel, upperLevel, relatFile,
								mergedDataFile, null, experimentName, null, true);
						// keep the results in a Map by experiment name
						result.addExperimentIntegrationResult(integrationResult, experimentName);
					} else {
						// add the replicate integration to the experiment
						// integration
						final IntegrationResultWrapper replicateIntegration = result
								.getReplicateIntegrationResultsByExperiment().get(experimentName).values().iterator()
								.next();
						result.addExperimentIntegrationResult(replicateIntegration, experimentName);
					}
				}
			}

			// now, merge results at experiment level if there is more than one
			// experiment
			// lowLevelPair = upperLevelPair;
			// lowLevel = lowLevelPair.getFirstelement();
			// upperLevelPair =
			// fileMappingResults.getNextAvailableLevel(upperLevel);
			// upperLevel = upperLevelPair.getFirstelement();

			// first integration without relation file (option -C)
			// for (String experimentName :
			// experimentIntegrationResults.keySet()) {
			// dataFile =
			// experimentIntegrationResults.get(experimentName).getHigherLevelDataFile();
			// final IntegrationResultWrapper integrationResult =
			// integrate(lowLevel, upperLevel, null, dataFile,
			// null, experimentName, null);
			// experimentIntegrationResults.put(experimentName,
			// integrationResult);
			// }
			if (experimentAndReplicateNames.size() > 1) {
				if (keepExperimentsSeparated) {
					firePropertyChange(END_ANALYSIS, null, result.getLastIntegrationResults());
					return;
				}
				// merge all higher data files
				final Set<File> higherLevelDataResults = new THashSet<File>();
				for (final IntegrationResultWrapper integrationResult : result.getExperimentIntegrationResults()
						.values()) {
					higherLevelDataResults.add(integrationResult.getHigherLevelDataFile());
				}
				final File mergedFile = new File(fileMappingResults.getWorkingFolder().getAbsolutePath()
						+ File.separator + "experiment_data_merged.tsv");
				FileUtils.mergeFiles(higherLevelDataResults, mergedFile, true);
				// integrate that in the next level
				lowLevelPair = upperLevelPair;
				lowLevel = lowLevelPair.getFirstelement();
				upperLevelPair = fileMappingResults.getNextAvailableLevel(upperLevel);
				if (upperLevelPair != null) {
					upperLevel = upperLevelPair.getFirstelement();
					final IntegrationResultWrapper integrationResult = integrate(lowLevel, upperLevel,
							lowLevelPair.getSecondElement(), mergedFile, null, "_TOTAL", null, true);
					result.addIntegrationResult(integrationResult);
				} else {
					upperLevel = lowLevel;
				}
				// make the last integration
				final IntegrationResultWrapper lastIntegrationResults = result.getLastIntegrationResults();
				final IntegrationResultWrapper integrationResult = integrate(lowLevel, upperLevel, null,
						lastIntegrationResults.getHigherLevelDataFile(), null, "", null, true);
				result.addIntegrationResult(integrationResult);
				firePropertyChange(END_ANALYSIS, null, integrationResult);
			} else {
				// make the last integration
				// integrate that in the next level
				lowLevelPair = upperLevelPair;
				lowLevel = lowLevelPair.getFirstelement();
				final IntegrationResultWrapper lastIntegrationResults = result.getLastIntegrationResults();
				final IntegrationResultWrapper integrationResult = integrate(lowLevel, upperLevel, null,
						lastIntegrationResults.getHigherLevelDataFile(), null, "", null, true);
				result.addIntegrationResult(integrationResult);
			}
		} catch (final NextLevelException e) {
			// make the last integration
			final IntegrationResultWrapper lastIntegrationResults = result.getLastIntegrationResults();
			final IntegrationResultWrapper integrationResult = integrate(lowLevel, upperLevel, null,
					lastIntegrationResults.getHigherLevelDataFile(), null, "", null, true);
			result.addIntegrationResult(integrationResult);
		}
	}

	private List<String> getDataSetNamesOLD(Map<String, List<String>> experimentAndReplicateNames,
			String experimentName) {
		final boolean onlyOneExperiment = experimentAndReplicateNames.size() == 1;
		final List<String> replicateNames = experimentAndReplicateNames.get(experimentName);
		final boolean onlyOneReplicate = replicateNames.size() == 1;

		final List<String> datasetNames = new ArrayList<String>();
		for (final String replicateName : replicateNames) {
			String datasetName = onlyOneExperiment ? "" : experimentName;
			if (!onlyOneReplicate) {
				if (!"".endsWith(datasetName))
					datasetName = "_" + datasetName;
				datasetName = replicateName + datasetName;
			}

			datasetNames.add(datasetName);
		}
		return datasetNames;
	}

	private List<String> getDataSetNames(Map<String, List<String>> experimentAndReplicateNames, String experimentName) {
		final List<String> replicateNames = experimentAndReplicateNames.get(experimentName);
		String experimentKey = "";
		if (experimentAndReplicateNames.size() > 1) {
			experimentKey = experimentName;
		}

		final List<String> datasetNames = new ArrayList<String>();
		if (replicateNames.size() == 1) {

			final String datasetName = experimentKey;
			datasetNames.add(datasetName);

		} else {
			for (final String replicateName : replicateNames) {
				String datasetName = replicateName;
				if (!"".equals(experimentKey)) {
					datasetName += "_" + experimentKey;
				}
				datasetNames.add(datasetName);
			}
		}
		return datasetNames;
	}

	private OutlierRemovalResultWrapper removeOutliers(int lowLevel, int upperLevel, File relatFile, File dataFile,
			File infoFile, QuantParameters quantParameters)
			throws IOException, InterruptedException, ExecutionException {
		final String msg = "Removing outliers data from level " + lowLevel + " to " + upperLevel + "...";
		log.info(msg);
		firePropertyChange(OUTLIER_REMOVAL, null, msg);
		final String prefix = OutlierRemovalResultWrapper.DEFAULT_OUTLIER_REMOVAL_PREFIX + lowLevel + "-" + upperLevel;
		final CommandLine removeOutlierCommandLine = getRemoveOutliersCommandLine(relatFile, prefix, dataFile, infoFile,
				fileMappingResults.getWorkingFolder(), quantParameters);

		final CommandLineRunner runner = runCommand(removeOutlierCommandLine, quantParameters.getTimeout());
		if (runner.getProcessExitCode().longValue() != 0) {
			throw new IllegalArgumentException(
					"Some error happen while outlier removal process: " + runner.getErrorMessage());
		}
		final OutlierRemovalResultWrapper outliersRemovalResults = new OutlierRemovalResultWrapper(
				fileMappingResults.getWorkingFolder(), prefix);

		log.info("Outlier removal performed. New relation file at:"
				+ outliersRemovalResults.getRelatFile().getAbsolutePath());

		firePropertyChange(OUTLIER_REMOVAL_DONE, null, outliersRemovalResults);
		return outliersRemovalResults;

	}

	private IntegrationResultWrapper integrate(int lowLevel, int upperLevel, File relatFile, File dataFile,
			File infoFile, String prefix, Double forzedVariance, boolean checkRelationshipValidity)
			throws IOException, InterruptedException, ExecutionException {
		if (checkRelationshipValidity && !SanXotInterfaze.checkDataValidity(relatFile, dataFile)) {
			SanXotInterfaze.checkDataValidity(relatFile, dataFile);
			throw new IllegalArgumentException("Combination of data file and relat file is not valid: "
					+ FilenameUtils.getName(relatFile.getAbsolutePath()) + " and "
					+ FilenameUtils.getName(dataFile.getAbsolutePath()));
		}

		final String msg = "Integrating data from level " + lowLevel + " to " + upperLevel + "...";
		log.info(msg);
		if (relatFile == null) {
			log.info("Not using relationship file. Using -C option to correct protein loading error");
		} else {
			log.info("Using relationship file " + FilenameUtils.getName(relatFile.getAbsolutePath()));
		}
		log.info("Using data file  " + FilenameUtils.getName(dataFile.getAbsolutePath()));
		firePropertyChange(INTEGRATING, null, msg);
		final String prefixString = lowLevel + "-" + upperLevel + "_" + prefix;
		CommandLine integratingCommandLine = getIntegrationCommandLine(relatFile, dataFile, infoFile, prefixString,
				forzedVariance, fileMappingResults.getWorkingFolder(), quantParameters);

		final CommandLineRunner runner = runCommand(integratingCommandLine, quantParameters.getTimeout());
		if (runner.getProcessExitCode().longValue() != 0) {
			if (runner.getProcessExitCode().longValue() == ProcessExecutor.TIMEOUT_ERROR_CODE) {
				final String message = "The process cound't finish before the timeout of "
						+ quantParameters.getTimeout() + " ms";
				log.warn(message);
				log.info("Trying to fix the problem by forzing variance to 0 (Using -f v0)");
				integratingCommandLine = getIntegrationCommandLine(relatFile, dataFile, infoFile, prefixString, 0.0,
						fileMappingResults.getWorkingFolder(), quantParameters);
				final CommandLineRunner newRunner = runCommand(integratingCommandLine, quantParameters.getTimeout());
				if (newRunner.getProcessExitCode().longValue() != 0) {
					if (newRunner.getProcessExitCode().longValue() == ProcessExecutor.TIMEOUT_ERROR_CODE) {
						log.warn(message);
						throw new IllegalArgumentException(message + ": " + newRunner.getErrorMessage());
					}
					throw new IllegalArgumentException(
							"Some error happen while integration process: " + newRunner.getErrorMessage());
				}
			} else {
				throw new IllegalArgumentException(
						"Some error happen while integration process: " + runner.getErrorMessage());
			}
		}
		final IntegrationResultWrapper integrationResults = new IntegrationResultWrapper(
				fileMappingResults.getWorkingFolder(), prefixString, lowLevel, upperLevel, fileMappingResults);

		log.info("Integration performed. Integration file at:"
				+ integrationResults.getHigherLevelDataFile().getAbsolutePath());
		log.info("Integration Variance=" + integrationResults.getIntegrationVariance());
		firePropertyChange(INTEGRATING_DONE, null, integrationResults);
		return integrationResults;
	}

	/**
	 * Checks whether there is at least two lower elements in dataFile that maps
	 * to the same upper level item
	 *
	 * @param relatFile
	 * @param dataFile
	 * @return
	 */
	public static boolean checkDataValidity(File relatFile, File dataFile) {
		if (relatFile == null) {
			// ignore the check
			return true;
		}
		final Map<String, SanxotQuantResult> sanXotQuantResultFromFile = IntegrationResultWrapper
				.getSanXotQuantResultFromDataFile(dataFile);
		final Map<String, Set<String>> relationShipsFromRelatFile = IntegrationResultWrapper
				.getRelationShipsFromRelatFile(relatFile);
		final Set<String> lowerLevelFromRelat = new THashSet<String>();
		for (final String upperLevel : relationShipsFromRelatFile.keySet()) {
			lowerLevelFromRelat.addAll(relationShipsFromRelatFile.get(upperLevel));
		}
		// Set<String> upperLevelsMapped = new THashSet<String>();
		for (final String dataFileKey : sanXotQuantResultFromFile.keySet()) {
			if (!lowerLevelFromRelat.contains(dataFileKey)) {
				log.info(dataFileKey + " is  not found as lower level item in the relationship file "
						+ relatFile.getAbsolutePath());
				return false;
			}
			// boolean found = false;
			// for (String upperLevel : relationShipsFromRelatFile.keySet()) {
			// final Set<String> lowerLevels =
			// relationShipsFromRelatFile.get(upperLevel);
			// if (lowerLevels.contains(dataFileKey)) {
			// found = true;
			// // if (upperLevelsMapped.contains(upperLevel)) {
			// // return true;
			// // } else {
			// // upperLevelsMapped.add(upperLevel);
			// // }
			// }
			// }
			// if (!found) {
			// log.info(dataFileKey + " is not found as lower level item in the
			// relationship file "
			// + relatFile.getAbsolutePath());
			// return false;
			// }
		}
		return true;
	}

	private KalibrateResultWrapper calibrate(int lowLevel, int upperLevel, File relatFile, File dataFile, String key,
			long timeout) throws IOException, InterruptedException, ExecutionException {
		final String msg = "Calibrating data " + lowLevel + " - " + upperLevel + "...";

		log.info(msg);
		firePropertyChange(CALIBRATING, null, msg);
		final String prefix = KalibrateResultWrapper.DEFAULT_CALIBRATED_PREFIX + lowLevel + "-" + upperLevel + key;
		final CommandLine calibratingCommandLine = getCalibratingCommandLine(relatFile, dataFile, prefix,
				fileMappingResults.getWorkingFolder(), quantParameters);

		final CommandLineRunner runner = runCommand(calibratingCommandLine, timeout);
		if (runner.getProcessExitCode().longValue() != 0) {
			throw new SanxotErrorException("Some error happen while calibration process:\n" + runner.getErrorMessage()
					+ "\nYou may consider the following actions:\n"
					+ "-If you are running the software in a remote terminal, make sure that you have X11 activated.\n"
					+ "-For developers only: Increase the timeout time by saxotInterface.setTimeout(long timeout) method");
		}
		final KalibrateResultWrapper calibrationResults = new KalibrateResultWrapper(
				fileMappingResults.getWorkingFolder(), prefix);
		if (calibrationResults.getCalibratedDataFile() != null) {
			log.info("Calibration performed. Calibrated file at:"
					+ calibrationResults.getCalibratedDataFile().getAbsolutePath());
			log.info("Calibration K=" + calibrationResults.getCalibrationKConstant());
			log.info("Calibration Variance=" + calibrationResults.getCalibrationVariance());
			firePropertyChange(CALIBRATING_DONE, null, calibrationResults);
		} else {
			log.warn("Calibration of data file " + dataFile.getAbsolutePath() + " failed");
			firePropertyChange(CALIBRATING_ERROR, null, null);
			throw new SanxotErrorException("Error while calibrating. Command: " + calibratingCommandLine.toString());
		}
		return calibrationResults;
	}

	public static CommandLine getCalibratingCommandLine(File relatFile, File dataFile, String prefix,
			File workingFolder, QuantParameters quantParameters) {
		// klibrate -p"%directory%" -d%datafile% -r%level_1_2_file%
		// -a%level_1_2_prefix% -g
		CommandLine commandline = null;
		if (quantParameters.isUsePython()) {
			commandline = new CommandLine(PYTHON);

		} else {
			commandline = new CommandLine(quantParameters.getSanxotScriptsFolder() + File.separator + KLIBRATE_EXE);
		}
		// String dataFileName = FilenameUtils.getName(fileMappingResults
		// .getDataFile().getAbsolutePath());
		// final String dataFileName =
		// FilenameUtils.getName(dataFile.getAbsolutePath());
		final String dataFileName = dataFile.getAbsolutePath();

		// final String level1File =
		// FilenameUtils.getName(relatFile.getAbsolutePath());
		final String level1File = relatFile.getAbsolutePath();

		if (quantParameters.isUsePython()) {
			commandline.addArgument(quantParameters.getSanxotScriptsFolder() + File.separator + KLIBRATE_PY);
		}
		commandline.addArgument("-d" + dataFileName);
		commandline.addArgument("-r" + level1File);
		commandline.addArgument("-a" + prefix);
		commandline.addArgument("-p" + workingFolder.getAbsolutePath());
		commandline.addArgument("-m" + quantParameters.getMaxIterations());
		commandline.addArgument("-g");
		return commandline;
	}

	public static CommandLine getIntegrationCommandLine(File relatFile, File dataFile, File infoFile, String prefix,
			Double forcedVariance, File workingFolder, QuantParameters quantParameters) {
		// sanxot -d%level_1_2_prefix%_calibrated.xls -p"%directory%"
		// -r%level_1_2_file% -a%level_1_2_prefix%%outliers_sufix% -g
		CommandLine commandline = null;
		if (quantParameters.isUsePython()) {
			commandline = new CommandLine(PYTHON);
		} else {
			commandline = new CommandLine(quantParameters.getSanxotScriptsFolder() + File.separator + SANXOT_EXE);
		}

		final String dataFileName = FilenameUtils.getName(dataFile.getAbsolutePath());

		String infoFileName = null;
		if (infoFile != null)
			infoFileName = FilenameUtils.getName(infoFile.getAbsolutePath());
		if (quantParameters.isUsePython()) {
			commandline.addArgument(quantParameters.getSanxotScriptsFolder() + File.separator + SANXOT_PY);
		}
		commandline.addArgument("-d" + dataFileName);
		if (relatFile != null) {
			final String level1FileName = FilenameUtils.getName(relatFile.getAbsolutePath());
			commandline.addArgument("-r" + level1FileName);
		} else {
			commandline.addArgument("-C");
		}
		if (infoFile != null) {
			commandline.addArgument("-V" + infoFileName);
			commandline.addArgument("-f");
		}
		if (forcedVariance != null) {
			commandline.addArgument("-v" + String.valueOf(forcedVariance));
			commandline.addArgument("-f");
		}
		commandline.addArgument("-a" + prefix);
		commandline.addArgument("-p" + workingFolder.getAbsolutePath());
		commandline.addArgument("-g");
		return commandline;
	}

	public static CommandLine getRemoveOutliersCommandLine(File relatFile, String prefix, File dataFile, File infoFile,
			File workingFolder, QuantParameters quantParameters) {
		// sanxotsieve -d%level_1_2_prefix%_calibrated.xls -p"%directory%"
		// -r%level_1_2_file% -V%level_1_2_prefix%%outliers_sufix%_infoFile.txt
		// -a%level_1_2_prefix%-sanxotsieve-results -f%fdr_outliers_removal%
		CommandLine commandline = null;
		if (quantParameters.isUsePython()) {
			commandline = new CommandLine(PYTHON);
		} else {
			commandline = new CommandLine(quantParameters.getSanxotScriptsFolder() + File.separator + SANXOT_SIEVE_EXE);
		}
		final String dataFileName = FilenameUtils.getName(dataFile.getAbsolutePath());

		final String level1FileName = FilenameUtils.getName(relatFile.getAbsolutePath());
		final String infoFileName = FilenameUtils.getName(infoFile.getAbsolutePath());
		if (quantParameters.isUsePython()) {
			commandline.addArgument(quantParameters.getSanxotScriptsFolder() + File.separator + SANXOT_SIEVE_PY);
		}
		commandline.addArgument("-d" + dataFileName);
		commandline.addArgument("-p" + workingFolder.getAbsolutePath());
		commandline.addArgument("-r" + level1FileName);
		commandline.addArgument("-V" + infoFileName);
		commandline.addArgument("-a" + prefix);
		commandline.addArgument("-f" + quantParameters.getOutlierRemovalFDR().toString());
		return commandline;
	}

	private CommandLineRunner runCommand(CommandLine commandLine, long timeout)
			throws IOException, InterruptedException, ExecutionException {

		final CommandLineRunner runner = new CommandLineRunner();
		runner.runCommand(commandLine, timeout);
		final Long processExitCode = runner.getProcessExitCode();
		//
		// final String commandString = commandLine.toString();
		// log.info("Running: " + commandString);
		// logFileWriter.write(commandString + "\n");
		// firePropertyChange(STARTING_COMMAND, null, commandString);
		//
		// final ProcessExecutorHandler handler = new ProcessExecutorHandler() {
		//
		// @Override
		// public void onStandardOutput(String msg) {
		// log.debug("OUTPUT:" + msg);
		//
		// }
		//
		// @Override
		// public void onStandardError(String msg) {
		// log.error("ERROR:" + msg);
		//
		// }
		// };
		// final Future<Long> runProcess =
		// ProcessExecutor.runProcess(commandLine, handler, timeout);
		// while (!runProcess.isDone() && !runProcess.isCancelled()) {
		// Thread.sleep(1000);
		// }
		// final Long processExitCode = runProcess.get();
		log.info("Process exitValue: " + processExitCode);
		firePropertyChange(END_COMMAND, null, processExitCode);
		return runner;
	}

	/**
	 * @return the result
	 */
	public SanXotAnalysisResult getResult() {
		return result;
	}

	public static boolean checkAnyDifferentRelationShip(File relatFile) {
		if (relatFile == null) {
			// ignore the check
			return true;
		}
		final Map<String, Set<String>> relationShipsFromRelatFile = IntegrationResultWrapper
				.getRelationShipsFromRelatFile(relatFile);
		for (final String upperLevel : relationShipsFromRelatFile.keySet()) {
			if (relationShipsFromRelatFile.get(upperLevel).size() > 1) {
				return true;
			}
			if (relationShipsFromRelatFile.get(upperLevel).size() == 1) {
				if (!relationShipsFromRelatFile.get(upperLevel).iterator().next().equals(upperLevel)) {
					return true;
				}
			}
		}
		return false;
	}

}
