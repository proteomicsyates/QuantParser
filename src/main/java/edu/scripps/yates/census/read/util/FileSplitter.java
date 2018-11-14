package edu.scripps.yates.census.read.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Map;

import org.apache.commons.io.FilenameUtils;

import gnu.trove.map.hash.THashMap;

public class FileSplitter {
	/**
	 * Reads the relatFile and split the file into as many dataSetNames as
	 * exist, each one containing the data corresponding to each dataSet
	 *
	 * @param inputFile
	 * @param dataSetNames
	 * @return
	 * @throws IOException
	 */
	public static Map<String, File> splitFiles(File inputFile, Collection<String> dataSetNames) throws IOException {

		// if (dataSetNames.size() == 1) {
		// Map<String, File> ret = new THashMap<String, File>();
		// ret.put(dataSetNames.iterator().next(), inputFile);
		// return ret;
		// }

		final Map<String, BufferedWriter> mapOfFiles = new THashMap<String, BufferedWriter>();
		final Map<String, File> listOfFiles = new THashMap<String, File>();
		// create an output file per each dataSetName
		for (final String dataSetName : dataSetNames) {
			final File file = new File(inputFile.getParentFile().getAbsoluteFile() + File.separator
					+ FilenameUtils.getBaseName(inputFile.getAbsolutePath()) + "_" + dataSetName + "."
					+ FilenameUtils.getExtension(inputFile.getAbsolutePath()));
			listOfFiles.put(dataSetName, file);
			final FileWriter fstream = new FileWriter(file);
			final BufferedWriter out = new BufferedWriter(fstream);
			mapOfFiles.put(dataSetName, out);
		}

		FileInputStream fis;
		try {
			fis = new FileInputStream(inputFile);
			final BufferedReader in = new BufferedReader(new InputStreamReader(fis));

			String aLine;
			String firstLine = null;
			while ((aLine = in.readLine()) != null) {
				if (firstLine == null) {
					firstLine = aLine;
					// write the first line in all the writters
					final Collection<BufferedWriter> outs = mapOfFiles.values();
					for (final BufferedWriter out : outs) {
						out.write(firstLine);
						out.newLine();
					}
					continue;
				}
				final String[] split = aLine.split("\t");
				for (final String dataSetName : dataSetNames) {
					if (split[0].contains(dataSetName) || split[1].contains(dataSetName)) {
						// get the file writter corresponding to the dataset
						// detected as sufix
						final BufferedWriter out = mapOfFiles.get(dataSetName);
						out.write(aLine);
						out.newLine();
					}
				}
			}

			in.close();
			// close all file writters
			for (final BufferedWriter out : mapOfFiles.values()) {
				out.close();
			}
		} catch (final IOException e) {
			e.printStackTrace();
		}
		return listOfFiles;
	}

}
