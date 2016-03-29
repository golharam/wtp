package com.lifetechnologies.solid.wt.config;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.text.ParseException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.PropertyResourceBundle;
import java.util.ResourceBundle;
import java.util.Set;

import com.lifetechnologies.solid.wt.Utilities;
import com.lifetechnologies.solid.wt.WholeTranscriptomeAnalyzer;
import com.lifetechnologies.util.TextUtils;

/**
 * A set of utils for managing WT configuration.
 * WT config is represented as a Map of String to Set of Strings.
 * Configs can be specified in the global location, user's location and supplied config file.
 * @author mullermw
 *
 */
public class WTConfigUtils {

	public static final String CONFIG_FILE_NAME = "config.ini";
	
	/**
	 * Construct a config map from the global, user's home and the supplied file.
	 * @param map to store the config.
	 * @param f client specified config file.
	 * @throws IOException if an error occurs in reading the config sources.
	 * @throws ParseException if one of the config sources cannot be parsed.
	 */
	public static void mapConfig(Map<ConfigKey, Set<String>> map, File f) throws IOException, ParseException {
		mapGlobalConfig(map);
		mapHomeConfig(map);
		mapConfigFile(f, map);
	}
	
	/**
	 * Convenience method for mapping config from global config and user's home config.
	 * @param map
	 * @throws IOException
	 * @throws ParseException
	 */
	public static void mapConfig(Map<ConfigKey, Set<String>> map) throws IOException, ParseException {
		mapConfig(map, null);
	}
	
	/**
	 * @return a new empty config.
	 */
	public static Map<ConfigKey, Set<String>> newConfig() {
		return new HashMap<ConfigKey, Set<String>>();
	}
	
	/**
	 * Parse the global config file located in $WT_HOME/etc/
	 * @param map to store the config.
	 * @throws IOException
	 * @throws ParseException
	 */
	public static void mapGlobalConfig(Map<ConfigKey, Set<String>> map) throws IOException, ParseException {
		if (map == null) return;
		String wtHome = System.getProperty(WholeTranscriptomeAnalyzer.WT_HOME_SYSTEM_PROPERTY);
		if (Utilities.isBlank(wtHome)) return;
		File file = new File(wtHome);
		file = new File(file, "etc");
		file = new File(file, CONFIG_FILE_NAME);
		if (!file.exists()) return;
		if (!file.isFile()) return;
		if (!file.canRead()) return;
		mapConfigFile(file, map);
	}
	
	/**
	 * Parse the user's home config, located in $user.home/.ab_wtp/
	 * @param map
	 * @throws IOException
	 * @throws ParseException
	 */
	public static void mapHomeConfig(Map<ConfigKey, Set<String>> map) throws IOException, ParseException {
		if (map == null) return;
		if (map == null) return;
		String userHome = System.getProperty("user.home");
		if (Utilities.isBlank(userHome)) return;
		File file = new File(userHome);
		file = new File(file, ".ab_wtp");
		file = new File(file, CONFIG_FILE_NAME);
		if (!file.exists()) return;
		if (!file.isFile()) return;
		if (!file.canRead()) return;
		mapConfigFile(file, map);
	}
	
	/**
	 * Parse the supplied config file.
	 * @param f the file to parse
	 * @param map to hold the config.
	 * @throws IOException
	 * @throws ParseException
	 */
	private static void mapConfigFile(File f, Map<ConfigKey, Set<String>> map) throws IOException, ParseException {
		if (f == null) return;
		if (!f.exists()) return;
		if (!f.isFile()) return;
		if (!f.canRead()) return;
		InputStream stream = new FileInputStream(f);
		try {
			mapConfigStream(stream, map);
		}catch (ParseException e) {
			throw new ParseException("Failed to parse config file: " + f + " at line " + e.getErrorOffset()+"\n"+e.getMessage(), e.getErrorOffset());
		} finally {
			stream.close();
		}
	}
	
	/**
	 * Parse the supplied config stream.
	 * @param s the stream to parse.
	 * @param map to hold the config.
	 * @throws IOException
	 * @throws ParseException
	 */
	private static void mapConfigStream(InputStream s, Map<ConfigKey, Set<String>> map) throws IOException, ParseException {
		Map<ConfigKey, Set<String>> myMap = new HashMap<ConfigKey, Set<String>>();
		Properties properties = new Properties();
		properties.load(s);
		for (Object obj : properties.keySet()) {
			String qname = (String)obj;
			String value = properties.getProperty(qname).trim();
			List<String> values = TextUtils.splitCSV(value);
			ConfigKey key = ConfigKey.getKeyByQName(qname);
			if (key == null) continue;
			myMap.put(key, new HashSet<String>(values));
		}
		
		for (Map.Entry<ConfigKey, Set<String>> entry : myMap.entrySet())
			map.put(entry.getKey(), entry.getValue());
	}
	
	public static ResourceBundle getApplicationProperties() {
		return (PropertyResourceBundle) ResourceBundle.getBundle(WholeTranscriptomeAnalyzer.class.getPackage().getName()+".application");
	}
	
	public static Set<String> asSingleton(String s) {
		Set<String> set = new HashSet<String>();
		set.add(s);
		return set;
	}
	
}
