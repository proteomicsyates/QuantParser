package edu.scripps.yates.census.read.model.interfaces;

import java.util.Map;

import gnu.trove.map.hash.THashMap;

public class StaticItemStorage<T extends HasKey> {
	private final Map<String, T> map = new THashMap<String, T>();

	public boolean contains(T hasKeyObj) {
		return map.containsKey(hasKeyObj.getKey());
	}

	public boolean containsKey(String key) {
		return map.containsKey(key);
	}

	public T addItem(T hasKeyObj) {

		return map.put(hasKeyObj.getKey(), hasKeyObj);

	}

	public T getItem(String key) {
		return map.get(key);
	}

	public int size() {
		return map.size();
	}

	public void clear() {
		map.clear();
	}

	public T remove(T hasKeyObj) {
		return map.remove(hasKeyObj.getKey());
	}

	public T remove(String key) {
		return map.remove(key);

	}

	public boolean isEmpty() {
		return map.isEmpty();
	}
}
