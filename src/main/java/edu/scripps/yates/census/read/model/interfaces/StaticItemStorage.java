package edu.scripps.yates.census.read.model.interfaces;

import java.util.Map;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock.ReadLock;
import java.util.concurrent.locks.ReentrantReadWriteLock.WriteLock;

import edu.scripps.yates.utilities.proteomicsmodel.HasKey;
import gnu.trove.map.hash.THashMap;

public class StaticItemStorage<T extends HasKey> {
	private final Map<String, T> map = new THashMap<String, T>();
	private final ReentrantReadWriteLock lock = new ReentrantReadWriteLock();

	public boolean contains(T hasKeyObj) {
		ReadLock readLock = lock.readLock();
		try {
			readLock.lock();
			return map.containsKey(hasKeyObj.getKey());
		} finally {
			readLock.unlock();
		}
	}

	public boolean containsKey(String key) {
		ReadLock readLock = lock.readLock();
		try {
			readLock.lock();
			return map.containsKey(key);
		} finally {
			readLock.unlock();
		}
	}

	public T addItem(T hasKeyObj) {
		WriteLock writeLock = lock.writeLock();
		try {
			writeLock.lock();
			return map.put(hasKeyObj.getKey(), hasKeyObj);
		} finally {
			writeLock.unlock();
		}
	}

	public T getItem(String key) {
		ReadLock readLock = lock.readLock();
		try {
			readLock.lock();
			return map.get(key);
		} finally {
			readLock.unlock();
		}
	}

	public int size() {
		ReadLock readLock = lock.readLock();
		try {
			readLock.lock();
			return map.size();
		} finally {
			readLock.unlock();
		}
	}

	public void clear() {
		WriteLock writeLock = lock.writeLock();
		try {
			writeLock.lock();
			map.clear();
		} finally {
			writeLock.unlock();
		}
	}

	public T remove(T hasKeyObj) {
		WriteLock writeLock = lock.writeLock();
		try {
			writeLock.lock();
			return map.remove(hasKeyObj.getKey());
		} finally {
			writeLock.unlock();
		}
	}

	public T remove(String key) {
		WriteLock writeLock = lock.writeLock();
		try {
			writeLock.lock();
			return map.remove(key);
		} finally {
			writeLock.unlock();
		}

	}

	public boolean isEmpty() {
		ReadLock readLock = lock.readLock();
		try {
			readLock.lock();
			return map.isEmpty();
		} finally {
			readLock.unlock();
		}
	}
}
