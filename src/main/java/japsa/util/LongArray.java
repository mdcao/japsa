/*
 * Copyright (c) 1998 - 2005 Versant Corporation
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 * Versant Corporation - initial API and implementation
 */

/**
 * Growable int[]. This is based com.sosnoski.util.array.IntArray from
 * Sosnoski Software Solutions, Inc.
 */


/**************************     REVISION HISTORY    **************************
 * 14/03/2014 - Minh Duc Cao: added                                        
 *  
 ****************************************************************************/

package japsa.util;
/**
 * A resizeable array of long. This array bypass safe checks
 * @author minhduc
 *
 */
public final class LongArray {

	private long[] buf;
	private int size;

	public LongArray() {
		this(64);
	}

	public LongArray(int capacity) {
		buf = new long[capacity];
	}

	public int size() {
		return size;
	}

	private void ensureCapacity(int len) {
		if (size + len > buf.length) {
			int n = buf.length * 3 / 2 + 1;
			if (size + len > n) {
				n = size + len;
			}
			long[] a = new long[n];
			System.arraycopy(buf, 0, a, 0, size);
			buf = a;
		}
	}

	public void add(long v) {
		ensureCapacity(size + 1);
		buf[size++] = v;
	}

	/**
	 * Add a value at a specified index in the array.
	 */
	public void add(int index, long value) {
		ensureCapacity(size + 1);
		if (index == size) {
			buf[size++] = value;
		} else {
			System.arraycopy(buf, index, buf, index + 1, size - index);
			buf[index] = value;
		}
	}

	/**
	 * Constructs and returns a simple array containing the same data as held
	 * in this growable array.
	 */
	public long[] toArray() {
		long[] a = new long[size];
		System.arraycopy(buf, 0, a, 0, size);
		return a;
	}

	public void clear() {
		size = 0;
	}

	/**
	 * Retrieve the value present at an index position in the array.
	 */
	public long get(int index) {
		return buf[index];
	}

}

