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
 * A resizeable array of double. This array bypass safe checks
 * @author minhduc
 *
 */
public final class DoubleArray {

    private double[] buf;
    private int size;

    public DoubleArray() {
        this(64);
    }

    public DoubleArray(int capacity) {
        buf = new double[capacity];
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
            double[] a = new double[n];
            System.arraycopy(buf, 0, a, 0, size);
            buf = a;
        }
    }

    public void add(double v) {
        ensureCapacity(size + 1);
        buf[size++] = v;
    }

    /**
     * Add a value at a specified index in the array.
     */
    public void add(int index, double value) {
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
    public double[] toArray() {
    	double[] a = new double[size];
        System.arraycopy(buf, 0, a, 0, size);
        return a;
    }

    public void clear() {
        size = 0;
    }

    /**
     * Retrieve the value present at an index position in the array.
     */
    public double get(int index) {
        return buf[index];
    }

}

   