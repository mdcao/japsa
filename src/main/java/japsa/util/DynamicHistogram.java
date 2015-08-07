/*****************************************************************************
 * Copyright (c) Minh Duc Cao, Monash Uni & UQ, All rights reserved.         *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  * 
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 * 3. Neither the names of the institutions nor the names of the contributors*
 *    may be used to endorse or promote products derived from this software  *
 *    without specific prior written permission.                             *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   *
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, *
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    *
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 ****************************************************************************/

/**************************     REVISION HISTORY    **************************
 * 20 Apr 2015 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/
package japsa.util;

/**
 * @author minhduc
 * Adapt from Histogram
 *
 */


import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.util.ParamChecks;
import org.jfree.data.general.DatasetChangeEvent;
import org.jfree.data.statistics.HistogramBin;
import org.jfree.data.statistics.HistogramType;
import org.jfree.data.xy.AbstractIntervalXYDataset;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.util.ObjectUtilities;
import org.jfree.util.PublicCloneable;

/**
 * A dataset that can be used for creating histograms.
 *
 * @see SimpleHistogramDataset
 */
public class DynamicHistogram extends AbstractIntervalXYDataset
        implements IntervalXYDataset, Cloneable, PublicCloneable,
                   Serializable {

    /** For serialization. */
    private static final long serialVersionUID = -6341668077370231153L;

    /** A list of maps. */
    //private List seriesList;
    
    private ArrayList<MySeries>  mySeriesList;
    

    /** The histogram type. */
    private HistogramType type;

    /**
     * Creates a new (empty) dataset with a default type of
     * {@link HistogramType}.FREQUENCY.
     */
    public DynamicHistogram() {
        //this.seriesList = new ArrayList();
        
        mySeriesList = new ArrayList<MySeries>();
        this.type = HistogramType.FREQUENCY;
    }

    /**
     * Returns the histogram type.
     *
     * @return The type (never <code>null</code>).
     */
    public HistogramType getType() {
        return this.type;
    }

    /**
     * Sets the histogram type and sends a {@link DatasetChangeEvent} to all
     * registered listeners.
     *
     * @param type  the type (<code>null</code> not permitted).
     */
    public void setType(HistogramType type) {
        ParamChecks.nullNotPermitted(type, "type");
        this.type = type;
        fireDatasetChanged();
    }

      
    public void prepareSeries(String key, int bins,
            double minimum, double maximum) {
    	//FIXME: Make sure the key is not in
    	ParamChecks.nullNotPermitted(key, "key");
    	
    	 double lower = minimum;
         double upper;
    	
    	if (bins < 1) {
            throw new IllegalArgumentException(
                    "The 'bins' value must be at least 1.");
        }
    	
    	MySeries newSeries = new MySeries();    	
        double binWidth = (maximum - minimum) / bins;
        newSeries.seriesID = key;
        newSeries.binWidth = binWidth;
        newSeries.maximum = maximum;
        newSeries.minimum = minimum;
        
        ArrayList<HistogramBin>  binList = new ArrayList<HistogramBin>(bins);
        for (int i = 0; i < bins; i++) {
            HistogramBin bin;
            // make sure bins[bins.length]'s upper boundary ends at maximum
            // to avoid the rounding issue. the bins[0] lower boundary is
            // guaranteed start from min
            if (i == bins - 1) {
                bin = new HistogramBin(lower, maximum);
            }
            else {
                upper = minimum + (i + 1) * binWidth;
                bin = new HistogramBin(lower, upper);
                lower = upper;
            }
            binList.add(bin);
        }
        
        newSeries.binList = binList;
        mySeriesList.add(newSeries);
    }
    
    public void notifyChanged(){
    	this.fireDatasetChanged();
    }

        
    public int getSeriesIndex(String serieKey){    	
    	for (int index = 0;index < mySeriesList.size();index ++){
    		if (serieKey.equals(mySeriesList.get(index).seriesID))
    				return index;
    	}
    	return -1;
    }
    
    public void addSeries(int serieIndex, double value){
    	MySeries series = mySeriesList.get(serieIndex);    	
    	 int binIndex = series.binList.size() - 1;
    	 int bins = series.binList.size();
    	 
         if (value < series.maximum) {
             double fraction = (value - series.minimum) / (series.maximum - series.minimum);
             if (fraction < 0.0) {
                 fraction = 0.0;
             }
             binIndex = (int) (fraction * bins);
             // rounding could result in binIndex being equal to bins
             // which will cause an IndexOutOfBoundsException - see bug
             // report 1553088
             if (binIndex >= bins) {
                 binIndex = bins - 1;
             }
         }
         HistogramBin bin = series.binList.get(binIndex);
         bin.incrementCount();
         series.totalValues ++;
    }
    		
   
  
    /**
     * Returns the bins for a series.
     *
     * @param series  the series index (in the range <code>0</code> to
     *     <code>getSeriesCount() - 1</code>).
     *
     * @return A list of bins.
     *
     * @throws IndexOutOfBoundsException if <code>series</code> is outside the
     *     specified range.
     */
    List<HistogramBin> getBins(int series) {
    	return this.mySeriesList.get(series).binList;
    }

    /**
     * Returns the total number of observations for a series.
     *
     * @param series  the series index.
     *
     * @return The total.
     */
    private int getTotal(int series) {
    	return this.mySeriesList.get(series).totalValues;
    }

    /**
     * Returns the bin width for a series.
     *
     * @param series  the series index (zero based).
     *
     * @return The bin width.
     */
    private double getBinWidth(int series) {
    	return this.mySeriesList.get(series).binWidth;
    }

    /**
     * Returns the number of series in the dataset.
     *
     * @return The series count.
     */
    @Override
    public int getSeriesCount() {
        return this.mySeriesList.size();
    }

    /**
     * Returns the key for a series.
     *
     * @param series  the series index (in the range <code>0</code> to
     *     <code>getSeriesCount() - 1</code>).
     *
     * @return The series key.
     *
     * @throws IndexOutOfBoundsException if <code>series</code> is outside the
     *     specified range.
     */
    @Override
    public Comparable getSeriesKey(int series) {
    	return this.mySeriesList.get(series).seriesID;
    }

    /**
     * Returns the number of data items for a series.
     *
     * @param series  the series index (in the range <code>0</code> to
     *     <code>getSeriesCount() - 1</code>).
     *
     * @return The item count.
     *
     * @throws IndexOutOfBoundsException if <code>series</code> is outside the
     *     specified range.
     */
    @Override
    public int getItemCount(int series) {
        return getBins(series).size();
    }

    /**
     * Returns the X value for a bin.  This value won't be used for plotting
     * histograms, since the renderer will ignore it.  But other renderers can
     * use it (for example, you could use the dataset to create a line
     * chart).
     *
     * @param series  the series index (in the range <code>0</code> to
     *     <code>getSeriesCount() - 1</code>).
     * @param item  the item index (zero based).
     *
     * @return The start value.
     *
     * @throws IndexOutOfBoundsException if <code>series</code> is outside the
     *     specified range.
     */
    @Override
    public Number getX(int series, int item) {
        List<HistogramBin> bins = getBins(series);
        HistogramBin bin =  bins.get(item);
        double x = (bin.getStartBoundary() + bin.getEndBoundary()) / 2.;
        return new Double(x);
    }

    /**
     * Returns the y-value for a bin (calculated to take into account the
     * histogram type).
     *
     * @param series  the series index (in the range <code>0</code> to
     *     <code>getSeriesCount() - 1</code>).
     * @param item  the item index (zero based).
     *
     * @return The y-value.
     *
     * @throws IndexOutOfBoundsException if <code>series</code> is outside the
     *     specified range.
     */
    @Override
    public Number getY(int series, int item) {
        List<HistogramBin> bins = getBins(series);
        HistogramBin bin =  bins.get(item);
        double total = getTotal(series);
        double binWidth = getBinWidth(series);

        if (this.type == HistogramType.FREQUENCY) {
            return new Double(bin.getCount());
        }
        else if (this.type == HistogramType.RELATIVE_FREQUENCY) {
            return new Double(bin.getCount() / total);
        }
        else if (this.type == HistogramType.SCALE_AREA_TO_1) {
            return new Double(bin.getCount() / (binWidth * total));
        }
        else { // pretty sure this shouldn't ever happen
            throw new IllegalStateException();
        }
    }

    /**
     * Returns the start value for a bin.
     *
     * @param series  the series index (in the range <code>0</code> to
     *     <code>getSeriesCount() - 1</code>).
     * @param item  the item index (zero based).
     *
     * @return The start value.
     *
     * @throws IndexOutOfBoundsException if <code>series</code> is outside the
     *     specified range.
     */
    @Override
    public Number getStartX(int series, int item) {
        List<HistogramBin> bins = getBins(series);
        HistogramBin bin = bins.get(item);
        return new Double(bin.getStartBoundary());
    }

    /**
     * Returns the end value for a bin.
     *
     * @param series  the series index (in the range <code>0</code> to
     *     <code>getSeriesCount() - 1</code>).
     * @param item  the item index (zero based).
     *
     * @return The end value.
     *
     * @throws IndexOutOfBoundsException if <code>series</code> is outside the
     *     specified range.
     */
    @Override
    public Number getEndX(int series, int item) {
        List<HistogramBin> bins = getBins(series);
        HistogramBin bin = bins.get(item);
        return new Double(bin.getEndBoundary());
    }

    /**
     * Returns the start y-value for a bin (which is the same as the y-value,
     * this method exists only to support the general form of the
     * {@link IntervalXYDataset} interface).
     *
     * @param series  the series index (in the range <code>0</code> to
     *     <code>getSeriesCount() - 1</code>).
     * @param item  the item index (zero based).
     *
     * @return The y-value.
     *
     * @throws IndexOutOfBoundsException if <code>series</code> is outside the
     *     specified range.
     */
    @Override
    public Number getStartY(int series, int item) {
        return getY(series, item);
    }

    /**
     * Returns the end y-value for a bin (which is the same as the y-value,
     * this method exists only to support the general form of the
     * {@link IntervalXYDataset} interface).
     *
     * @param series  the series index (in the range <code>0</code> to
     *     <code>getSeriesCount() - 1</code>).
     * @param item  the item index (zero based).
     *
     * @return The Y value.
     *
     * @throws IndexOutOfBoundsException if <code>series</code> is outside the
     *     specified range.
     */
    @Override
    public Number getEndY(int series, int item) {
        return getY(series, item);
    }

    /**
     * Tests this dataset for equality with an arbitrary object.
     *
     * @param obj  the object to test against (<code>null</code> permitted).
     *
     * @return A boolean.
     */
    @Override
    public boolean equals(Object obj) {
        if (obj == this) {
            return true;
        }
        if (!(obj instanceof DynamicHistogram)) {
            return false;
        }
        DynamicHistogram that = (DynamicHistogram) obj;
        if (!ObjectUtilities.equal(this.type, that.type)) {
            return false;
        }
        if (!ObjectUtilities.equal(this.mySeriesList, that.mySeriesList)) {
            return false;
        }
        return true;
    }

    /**
     * Returns a clone of the dataset.
     *
     * @return A clone of the dataset.
     *
     * @throws CloneNotSupportedException if the object cannot be cloned.
     */
    @Override
    public Object clone() throws CloneNotSupportedException {
        //DynamicHistogram clone = (DynamicHistogram) super.clone();
        //int seriesCount = getSeriesCount();
        //clone.seriesList = new java.util.ArrayList(seriesCount);
        //for (int i = 0; i < seriesCount; i++) {
         //   clone.seriesList.add(new HashMap((Map) this.seriesList.get(i)));
        //}
        return null;
    }
    
    public static class MySeries{
    	String seriesID = "";
    	ArrayList<HistogramBin> binList;
    	int totalValues = 0;//number of values
    	double binWidth;
    	double minimum, maximum;
    }
}
