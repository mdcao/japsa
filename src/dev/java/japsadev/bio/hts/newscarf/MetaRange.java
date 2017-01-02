package japsadev.bio.hts.newscarf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * A meta range is composed of sub ranges. Its left endpoint is
 * the leftmost point of its subranges and its right endpoint is 
 * the rightmost point of its subranges. 
 */
public class MetaRange extends Range{

    /**
     * meta Ranges are composed of lists of subRanges.
     */
    private List<Range> subranges;

    public MetaRange(Range initialSubRange){
        super(initialSubRange.getLeft(),initialSubRange.getRight());
        this.subranges = new ArrayList<Range>();
        this.subranges.add(initialSubRange);
    }

    public void join(MetaRange other){
        verifyHomo(other);

        this.setLeft(Math.min(this.getLeft(),other.getLeft()));
        this.setRight(Math.max(this.getRight(),other.getRight()));
        this.subranges.addAll(other.getSubRanges());

        other.invalidate();
    }

    public void setSubRanges(List<Range> subRanges){
        this.subranges = subRanges;
    }

    public List<Range> getSubRanges(){
        return this.subranges;
    }

    private void invalidate(){
        this.setSubRanges(null);
        this.setLeft(0);
        this.setRight(0);
    }

    public boolean isInvalid(){
        return this.getSubRanges()==null;
    }

    /**
     * This code is retrieved from a similar one in StackOverFlow: credits to augray :). 
     * From here on, the term "overlap" actually means ranges that 
     * belong to the same homo-group 
     * 
     * Given a list of Ranges, returns a list of Range groups where
     * the Ranges in the groups all overlap. So if A overlaps with B and
     * B with C, then A,B,and C will be returned in a group. Supposing Ranges
     * D  and E share nothing with A, B, or C, but do share with each other, D and E
     * will be returned as a separate group from A,B,C.
     * @param baseRanges
     * @return
     */
    public static List<List<Range>> getOverlappingGroups(List<Range> baseRanges){
        List<MetaRange> baseMetaRanges = toMetaRanges(baseRanges);

        List<MetaRange> mergedMetaRanges = getMergedMetaRanges(baseMetaRanges);

        List<List<Range>> RangeGroups = metaRangesToGroups(mergedMetaRanges);
        return RangeGroups;
    }



    private static List<MetaRange> getMergedMetaRanges(
            List<MetaRange> metaRanges) {
        if(metaRanges.isEmpty()){
            return metaRanges;
        }
        //order the MetaRanges by their starting point.
        Collections.sort(metaRanges);

        //go through and merge the overlapping meta Ranges.
        //This relies on the logic that if Range i overlaps with
        //an Range that started before it, then once all the Ranges
        //before i have been merged, Range i will have a starting point
        //consecutive to the starting point of the the preceeding Range.
        for(int i=0; i< metaRanges.size()-1; i++){
            MetaRange thisRange = metaRanges.get(i);
            MetaRange nextRange = metaRanges.get(i+1);

            if(thisRange.isHomo(nextRange)){
                nextRange.join(thisRange);
            }
        }

        List<MetaRange> resultRanges = new ArrayList<MetaRange>();

        //All Ranges from the original list either:
        //(a) didn't overlap with any others
        //(b) overlapped with others and were chosen to represent the merged group or
        //(c) overlapped with others, were represented in the group in another 
        // MetaRange, and then marked as invalid.
        //Go through and only add the valid Ranges to be returned.

        for(MetaRange i : metaRanges){
            if(!i.isInvalid()){
                resultRanges.add(i);
            }
        }
        return resultRanges;
    }

    /**
     * Convert a list of MetaRanges into groups of Ranges.
     * @param mergedMetaRanges
     * @return
     */
    private static List<List<Range>> metaRangesToGroups(
            List<MetaRange> mergedMetaRanges) {
        List<List<Range>> groups = new ArrayList<>(mergedMetaRanges.size());
        for(MetaRange metaRange : mergedMetaRanges){
            groups.add(metaRange.getSubRanges());
        }
        return groups;
    }

    private static List<MetaRange> toMetaRanges(
            List<Range> baseRanges) {
        ArrayList<MetaRange> metaRanges = new ArrayList<MetaRange>(baseRanges.size());
        for(Range i : baseRanges){
            metaRanges.add(new MetaRange(i));
        }

        return metaRanges;
    }

}
