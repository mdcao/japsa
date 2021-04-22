package japsa.bio.np;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import japsa.bio.np.RealtimeSpeciesTyping.Coverage;
import japsa.tools.seq.SequenceUtils;
import japsa.tools.seq.SparseVector;
import japsa.tools.seq.SparseVectorCollection;

public  class AllRecords{
	public static String tag = null;//	String tag = "BS"; //AS PI
		String readnme = null;
		List<SAMRecord> records = new ArrayList<SAMRecord>();
		List<String> refs = new ArrayList<String>();
		List<Integer> species = new ArrayList<Integer>(); //specesIndex
		SparseVector all_species = new SparseVector();
		//SparseVector all_speciesLen = new SparseVector(); // for capturing number of bases

		//List<Integer> quality = new ArrayList<Integer>();
		public void clear(){
			readnme=null;
			records.clear();
			refs.clear();
			all_species.clear();
			species.clear();
		//	quality.clear();
		}
		
		//return 
		public int  getAll(Integer species2,List<SAMRecord> out) {
			int q = -1;
			int best =-1;
			for(int i=0; i<records.size(); i++){
				if(species.get(i).equals(species2)){
					int q1 = records.get(i).getMappingQuality();
					if(q1>q){
						best = out.size();
						q = q1;
					}
					out.add(records.get(i));
				}
			}
			return best;
			//return out;
		}
		
		public void add(SAMRecord sam, int spec){
			boolean secondary = sam.isSecondaryOrSupplementary();
			if(!SequenceUtils.secondary  && secondary){
				//only include supplementary alignments to same species
				if(species.size()==0 || species.get(0)!=spec){
					return;
				}
			}
			if(readnme==null){
				//if(sam.isSecondaryOrSupplementary()){
				//	throw new RuntimeException("@!!");
				//}
				readnme=sam.getReadName();
			}
			else if(!readnme.equals(sam.getReadName())) {
				throw new RuntimeException("!!");
			}
			this.records.add(sam);
			this.refs.add(sam.getReferenceName());
			this.species.add(spec);
			this.all_species.addToEntry(spec, 0); // set as placemarker
		
			
//			this.quality.add(sam.getMappingQuality());
		}
		public int size() {
			// TODO Auto-generated method stub
			return records.size();
		}
		public void transferReads(List<Coverage>species2ReadList ,SparseVectorCollection all_reads) {
			if(this.size()>0) {
				;
			//	if(this.all_species.size()>1){
				//	System.err.println("h");
				//}
				Iterator<Integer> specs = this.all_species.keySet().iterator();
				List<SAMRecord> sams= new ArrayList<SAMRecord>();
				List<SAMRecord> filtered  = new ArrayList<SAMRecord>();
				while(specs.hasNext()){
					Integer spec = specs.next();
				//	sv.addToEntry(spec,all_species.get(spec));
				//			-1*Math.pow(10, -1*all_species.get(spec).doubleValue()));
					Coverage coverage = species2ReadList.get(spec);
					
					
					
					
					sams.clear();filtered.clear();
					int besti = getAll(spec,sams);
					SAMRecord sam = sams.get(besti);
					int q = sam.getMappingQuality();
					
					boolean primary = !sam.isSecondaryOrSupplementary();
					RealtimeSpeciesTyping.filter(sams, filtered, besti);
					if(primary) coverage.addRead(filtered);
					// we use 0.1 as minimum to avoid zero probability for mq=0 reads
					double v; 
				
				//	List<SAMTagAndValue> attr = sams.get(besti).getAttributes();
					if(tag!=null ){
						Number attr1 = (Number) sam.getAttribute(tag);
						if(attr1==null) throw new RuntimeException ("tag not available");
						v = attr1.doubleValue();///100.0;
					}else{
						v= Math.max(0.00001, 1-Math.pow(RealtimeSpeciesTyping.base_,-1.0*(double)q));
					}
				/*	for(int i=0; i<attr.size(); i++){
						System.err.println(attr.get(i).tag+" "+attr.get(i).value);
					}*/
					
					
					this.all_species.update(spec,v);
				//	this.all_speciesLen.update(spec,getBases(filtered));
					
				}
				if(all_reads !=null) all_reads.add(this.all_species);

				}
			all_species = new SparseVector();
			this.clear();
		}

		
	}