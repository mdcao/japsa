
package test;

import japsa.seq.SequenceOutputStream;
import japsa.seq.nanopore.NanoporeReader;

public class TestNanopore {
	
	public static void main(String[] args) throws Exception{
		//Create a buffered stream from the 
        SequenceOutputStream sos = SequenceOutputStream.makeOutputStream("-");
        NanoporeReader reader = new NanoporeReader(args[0]);        
        
        if (reader.getSeqTemplate() != null)
        	reader.getSeqTemplate().print(sos);
        
        if (reader.getSeqComplement() != null)
        	reader.getSeqComplement().print(sos);
                
        if (reader.getSeq2D() != null)
        	reader.getSeq2D().print(sos);
        
        sos.flush();
	}
}
