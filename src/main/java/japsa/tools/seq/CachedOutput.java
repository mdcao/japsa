package japsa.tools.seq;

import htsjdk.samtools.SAMRecord;

public interface CachedOutput {
	 public void write(SAMRecord sam);
	 public void close();
}
