/*
 * Copyright (c) 2017  Minh Duc Cao (minhduc.cao@gmail.com).
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the names of the institutions nor the names of the contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

package japsadev.bio.np.phage;

/**
 * Created by minhduc on 22/04/17.
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;


import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.HTSUtilities;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class VectorSequenceExtraction {

    private static final Logger LOG = LoggerFactory.getLogger(VectorSequenceExtraction.class);

    static final int FLANKING=100;
    boolean trim=false;
    //Sequence plasmid;
    String plasmidFile; //plasmid fasta/q file that are already indexed by bwa
    int s5, e5,
            s3, e3;
    String bwaExe = "bwa";

    Process bwaProcess = null;

    public VectorSequenceExtraction(String seqFile, String bwaExe, boolean trim, int e5, int s3) throws IOException{
        plasmidFile=seqFile;
        this.e5=e5;
        this.s5=this.e5-FLANKING;
        this.s3=s3;
        this.e3=this.s3+FLANKING;
        this.bwaExe = bwaExe;
        this.trim = trim;

        SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
    }


    private SamReader getSamStream(String inFile, String format, int bwaThread) throws IOException, InterruptedException{
        SamReader reader = null;

        if (format.endsWith("am")){//bam or sam
            if ("-".equals(inFile))
                reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
            else
                reader = SamReaderFactory.makeDefault().open(new File(inFile));
        }else{
            LOG.info("Starting bwa  at " + new Date());

            ProcessBuilder pb = null;
            if ("-".equals(inFile)){
                pb = new ProcessBuilder(bwaExe,
                        "mem",
                        "-t",
                        "" + bwaThread,
                        "-k11",
                        "-W20",
                        "-r10",
                        "-A1",
                        "-B1",
                        "-O1",
                        "-E1",
                        "-L0",
                        "-a",
                        "-Y",
                        "-K",
                        "20000",
                        plasmidFile,
                        "-"
                ).redirectInput(ProcessBuilder.Redirect.INHERIT);
            }else{
                pb = new ProcessBuilder(bwaExe,
                        "mem",
                        "-t",
                        "" + bwaThread,
                        "-k11",
                        "-W20",
                        "-r10",
                        "-A1",
                        "-B1",
                        "-O1",
                        "-E1",
                        "-L0",
                        "-a",
                        "-Y",
                        "-K",
                        "20000",
                        plasmidFile,
                        inFile
                );
            }
            bwaProcess  = pb.redirectError(ProcessBuilder.Redirect.to(new File("/dev/null"))).start();
            reader = SamReaderFactory.makeDefault().open(SamInputResource.of(bwaProcess.getInputStream()));
        }
        return reader;
    }

    public void extractInsertSequence(String inFile, int qual, String format, int bwaThread, String output) throws IOException, InterruptedException{
        int [] refP5 = {s5,e5};
        int [] refP3 = {s3,e3};

        SamReader reader = getSamStream(inFile, format, bwaThread);
        SAMRecordIterator iter = reader.iterator();

        String currentReadName = "";

        SequenceOutputStream outFile = SequenceOutputStream.makeOutputStream(output);
        int count=0;
        SAMRecord currentRecord= null, prevRecord = null;

        boolean extracted = false;
        boolean firstDirection = false;
        int startNeg = 0,  endNeg = 0, startPos = 0, endPos = 0;
        String sequenceStr = "";
        while (iter.hasNext()) {
            currentRecord = iter.next();

            if (currentRecord.getReadUnmappedFlag())
                continue;
            if (currentRecord.getMappingQuality() <= qual)
                continue;
            count++;

            String currentRecordName = currentRecord.getReadName();

            if (!currentReadName.equals(currentRecordName)) {
                //start of a new record
                extracted = false;
                startNeg = endNeg = startPos = endPos = 0;
                currentReadName = currentRecordName;
                firstDirection = currentRecord.getReadNegativeStrandFlag();
                sequenceStr = currentRecord.getReadString();
            }else if (extracted){
                continue;//I have extracted from this read
            }else if (firstDirection != currentRecord.getReadNegativeStrandFlag()) {
                LOG.warn("Read " + currentReadName + " not match first direction " + firstDirection);
                continue;
            }

            //Neg
            //if (currentRecord.getReadNegativeStrandFlag()) {
            //Check left
            if (currentRecord.getAlignmentStart() <= s5 && currentRecord.getAlignmentEnd() >= e5) {
                int[] pos = HTSUtilities.positionsInRead(currentRecord, refP5);
                //if(pos[0] > s5*0.8 && pos[0] < s5*1.2){
                if (pos[0] > 0) {
                    startNeg = pos[0];
                }
            }

            //Check right
            if (currentRecord.getAlignmentStart() <= s3 && currentRecord.getAlignmentEnd() >= e3) {
                int[] pos = HTSUtilities.positionsInRead(currentRecord, refP3);
                //if(pos[0] > s5*0.8 && pos[0] < s5*1.2){
                if (pos[1] > 0) {
                    endNeg = pos[1];
                }
            }

            //Check if both startNeg and endNeg were anchored
            if (endNeg > 0 && startNeg > 0) {
                if (endNeg < startNeg) {
                    LOG.warn("Read " + currentReadName + ": find end (" + endNeg + ") < start (" + startNeg + ") of " + firstDirection);
                    continue;
                }

                if (endNeg > sequenceStr.length()) {
                    LOG.warn("Read " + currentReadName + ": find end (" + endNeg + ") > length (" + sequenceStr.length() + ")");
                    continue;
                }

                String readSub = sequenceStr.substring(startNeg,endNeg);
                Sequence rs = new Sequence(Alphabet.DNA16(), readSub, currentReadName);
                rs.writeFasta(outFile);
                process(rs);

                extracted = true;
                //TODO: yay
            }
        }
        iter.close();
        outFile.close();
        reader.close();
        if (bwaProcess != null){
            bwaProcess.waitFor();
        }
    }

    //HashMap<String, Sequence> seqs = new HashMap<String, Sequence>();
    HashMap<String, Group> groups = new HashMap<String, Group>();

    static class Group implements Comparable<Group>{
        Sequence repSequence;
        int countSequence = 1;
        Group(Sequence seq){
            repSequence = seq;
        }

        @Override
        public int compareTo(Group o) {
            return Integer.compare(o.countSequence, countSequence);
        }
    }

    private void checkRatio(int currentCount) throws IOException {
        System.out.println("====================== " + insertCount + " =======================");
        ArrayList<Group> list = new ArrayList<Group>(groups.values());
        SequenceOutputStream out = SequenceOutputStream.makeOutputStream("top" + currentCount + ".japsa");
        Collections.sort(list);
        for (int i = 0; i < 20 && i < list.size(); i++) {
            Sequence seq = list.get(i).repSequence;
            System.out.println(seq.getName() + " " + list.get(i).countSequence);
            seq.writeFasta(out);
        }
        out.close();
    }
    int insertCount = 0;
    private void process(Sequence seq) throws IOException, InterruptedException {
        insertCount ++;
        if (insertCount % 100 == 0){
            checkRatio(insertCount);
        }
        if (groups.size() == 0)
            groups.put(seq.getName(), new Group(seq));
        else{

            String subjectFileName = "subject" + groups.size() + ".fasta";
            File subjectFile = new File(subjectFileName);
            if (!subjectFile.exists()){
                SequenceOutputStream sos = SequenceOutputStream.makeOutputStream("subject" + groups.size() + ".fasta");
                for (Group g:groups.values())
                    g.repSequence.writeFasta(sos);
                sos.close();
            }

            SequenceOutputStream sos = SequenceOutputStream.makeOutputStream("query.fasta");
            seq.writeFasta(sos);
            sos.close();
            //run blastn

            ProcessBuilder pb = new ProcessBuilder("blastn",
                    "-outfmt",
                    "6 qseqid qlen qstart qend sseqid slen sstart send length frames pident nident gaps mismatch score bitscore",
                    "-subject",
                    subjectFileName,
                    "-query",
                    "query.fasta"
                    );

            Process process  = pb.redirectError(ProcessBuilder.Redirect.to(new File("/dev/null"))).start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            String line = "";
            String match = null;
            while ( (line = reader.readLine()) != null){
                String [] toks = line.trim().split("\t");
                int qLen = Integer.parseInt(toks[1]);
                int qStart = Integer.parseInt(toks[2]);
                int qEnd = Integer.parseInt(toks[3]);
                if ((0.0 + qEnd - qStart)/qLen < 0.85)
                    continue;

                int sLen = Integer.parseInt(toks[5]);
                int sStart = Integer.parseInt(toks[6]);
                int sEnd = Integer.parseInt(toks[7]);
                if ((0.0 + sEnd - sStart)/sLen < 0.85)
                    continue;

                double identity = Double.parseDouble(toks[10]);
                if (identity < 0.85)
                    continue;
                //great
                match = toks[4];
                //break;
            }
            process.waitFor();
            if (match == null){
                groups.put(seq.getName(), new Group(seq));
                LOG.info("Add " + seq.getName() + "    " + groups.size() + " from " + insertCount);
            }else {
                groups.get(match).countSequence = groups.get(match).countSequence + 1;
                LOG.info("Pass " + seq.getName() + " -> " + match + " " +  groups.get(match).countSequence);
            }
        }
    }

}
