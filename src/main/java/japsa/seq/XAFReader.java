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

/*                           Revision History                                
 * 05/06/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package japsa.seq;

import japsa.util.Logging;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * eXtendable Annotation Format (xaf) is a format to store information at
 * certain positions along a genome or a group of sequences. It is a text format
 * where a record is stored in a row, and fields are separated by a separator (a
 * tab by default, but is changable). The format is inspired by Variant Call
 * Format where the fields (columns) and meta data can be encapsulated.
 * 
 * @author Minh Duc Cao
 */

public class XAFReader implements Closeable {
	String sep = "\t";
	public static final String XAF_HEADER = "#XAF";
	public static final String ROW_HEADER = "#H:";
	public static final String COMMENT = "##";
	public static final String INFORMATION_HEADER = "#I:";
	// #I:ID=NS,Type=int,Description="Number of Samples With Data"
	public static final String I_ID = "ID=", I_TYPE = "Type=",
			I_DESC = "Description=";

	/**
	 * Some default column headers
	 */
	public static final FieldHeader FH_ID = new FieldHeader("ID", "String",
			"ID of the record"), FH_CHROM = new FieldHeader("chrom", "String",
			"Chromosome the record"), FH_START = new FieldHeader("start",
			"int", "The start position (1-index, inclusive) of the record"),
			FH_END = new FieldHeader("end", "int",
					"The end position (1-index, inclusive) of the record");

	// These columns would have been defined as:
	// #I:ID=ID,Type=String,Description="ID of the record"
	// #I:ID=chrom,Type=String,Description="Chromosome the record"
	// #I:ID=start,Type=int,Description="The start position (1-index, inclusive) of the record"
	// #I:ID=end,Type=int,Description="The end position (1-index, inclusive) of the record"

	private String currentRecord = null;
	private String[] fields = null;
	private int lineNo = 0;
	private int recordNo = 0;

	private BufferedReader bf;
	private HashMap<String, FieldHeader> headerPool;
	private ArrayList<FieldHeader> headerList;

	public XAFReader(String fileName) throws IOException {
		bf = SequenceReader.openFile(fileName);

		// Set up headers, start by the default fields
		headerPool = new HashMap<String, FieldHeader>();
		headerPool.put(FH_ID.headerStr, FH_ID);
		headerPool.put(FH_CHROM.headerStr, FH_CHROM);
		headerPool.put(FH_START.headerStr, FH_START);
		headerPool.put(FH_END.headerStr, FH_END);

		headerList = new ArrayList<FieldHeader>();
		// set up default headerList
		headerList.add(FH_ID);
		FH_ID.index = 0;

		headerList.add(FH_CHROM);
		FH_CHROM.index = 1;

		headerList.add(FH_START);
		FH_START.index = 2;

		headerList.add(FH_END);
		FH_END.index = 3;
	}

	public void close() throws IOException {
		bf.close();
	}

	/**
	 * Get the string of the current record
	 * 
	 * @return
	 */
	public String getCurrentRecord() {
		return currentRecord;
	}

	/**
	 * Get the current list of headers
	 * 
	 * @return
	 */

	public ArrayList<FieldHeader> getHeaderList() {
		return headerList;
	}

	/**
	 * Get field of the current record specified by fieldName
	 * @param fieldName: the name of the column
	 * @return the null if 1. fieldName invalid or there are fewer fields in the record;
	 *         else return the field
	 */
	public String getField(String fieldName) {
		FieldHeader header = headerPool.get(fieldName);
		if (header == null || header.index < 0) {
			//throw new RuntimeException("Field " + fieldName + " not found");
			//Logging.warn("Field " + fieldName + " not found");
			return null;
		}

		if (header.index >= fields.length) {
			//throw new RuntimeException("Only " + fields.length	+ " fields at line " + lineNo);
			//Logging.warn("Only " + fields.length	+ " fields at line " + lineNo);
			return null;
		}
		return fields[header.index];
	}
	/**
	 * Get a field based on the index in line (first field = 0)
	 * @param fieldNo
	 * @return
	 */
	public String getField(int fieldNo) {
		if (fieldNo >= fields.length) {
			//throw new RuntimeException("Only " + fields.length	+ " fields at line " + lineNo);
			Logging.warn("Only " + fields.length	+ " fields at line " + lineNo);
			return null;
		}
		return fields[fieldNo];
	}
	

	public int lineNo() {
		return lineNo;
	}
	
	/**
	 * Get record number (first one has number 1)
	 * @return
	 */
	public int recordNo() {
		return recordNo;
	}

	public String next() throws IOException {
		while (true) {
			if ((currentRecord = bf.readLine()) == null) {
				fields = null;
				return null;
			}
			lineNo++;
			currentRecord = currentRecord.trim();
			String[] toks;
			if (currentRecord.isEmpty())
				continue;
			
			// Get list of headers
			if (currentRecord.startsWith(INFORMATION_HEADER)) {
				// Get header information
				FieldHeader header = new FieldHeader();
				toks = currentRecord.substring(INFORMATION_HEADER.length())
						.split(",", 3);
				for (String tok : toks) {
					if (tok.startsWith("I_ID"))
						header.headerStr = tok.substring(I_ID.length());
					else if (tok.startsWith(I_TYPE)) {
						header.headerType = tok.substring(I_TYPE.length());
					} else if (tok.startsWith(I_DESC)) {
						header.headerDesc = tok.substring(I_DESC.length());
					}
				}

				if (header.headerStr != null && header.headerType != null) {
					headerPool.put(header.headerStr, header);
				}
			}// if
			else if (currentRecord.startsWith(ROW_HEADER)) {
				// Reset the headerList
				toks = currentRecord.substring(ROW_HEADER.length()).split(sep);
				if (toks.length > 0) {
					headerList.clear();
					for (FieldHeader header : headerPool.values())
						header.index = -1;

					for (int i = 0; i < toks.length; i++) {
						FieldHeader header = headerPool.get(toks[i]);
						if (header == null) {
							Logging.warn("Header " + toks[i] + " not defined");
							header = new FieldHeader();
							header.headerStr = toks[i];
							headerPool.put(toks[i], header);
						}
						header.index = i;
						headerList.add(header);
					}
				}
			} else if (currentRecord.startsWith("#")) {
			} else {
				fields = currentRecord.split(sep);
				recordNo ++;
				return currentRecord;
			}
		}

	}

	public static class FieldHeader {
		private String headerStr;
		private String headerType;
		private String headerDesc;
		private int index = -1;

		FieldHeader(String str, String type, String desc) {
			headerStr = str;
			headerType = type;
			headerDesc = desc;
		}

		FieldHeader() {
			this(null, null, null);
		}
		public String toString(){
			return headerStr;
		}
		
		public String getHeaderDesc(){
			return headerDesc;
		}
	}
}