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

/*****************************************************************************
 *                           Revision History                                
 * 3 Jul 2015 - Minh Duc Cao: Created                                        
 * 
 ****************************************************************************/
package japsadev.tools.work;

import java.io.File;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * @author minhduc
 *
 */
public class VariableVariantVCF {
	public static void main(String[] args){
		VCFFileReader vcf = new  VCFFileReader(new File(args[0]));
		CloseableIterator<VariantContext> iter = vcf.iterator();
		
		System.out.print("#CHROM start end qual");
		
		for (String sample:vcf.getFileHeader().getSampleNamesInOrder()){
			System.out.print(" " + sample);
		}
		System.out.println();
		while (iter.hasNext()){
			VariantContext var = iter.next();
			var.getChr();
			GenotypesContext gtypes = var.getGenotypes();
			boolean dongbo = true;
			Allele firstAllele = gtypes.get(0).getAllele(0);			
			for (Genotype genotype:gtypes){
				Allele allele = genotype.getAllele(0);
				if (allele.compareTo(firstAllele) != 0){
					dongbo = false;
					break;//for
				}				
			}//for
			if (dongbo){
				continue;
			}
			//!dongbo
			System.out.print(var.getChr() + " " + var.getStart() + " " + var.getEnd() + " " + var.getPhredScaledQual());
			for (Genotype genotype:gtypes){
				System.out.print(" " + genotype.getAllele(0));
			}
			System.out.println();
			
			
		}
		
		vcf.close();
		 
	}
}
