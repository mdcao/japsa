/*****************************************************************************
 * Copyright (c) 2010 Minh Duc Cao, Monash University.  All rights reserved. *
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
 * 3. Neither the name of Monash University nor the names of its contributors*
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

//This class is written by Julie Bernal and subsequently modified and maintained
//by Minh Duc Cao

package japsadev.misc.dnaPlatform.function;

import java.io.IOException;

import japsadev.misc.dnaPlatform.OptionsHandle;
import japsadev.misc.dnaPlatform.sequence.*;

/**
 * <p>
 * Title: Function
 * </p>
 * 
 * <p>
 * Description: This interface defines the methods all functions within the
 * DNAPlatform must implement
 * </p>
 * 
 * @author Julie Bernal
 * @version 1.0
 */
@SuppressWarnings("rawtypes")
public interface Function {

	/**
	 * Returns an OptionsHandle object for a particular function which contains
	 * function options.
	 * 
	 * @return OptionsHandle
	 */
	public OptionsHandle getOptionsHandle();

	/**
	 * Returns an array of classes a model can be applied to
	 * 
	 * @return Class[]
	 */
	public Class[] getTypeSequenceData();

	/**
	 * This method is used to execute funcions, which map a SequenceData object
	 * to another SequenceData object
	 * 
	 * @param data
	 *            SequenceData
	 * @return SequenceData
	 * @throws IOException
	 * @throws RuntimeException
	 */
	public SequenceData execute(OptionsHandle options, SequenceData data)
			throws IOException, RuntimeException;

}
