#############################################################################
# Copyright (c) Minh Duc Cao, Monash Uni & UQ, All rights reserved.         #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions        #
# are met:                                                                  # 
#                                                                           #
# 1. Redistributions of source code must retain the above copyright notice, #
#    this list of conditions and the following disclaimer.                  #
# 2. Redistributions in binary form must reproduce the above copyright      #
#    notice, this list of conditions and the following disclaimer in the    #
#    documentation and/or other materials provided with the distribution.   #
# 3. Neither the names of the institutions nor the names of the contributors#
#    may be used to endorse or promote products derived from this software  #
#    without specific prior written permission.                             #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   #
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, #
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    #
# PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         #
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     #
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       #
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
############################################################################/


############################################################################
# A generic makefile for java
#  
############################################################################


###########################################################################
# Modify the following parameter to suit your setting

#Class directory
CLASS_DIR=target/main-classes
#Source directory
SRC_DIR=src/main/java

#Java compiler
JCC=javac

#target jar file
JAR_FILE=japsa.jar

#External library directory
LIB_DIR=libs

#List of external libraries
EXT_LIBS=colloquial.jar commons-math3-3.0.jar jhdf5.jar jhdfobj.jar htsjdk-1.126.jar guava-18.0.jar jcommon-1.0.23.jar jfreechart-1.0.19.jar JRIEngine.jar JRI.jar jaligner.jar 

###########################################################################
##What this scripts does:
# 1. Find all packages by searching SRC_DIR subdirectories
# 2. Get all the java file, and compile them to the CLASS_DIR directory
# 5. Clean back up files (*.bak, ~, etc)
# To be implemented:
# 3. Copy all the resources to class directory
# 4. Make java file
##########################################################################

COMMA:= ,
EMPTY:=
SPACE:= $(SPACE) $(SPACE)
LIBS:=$(subst $(SPACE),:, $(addprefix $(LIB_DIR)/, $(EXT_LIBS)))

##Get all the packages in $(SRC_DIR). 
PACKAGE_DIRS := $(shell echo `cd $(SRC_DIR);find . -type d`)
SRC_DIRS := $(addprefix $(SRC_DIR)/, $(PACKAGE_DIRS))

##For each package directory, get its equivalent class directory
CLASS_DIRS := $(addprefix $(CLASS_DIR)/, $(PACKAGE_DIRS))

##Get the list of java files
JAVA_FILES  := $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.java))
CLASS_FILES := $(subst $(SRC_DIR)/,$(CLASS_DIR)/,$(JAVA_FILES:.java=.class))


#####################Make targets 

VPATH=$(subst ' ',':',$(PACKAGE_DIRS)) 
$(CLASS_DIR)/%.class: $(SRC_DIR)/%.java $(CLASS_DIR)
	$(JCC) -sourcepath $(SRC_DIR) -cp $(CLASS_DIR):$(LIBS) -nowarn -d $(CLASS_DIR) $(JDEBUGFLAGS) $<


all: jar

classes:  $(CLASS_FILES)


###Create the class directory if neccessary
$(CLASS_DIR):
	mkdir -p $(CLASS_DIR)

images:
	cp -r  $(SRC_DIR)/japsa/bio/misc/dnaPlatform/gui/images $(CLASS_DIR)/japsa/bio/misc/dnaPlatform/gui/
	cp -r  $(SRC_DIR)/japsa/seq/nanopore/icons $(CLASS_DIR)/japsa/seq/nanopore/  

$(JAR_FILE):classes images
	jar cf $(JAR_FILE) -C $(CLASS_DIR) . 

jar: $(JAR_FILE)

clean:
	@@for i in $(PACKAGE_DIRS); do \
		echo "rm -f $(SRC_DIR)/$$i/*.bak $(SRC_DIR)/$$i/*.class $(SRC_DIR)/$$i/*~ $(CLASS_DIR)/$$i/*.class"; \
		rm -f $(SRC_DIR)/$$i/*.bak $(SRC_DIR)/$$i/*.class $(SRC_DIR)/$$i/*~ $(CLASS_DIR)/$$i/*.class; \
	done
	rm -f install.sh install.bat

##############################################################################
#                            INSTALLATION                   
##############################################################################

ifdef INSTALL_DIR
O_INS_DIR=-installDir=\"${INSTALL_DIR}\"
else
O_INS_DIR=
endif


ifdef JLP
O_JLP=-jlp=\"${JLP}\"
else
O_JLP=
endif

ifdef MXMEM
O_MXMEM=-xmx=${MXMEM}
else
O_MXMEM=
endif


ifdef SERVER
O_SERVER=--server=${SERVER}
else
O_SERVER=
endif


RELEASE=JapsaRelease

pre-install: jar
	@@echo "java -cp $(JAR_FILE):$(LIB_DIR)/guava-18.0.jar japsa.util.deploy.Deploy --mode install --libs $(subst $(SPACE),:, $(EXT_LIBS)) $(O_INS_DIR) $(O_JLP) $(O_MXMEM) $(O_SERVER) --compiler \"`$(JCC) -version 2>&1`\"" > install.sh && \
	chmod u+x install.sh && \
	echo "java -cp $(JAR_FILE);$(LIB_DIR)\guava-18.0.jar japsa.util.deploy.Deploy --mode install --libs $(subst $(SPACE),:, $(EXT_LIBS)) $(O_INS_DIR) $(O_JLP) $(O_MXMEM) $(O_SERVER) --compiler \"`$(JCC) -version 2>&1`\"" > install.bat && \
	echo "Installation scripts created"

install: pre-install
	./install.sh

uninstall:
	@@java -cp $(JAR_FILE):$(LIB_DIR)/guava-18.0.jar japsa.util.deploy.Deploy --mode uninstall --libs $(subst $(SPACE),:, $(EXT_LIBS)) ${O_INS_DIR} && echo "Japsa uninstalled!"

release: pre-install
	@@mkdir -p $(RELEASE)/libs/ && \
	cp -f $(addprefix $(LIB_DIR)/, $(EXT_LIBS)) $(JAR_FILE)  $(RELEASE)/libs && \
	cp -f $(JAR_FILE) install.sh install.bat LICENSE.md README.md  $(RELEASE)/ && \
	sleep 5 && \
	zip -r  $(RELEASE).zip  $(RELEASE)/ && \
	tar zcvf  $(RELEASE).tar.gz  $(RELEASE)/
#############################################################################

