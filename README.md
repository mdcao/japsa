### Welcome to Japsa

Japsa is a Java Package for Sequence Analysis. As the name implies, Japsa was 
written primarily in Java. It contains a large number of ready to run programs 
as well as a Java API. Please note that from our newest versions, it is required 
Java 8 to compile and run from the source code. Prebuilt releases would need 
JVM 1.8.0_144 or newer to run properly.

### Installation and Usage

Quick installation guide:
```
git clone https://github.com/mdcao/japsa.git
cd japsa
make install \
  [INSTALL_DIR=~/.usr/local \] 
  [MXMEM=7000m \] 
  [SERVER=true \] 
  [JLP=/usr/lib/jni]
```

Details of installation (including for Windows) and usage of Japsa can be found 
in its documentation hosted on [ReadTheDocs](http://japsa.readthedocs.org/en/latest/index.html) 

Alternatively, build with maven (experimental). 
First you need to manually install a JAR file that is not available on Maven
public repo to your local repo
```
mvn install:install-file \
 	-Dfile=./libs/colloquial.jar \
	-DgroupId=com.colloquial \
	-DartifactId=arithcode \
	-Dversion=1.1 \
	-Dpackaging=jar
```
then 
```
mvn clean package install
```
you might need to try packaging Japsa again if failed. A SNAPSHOT is then created and you can invoke
tool, e.g. species typer, by
```
java -cp ./target/japsa-1.0-SNAPSHOT.jar japsa.tools.bio.np.RealtimeSpeciesTypingCmd --bam <bam> --index <index>
```
Convenient scripts (like install by make) would be added soon.

### Authors and Contributors
Japsa is currently maintained by  ~~Minh Duc Cao (@mdcao)~~ [Son Hoang Nguyen](https://github.com/hsnguyen). The following 
people (in alphatical order) have contributed to the development of Japsa, including ideas, 
algorithms, implementation, documentation and feedback:

* Allen Day
* Bhuvan Sankar
* David Powell
* Devika Ganesamoorthy
* Hoang Anh Nguyen
* Julia Bernal
* [Lachlan Coin](http://www.imb.uq.edu.au/lachlan-coin)
* [Lloyd Allison](http://www.allisons.org/ll/)
* [Michael Hall](https://github.com/mbhall88)
* Mikael Boden
* Minh Duc Cao
* [Son Hoang Nguyen](https://github.com/hsnguyen)
* Trevor I Dix


### Other projects based on Japsa


* [eXpert Model](https://github.com/mdcao/xm): The expert model compression model
* [XMas](https://github.com/mdcao/XMas): Phylogenetic distance method using information theory
* [capsim](https://github.com/mdcao/capsim): Simulation of capture sequencing
* [npScarf](https://github.com/mdcao/npScarf): Scaffold and Complete assemblies in real-time fashion
* [npAnalysis](https://github.com/mdcao/npAnalysis): Realtime identification of bacterial sample
* [npReader](https://github.com/mdcao/npReader): Real-time extraction and analysis Oxford Nanopore sequencing data
* [npBarcode](https://github.com/hsnguyen/npBarcode): Demultiplex barcoded Oxford Nanopore sequencing 
* [PhageXpress](https://github.com/mdcao/phagexpress)

and more to come.


### License
Japsa is released under the accompanying BSD-like license.

