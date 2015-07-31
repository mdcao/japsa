
Dependencies
------------

Japsa has the following dependencies, some of which are included in the
package.

* Java Runtime Environment (java) >=1.6
* Jfreechart >= 1.0.19 (http://www.jfree.org/jfreechart/) and Jcommon (http://www.jfree.org/jcommon/)-- included
* colloquial.jar, now java-arithcode (https://github.com/bob-carpenter/java-arithcode)  -- included
* common-math.jar >=3.3.0 (http://commons.apache.org/proper/commons-math/) --included
* htsjdk >=1.126 -- included
* sam >=1.84 -- included. This will be replaced by htsjdk above
* guava >=16.0 -- included
* jhdf5 >= 18 -- included
* JRI -- included. Dependency on JRI will be removed in the future
* HDF-Java https://www.hdfgroup.org/products/java/release/download.html -- needed for npReader
* R and rJava -- needed only for speciesTyping

Japsa is provided with a ready to run package at every stable release.
These pre-built releases are compiled with javac 1.6 to ensure compatibility.
If you wish to use the latest version of Japsa, or use japsa compiled with your
version of Java, you will need to build from source, which requires:

* Java Development Kit (javac) >= 1.6
* make

