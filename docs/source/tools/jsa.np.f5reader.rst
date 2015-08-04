-----------------------------------------------------------------------
npReader: real-time conversion and analysis of Nanopore sequencing data
-----------------------------------------------------------------------

npReader (jsa.np.f5reader) is a program that extract Oxford Nanopore sequencing
data from FAST5 files, perform an initial analysis of the date and stream them
to real-time analysis pipelines. These pipelines can run on the same computer
or on computing clouds/high performance clusters.

npReader is included in the Japsa package. It requires
`JAVA HDF5 INTERFACE (JHI5) library <https://www.hdfgroup.org/products/java/JNI/jhi5/index.html>`_
to be installed prior to setting up Japsa. Details of installation as follows:

**On Windows/Mac**

#. Download and install HDF-View from https://www.hdfgroup.org/products/java/release/download.html.
Note the folder that the JHI library is installed, e.g.,
*C:\\Program Files\\HDF_Group\\HDFView\\2.11.0\\lib*

#. Follow the instructions to install Japsa on (install_). Upon prompting for
"Path to HDF library", enter the above path.

**On Linux**

You can either install the JHI5 library by downloading the software from
*https://www.hdfgroup.org/products/java/JNI/jhi5/index.html* or from your
Linux distribution software repository, such as::
   
   sudo apt-get install libjhdf5-jni
   
The library is typically installed to */usr/lib/jni*. Enter this path when
prompt for "Path to HDF library" during installation of Japsa.

Synopsis
********
jsa.np.f5reader: Extract Oxford Nanopore sequencing data from FAST5 files,
perform an initial analysis of the date and stream them to realtime analysis
pipelines

Usage
*****
::

   jsa.np.f5reader [options]

Options
*******
  --GUI           Run with a Graphical User Interface
                  (default='false')
  --realtime      Run the program in real-time mode, i.e., keep waiting for new data from Metrichon agent
                  (default='false')
  --folder=s      The folder containing base-called reads
                  (default='null')
  --fail          Get sequence reads from fail folder
                  (default='false')
  --output=s      Name of the output file, - for stdout
                  (default='-')
  --streams=s     Stream output to some servers, format "IP:port,IP:port" (no spaces)
                  (default='null')
  --format=s      Format of sequence reads (fastq or fasta)
                  (default='fastq')
  --minLength=i   Minimum read length
                  (default='0')
  --number        Add a unique number to read name
                  (default='false')
  --stats         Generate a report of read statistics
                  (default='false')
  --time          Extract the sequencing time of each read -- only work with Metrichon > 1.12
                  (default='false')
  --help          Display this usage and exit
                  (default='false')


Usage examples
**************




