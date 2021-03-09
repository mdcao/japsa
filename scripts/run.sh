##general script for running analysis
##ASSUMES THE PROJECT IS STORED IN $HOME/github/japsa_coverage
if [ ! $JSA_MEM ] ; then
	JSA_MEM=8000m
fi


#npTranscript=$HOME/github/npTranscript
#classp=$(ls ${npTranscript}/libs | xargs -I {} echo ${npTranscript}/libs/{} )
#classpath=$(echo $classp | sed 's/ /:/g')
JSA_CP=${japsa_coverage}/target/japsacov-1.9.5e.jar
#echo $JSA_CP



#params="--resDB=/home/lachlan/ResistanceTyping/resFinder --bamFile=rt_resType.sam --log=true    --dbpath=/home/lachlan/SpeciesTyping  --mm2_path=/home/lachlan/github/minimap2/minimap2 --readList=readList_fos1.txt --time=1 --dbs=clinicaldb"


mainclass=$1
shift;



str="java -Xmx${JSA_MEM} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8 -classpath ${JSA_CP} ${mainclass} $@"
echo "running .."
echo $str
$str
echo "finished"


