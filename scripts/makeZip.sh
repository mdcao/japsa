# while read line; do bash ../makeZip.sh ${line} 2; done< groups.txt
# bash ../makeZip.sh blast 3
# bash ../makeZip.sh blast 3 results_combined.krkn
tofind=results.krkn
#rm -rf ${1}.totransf
if [ $2 ]; then
	pos=$2
else
	pos=3
fi
if [ $3 ] ; then
	tofind=$3
fi
echo $tofind
find . -name ${tofind} -size +1b | grep $1 > totransf.${1}.txt
cmd="find . -name ${tofind} -size +1b | grep blast | grep $1 > totransf.${1}.txt"
totransf=${1}.totransf
 mkdir ${totransf}
 outnme=$( pwd | rev | cut -f 1 -d /  | rev)
while read line; do
	echo $line
	nme=$(echo $line | cut -f $pos -d /)
	cp $line ${totransf}/${nme}.krkn

done < totransf.${1}.txt
cd ${totransf}
ls  | grep krkn | xargs zip out.zip
ls | grep krkn | xargs rm 
mv out.zip ../${outnme}.${1}.zip
cd ..
rmdir ${totransf}
