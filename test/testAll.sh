#!/bin/bash
#test script for metadamage, 
#

PRG=../metaDMG-cpp
BAM1=./data/f570b1db7c.dedup.filtered.bam


echo "--------------------"
echo "Using PRG: '${PRG}' and BAMFILE: '${BAM1}'"
echo "--------------------"

LOG=${0}.log
echo Using Logfile: ${LOG}
rm -f ${LOG} *.bin #remove old logfile and binary tempfile

RVAL=0

echo "Testing Existence of ${PRG}"
if [[ ! -f "${PRG}" ]]; then
    echo "Problem finding program: ${PRG}"
    RVAL=1
fi

echo "Testing Existence of ${BAM1}"
if [[ ! -f "${BAM1}" ]]; then
    echo "Problem finding file: ${PRG}"
    RVAL=$((2+RVAL))
fi

echo "Sorting bamfile"
BAM="$(dirname ${BAM1})/$(basename ${BAM1} .bam).rname.bam"
CMD="samtools sort -n ${BAM1} -o ${BAM}"
${CMD}
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((4+RVAL))
fi

mkdir -p output

echo "Validating getdamage local"
CMD="${PRG} getdamage -r 0 ${BAM} -o output/test1"
${CMD} &>>${LOG}
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((8+RVAL))
fi

echo "Validating getdamage global"
CMD="${PRG} getdamage -r 1 ${BAM} -o output/test2"
${CMD} &>>${LOG}
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((16+RVAL))
fi

echo "Validating lca"
CMD="${PRG} lca -bam ${BAM} -names data/names.dmp.gz -nodes data/nodes.dmp.gz -acc2tax data/acc2taxid.map.gz -simscorelow 0.9 -simscorehigh 1.0 -editdistmin 0 -editdistmax 10 -minmapq 0 -howmany 35 -weighttype 1 -fix_ncbi 0 -out output/test3"
${CMD} &>>${LOG}
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((32+RVAL))
fi
## below is for validating print_ugly with different parameters, should be updated with a print_ugly --out argument so it is not nescessary
CMD="${PRG} lca -bam ${BAM} -names data/names.dmp.gz -nodes data/nodes.dmp.gz -acc2tax data/acc2taxid.map.gz -simscorelow 0.9 -simscorehigh 1.0 -editdistmin 0 -editdistmax 10 -minmapq 0 -howmany 35 -weighttype 1 -fix_ncbi 0 -out output/test4"
${CMD} &>>${LOG}
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((32+RVAL))
fi



echo "Validating printoptions"

CMD="${PRG} print output/test1.bdamage.gz "
${CMD} 1>output/test1.bdamage.gz.txt 2>>${LOG}
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((64+RVAL))
fi

CMD="${PRG} print output/test2.bdamage.gz"
${CMD} 1>output/test2.bdamage.gz.txt 2>>${LOG}
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((128+RVAL))
fi

CMD="${PRG} print_ugly output/test3.bdamage.gz"
${CMD} &>>${LOG}
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((256+RVAL))
fi

CMD="${PRG} print_ugly output/test4.bdamage.gz -names data/names.dmp.gz -nodes data/nodes.dmp.gz -lcastat output/test4.stat "
${CMD} &>>${LOG}
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((512+RVAL))
fi


echo "Copying logfile and validating checksum"
echo "========================"
grep -v version ${LOG}|grep -v VERSION|grep -v tempfolder |grep -v walltime|grep -v taken|grep -v thread1 >output/logfile


gunzip -c output/test3.lca.gz|sed 1d |md5sum -c files2.md5
if [[ $? -ne 0 ]]; then
    echo "Problem with md5sum for lca file"
    RVAL=$((1024+${RVAL}))
fi

#md5sum output/*|grep -v lca.gz >files.md5

md5sum -c files.md5 
if [[ $? -ne 0 ]]; then
    echo "Problem with md5sums"
    RVAL=$((2048+${RVAL}))
fi
echo "=====RVAL:${RVAL}======="


if [[ ${RVAL} -ne 0 ]];then
    echo "====Prunedlog==="
    cat output/logfile
    echo "====Purelog==="
    cat ${LOG}
    echo "====EndOfLog==="  
fi

exit ${RVAL}
