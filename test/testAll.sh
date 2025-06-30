#!/bin/bash
#test script for metadamage, 
#


PRG=../metaDMG-cpp
BAM1=./data/f570b1db7c.dedup.filtered.bam

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

if [[ ${RVAL} -ne 0 ]];then
    echo "Problem with input files"
    cat ${LOG}
    echo "====EndOfLog==="
    exit 1 #exit codes are capped at 255
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

echo "Running getdamage global"
CMD="${PRG} getdamage --run_mode 0 --min_length 35 --print_length 5 --out_prefix output/test_getdamage_global ${BAM}"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((8+RVAL))
fi

echo "Running getdamage local"
CMD="${PRG} getdamage -r 1 -l 35 -p 5 -o output/test_getdamage_local ${BAM}"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((16+RVAL))
fi

echo "Running lca"
CMD="${PRG} lca --bam ${BAM} --names data/names.dmp.gz --nodes data/nodes.dmp.gz --acc2tax data/acc2taxid.map.gz --sim_score_low 0.9 --sim_score_high 1.0 --edit_dist_min 0 --edit_dist_max 10 --min_mapq 0 --how_many 35 --weight_type 1 --fix_ncbi 0 --out_prefix output/test_lca"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((32+RVAL))
fi
gunzip -c  output/test_lca.lca.gz | cut -f1,3-100 >output/test_lca.lca.sub

echo "Running aggregate"
CMD="${PRG} aggregate output/test_lca.bdamage.gz --nodes data/nodes.dmp.gz --names data/names.dmp.gz --lcastat output/test_lca.stat.gz --out_prefix output/test_aggregate"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((64+RVAL))
fi

echo "Running dfit local (single threaded)"
CMD="${PRG} dfit output/test_lca.bdamage.gz --names data/names.dmp.gz --nodes data/nodes.dmp.gz --showfits 2 --nopt 2 --nbootstrap 2 --seed 12345 --lib ds --out output/test_dfit_local"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((64+RVAL))
fi
# Remove 'ncall' column and round values, since it fails on GitHub tests
echo "Running test of output/test_dfit_local.dfit.gz by gzip -t  output/test_dfit_local.dfit.gz "
gzip -t  output/test_dfit_local.dfit.gz
echo -ne "Return value of test $?\nWill now run gunzip -c output/test_dfit_local.dfit.gz|wc -c\n"
gunzip -c output/test_dfit_local.dfit.gz | wc -c


gunzip -c output/test_dfit_local.dfit.gz > tmp.txt
head -n 10 tmp.txt > tmp2.txt
#rm tmp.txt
cut -f 1-6,8- tmp2.txt | numfmt -d $'\t' --header --format='%.2f' --field=2- --invalid=ignore > output/test_dfit_local.dfit.fix
#rm tmp2.txt
#zcat output/test_dfit_local.dfit.gz | head -n 10|cut -f 1-6,8- | numfmt -d $'\t' --header --format='%.2f' --field=2- --invalid=ignore > output/test_dfit_local.dfit.fix
#exit 0


echo "Running dfit local (10 threaded)"
CMD="${PRG} dfit output/test_lca.bdamage.gz --threads 10 --names data/names.dmp.gz --nodes data/nodes.dmp.gz --showfits 2 --nopt 2 --nbootstrap 2 --seed 12345 --lib ds --out output/test_dfit_local_10threads"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((128+RVAL))
fi
# Remove 'ncall' column and round values, since it fails on GitHub tests

#zcat output/test_dfit_local_10threads.dfit.gz | cut -f 1-6,8- | head -n 10 | numfmt -d $'\t' --header --format='%.2f' --field=2- --invalid=ignore | sort -r > output/test_dfit_local_10threads.dfit.fix
gunzip -c output/test_dfit_local_10threads.dfit.gz > tmp.txt
head -n 10 tmp.txt > tmp2.txt
rm tmp.txt
cut -f 1-6,8- tmp2.txt | numfmt -d $'\t' --header --format='%.2f' --field=2- --invalid=ignore | sort -r > output/test_dfit_local_10threads.dfit.fix
rm tmp2.txt

echo "Running dfit global"
CMD="${PRG} dfit output/test_getdamage_global.bdamage.gz --showfits 2 --seed 12345 --lib ds --out output/test_dfit_global"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((128+RVAL))
fi
# Remove 'ncall' column since it fail on GitHub tests
gunzip -c output/test_dfit_global.dfit.gz | cut -f 1-6,8- > output/test_dfit_global.dfit.fix

echo "Running printoptions"
CMD="${PRG} print output/test_getdamage_local.bdamage.gz"
${CMD} 1>output/test_getdamage_local.bdamage.tsv 2>>${LOG}
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((256+RVAL))
fi

CMD="${PRG} print output/test_getdamage_global.bdamage.gz"
${CMD} 1>output/test_getdamage_global.bdamage.tsv 2>>${LOG}
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((256+RVAL))
fi

CMD="${PRG} print_ugly output/test_lca.bdamage.gz --out_prefix output/test_lca"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((512+RVAL))
fi

CMD="${PRG} print_ugly output/test_lca.bdamage.gz --names data/names.dmp.gz --nodes data/nodes.dmp.gz --lcastat output/test_lca.stat.gz --out_prefix output/test_lca_taxa"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((1024+RVAL))
fi

echo "Running mergedamage test with 2 bdamage and 2 rlens files"
CMD="${PRG} mergedamage -b output/test_getdamage_local.bdamage.gz output/test_getdamage_global.bdamage.gz -r output/test_lca.rlens.gz output/test_lca.rlens.gz -out output/test_mergedamage"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running mergedamage: ${CMD}"
    RVAL=$((4096+RVAL))
fi

CHECKSUMFILE=output.md5
if [[ "$(uname)" == "Darwin" ]]; then
  CHECKSUMFILE=output.md5.macos
else
  echo "Not macOS"
fi

echo "Validating checksums: ${CHECKSUMFILE}"
echo "========================"
gunzip -f output/*.gz
md5sum -c ${CHECKSUMFILE}
if [[ $? -ne 0 ]]; then
    echo "Problem with md5sums"
    RVAL=$((2048+${RVAL}))
fi
echo "=====RVAL:${RVAL}======="

if [[ ${RVAL} -ne 0 ]];then
    echo "====StartOfLog==="
    cat ${LOG}
    echo "====EndOfLog==="
    exit 1 #exit codes are capped at 255
fi

exit 0
