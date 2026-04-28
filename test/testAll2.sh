#!/usr/bin/env bash
# test script for metaDMG

set -u
set -o pipefail

if [[ -x ../metaDMG-cpp ]]; then
    PRG=../metaDMG-cpp
else
    PRG=metaDMG-cpp
fi

BAM1=./data/f570b1db7c.dedup.filtered.bam
BAM="$(dirname "${BAM1}")/$(basename "${BAM1}" .bam).rname.bam"
LOG="$(basename "$0").log"
RVAL=0

note() {
    echo "$*"
}

mark_fail() {
    echo "$*"
    RVAL=1
}

run_logged() {
    local label="$1"
    shift

    note "${label}"
    "$@" >>"${LOG}" 2>&1
    local rc=$?
    if [[ ${rc} -ne 0 ]]; then
        mark_fail "Problem running ${label}"
    fi
    return 0
}

require_file() {
    local path="$1"
    if [[ ! -f "${path}" ]]; then
        mark_fail "Problem finding file: ${path}"
        return 1
    fi
    return 0
}

require_cmd() {
    local cmd="$1"
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        mark_fail "Problem finding program: ${cmd}"
        return 1
    fi
    return 0
}

prepare() {
    note "Using Logfile: ${LOG}"
    rm -f "${LOG}" ./*.bin
    mkdir -p output
    mkdir -p output_data2
    mkdir -p output_data2_gd
    mkdir -p output_compressbam
    mkdir -p output_extract_reads
}

assert_file_contains() {
    local file="$1"
    local pattern="$2"

    if ! grep -F -- "${pattern}" "${file}" >/dev/null; then
        mark_fail "Expected pattern not found in ${file}: ${pattern}"
    fi
}

assert_file_not_contains() {
    local file="$1"
    local pattern="$2"

    if grep -F -- "${pattern}" "${file}" >/dev/null; then
        mark_fail "Unexpected pattern found in ${file}: ${pattern}"
    fi
}

assert_file_line_count() {
    local file="$1"
    local expected="$2"
    local actual

    actual="$(wc -l < "${file}" | tr -d '[:space:]')"
    if [[ "${actual}" != "${expected}" ]]; then
        mark_fail "Unexpected line count in ${file}: expected ${expected}, got ${actual}"
    fi
}

assert_gzip_contains() {
    local file="$1"
    local pattern="$2"

    if ! gunzip -c "${file}" | grep -F -- "${pattern}" >/dev/null; then
        mark_fail "Expected pattern not found in ${file}: ${pattern}"
    fi
}

assert_gzip_not_contains() {
    local file="$1"
    local pattern="$2"

    if gunzip -c "${file}" | grep -F -- "${pattern}" >/dev/null; then
        mark_fail "Unexpected pattern found in ${file}: ${pattern}"
    fi
}

test_dust_unit() {
    local unit_src="test_dust_score.cpp"
    local unit_bin="output/test_dust_score_unit"
    local hts_inc=""
    local hts_link=""

    if [[ -f ../htslib/libhts.a ]]; then
        hts_inc="-I../htslib"
        hts_link="../htslib/libhts.a -lbz2 -llzma -lcurl -lz -lm -lpthread"
    else
        hts_inc=""
        hts_link="-lhts -lz -lm -lpthread"
    fi

    run_logged "Compiling dust unit test" \
        c++ -O2 -Wall -Wextra -I.. ${hts_inc} "${unit_src}" ../shared.o \
        -o "${unit_bin}" ${hts_link}

    if [[ -x "${unit_bin}" ]]; then
        run_logged "Running dust unit test" "${unit_bin}"
    else
        mark_fail "Problem building dust unit test binary: ${unit_bin}"
    fi
}

test_getdamage() {
    run_logged "Running getdamage global" \
        "${PRG}" getdamage --run_mode 0 --min_length 35 --print_length 5 \
        --out_prefix output/test_getdamage_global "${BAM}"

    run_logged "Running getdamage local" \
        "${PRG}" getdamage -r 1 -l 35 -p 5 \
        -o output/test_getdamage_local "${BAM}"

    run_logged "Running getdamage taxid" \
        "${PRG}" getdamage --run_mode 2 --acc2tax data/acc2taxid.map.gz \
        --min_length 35 --print_length 5 \
        --out_prefix output/test_getdamage_taxid "${BAM}"

    assert_gzip_contains output/test_getdamage_taxid.stat.gz $'taxid\tnreads\tmean_len\tvar_len\tmean_gc\tvar_gc\tlca\trank'
    assert_gzip_not_contains output/test_getdamage_taxid.stat.gz $'global\t'
}

test_lca_aggregate() {
    run_logged "Running lca" \
        "${PRG}" lca --bam "${BAM}" --names data/names.dmp.gz \
        --nodes data/nodes.dmp.gz --acc2tax data/acc2taxid.map.gz \
        --sim_score_low 0.9 --sim_score_high 1.0 --edit_dist_min 0 \
        --edit_dist_max 10 --min_mapq 0 --how_many 35 --weight_type 1 \
        --fix_ncbi 0 --out_prefix output/test_lca

    if ! gunzip -c output/test_lca.lca.gz | cut -f1,3-100 > output/test_lca.lca.sub; then
        mark_fail "Problem preparing output/test_lca.lca.sub"
    fi

    run_logged "Running aggregate" \
        "${PRG}" aggregate output/test_lca.bdamage.gz --nodes data/nodes.dmp.gz \
        --names data/names.dmp.gz --lcastat output/test_lca.stat.gz \
        --out_prefix output/test_aggregate
}

test_dfit() {
    run_logged "Running dfit local (single threaded)" \
        "${PRG}" dfit output/test_lca.bdamage.gz --names data/names.dmp.gz \
        --nodes data/nodes.dmp.gz --showfits 2 --nopt 2 --nbootstrap 2 \
        --seed 12345 --lib ds --out output/test_dfit_local

    note "Running test of output/test_dfit_local.dfit.gz by gzip -t"
    if ! gzip -t output/test_dfit_local.dfit.gz; then
        mark_fail "Problem validating output/test_dfit_local.dfit.gz with gzip -t"
    fi

    note "Will now run gunzip -c output/test_dfit_local.dfit.gz | wc -c"
    if ! gunzip -c output/test_dfit_local.dfit.gz | wc -c; then
        mark_fail "Problem reading output/test_dfit_local.dfit.gz"
    fi

    if ! gunzip -c output/test_dfit_local.dfit.gz | \
        cut -f1-6,8- | \
        awk 'BEGIN{OFS=FS="\t"} NR==1 {print; next} {for(i=2;i<=NF;i++) if($i ~ /^[0-9.eE+-]+$/) $i = sprintf("%.2f", $i); print}' \
        > output/test_dfit_local.dfit.fix; then
        mark_fail "Problem creating output/test_dfit_local.dfit.fix"
    fi

    run_logged "Running dfit local (4 threaded)" \
        "${PRG}" dfit output/test_lca.bdamage.gz --threads 4 \
        --names data/names.dmp.gz --nodes data/nodes.dmp.gz --showfits 2 \
        --nopt 2 --nbootstrap 2 --seed 12345 --lib ds \
        --out output/test_dfit_local_10threads

    if ! gunzip -c output/test_dfit_local_10threads.dfit.gz | \
        cut -f1-6,8- | ./round_file.sh | sort -k1,1n \
        > output/test_dfit_local_10threads.dfit.fix; then
        mark_fail "Problem creating output/test_dfit_local_10threads.dfit.fix"
    fi

    run_logged "Running dfit global" \
        "${PRG}" dfit output/test_getdamage_global.bdamage.gz --showfits 2 \
        --seed 12345 --lib ds --out output/test_dfit_global

    if ! gunzip -c output/test_dfit_global.dfit.gz | cut -f1-6,8- | \
        ./round_file.sh > output/test_dfit_global.dfit.fix; then
        mark_fail "Problem creating output/test_dfit_global.dfit.fix"
    fi
}

test_prints() {
    note "Running printoptions"

    if ! "${PRG}" print output/test_getdamage_local.bdamage.gz \
        1>output/test_getdamage_local.bdamage.tsv 2>>"${LOG}"; then
        mark_fail "Problem running print on output/test_getdamage_local.bdamage.gz"
    fi

    if ! "${PRG}" print output/test_getdamage_global.bdamage.gz \
        1>output/test_getdamage_global.bdamage.tsv 2>>"${LOG}"; then
        mark_fail "Problem running print on output/test_getdamage_global.bdamage.gz"
    fi

    if ! "${PRG}" print output/test_getdamage_taxid.bdamage.gz \
        1>output/test_getdamage_taxid.bdamage.tsv 2>>"${LOG}"; then
        mark_fail "Problem running print on output/test_getdamage_taxid.bdamage.gz"
    fi
    assert_file_contains output/test_getdamage_taxid.bdamage.tsv $'taxid\tNalignments\tDirection\tPos'

    run_logged "Running print_ugly basic" \
        "${PRG}" print_ugly output/test_lca.bdamage.gz --out_prefix output/test_lca

    run_logged "Running print_ugly taxa" \
        "${PRG}" print_ugly output/test_lca.bdamage.gz --names data/names.dmp.gz \
        --nodes data/nodes.dmp.gz --lcastat output/test_lca.stat.gz \
        --out_prefix output/test_lca_taxa
}

test_data2() {
    run_logged "Running data2 lca sam2" \
        "${PRG}" lca --bam data2/sam2.sam --names data2/names.dmp \
        --nodes data2/nodes.dmp --acc2tax data2/acc2taxid.map \
        --sim_score_low 0.9 --sim_score_high 1.0 --edit_dist_min 0 \
        --edit_dist_max 10 --min_mapq 0 --how_many 35 --weight_type 1 \
        --fix_ncbi 0 --out_prefix output_data2/sam2

    run_logged "Running data2 lca sam3" \
        "${PRG}" lca --bam data2/sam3.sam --names data2/names.dmp \
        --nodes data2/nodes.dmp --acc2tax data2/acc2taxid.map \
        --sim_score_low 0.9 --sim_score_high 1.0 --edit_dist_min 0 \
        --edit_dist_max 10 --min_mapq 0 --how_many 35 --weight_type 1 \
        --fix_ncbi 0 --out_prefix output_data2/sam3

    run_logged "Running data2 lca sam5" \
        "${PRG}" lca --bam data2/sam5.sam --names data2/names.dmp \
        --nodes data2/nodes.dmp --acc2tax data2/acc2taxid.map \
        --sim_score_low 0.9 --sim_score_high 1.0 --edit_dist_min 0 \
        --edit_dist_max 10 --min_mapq 0 --how_many 35 --weight_type 1 \
        --fix_ncbi 0 --out_prefix output_data2/sam5

    assert_gzip_contains output_data2/sam2.lca.gz $'queryid\tseq\tlen\tnaln\tnspec\tdustscore\tgc\tlca\ttaxa_path'
    assert_gzip_contains output_data2/sam2.lca.gz $'read1\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t30\t3\t1\t'
    assert_gzip_contains output_data2/sam2.lca.gz $'\t10:"l__Bacteria":"species"'
    assert_gzip_contains output_data2/sam2.stat.gz $'10\t1\t30.000000\t0.000000\t0.000000\t0.000000\t"l__Bacteria"\t"species"'

    assert_gzip_contains output_data2/sam3.lca.gz $'read1\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\t30\t2\t1\t'
    assert_gzip_contains output_data2/sam3.lca.gz $'read2\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t30\t1\t1\t'
    assert_gzip_contains output_data2/sam3.lca.gz $'\t10:"l__Bacteria":"species"'
    assert_gzip_contains output_data2/sam3.lca.gz $'\t11:"l__Bacteria":"subspecies"'
    assert_gzip_contains output_data2/sam3.stat.gz $'10\t1\t30.000000\t0.000000\t0.000000\t0.000000\t"l__Bacteria"\t"species"'
    assert_gzip_contains output_data2/sam3.stat.gz $'11\t1\t30.000000\t0.000000\t0.000000\t0.000000\t"l__Bacteria"\t"subspecies"'

    assert_gzip_contains output_data2/sam5.lca.gz $'read1\t'
    assert_gzip_contains output_data2/sam5.lca.gz $'\t607\t1\t1\t'
    assert_gzip_contains output_data2/sam5.lca.gz $'\t11:"l__Bacteria":"subspecies"'
    assert_gzip_contains output_data2/sam5.stat.gz $'11\t1\t607.000000\t0.000000\t0.000000\t0.000000\t"l__Bacteria"\t"subspecies"'
}

test_data2_getdamage() {
    run_logged "Running data2 getdamage sam2" \
        "${PRG}" getdamage --run_mode 1 --min_length 0 --print_length 5 \
        --out_prefix output_data2_gd/sam2 data2/sam2.sam

    run_logged "Running data2 getdamage sam3" \
        "${PRG}" getdamage --run_mode 1 --min_length 0 --print_length 5 \
        --out_prefix output_data2_gd/sam3 data2/sam3.sam

    run_logged "Running data2 getdamage sam5" \
        "${PRG}" getdamage --run_mode 1 --min_length 0 --print_length 5 \
        --out_prefix output_data2_gd/sam5 data2/sam5.sam

    run_logged "Running data2 getdamage taxid sam2 maxdust0" \
        "${PRG}" getdamage --run_mode 2 --acc2tax data2/acc2taxid.map \
        --maxdust 0 --min_length 0 --print_length 5 \
        --out_prefix output_data2_gd/sam2_taxid_maxdust0 data2/sam2.sam

    assert_gzip_contains output_data2_gd/sam2.stat.gz $'ref1\t1\t30.000000\t0.000000\t0.000000\t0.000000\tNA\tNA'
    assert_gzip_contains output_data2_gd/sam2.stat.gz $'ref2\t1\t30.000000\t0.000000\t0.000000\t0.000000\tNA\tNA'
    assert_gzip_contains output_data2_gd/sam2.stat.gz $'ref3\t1\t30.000000\t0.000000\t0.000000\t0.000000\tNA\tNA'
    assert_gzip_contains output_data2_gd/sam2.rlens.gz $'0\t30:1'
    assert_gzip_contains output_data2_gd/sam2.rlens.gz $'1\t30:1'
    assert_gzip_contains output_data2_gd/sam2.rlens.gz $'2\t30:1'

    assert_gzip_contains output_data2_gd/sam3.stat.gz $'ref1\t2\t30.000000\t0.000000\t0.000000\t0.000000\tNA\tNA'
    assert_gzip_contains output_data2_gd/sam3.stat.gz $'ref2\t1\t30.000000\t0.000000\t0.000000\t0.000000\tNA\tNA'
    assert_gzip_contains output_data2_gd/sam3.rlens.gz $'0\t30:2'
    assert_gzip_contains output_data2_gd/sam3.rlens.gz $'1\t30:1'

    assert_gzip_contains output_data2_gd/sam5.stat.gz $'ref1\t1\t607.000000\t0.000000\t0.000000\t0.000000\tNA\tNA'
    assert_gzip_contains output_data2_gd/sam5.rlens.gz $'0\t607:1'

    if ! gzip -dc output_data2_gd/sam2_taxid_maxdust0.stat.gz \
        > output_data2_gd/sam2_taxid_maxdust0.stat.tsv; then
        mark_fail "Problem decompressing output_data2_gd/sam2_taxid_maxdust0.stat.gz"
    fi
    assert_file_line_count output_data2_gd/sam2_taxid_maxdust0.stat.tsv 1

    if ! "${PRG}" print output_data2_gd/sam2_taxid_maxdust0.bdamage.gz \
        1>output_data2_gd/sam2_taxid_maxdust0.bdamage.tsv 2>>"${LOG}"; then
        mark_fail "Problem running print on output_data2_gd/sam2_taxid_maxdust0.bdamage.gz"
    fi
    assert_file_line_count output_data2_gd/sam2_taxid_maxdust0.bdamage.tsv 1
}


test_filter_bdamage() {
    run_logged "Running filter_bdamage nodes expansion (lca taxid)" \
        "${PRG}" filter_bdamage output_data2/sam3.bdamage.gz \
        --id 10 --nodes data2/nodes.dmp --out_prefix output_data2/filter_nodes

    require_file output_data2/filter_nodes.bdamage.gz || return 0
    require_file output_data2/filter_nodes.stat.gz || return 0
    require_file output_data2/filter_nodes.rlens.gz || return 0

    assert_gzip_contains output_data2/filter_nodes.stat.gz $'11\t1\t30.000000\t0.000000\t0.000000\t0.000000\t"l__Bacteria"\t"subspecies"'
    assert_gzip_not_contains output_data2/filter_nodes.stat.gz $'10\t1\t30.000000\t0.000000\t0.000000\t0.000000\t"l__Bacteria"\t"species"'

    if ! "${PRG}" print output_data2/filter_nodes.bdamage.gz \
        1>output_data2/filter_nodes.bdamage.tsv 2>>"${LOG}"; then
        mark_fail "Problem running print on output_data2/filter_nodes.bdamage.gz"
    fi
    assert_file_line_count output_data2/filter_nodes.bdamage.tsv 61

    run_logged "Running filter_bdamage resolve with BAM (getdamage local)" \
        "${PRG}" filter_bdamage output_data2_gd/sam2.bdamage.gz \
        --resolve ref2 --bam data2/sam2.sam --out_prefix output_data2_gd/filter_ref2

    require_file output_data2_gd/filter_ref2.bdamage.gz || return 0
    require_file output_data2_gd/filter_ref2.stat.gz || return 0
    require_file output_data2_gd/filter_ref2.rlens.gz || return 0

    assert_gzip_contains output_data2_gd/filter_ref2.stat.gz $'ref2\t1\t30.000000\t0.000000\t0.000000\t0.000000\tNA\tNA'
    assert_gzip_not_contains output_data2_gd/filter_ref2.stat.gz $'ref1\t1\t30.000000\t0.000000\t0.000000\t0.000000\tNA\tNA'
    assert_gzip_not_contains output_data2_gd/filter_ref2.stat.gz $'ref3\t1\t30.000000\t0.000000\t0.000000\t0.000000\tNA\tNA'

    assert_gzip_contains output_data2_gd/filter_ref2.rlens.gz $'1\t30:1'
    assert_gzip_not_contains output_data2_gd/filter_ref2.rlens.gz $'0\t30:1'
    assert_gzip_not_contains output_data2_gd/filter_ref2.rlens.gz $'2\t30:1'

    if ! "${PRG}" print output_data2_gd/filter_ref2.bdamage.gz \
        1>output_data2_gd/filter_ref2.bdamage.tsv 2>>"${LOG}"; then
        mark_fail "Problem running print on output_data2_gd/filter_ref2.bdamage.gz"
    fi
    assert_file_line_count output_data2_gd/filter_ref2.bdamage.tsv 11
}

test_compressbam() {
    local compressbam="../misc/compressbam"
    local header_all="output_compressbam/basic.all.header.sam"
    local reads_all="output_compressbam/basic.all.reads.sam"
    local header_min="output_compressbam/basic.min.header.sam"
    local reads_min="output_compressbam/basic.min.reads.sam"

    note "Testing Existence of ${compressbam}"
    if [[ ! -x "${compressbam}" ]]; then
        mark_fail "Problem finding program: ${compressbam}"
        return 0
    fi

    run_logged "Running compressbam basic" \
        "${compressbam}" --threads 1 --input data_compressbam/basic.sam \
        --output output_compressbam/basic.all.bam

    run_logged "Running compressbam basic with min-match" \
        "${compressbam}" --threads 1 --input data_compressbam/basic.sam \
        --output output_compressbam/basic.min.bam --min-match 10

    if ! samtools view -H output_compressbam/basic.all.bam > "${header_all}"; then
        mark_fail "Problem extracting header from output_compressbam/basic.all.bam"
    fi
    if ! samtools view output_compressbam/basic.all.bam > "${reads_all}"; then
        mark_fail "Problem extracting alignments from output_compressbam/basic.all.bam"
    fi
    if ! samtools view -H output_compressbam/basic.min.bam > "${header_min}"; then
        mark_fail "Problem extracting header from output_compressbam/basic.min.bam"
    fi
    if ! samtools view output_compressbam/basic.min.bam > "${reads_min}"; then
        mark_fail "Problem extracting alignments from output_compressbam/basic.min.bam"
    fi

    assert_file_contains "${header_all}" $'@SQ\tSN:refA\tLN:1000'
    assert_file_not_contains "${header_all}" $'@SQ\tSN:refB\tLN:1000'
    assert_file_not_contains "${header_all}" $'@SQ\tSN:refC\tLN:1000'
    assert_file_contains "${reads_all}" $'read_keep\t0\trefA'
    assert_file_not_contains "${reads_all}" $'read_mid\t0\trefC'

    assert_file_contains "${header_min}" $'@SQ\tSN:refA\tLN:1000'
    assert_file_contains "${header_min}" $'@SQ\tSN:refC\tLN:1000'
    assert_file_not_contains "${header_min}" $'@SQ\tSN:refB\tLN:1000'
    assert_file_contains "${reads_min}" $'read_keep\t0\trefA'
    assert_file_contains "${reads_min}" $'read_mid\t0\trefC'
}

test_extract_reads() {
    local extract_reads="../misc/extract_reads"
    local input_sam="data_extract_reads/basic.sam"
    local input_bam="output_extract_reads/basic.qname.bam"
    local byreadid_sam="output_extract_reads/byreadid.read3.sam"
    local byreadid_reads="output_extract_reads/byreadid.read3.reads.sam"
    local byreadid_complement_sam="output_extract_reads/byreadid.read3.complement.sam"
    local byreadid_complement_reads="output_extract_reads/byreadid.read3.complement.reads.sam"
    local byreadid_stdout_reads="output_extract_reads/byreadid.read3.stdout.reads.sam"
    local strict1_sam="output_extract_reads/byrefid.refC.strict1.sam"
    local strict1_reads="output_extract_reads/byrefid.refC.strict1.reads.sam"
    local strict0_sam="output_extract_reads/byrefid.refC.strict0.sam"
    local strict0_reads="output_extract_reads/byrefid.refC.strict0.reads.sam"
    local type_bam="output_extract_reads/byrefid.refD.strict1.bam"
    local type_bam_reads="output_extract_reads/byrefid.refD.strict1.reads.sam"

    note "Testing Existence of ${extract_reads}"
    if [[ ! -x "${extract_reads}" ]]; then
        mark_fail "Problem finding program: ${extract_reads}"
        return 0
    fi

    if ! samtools view -b -o "${input_bam}" "${input_sam}"; then
        mark_fail "Problem converting ${input_sam} to ${input_bam}"
        return 0
    fi

    run_logged "Running extract_reads byreadid" \
        "${extract_reads}" byreadid -i "${input_bam}" \
        -N data_extract_reads/read3.txt -o "${byreadid_sam}" --output-fmt sam

    if ! samtools view "${byreadid_sam}" > "${byreadid_reads}"; then
        mark_fail "Problem extracting alignments from ${byreadid_sam}"
    fi
    assert_file_line_count "${byreadid_reads}" 3
    assert_file_contains "${byreadid_reads}" $'read3\t0\trefA'
    assert_file_contains "${byreadid_reads}" $'read3\t0\trefB'
    assert_file_contains "${byreadid_reads}" $'read3\t0\trefC'
    assert_file_not_contains "${byreadid_reads}" $'read4\t0\trefA'

    if ! "${extract_reads}" byreadid -i "${input_bam}" -N data_extract_reads/read3.txt \
        > "${byreadid_stdout_reads}"; then
        mark_fail "Problem running extract_reads byreadid to stdout"
    fi
    assert_file_line_count "${byreadid_stdout_reads}" 9
    assert_file_contains "${byreadid_stdout_reads}" $'@SQ\tSN:refC\tLN:1000'
    assert_file_contains "${byreadid_stdout_reads}" $'read3\t0\trefA'
    assert_file_contains "${byreadid_stdout_reads}" $'read3\t0\trefC'
    assert_file_not_contains "${byreadid_stdout_reads}" $'read4\t0\trefA'

    run_logged "Running extract_reads byreadid complement" \
        "${extract_reads}" byreadid -i "${input_bam}" \
        -N data_extract_reads/read3.txt -o "${byreadid_complement_sam}" \
        --output-fmt sam --exclude

    if ! samtools view "${byreadid_complement_sam}" > "${byreadid_complement_reads}"; then
        mark_fail "Problem extracting alignments from ${byreadid_complement_sam}"
    fi
    assert_file_line_count "${byreadid_complement_reads}" 7
    assert_file_not_contains "${byreadid_complement_reads}" $'read3\t0\trefA'
    assert_file_contains "${byreadid_complement_reads}" $'read4\t0\trefD'

    run_logged "Running extract_reads byrefid strict=1" \
        "${extract_reads}" byrefid -i "${input_bam}" \
        -R data_extract_reads/refC.txt -o "${strict1_sam}" --output-fmt sam --strict 1

    if ! samtools view "${strict1_sam}" > "${strict1_reads}"; then
        mark_fail "Problem extracting alignments from ${strict1_sam}"
    fi
    assert_file_line_count "${strict1_reads}" 2
    assert_file_contains "${strict1_reads}" $'read3\t0\trefC'
    assert_file_contains "${strict1_reads}" $'read4\t0\trefC'
    assert_file_not_contains "${strict1_reads}" $'read3\t0\trefA'

    run_logged "Running extract_reads byrefid strict=0" \
        "${extract_reads}" byrefid -i "${input_bam}" \
        -R data_extract_reads/refC.txt -o "${strict0_sam}" --output-fmt sam --strict 0

    if ! samtools view "${strict0_sam}" > "${strict0_reads}"; then
        mark_fail "Problem extracting alignments from ${strict0_sam}"
    fi
    assert_file_line_count "${strict0_reads}" 7
    assert_file_contains "${strict0_reads}" $'read3\t0\trefA'
    assert_file_contains "${strict0_reads}" $'read4\t0\trefD'
    assert_file_not_contains "${strict0_reads}" $'read2\t0\trefA'

    run_logged "Running extract_reads byrefid BAM output" \
        "${extract_reads}" byrefid -i "${input_bam}" \
        -R data_extract_reads/refD.txt -o "${type_bam}" -b --strict 1

    if ! samtools view "${type_bam}" > "${type_bam_reads}"; then
        mark_fail "Problem extracting alignments from ${type_bam}"
    fi
    assert_file_line_count "${type_bam_reads}" 1
    assert_file_contains "${type_bam_reads}" $'read4\t0\trefD'
}

validate_checksums() {
    local checksum_file="output.md5"

    note "Validating checksums: ${checksum_file}"
    note "========================"

    if ! gunzip -f output/*.gz; then
        mark_fail "Problem decompressing output/*.gz before checksum validation"
        return 0
    fi

    while read -r f; do
        [[ -f "${f}" ]] || {
            mark_fail "Missing: ${f}"
            return 0
        }
    done < <(awk '!/^[[:space:]]*#/ {print $2}' "${checksum_file}")

    if ! grep -v '^[[:space:]]*#' "${checksum_file}" | md5sum -c -; then
        mark_fail "Problem with md5sums"
    fi
}

main() {
    prepare

    note "Testing Existence of ${PRG}"
    if ! command -v "${PRG}" >/dev/null 2>&1 && [[ ! -x "${PRG}" ]]; then
        mark_fail "Problem finding program: ${PRG}"
        exit 1
    fi

    note "Testing Existence of samtools"
    require_cmd samtools || exit 1

    note "Testing Existence of c++"
    require_cmd c++ || exit 1

    note "Testing Existence of ${BAM1}"
    require_file "${BAM1}" || exit 1

    note "Sorting bamfile"
    if ! samtools sort -n "${BAM1}" -o "${BAM}"; then
        mark_fail "Problem running samtools sort -n"
    fi

    test_dust_unit
    test_getdamage
    test_lca_aggregate
    test_dfit
    test_prints
    test_data2
    test_data2_getdamage
    test_filter_bdamage
    test_compressbam
    test_extract_reads
    validate_checksums

    note "=====RVAL:${RVAL}======="
    if [[ ${RVAL} -ne 0 ]]; then
        note "====StartOfLog==="
        cat "${LOG}"
        note "====EndOfLog==="
        exit 1
    fi
    exit 0
}

main "$@"
