#!/bin/bash

INPUT="/mnt/c/Users/USUARIO/Desktop/TFM VARIOS/genomas_entrenamiento/proteinas_kofam"
PROFILES="/home/melaniepolo/kofam_db/profiles_107"
KO_LIST="/home/melaniepolo/kofam_database/ko_list"
OUTPUT="/home/melaniepolo/kofam_db/resultados_53"
TMP="/home/melaniepolo/kofam_db/tmp_53"

mkdir -p "$OUTPUT"
mkdir -p "$TMP"

run_kofam() {

    faa="$1"
    accession=$(basename "$faa" .faa)

    outfile="$OUTPUT/${accession}_kofam107.tsv"
    tmpdir="$TMP/$accession"
    logfile="$OUTPUT/${accession}.log"

    if [ -f "$outfile" ]; then
        echo "YA EXISTE: $accession"
        return
    fi

    mkdir -p "$tmpdir"

    echo "========================================" | tee "$logfile"
    echo "Procesando: $accession" | tee -a "$logfile"
    echo "Inicio: $(date)" | tee -a "$logfile"
    echo "========================================" | tee -a "$logfile"

    exec_annotation \
        -p "$PROFILES" \
        -k "$KO_LIST" \
        --cpu 4 \
        --tmp-dir "$tmpdir" \
        -f detail-tsv \
        -o "$outfile" \
        "$faa" >> "$logfile" 2>&1

    status=$?

    if [ $status -eq 0 ]; then
        echo "OK: $accession" | tee -a "$logfile"
    else
        echo "ERROR: $accession (código $status)" | tee -a "$logfile"
    fi

    echo "Fin: $(date)" | tee -a "$logfile"
}

export -f run_kofam

export INPUT
export PROFILES
export KO_LIST
export OUTPUT
export TMP

find "$INPUT" -name "*.faa" -print0 | \
xargs -0 -n 1 -P 2 bash -c 'run_kofam "$0"'
