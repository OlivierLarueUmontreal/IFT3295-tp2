#!/bin/bash
# IFT3295 - TP2 - Script d'execution (fourni, NE PAS MODIFIER).
# Usage: ./run.sh [sequences.fasta] [BLOSUM62.txt]
set -uo pipefail

cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
FASTA="${1:-sequences.fasta}"
BLOSUM="${2:-BLOSUM62.txt}"

echo "=== Q1: matrice de similarite ==="
"$PY" main.py paires "$FASTA" "$BLOSUM" || exit 1

echo "=== Q2 a Q7: etoile, score SP, consensus, centre modifie, realignement ==="
"$PY" main.py etoile "$FASTA" "$BLOSUM" || exit 1

echo "=== Verification: score SP relu depuis alignement.fasta ==="
"$PY" main.py sp alignement.fasta "$BLOSUM" || exit 1

for f in similarite.csv alignement.fasta realignement.fasta; do
	if [ ! -s "$f" ]; then
		echo "ERREUR: fichier de sortie manquant ou vide: $f" >&2
		exit 1
	fi
done
echo "=== OK: tous les fichiers de sortie sont presents ==="
