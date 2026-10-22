"""Fonctions utilitaires fournies avec le TP — NE PAS MODIFIER CE FICHIER.

Lecture du FASTA et de BLOSUM62. Les retours sont entierement types.
"""

import re


def read_fasta_sequences(file_path: str) -> dict[str, str]:
    """Lit un fichier FASTA et retourne identifiants -> sequences.

    Args:
        file_path: Chemin vers le fichier FASTA.

    Returns:
        Dictionnaire dont chaque cle est l'identifiant (sans le ``>``) et
        chaque valeur la sequence correspondante, dans l'ordre du fichier.
        Les sequences peuvent etre enroulees sur plusieurs lignes.
    """
    sequences: dict[str, str] = {}
    current_id: str | None = None
    current_seq: list[str] = []
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                # Entete sans identifiant : l'enregistrement est ignore.
                morceaux: list[str] = line[1:].split()
                current_id = morceaux[0] if morceaux else None
                current_seq = []
            elif current_id is not None:
                current_seq.append(line)
        if current_id is not None:
            sequences[current_id] = "".join(current_seq)
    return sequences


def read_alignment(file_path: str) -> dict[str, str]:
    """Lit un alignement multiple au format FASTA ou Clustal.

    Args:
        file_path: Chemin vers l'alignement. Un fichier qui commence par
            ``CLUSTAL`` (sortie par defaut de Clustal Omega, avec ou sans
            numeros de position) est lu comme Clustal, sinon comme FASTA.

    Returns:
        Dictionnaire identifiant -> ligne alignee (gaps ``-`` inclus).
    """
    with open(file_path, "r") as file:
        lignes: list[str] = file.read().splitlines()
    if not lignes or not lignes[0].startswith("CLUSTAL"):
        return read_fasta_sequences(file_path)
    alignement: dict[str, str] = {}
    for ligne in lignes[1:]:
        # Les lignes de conservation (espaces puis *:) sont ignorees.
        # Format : ID SEQUENCE [compteur] — le compteur de colonne eventuel
        # (chiffres en fin de ligne) ne fait pas partie de la sequence.
        m = re.match(r"^(\S+)\s+([A-Z\-]+)\s*\d*\s*$", ligne)
        if m:
            identifiant, sequence = m.group(1), m.group(2)
            alignement[identifiant] = alignement.get(identifiant, "") + sequence
    return alignement


def read_blosum62(file_path: str) -> dict[tuple[str, str], int]:
    """Lit BLOSUM62 et retourne un dictionnaire symetrique (aa1, aa2) -> score.

    Args:
        file_path: Chemin vers ``BLOSUM62.txt``.

    Returns:
        Dictionnaire indexe par les couples d'acides amines, dans les deux
        sens. La premiere ligne du fichier donne l'ordre des colonnes.
    """
    with open(file_path, "r") as file:
        lignes: list[list[str]] = [ligne.split() for ligne in file if ligne.strip()]
    acides: list[str] = lignes[0]
    matrice: dict[tuple[str, str], int] = {}
    for ligne in lignes[1:]:
        acide: str = ligne[0]
        for indice, valeur in enumerate(ligne[1:]):
            matrice[(acide, acides[indice])] = int(valeur)
    return matrice
