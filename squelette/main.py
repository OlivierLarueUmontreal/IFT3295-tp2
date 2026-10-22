"""IFT3295 - TP2 - Programme principal (fourni, NE PAS MODIFIER).

Guide pour les etudiants
-----------------------
Ce fichier est entierement fourni : vous n'avez rien a y programmer.
Completez les 8 fonctions TODO de ``alignement.py`` (Q1) et ``etoile.py``
(Q2 a Q7). Il n'y a aucun parseur ni fichier de sortie a ecrire vous-memes :
* ``utils.read_fasta_sequences`` lit les identifiants et sequences FASTA ;
* ``utils.read_blosum62`` lit les scores BLOSUM62 ;
* ``utils.read_alignment`` lit un alignement FASTA ou Clustal ;
* ce programme ecrit la matrice CSV et les alignements FASTA.
Vos fonctions recoivent donc deja des chaines, listes et dictionnaires,
pas des chemins de fichiers. Par exemple, une substitution se consulte
avec ``matrice[("A", "W")]`` ; les gaps sont traites par vos algorithmes.

Utilisation
-----------
Depuis le dossier du squelette, ``bash run.sh`` lance tout le pipeline
avec les fichiers fournis, sans argument. Pour lancer une seule commande,
utilisez les exemples ci-dessous : les chemins FASTA et BLOSUM62 sont
obligatoires avec ``main.py``. Les chemins relatifs et les sorties de
``main.py`` sont relatifs au dossier depuis lequel vous lancez Python ;
``run.sh`` se place automatiquement dans le dossier du squelette.
Les fichiers de sortie existants sont remplaces.

Pour obtenir l'aide : ``python3 main.py --help`` ou, par exemple,
``python3 main.py paires --help``. Une erreur ``NotImplementedError``
est normale tant qu'une fonction appelee n'est pas completee.
La commande ``etoile`` necessite les fonctions de Q1 a Q7 ; ``sp``
necessite seulement ``score_sp`` (Q4), en plus des lecteurs fournis.
Les questions sur l'ARN, les illustrations des etapes, le bonus et BLAST
se traitent dans le rapport : ce programme ne les automatise pas.

Orchestration du TP: lecture des donnees, appels aux fonctions a completer
dans ``alignement.py`` et ``etoile.py``, ecriture des fichiers de sortie et
affichage des resultats. Les formats de sortie ci-dessous constituent le
contrat de correction: ils sont identiques pour toutes les copies.

Commandes:
    python3 main.py paires sequences.fasta BLOSUM62.txt
        1. ecrit ``similarite.csv`` (n lignes, n entiers separes par des
           virgules, sans entete; la ligne i correspond a la sequence i
           du FASTA; la diagonale vaut 0);
        2. affiche la matrice avec les identifiants FASTA.

    python3 main.py etoile sequences.fasta BLOSUM62.txt
        1. ecrit ``alignement.fasta`` (Q3) et ``realignement.fasta`` (Q7),
           un enregistrement par sequence, dans l'ordre du FASTA, gaps
           inclus;
        2. affiche:
               Centre: <identifiant FASTA>
               Score SP: <entier>
               Consensus: <chaine>
               Centre modifie: <chaine>
               Changements: <entier>
               Score SP realigne: <entier>

    python3 main.py sp alignement.fasta BLOSUM62.txt
        affiche le score SP d'un alignement multiple quelconque, au format
        FASTA ou Clustal (par exemple la sortie par defaut de Clustal Omega):
               Score SP: <entier>
"""

import argparse

from alignement import matrice_similarite
from etoile import (
    alignement_etoile,
    realigner_sur_centre,
    score_sp,
    sequence_centrale,
    sequence_consensus,
    sequence_modifiee,
)
from utils import read_alignment, read_blosum62, read_fasta_sequences

PAIRES: str = "paires"
ETOILE: str = "etoile"
SP: str = "sp"
GAP: str = "-"


def _lire(
    fasta: str, blosum: str
) -> tuple[list[str], list[str], dict[tuple[str, str], int]]:
    """Lit les deux fichiers avec les lecteurs fournis dans utils.py.

    Retourne les identifiants et sequences dans le meme ordre (celui du
    FASTA), puis le dictionnaire de scores indexe par couples de lettres.
    Exige au moins deux sequences pour les commandes paires et etoile.
    """
    sequences: dict[str, str] = read_fasta_sequences(fasta)
    if len(sequences) < 2:
        raise SystemExit(f"Erreur: {fasta} doit contenir au moins 2 sequences")
    identifiants: list[str] = list(sequences)
    return identifiants, [sequences[nom] for nom in identifiants], read_blosum62(blosum)


def _ecrire_matrice(scores: list[list[int]]) -> None:
    """Ecrit la matrice n x n dans ``similarite.csv`` (sans entete)."""
    with open("similarite.csv", "w") as fichier_csv:
        for ligne in scores:
            fichier_csv.write(",".join(str(valeur) for valeur in ligne) + "\n")


def _ecrire_alignement(
    identifiants: list[str], alignement: list[str], chemin: str
) -> None:
    """Ecrit un alignement multiple au format FASTA dans ``chemin``."""
    with open(chemin, "w") as fichier_fasta:
        for nom, ligne in zip(identifiants, alignement, strict=True):
            fichier_fasta.write(f">{nom}\n{ligne}\n")


def _verifier_matrice(scores: list[list[int]], n: int) -> None:
    """Verifie la forme de la matrice de similarite."""
    if len(scores) != n or any(len(ligne) != n for ligne in scores):
        raise SystemExit("Erreur: la matrice de similarite doit etre carree")
    if any(scores[i][i] != 0 for i in range(n)):
        raise SystemExit("Erreur: la diagonale de la matrice doit valoir 0")


def _verifier_alignement(sequences: list[str], alignement: list[str]) -> None:
    """Controle le format renvoye par les fonctions etudiantes.

    Il faut une ligne par sequence, dans l'ordre d'entree, de meme longueur
    et sans changement de lettres. Ce controle ne prouve pas l'optimalite
    ni la conservation des alignements centre-sequence de la question 3.
    """
    if len(alignement) != len(sequences):
        raise SystemExit("Erreur: l'alignement n'a pas une ligne par sequence")
    longueur: int = len(alignement[0])
    for sequence, ligne in zip(sequences, alignement, strict=True):
        if len(ligne) != longueur:
            raise SystemExit(
                "Erreur: les lignes de l'alignement n'ont pas la meme longueur"
            )
        if ligne.replace(GAP, "") != sequence:
            raise SystemExit("Erreur: une ligne ne redonne pas sa sequence")


def _commande_paires(args: argparse.Namespace) -> None:
    """Calcule et affiche la matrice de similarite."""
    identifiants, sequences, blosum = _lire(args.fasta, args.blosum)
    scores: list[list[int]] = matrice_similarite(sequences, blosum)
    _verifier_matrice(scores, len(sequences))
    _ecrire_matrice(scores)
    entete: str = "".join(f"{nom:>8}" for nom in identifiants)
    print(f"{'':8}{entete}")
    for nom, ligne in zip(identifiants, scores, strict=True):
        print(f"{nom:<8}" + "".join(f"{valeur:8d}" for valeur in ligne))


def _commande_etoile(args: argparse.Namespace) -> None:
    """Execute Q2 a Q7, en recalculant la matrice Q1.

    Le centre est un indice dans l'ordre FASTA, pas un identifiant.
    S*' est une copie sans gap : les sequences originales restent celles
    de l'entree. Le nombre de substitutions (Q6) est calcule ici, et non
    dans sequence_modifiee. Les deux alignements sont ecrits en FASTA.
    """
    identifiants, sequences, blosum = _lire(args.fasta, args.blosum)
    # Q1 : scores deux a deux ; Q2 : choix du centre.
    scores: list[list[int]] = matrice_similarite(sequences, blosum)
    _verifier_matrice(scores, len(sequences))
    centre: int = sequence_centrale(scores)
    if not 0 <= centre < len(sequences):
        raise SystemExit("Erreur: indice de sequence centrale hors limites")
    # Q3 : alignement etoile ; Q4 : score SP ; Q5 : consensus.
    alignement: list[str] = alignement_etoile(sequences, blosum)
    _verifier_alignement(sequences, alignement)
    print(f"Centre: {identifiants[centre]}")
    print(f"Score SP: {score_sp(alignement, blosum)}")
    print(f"Consensus: {sequence_consensus(alignement)}")
    _ecrire_alignement(identifiants, alignement, "alignement.fasta")

    # Q6 : copie modifiee du centre et comptage des substitutions.
    modifiee: str = sequence_modifiee(alignement, centre, blosum)
    if len(modifiee) != len(sequences[centre]) or GAP in modifiee:
        raise SystemExit("Erreur: S*' doit avoir la longueur du centre, sans gap")
    changements: int = sum(
        a != b for a, b in zip(sequences[centre], modifiee, strict=True)
    )
    # Q7 : realigner les sequences originales, puis calculer leur score SP.
    realignement: list[str] = realigner_sur_centre(modifiee, sequences, blosum)
    _verifier_alignement(sequences, realignement)
    print(f"Centre modifie: {modifiee}")
    print(f"Changements: {changements}")
    print(f"Score SP realigne: {score_sp(realignement, blosum)}")
    _ecrire_alignement(identifiants, realignement, "realignement.fasta")


def _commande_sp(args: argparse.Namespace) -> None:
    """Affiche le score SP d'un alignement multiple (FASTA ou Clustal)."""
    alignement: list[str] = list(read_alignment(args.alignement).values())
    if len(alignement) < 2 or len({len(ligne) for ligne in alignement}) != 1:
        raise SystemExit("Erreur: il faut au moins 2 lignes de meme longueur")
    print(f"Score SP: {score_sp(alignement, read_blosum62(args.blosum))}")


def _construire_parseur() -> argparse.ArgumentParser:
    """Declare les commandes et leurs arguments avec argparse (fourni).

    Ce parseur lit les arguments du terminal, pas les sequences. Les
    lecteurs de fichiers sont dans utils.py. Chaque sous-commande fixe
    args.action a la fonction _commande_* qui sera appelee par main().
    """
    parseur: argparse.ArgumentParser = argparse.ArgumentParser(
        prog="IFT3295-TP2",
        description="Alignement multiple par arbre etoile (voir docstring du module).",
    )
    sous_commandes = parseur.add_subparsers(dest="commande", required=True)

    paires: argparse.ArgumentParser = sous_commandes.add_parser(
        PAIRES, help="matrice des scores de similarite"
    )
    paires.add_argument("fasta", help="fichier FASTA des sequences")
    paires.add_argument("blosum", help="fichier BLOSUM62")
    paires.set_defaults(action=_commande_paires)

    etoile: argparse.ArgumentParser = sous_commandes.add_parser(
        ETOILE, help="alignement par sequence centrale"
    )
    etoile.add_argument("fasta", help="fichier FASTA des sequences")
    etoile.add_argument("blosum", help="fichier BLOSUM62")
    etoile.set_defaults(action=_commande_etoile)

    sp: argparse.ArgumentParser = sous_commandes.add_parser(
        SP, help="score SP d'un alignement multiple (FASTA ou Clustal)"
    )
    sp.add_argument("alignement", help="alignement FASTA ou Clustal (gaps '-')")
    sp.add_argument("blosum", help="fichier BLOSUM62")
    sp.set_defaults(action=_commande_sp)

    return parseur


def main() -> None:
    """Lit les arguments puis execute la commande choisie, sans code a completer."""
    parseur: argparse.ArgumentParser = _construire_parseur()
    args: argparse.Namespace = parseur.parse_args()
    args.action(args)


if __name__ == "__main__":
    main()
