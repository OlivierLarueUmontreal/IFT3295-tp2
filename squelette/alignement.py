"""IFT3295 - TP2 - Alignement global et matrice de similarite (Q1).

Ce fichier est A COMPLETER. Regles:
* Completez les fonctions marquees ``NotImplementedError # TODO``.
* Vous pouvez ajouter des imports de la bibliotheque standard et definir
  des fonctions auxiliaires dans ``alignement.py`` et ``etoile.py``.
* Ne modifiez ni les noms, ni les signatures, ni les annotations des
  fonctions fournies, ni les constantes.
* Annotez les parametres et retours de vos fonctions auxiliaires.
* Ne modifiez pas ``main.py``, ``utils.py`` ou ``run.sh``.
* Le programme principal est ``main.py`` (fourni): il appelle les fonctions
  de ce module, il ne faut donc rien executer ici directement.
"""

# Penalite affine d'un run de k gaps : -(OUVERTURE + k * EXTENSION).
OUVERTURE: int = 10
EXTENSION: int = 1


def alignement_global(
    seq1: str, seq2: str, matrice: dict[tuple[str, str], int]
) -> tuple[int, str, str]:
    """Calcule un alignement global optimal de deux sequences.

    Algorithme de programmation dynamique a penalites affines. Le score a maximiser utilise ``matrice`` pour une
    colonne acide amine/acide amine. Un run de ``k`` gaps coute
    ``-(OUVERTURE + k * EXTENSION)``, y compris en debut et en fin
    d'alignement. Les deux lignes retournees ont la meme longueur ; un gap
    est le caractere ``-``.

    Plusieurs alignements peuvent avoir le score optimal. Pour que toutes
    les copies retournent le meme, trancher les egalites ainsi pendant le
    retour arriere :
    * si plusieurs etats donnent le meme score, preferer, dans cet ordre :
      un acide amine de chaque sequence ; un acide amine de ``seq1`` face
      a un gap dans ``seq2`` ; un gap dans ``seq1`` face a un acide amine
      de ``seq2`` ;
    * dans un etat de gap, si prolonger un gap existant et ouvrir un nouveau
      gap donnent le meme score, preferer prolonger le gap existant.

    Args:
        seq1 (str): Premiere sequence, sans gap.
        seq2 (str): Deuxieme sequence, sans gap.
        matrice (dict[tuple[str, str], int]): Scores BLOSUM62.

    Returns:
        tuple[int, str, str]: Score optimal, ``seq1`` alignee, ``seq2``
        alignee.

    Examples:
        Avec BLOSUM62, ``alignement_global("AC", "A", matrice)`` a le score
        -7 : la colonne A/A vaut 4 et le gap de longueur 1 coute 11.
    """
    raise NotImplementedError  # TODO


def matrice_similarite(
    sequences: list[str], matrice: dict[tuple[str, str], int]
) -> list[list[int]]:
    """Construit la matrice des scores d'alignement global de toutes les paires.

    Args:
        sequences (list[str]): Sequences dans l'ordre du fichier FASTA.
        matrice (dict[tuple[str, str], int]): Scores BLOSUM62.

    Returns:
        list[list[int]]: Matrice carree. L'element ``[i][j]`` est le score
        de ``alignement_global(sequences[i], sequences[j], matrice)``.
        La diagonale contient des zeros.
    """
    raise NotImplementedError  # TODO: remplir la matrice
