"""IFT3295 - TP2 - Alignement par arbre etoile, score SP et consensus (Q2 a Q7).

Ce fichier est A COMPLETER. Regles:
* Completez les fonctions marquees ``NotImplementedError # TODO``.
* Vous pouvez ajouter des imports de la bibliotheque standard et definir
  des fonctions auxiliaires dans ``alignement.py`` et ``etoile.py``.
* Ne modifiez ni les noms, ni les signatures, ni les annotations des
  fonctions fournies, ni les constantes.
* Annotez les parametres et retours de vos fonctions auxiliaires.
* Ne modifiez pas ``main.py``, ``utils.py`` ou ``run.sh``.
* Reutilisez ``alignement_global`` de ``alignement.py`` pour les paires.
* Le programme principal est ``main.py`` (fourni): il appelle les fonctions
  de ce module, il ne faut donc rien executer ici directement.
"""

# Meme penalite affine que dans alignement.py. Ne pas modifier.
OUVERTURE: int = 10
EXTENSION: int = 1
# Les 20 acides amines standards, candidats de ``sequence_modifiee``.
ACIDES_AMINES: str = "ACDEFGHIKLMNPQRSTVWY"


def sequence_centrale(matrice: list[list[int]]) -> int:
    """Choisit la sequence centrale a partir de la matrice de similarite.

    Args:
        matrice (list[list[int]]): Matrice carree produite par
            ``matrice_similarite``. La diagonale vaut 0.

    Returns:
        int: Indice de la sequence dont la somme des scores avec les autres
        est maximale. En cas d'egalite, le plus petit indice.
    """
    raise NotImplementedError  # TODO


def alignement_etoile(
    sequences: list[str], matrice: dict[tuple[str, str], int]
) -> list[str]:
    """Construit l'alignement multiple par la methode de la sequence centrale.

    Aligner chaque sequence autre que le centre avec le centre, par
    ``alignement_global(centre, autre, matrice)`` : le centre est toujours
    le premier argument. Ajouter les sequences dans l'ordre de
    ``sequences``, en sautant le centre. Dans le resultat, la paire
    (centre, autre), une fois retirees les colonnes gap/gap, doit etre
    exactement l'alignement retourne par ``alignement_global``.

    Args:
        sequences (list[str]): Sequences dans l'ordre du fichier FASTA.
        matrice (dict[tuple[str, str], int]): Scores BLOSUM62.

    Returns:
        list[str]: Une ligne par sequence, dans le meme ordre que
        ``sequences``, gaps ``-`` inclus. Toutes les lignes ont la meme
        longueur. Retirer les gaps d'une ligne redonne la sequence d'origine.
    """
    raise NotImplementedError  # TODO: fusion des alignements avec le centre


def score_sp(alignement: list[str], matrice: dict[tuple[str, str], int]) -> int:
    """Calcule le score sum-of-pairs d'un alignement multiple.

    Somme, sur chaque paire de lignes, du score de l'alignement induit.
    Une colonne ou les deux symboles sont des gaps est ignoree. Une colonne
    acide amine/acide amine vaut ``matrice[(a, b)]``. Un run de ``k`` gaps
    face a des acides amines coute ``-(OUVERTURE + k * EXTENSION)``, extremites comprises.
    Un gap/gap au milieu d'un run ne le coupe pas et n'augmente pas ``k``.
    Un run est dans une seule ligne : un gap dans une ligne suivi d'un gap
    dans l'autre forme deux runs. Exemple : ``"A---CG-"`` / ``"AW-Y-GK"``
    vaut 4 - 12 - 11 + 6 - 11 = -24.

    Args:
        alignement (list[str]): Lignes de meme longueur, gaps ``-``.
        matrice (dict[tuple[str, str], int]): Scores BLOSUM62.

    Returns:
        int: Score SP.
    """
    raise NotImplementedError  # TODO: somme des paires


def sequence_consensus(alignement: list[str]) -> str:
    """Construit la sequence consensus d'un alignement multiple.

    Args:
        alignement (list[str]): Lignes de meme longueur, gaps ``-``.

    Returns:
        str: Une chaine de la longueur de l'alignement. Chaque position est
        l'acide amine le plus frequent de la colonne, gaps exclus (un gap
        n'est pas compte). En cas d'egalite, le premier acide amine dans
        l'ordre de la matrice BLOSUM62 (ligne d'en-tete de
        ``BLOSUM62.txt``). Une colonne sans acide amine donne ``-``.
    """
    raise NotImplementedError  # TODO: majorite par colonne


def sequence_modifiee(
    alignement: list[str], centre: int, matrice: dict[tuple[str, str], int]
) -> str:
    """Construit S*', une copie modifiee de la sequence centrale (Q6).

    La sequence centrale elle-meme n'est pas modifiee.

    Pour chaque colonne ou la ligne ``centre`` porte un acide amine,
    choisir l'acide amine ``x`` de ``ACIDES_AMINES`` qui maximise la somme
    des ``matrice[(x, c)]``, ou ``c`` parcourt les acides amines de la
    colonne (y compris celui du centre, gaps ignores). En cas d'egalite, la
    plus petite lettre. Les colonnes ou le centre a un gap sont ignorees.

    Args:
        alignement (list[str]): Alignement produit par ``alignement_etoile``.
        centre (int): Indice de la sequence centrale.
        matrice (dict[tuple[str, str], int]): Scores BLOSUM62.

    Returns:
        str: S*', sans gap, de la meme longueur que la sequence centrale.
    """
    raise NotImplementedError  # TODO: meilleur acide amine par colonne


def realigner_sur_centre(
    centre_impose: str, sequences: list[str], matrice: dict[tuple[str, str], int]
) -> list[str]:
    """Refait l'alignement etoile autour d'un centre impose (Q7).

    Aligner chacune des ``sequences``, dans l'ordre, avec
    ``alignement_global(centre_impose, sequence, matrice)`` et fusionner
    comme dans ``alignement_etoile``. Retirer ensuite la ligne du centre
    impose, puis les colonnes qui ne contiennent plus que des gaps.

    Args:
        centre_impose (str): Sequence centrale imposee, sans gap (S*').
        sequences (list[str]): Sequences dans l'ordre du fichier FASTA.
        matrice (dict[tuple[str, str], int]): Scores BLOSUM62.

    Returns:
        list[str]: Une ligne par sequence de ``sequences``, dans le meme
        ordre, gaps ``-`` inclus. Toutes les lignes ont la meme longueur.
    """
    raise NotImplementedError  # TODO: etoile avec centre impose, puis retrait
