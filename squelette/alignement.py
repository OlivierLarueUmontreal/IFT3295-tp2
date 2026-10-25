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
    m: int = len(seq1)
    n: int = len(seq2)

    # inspirée notes de cours (2.1-Alignement-Prog-Dyn)
    # E(i,j) : meilleur alignement qui se termine par (-, wj)
    # F(i,j) : meilleur alignement qui se termine par (vi, -)
    # G(i,j) : meilleur alignement qui se termine par (vi, wj)
    # V(i,j) = max[E(i,j), F(i,j), G(i,j)].
    # init tout a -infini
    V: list[list[float]] = [[float("-inf")] * (n + 1) for _ in range(m + 1)]
    E: list[list[float]] = [[float("-inf")] * (n + 1) for _ in range(m + 1)]
    F: list[list[float]] = [[float("-inf")] * (n + 1) for _ in range(m + 1)]
    G: list[list[float]] = [[float("-inf")] * (n + 1) for _ in range(m + 1)]

    V[0][0] = 0
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0 and j == 0: continue

            if i > 0 and j > 0:
                G[i][j] = V[i - 1][j - 1] + matrice[(seq1[i - 1], seq2[j - 1])]

            if j > 0:
                E[i][j] = max(E[i][j - 1], V[i][j - 1] - OUVERTURE) - EXTENSION

            if i > 0:
                F[i][j] = max(F[i - 1][j], V[i - 1][j] - OUVERTURE) - EXTENSION

            V[i][j] = max(E[i][j], F[i][j], G[i][j])

   
   # backtracking, on part de la fin en bas a droite du tableau et on remonte (V[m][n])
    ligne1: list[str] = []
    ligne2: list[str] = []
    i, j = m, n # pour garder les valeurs a m et n lors du retour

    while i > 0 or j > 0:
        if V[i][j] == G[i][j]:
            # G(i,j) = V(i-1,j-1) + delta : colonne (vi, wj)
            ligne1.append(seq1[i - 1])
            ligne2.append(seq2[j - 1])
            i = i - 1
            j = j - 1
        elif V[i][j] == F[i][j]:
            # Il y a un gap, on reste dans cette table pour tout la longueur du Gap
            prolonger_gap: bool = True
            while prolonger_gap:
                ligne1.append(seq1[i - 1])
                ligne2.append("-")
                prolonger_gap = F[i][j] == F[i - 1][j] - EXTENSION
                i -= 1
        else:
            # On continue dans le Gap
            prolonger_gap = True
            while prolonger_gap:
                ligne1.append("-")
                ligne2.append(seq2[j - 1])
                prolonger_gap = E[i][j] == E[i][j - 1] - EXTENSION
                j -= 1

    # Les colonnes ont ete construites de la fin vers le debut.
    return int(V[m][n]), "".join(reversed(ligne1)), "".join(reversed(ligne2))

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

    nombre: int = len(sequences)
    similarite: list[list[int]] = [[0] * nombre for _ in range(nombre)]
    for i in range(nombre):
        for j in range(nombre):
            if i != j:
                similarite[i][j] = alignement_global(sequences[i], sequences[j], matrice)[0]

    return similarite
