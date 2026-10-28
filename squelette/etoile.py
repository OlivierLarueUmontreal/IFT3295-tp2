import alignement
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
    scores = []
    max_score = float("-inf")
    index = 0 # une seul seq, elle est la seq etoile
    for i in range(len(matrice)): 
        row = matrice[i] # chq row est une seq
        score = sum(row)
        if score > max_score: 
            max_score = score # cas d'égalité, on commence par l'indice le plus petit
            index = i

    return index #attention le nom de la seq serait donc i + 1

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
    
    # ceci est le nom de la seq 
    index_centrale = sequence_centrale(
        alignement.matrice_similarite(sequences, matrice))

    centre = sequences[index_centrale]
    new_centre = ""
    tmp_centre = centre # tmp
    aln_mult = []
    for i in range(len(sequences)):
        autre = sequences[i]
        
        aln = ""
        if i == index_centrale: 
            aln_mult.append(tmp_centre)
            continue
        else : 
            res = alignement.alignement_global(centre, autre, matrice) # centre reste tjrs le même
            aln = res[2]
            new_centre = res[1]
        
        # fusioner avec l'alignement 
        # comparer les centres : 
        # si centre == _ and new_centre != _ : cor_centre += _; cor_aln += _;
        # si centre != _ and new_centre == _ : cor_centre += _; for i in aln_mult += 1;
        ci = 0 #cor centre index
        ni = 0 # new centre index        

        cor_centre = ""
        cor_aln = aln
        while ci < len(tmp_centre) or ni < len(new_centre): 
            c = tmp_centre[ci] if ci < len(tmp_centre) else None
            n = new_centre[ni] if ni < len(new_centre) else None
            pos = len(cor_centre) # pos de centre qu'on construit
            if c == "-" and n != "-" : 
                cor_centre += "-"
                cor_aln = cor_aln[:pos] + "-" + cor_aln[pos:]
                ci += 1

            elif c != "-" and n == "-" :
                cor_centre += "-"
                for iseq in range(len(aln_mult)):
                    seq = aln_mult[iseq] 
                    aln_mult[iseq] = seq[:pos] + "-" + seq[pos:]
                ni += 1

            else: 
                cor_centre += c
                ci += 1
                ni +=1
        tmp_centre = cor_centre #corriger le centre
        aln_mult.append(cor_aln)
    
    aln_mult[index_centrale] = tmp_centre

    return aln_mult


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
    def kgaps(k):
        return -(OUVERTURE + k * EXTENSION)
    # on vx faire des combinaisons de paires, garder la somme
    scores = []

    for i in range(len(alignement)):
        for u in range(len(alignement)):
            if i == u: continue
            if i < u : continue
            seq1 = alignement[i]
            seq2 = alignement[u]

            seqind = 0
            while seqind < len(seq1):
            # enlever les matching gaps
                if seq1[seqind] == "-" and seq2[seqind] == "-":
                    seq1 = seq1[:seqind] + seq1[seqind+1:]
                    seq2 = seq2[:seqind] + seq2[seqind+1:]
                else : seqind += 1
            
            score = 0
            gap = False
            gap_count = 0 # accumuler
            for seqind in range(len(seq1)):
                if seq1[seqind] != "-" and seq2[seqind] != "-":
                    if gap: 
                        score += kgaps(gap_count)
                        gap_count = 0
                        gap = False
                    score += matrice[seq1[seqind],seq2[seqind]]
                if seq1[seqind] == "-" or seq2[seqind] == "-":
                    # si on avait un trou ouvert
                    if gap == False:
                        gap = True
                        gap_count += 1
                    # verifier le gap continue est dans le même strand, sinon recommencer
                    else:
                        # A-
                        # -S 
                        if (seq1[seqind] == "-" and seq1[seqind-1] == "-") or (seq2[seqind] == "-" and seq2[seqind-1] == "-"): 
                            gap_count += 1
                            continue
                        else: 
                            # changer de brin 
                            score += kgaps(gap_count)
                            gap_count = 1
            if gap: score += kgaps(gap_count)
                            
            scores.append(score)
    return(sum(scores))


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
    def get_aa_freq(aa_table):
        # freq = naa/N
        if len(aa_table) == aa_table.count("-"): return "-"

        N = len("".join(("-").split("".join(aa_table))))
        max_freq = 0
        max_aa = ""
        for aa in ACIDES_AMINES:
            f = aa_table.count(aa)/N
            if f > max_freq:
                max_freq = f
                max_aa = aa
        return max_aa

    col_dict = {}
    seq_consensus = ""
    for i in range(len(alignement[0])): # longeur de tt sequence
        col_dict[i] = []
        for seq in alignement:
            col = col_dict[i]
            col.append(seq[i])
            col_dict[i] = col
        seq_consensus += get_aa_freq(col_dict[i])

    return seq_consensus


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
    S_start_prime = alignement[centre]
    col_dict = {}

    for i in range(len(S_start_prime)): # longeur de tt sequence
        col_dict[i] = []
        for seq in alignement:
            col = col_dict[i]
            col.append(seq[i])
            col_dict[i] = col

        x_max = ""
        max_score = float("-inf")        
        # determiner l'aa à l'index i de S *'
        for x in col_dict[i]: #tt aa de la colonne
            if x == "-": continue # skip

            x_score = 0
            for aa in ACIDES_AMINES:
                x_score += matrice[x,aa]
            # la somme est plus grande que le dernier x evalué?
            if x_score > max_score :
                max_score = x_score
                x_max = x
        
        # remplacer par le nouveau x
        S_start_prime = S_start_prime[:i] + x_max + S_start_prime[i+1:] 
    S_start_prime = S_start_prime.replace("-","")
    print(alignement[centre])
    
    print(S_start_prime.replace("-",""))
    print(len(alignement[centre]))
    
    print(len(S_start_prime.replace("-","")))

    return S_start_prime.replace("-","")



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
