# dataset_entrez.py
"""
dataset_entrez.py
=================

Baixa sequências do NCBI via Bio.Entrez.
Pergunta tudo via prompt: banco, termo, tamanho, quantidade.
"""

from typing import List, Tuple, Dict, Any
from Bio import Entrez, SeqIO
import random
import logging
from utils.config import ENTREZ_DEFAULTS

logger = logging.getLogger(__name__)

def fetch_dataset() -> Tuple[List[str], Dict[str, Any]]:
    """Baixa um dataset do NCBI e retorna as sequências e os parâmetros."""
    defaults = ENTREZ_DEFAULTS

    # Carregar email do config se não estiver definido
    if not getattr(Entrez, "email", None):
        email_input = input(f"Informe seu e-mail para o NCBI [{defaults['email']}]: ").strip()
        Entrez.email = email_input or defaults['email']
        if not Entrez.email:
            raise ValueError("É necessário fornecer um e-mail para acessar o NCBI")
        logger.debug(f"Entrez.email definido como '{Entrez.email}'")
    else:
        # Se já estiver definido, não faz nada
        pass

    # Carregar API key do config se não estiver definida
    if not getattr(Entrez, "api_key", None):
        api_key_input = input(f"Informe sua NCBI API key (opcional) [{defaults.get('api_key','')}]: ").strip()
        Entrez.api_key = api_key_input or defaults.get('api_key', '')
        if Entrez.api_key:
            logger.debug("Entrez.api_key definida.")
    else:
        # Se já estiver definida, não faz nada
        pass

    db_input = input(f"Base (nucleotide / protein) [{defaults['db']}]: ").strip()
    db = db_input or defaults['db']
    
    term_input = input(f"Termo de busca Entrez [{defaults['term']}]: ").strip()
    term = term_input or defaults['term']
    
    n_input = input(f"Quantos registros deseja baixar? [{defaults['n']}]: ").strip()
    n = int(n_input) if n_input else defaults['n']
    
    params = {'db': db, 'term': term, 'n': n}
    logger.debug(f"fetch_dataset chamado com db={db}, term='{term}', n={n}")

    rng = random.Random(0)

    logger.debug("Iniciando ESearch...")
    handle = Entrez.esearch(db=db, term=term, retmax=50000)
    search_result = Entrez.read(handle)
    
    if not isinstance(search_result, dict):
        raise TypeError(f"Resultado da busca Entrez não é um dicionário: {type(search_result)}")
        
    logger.debug(f"ESearch retornou {len(search_result.get('IdList', []))} IDs")

    ids = search_result.get("IdList", [])
    if not ids:
        raise RuntimeError("A busca não retornou nenhum registro. Verifique seu termo de busca.")

    if len(ids) < n:
        logger.warning(f"A busca retornou {len(ids)} IDs, menos que os {n} solicitados. Usando todos os IDs encontrados.")
        params['n'] = len(ids)
        sample_ids = ids
    else:
        sample_ids = rng.sample(ids, n)
    logger.debug(f"IDs amostrados (primeiros 5): {sample_ids[:5]} ...")

    logger.debug("Iniciando EFetch...")
    fetch_handle = Entrez.efetch(db=db, id=",".join(sample_ids),
                                rettype="fasta", retmode="text")
    records = list(SeqIO.parse(fetch_handle, "fasta"))
    logger.debug(f"{len(records)} seqüências baixadas")

    seqs = [str(rec.seq).upper() for rec in records]
    logger.info(f"Dataset Entrez obtido: n={len(seqs)}, L={len(seqs[0])}")
    return seqs, params
    seqs = [str(rec.seq).upper() for rec in records]
    print(f"Dataset Entrez obtido: n={len(seqs)}, L={len(seqs[0])}")
    return seqs, params
