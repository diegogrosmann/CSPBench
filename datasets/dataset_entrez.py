"""
Download de sequências do NCBI via Bio.Entrez.

Funções:
    fetch_dataset(): Baixa dataset do NCBI e retorna sequências e parâmetros.
"""

# dataset_entrez.py
"""
dataset_entrez.py
=================

Baixa sequências do NCBI via Bio.Entrez.
Pergunta tudo via prompt: banco, termo, tamanho, quantidade.
"""

from typing import List, Tuple, Dict, Any
import logging
import random
from utils.config import ENTREZ_DEFAULTS, safe_input
from src.console_manager import console

try:
    from Bio import Entrez, SeqIO
except ImportError:
    raise ImportError("Biopython não encontrado. Instale com: pip install biopython")

logger = logging.getLogger(__name__)

def fetch_dataset() -> Tuple[List[str], Dict[str, Any]]:
    """Baixa um dataset do NCBI e retorna as sequências e os parâmetros."""
    defaults = ENTREZ_DEFAULTS

    # Carregar email do config se não estiver definido
    if not getattr(Entrez, "email", None):
        email_input = safe_input(f"Informe seu e-mail para o NCBI [{defaults['email']}]: ")
        Entrez.email = email_input or defaults['email']
        if not Entrez.email:
            raise ValueError("É necessário fornecer um e-mail para acessar o NCBI")
        logger.debug(f"Entrez.email definido como '{Entrez.email}'")

    # Carregar API key do config se não estiver definida
    if not getattr(Entrez, "api_key", None):
        api_key_input = safe_input(f"Informe sua NCBI API key (opcional) [{defaults.get('api_key','')}]: ")
        Entrez.api_key = api_key_input or defaults.get('api_key', '')
        if Entrez.api_key:
            logger.debug("Entrez.api_key definida.")

    db_input = safe_input(f"Base (nucleotide / protein) [{defaults['db']}]: ")
    db = db_input or defaults['db']
    
    term_input = safe_input(f"Termo de busca Entrez [{defaults['term']}]: ")
    term = term_input or defaults['term']
    
    n_input = safe_input(f"Quantos registros deseja baixar? [{defaults['n']}]: ")
    n = int(n_input) if n_input else defaults['n']
    
    params = {'db': db, 'term': term, 'n': n}
    logger.debug(f"fetch_dataset chamado com db={db}, term='{term}', n={n}")

    rng = random.Random(0)

    logger.debug("Iniciando ESearch...")
    handle = Entrez.esearch(db=db, term=term, retmax=50000)
    search_result = Entrez.read(handle)
    handle.close()
    
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
    fetch_handle.close()
    logger.debug(f"{len(records)} sequências baixadas")

    seqs = [str(rec.seq).upper() for rec in records]
    logger.info(f"Dataset Entrez obtido: n={len(seqs)}, L={len(seqs[0])}")
    return seqs, params

def fetch_dataset_with_params(params: Dict[str, Any]) -> Tuple[List[str], Dict[str, Any]]:
    """
    Busca dataset do NCBI com parâmetros específicos.
    
    Args:
        params: Dicionário com parâmetros (email, db, term, n, api_key)
        
    Returns:
        Tupla (sequências, parâmetros_usados)
    """
    # Merge com defaults
    merged_params = {**ENTREZ_DEFAULTS}
    merged_params.update(params)
    
    email = merged_params['email']
    db = merged_params['db']
    term = merged_params['term']
    n = merged_params['n']
    api_key = merged_params.get('api_key')
    
    console.print(f"Buscando no NCBI: db={db}, term='{term}', n={n}")
    
    # Configurar Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    
    try:
        # Buscar IDs
        search_handle = Entrez.esearch(db=db, term=term, retmax=n*2)  # Busca extra para filtrar
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        # Verificar se search_results é válido e converter para dict se necessário
        if not search_results:
            raise ValueError(f"Resposta vazia do NCBI para o termo: {term}")
        
        # Converter para dict se não for
        if not isinstance(search_results, dict):
            raise ValueError(f"Formato de resposta inesperado do NCBI: {type(search_results)}")
        
        # Verificar se tem IdList usando get() que é mais seguro
        ids_list = search_results.get('IdList')
        if not ids_list:
            raise ValueError(f"Nenhum resultado encontrado para o termo: {term}")
        
        # Converter para lista de strings se necessário
        ids = [str(id_item) for id_item in ids_list]
        
        console.print(f"Encontrados {len(ids)} registros, baixando sequências...")
        
        # Buscar sequências
        fetch_handle = Entrez.efetch(db=db, id=ids, rettype="fasta", retmode="text")
        fasta_data = fetch_handle.read()
        fetch_handle.close()
        
    except Exception as e:
        raise ValueError(f"Erro ao acessar NCBI: {e}")
    
    # Processar FASTA
    sequences = []
    current_seq = ""
    
    for line in fasta_data.split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_seq:
                clean_seq = ''.join(c.upper() for c in current_seq if c.isalpha())
                if clean_seq:
                    sequences.append(clean_seq)
                    if len(sequences) >= n:
                        break
            current_seq = ""
        else:
            current_seq += line
    
    # Adicionar última sequência
    if current_seq and len(sequences) < n:
        clean_seq = ''.join(c.upper() for c in current_seq if c.isalpha())
        if clean_seq:
            sequences.append(clean_seq)
    
    if not sequences:
        raise ValueError("Nenhuma sequência válida encontrada")
    
    # Filtrar por comprimento uniforme
    L = len(sequences[0])
    filtered_sequences = []
    for seq in sequences:
        if len(seq) == L:
            filtered_sequences.append(seq)
    
    if len(filtered_sequences) < min(n, len(sequences) // 2):
        raise ValueError(f"Muito poucas sequências de comprimento uniforme: {len(filtered_sequences)}")
    
    # Limitar ao número solicitado
    final_sequences = filtered_sequences[:n]
    
    used_params = {
        'email': email,
        'db': db,
        'term': term,
        'n_requested': n,
        'n_obtained': len(final_sequences),
        'L': L,
        'api_key_used': bool(api_key),
        'total_found': len(sequences),
        'uniform_length_count': len(filtered_sequences)
    }
    
    console.print(f"✓ {len(final_sequences)} sequências obtidas (L={L})")
    
    return final_sequences, used_params
