"""
Download de sequências do NCBI via Bio.Entrez.

Funções:
    fetch_dataset(): Baixa dataset do NCBI e retorna sequências e parâmetros.
"""

import logging
import random
from typing import Any

from src.core.io import validate_sequences
from src.ui.cli.console_manager import console
from src.utils.config import ENTREZ_DEFAULTS, safe_input

try:
    from Bio import Entrez, SeqIO
except ImportError:
    raise ImportError("Biopython não encontrado. Instale com: pip install biopython")

logger = logging.getLogger(__name__)


def fetch_dataset() -> tuple[list[str], dict[str, Any]]:
    """Baixa um dataset do NCBI e retorna as sequências e os parâmetros."""
    defaults = ENTREZ_DEFAULTS

    # Carregar email do config se não estiver definido
    if not getattr(Entrez, "email", None):
        email_input = safe_input(f"Informe seu e-mail para o NCBI [{defaults['email']}]: ")
        Entrez.email = email_input or defaults["email"]
        if not Entrez.email:
            raise ValueError("É necessário fornecer um e-mail para acessar o NCBI")
        logger.debug(f"Entrez.email definido como '{Entrez.email}'")

    # Carregar API key do config se não estiver definida
    if not getattr(Entrez, "api_key", None):
        api_key_input = safe_input(f"Informe sua NCBI API key (opcional) [{defaults.get('api_key','')}]: ")
        Entrez.api_key = api_key_input or defaults.get("api_key", "")
        if Entrez.api_key:
            logger.debug("Entrez.api_key definida.")

    db_input = safe_input(f"Base (nucleotide / protein) [{defaults['db']}]: ")
    db = db_input or defaults["db"]

    term_input = safe_input(f"Termo de busca Entrez [{defaults['term']}]: ")
    term = term_input or defaults["term"]

    n_input = safe_input(f"Quantos registros deseja baixar? [{defaults['n']}]: ")
    n = int(n_input) if n_input else defaults["n"]

    params = {"db": db, "term": term, "n": n}
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
        logger.warning(
            f"A busca retornou {len(ids)} IDs, menos que os {n} solicitados. Usando todos os IDs encontrados."
        )
        params["n"] = len(ids)
        sample_ids = ids
    else:
        sample_ids = rng.sample(ids, n)
    logger.debug(f"IDs amostrados (primeiros 5): {sample_ids[:5]} ...")

    logger.debug("Iniciando EFetch...")
    fetch_handle = Entrez.efetch(db=db, id=",".join(sample_ids), rettype="fasta", retmode="text")
    records = list(SeqIO.parse(fetch_handle, "fasta"))
    fetch_handle.close()
    logger.debug(f"{len(records)} sequências baixadas")

    seqs = [str(rec.seq).upper() for rec in records]
    logger.info(f"Dataset Entrez obtido: n={len(seqs)}, L={len(seqs[0])}")
    return seqs, params


def fetch_dataset_with_params(
    params: dict[str, Any],
) -> tuple[list[str], dict[str, Any]]:
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

    email = merged_params["email"]
    db = merged_params["db"]
    term = merged_params["term"]
    n = merged_params["n"]
    api_key = merged_params.get("api_key")

    console.print(f"Buscando no NCBI: db={db}, term='{term}', n={n}")

    # Configurar Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    try:
        # Buscar IDs
        search_handle = Entrez.esearch(db=db, term=term, retmax=n * 2)  # Busca extra para filtrar
        search_results = Entrez.read(search_handle)
        search_handle.close()

        # Verificar se search_results é válido e converter para dict se necessário
        if not search_results:
            raise ValueError(f"Resposta vazia do NCBI para o termo: {term}")

        # Converter para dict se não for
        if not isinstance(search_results, dict):
            raise ValueError(f"Formato de resposta inesperado do NCBI: {type(search_results)}")

        # Verificar se tem IdList usando get() que é mais seguro
        ids_list = search_results.get("IdList")
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

    # Processar FASTA usando método unificado
    sequences = _parse_fasta_data(fasta_data)

    if not sequences:
        raise ValueError("Nenhuma sequência válida encontrada")

    # Usar validação unificada
    validation = validate_sequences(sequences)

    if not validation["valid"]:
        for error in validation["errors"]:
            logger.warning(f"Validação NCBI: {error}")

    # Filtrar por comprimento uniforme se necessário
    if not validation["uniform_length"]:
        L = len(sequences[0])
        filtered_sequences = []
        for seq in sequences:
            if len(seq) == L:
                filtered_sequences.append(seq)

        if len(filtered_sequences) < min(n, len(sequences) // 2):
            raise ValueError(f"Muito poucas sequências de comprimento uniforme: {len(filtered_sequences)}")
        sequences = filtered_sequences

    # Limitar ao número solicitado
    final_sequences = sequences[:n]

    used_params = {
        "email": email,
        "db": db,
        "term": term,
        "n_requested": n,
        "n_obtained": len(final_sequences),
        "L": len(final_sequences[0]) if final_sequences else 0,
        "api_key_used": bool(api_key),
        "total_found": validation["count"],
        "uniform_length_count": len(sequences),
        "alphabet": validation["alphabet"],
    }

    console.print(
        f"✓ {len(final_sequences)} sequências obtidas (L={len(final_sequences[0]) if final_sequences else 0})"
    )

    return final_sequences, used_params


def _parse_fasta_data(fasta_data: str) -> list[str]:
    """
    Processa dados FASTA baixados do NCBI.

    Args:
        fasta_data: String com dados FASTA

    Returns:
        Lista de sequências limpas
    """
    sequences = []
    current_seq = ""

    for line in fasta_data.split("\n"):
        line = line.strip()
        if line.startswith(">"):
            if current_seq:
                clean_seq = "".join(c.upper() for c in current_seq if c.isalpha())
                if clean_seq:
                    sequences.append(clean_seq)
            current_seq = ""
        else:
            current_seq += line

    # Adicionar última sequência
    if current_seq:
        clean_seq = "".join(c.upper() for c in current_seq if c.isalpha())
        if clean_seq:
            sequences.append(clean_seq)

    return sequences
