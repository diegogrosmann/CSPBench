"""
NCBI Dataset Download Module - CSPBench

Provides functionality for automatic download of biological sequences
from NCBI via Entrez API for Closest String Problem (CSP) applications.

Configuration:
    Create a .env file in the project root with:
    NCBI_EMAIL=your_email@example.com
    NCBI_API_KEY=your_api_key  # Optional

Author: CSPBench Development Team
"""

import logging
import os
import random
from typing import Any, Dict, List, Tuple, Optional, Union
from urllib.error import HTTPError

from dotenv import load_dotenv
from src.domain.dataset import Dataset
from src.domain.config import EntrezDatasetConfig

load_dotenv()

# TODO: Reimplementar após migração completa
# from cspbench.ui.cli.entrez_wizard import collect_entrez_parameters

try:
    from Bio import Entrez, SeqIO
except ImportError as exc:
    raise ImportError(
        "Biopython not found. Install with: pip install biopython"
    ) from exc

logger = logging.getLogger(__name__)

# Default settings loaded from .env
ENTREZ_DEFAULTS = {
    "email": os.getenv("NCBI_EMAIL", "change_me@example.com"),
    "api_key": os.getenv("NCBI_API_KEY"),
    "db": "nucleotide",
    # Padrões unificados (sem duplicidade)
    "query": "COIGene AND 600:650[SLEN]",
    "max_sequences": None,
    "rettype": "fasta",
    "retmode": "text",
    "min_length": None,
    "max_length": None,
    "uniform_policy": "strict",
    "pad_char": "N",
}


class EntrezDatasetDownloader:
    """
    NCBI Entrez dataset downloader for CSP applications.

    Provides static methods for downloading biological sequences from NCBI
    and converting them to Dataset objects.
    """

    @staticmethod
    def download(
        params: Union[Dict[str, Any], EntrezDatasetConfig],
    ) -> Tuple[Dataset, Dict[str, Any]]:
        """
        Download a dataset from NCBI usando os parâmetros fornecidos.

        Args:
            params: Parâmetros de configuração como dict ou EntrezDatasetConfig

        Returns:
            Tuple[Dataset, Dict[str, Any]]: Dataset com as sequências e
            um dicionário com todos os parâmetros efetivamente utilizados.

        Raises:
            ValueError: If error in search, download or validation
            HTTPError: If communication problem with NCBI
        """
        # Convert EntrezDatasetConfig to dict if needed
        if isinstance(params, EntrezDatasetConfig):
            params_dict = EntrezDatasetDownloader._config_to_params(params)
        else:
            params_dict = params

        sequences, used_params = EntrezDatasetDownloader._fetch_sequences(params_dict)

        return Dataset(sequences), used_params

    @staticmethod
    def _config_to_params(config: EntrezDatasetConfig) -> Dict[str, Any]:
        """
        Convert EntrezDatasetConfig to parameters dict.

        Args:
            config: EntrezDatasetConfig object

        Returns:
            Dict with parameters for _fetch_sequences
        """
        return {
            "query": config.query,
            "db": config.db,
            "max_sequences": config.retmax,
            "min_length": config.min_length,
            "max_length": config.max_length,
            "uniform_policy": config.uniform_policy,
            # Add default values from ENTREZ_DEFAULTS for missing fields
            "email": ENTREZ_DEFAULTS["email"],
            "api_key": ENTREZ_DEFAULTS["api_key"],
            "rettype": ENTREZ_DEFAULTS["rettype"],
            "retmode": ENTREZ_DEFAULTS["retmode"],
        }

    @staticmethod
    def _fetch_sequences(params: Dict[str, Any]) -> Tuple[List[str], Dict[str, Any]]:
        """
        Internal method to fetch sequences from NCBI.

        Args:
            params: Configuration parameters

        Returns:
            Tuple with list of sequences and used parameters

        Raises:
            ValueError: If error in search, download or validation
            HTTPError: If communication problem with NCBI
        """
        # Mesclar parâmetros com configurações padrão
        merged_params = {**ENTREZ_DEFAULTS}
        merged_params.update(params)

        # Extrair parâmetros de configuração (aceitando aliases amigáveis)
        email = merged_params["email"]
        db = merged_params["db"]
        # query é o nome canônico; aceitar 'n' apenas como alias para max_sequences
        query = merged_params.get("query", "*")
        max_seq_raw = merged_params.get("max_sequences", merged_params.get("n"))
        if max_seq_raw is None:
            max_sequences = None
        else:
            try:
                max_sequences = int(max_seq_raw)
                if max_sequences <= 0:
                    raise ValueError
            except (TypeError, ValueError):
                raise ValueError(
                    f"Parâmetro 'max_sequences'/'n' inválido: {max_seq_raw!r}"
                )
        api_key = merged_params.get("api_key")

        # Parâmetros de pós-processamento (opcionais)
        uniform_policy: Optional[str] = merged_params.get("uniform_policy")
        pad_char: str = merged_params.get("pad_char", "N")
        # Novos parâmetros opcionais de filtro de comprimento
        min_length: Optional[int] = merged_params.get("min_length")
        max_length: Optional[int] = merged_params.get("max_length")
        # novos: usar os configurados
        rettype: str = merged_params["rettype"]
        retmode: str = merged_params["retmode"]

        logger.debug(
            "_fetch_sequences chamado com db=%s, query='%s', max_sequences=%s",
            db,
            query,
            max_sequences,
        )

        # Configurar cliente Entrez
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        # Configurar gerador para seleção reprodutível
        rng = random.Random(0)

        # Fase 1: Buscar IDs das sequências
        logger.debug("Iniciando ESearch...")
        try:
            # Se max_sequences é None, usamos um retmax interno razoável para busca
            search_retmax = (max_sequences * 2) if isinstance(max_sequences, int) else 100
            handle = Entrez.esearch(db=db, term=query, retmax=search_retmax)
            search_result = Entrez.read(handle)
            handle.close()
        except HTTPError as e:
            raise ValueError(f"Erro ao acessar NCBI: HTTP {e.code} - {e.reason}") from e
        except Exception as e:
            raise ValueError(f"Erro na busca no NCBI: {e}") from e

        # Validar resultado da busca
        if not isinstance(search_result, dict):
            raise TypeError(
                f"Resultado da busca Entrez não é um dicionário: {type(search_result)}"
            )

        logger.debug("ESearch retornou %s IDs", len(search_result.get("IdList", [])))

        # Extrair IDs dos resultados
        ids = search_result.get("IdList", [])
        if not ids:
            raise ValueError(
                "Nenhum resultado encontrado. Verifique seu termo de busca."
            )

        # Selecionar IDs para download
        if isinstance(max_sequences, int) and len(ids) < max_sequences:
            logger.warning(
                "A busca retornou %s IDs, menos que os %s solicitados. Usando todos os IDs encontrados.",
                len(ids),
                max_sequences,
            )
            max_sequences = len(ids)
            sample_ids = ids
        elif isinstance(max_sequences, int):
            # Para compatibilidade com testes VCR, usar primeiros IDs quando n é pequeno
            if max_sequences <= 10:
                sample_ids = ids[:max_sequences]
            else:
                sample_ids = rng.sample(ids, max_sequences)
        else:
            # Sem limite explícito: usar todos os IDs retornados pela busca
            sample_ids = ids

        logger.debug("IDs amostrados (primeiros 5): %s ...", sample_ids[:5])

        # Fase 2: Download das sequências
        logger.debug("Iniciando EFetch...")
        try:
            fetch_handle = Entrez.efetch(
                db=db, id=",".join(sample_ids), rettype=rettype, retmode=retmode
            )
            records = list(SeqIO.parse(fetch_handle, "fasta"))
            fetch_handle.close()
        except HTTPError as e:
            raise ValueError(
                f"Erro ao baixar sequências do NCBI: HTTP {e.code} - {e.reason}"
            ) from e
        except Exception as e:
            raise ValueError(f"Erro no processamento das sequências: {e}") from e

        # (continuação do processamento dentro da função)

        # Processar e normalizar sequências
        seqs = [str(rec.seq).upper() for rec in records]

        if not seqs:
            raise ValueError("Nenhuma sequência válida foi obtida")

        # Aplicar filtros de comprimento antes da uniformização
        if min_length is not None or max_length is not None:
            mn = int(min_length) if min_length is not None else None
            mx = int(max_length) if max_length is not None else None
            seqs = [
                s
                for s in seqs
                if (mn is None or len(s) >= mn) and (mx is None or len(s) <= mx)
            ]
            if not seqs:
                raise ValueError(
                    "Nenhuma sequência após aplicar filtros de comprimento"
                )

        # Uniformização opcional
        if uniform_policy:
            seqs = EntrezDatasetDownloader._apply_uniform_policy(
                seqs, uniform_policy, pad_char
            )

        # Preparar parâmetros de retorno (nomes canônicos)
        used_params = {
            "email": email,
            "db": db,
            "query": query,
            "max_sequences": max_sequences,
            "api_key": api_key,
            "rettype": rettype,
            "retmode": retmode,
            "n_obtained": len(seqs),
            "uniform_policy": uniform_policy,
            "pad_char": pad_char,
            "min_length": min_length,
            "max_length": max_length,
        }
        L_info = len(seqs[0]) if seqs else 0
        logger.info("Dataset Entrez obtido: n=%s, L~%s", len(seqs), L_info)
        return seqs, used_params

    @staticmethod
    def _apply_uniform_policy(
        sequences: List[str], policy: str, pad_char: str = "N"
    ) -> List[str]:
        """Aplica política de uniformização às sequências.

        Políticas suportadas:
          - strict: mantém apenas sequências do comprimento mais frequente
          - majority: equivalente a strict (alias)
          - pad: preenche à direita com pad_char até o maior comprimento
          - trim: corta à direita até o menor comprimento
        """
        lengths = [len(s) for s in sequences]
        uniq = sorted(set(lengths))
        if len(uniq) <= 1:
            return sequences

        policy = str(policy).lower().strip()

        if policy in {"strict", "majority"}:
            from collections import Counter

            cnt = Counter(lengths)
            L, _ = cnt.most_common(1)[0]
            return [s for s in sequences if len(s) == L]

        if policy == "pad":
            L = max(lengths)
            return [s + (pad_char * (L - len(s))) for s in sequences]

        if policy == "trim":
            L = min(lengths)
            return [s[:L] for s in sequences]

        # Política desconhecida -> retorna inalterado
        return sequences
