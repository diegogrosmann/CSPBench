"""
Entrez Dataset Repository Implementation

Implements EntrezDatasetRepository port for fetching datasets from NCBI using Entrez API.
Part of the infrastructure layer - handles external API communication.
"""

import logging
import os
from typing import Any, Dict, List, Tuple

from dotenv import load_dotenv

from src.application.ports import EntrezDatasetRepository
from src.domain.errors import DatasetNotFoundError, InfrastructureError

# Load environment variables
load_dotenv()

logger = logging.getLogger(__name__)


class NCBIEntrezDatasetRepository:
    """
    Implementation of EntrezDatasetRepository for NCBI datasets.
    
    Uses Biopython's Entrez module to fetch biological sequences from NCBI databases.
    Requires NCBI_EMAIL and optionally NCBI_API_KEY in environment variables.
    """

    def __init__(self, email: str = None, api_key: str = None):
        """
        Initialize Entrez repository.
        
        Args:
            email: NCBI email (required by NCBI policy)
            api_key: NCBI API key (optional, for higher rate limits)
        """
        # Load credentials from environment or parameters
        self.email = email or os.getenv("NCBI_EMAIL")
        self.api_key = api_key or os.getenv("NCBI_API_KEY")
        
        if not self.email:
            raise InfrastructureError(
                "NCBI email is required. Set NCBI_EMAIL environment variable "
                "or provide email parameter."
            )
        
        # Check if Biopython is available
        try:
            from Bio import Entrez
            self._entrez = Entrez
            self._entrez.email = self.email
            if self.api_key:
                self._entrez.api_key = self.api_key
        except ImportError as e:
            raise InfrastructureError(
                "Biopython is required for Entrez datasets. "
                "Install with: pip install biopython"
            ) from e
        
        logger.info("Entrez repository initialized with email: %s", self.email)

    def fetch_dataset(
        self, query: str, db: str = "nucleotide", retmax: int = 20, **kwargs
    ) -> Tuple[List[str], Dict[str, Any]]:
        """
        Fetch dataset from NCBI using Entrez API.
        
        Args:
            query: Search query for NCBI
            db: Database to search (nucleotide, protein, etc.)
            retmax: Maximum number of sequences to fetch
            **kwargs: Additional parameters (seed, etc.)
            
        Returns:
            Tuple[List[str], Dict[str, Any]]: (sequences, metadata)
            
        Raises:
            DatasetNotFoundError: If no sequences found
            InfrastructureError: If API call fails
        """
        try:
            # Import the fetch function from the existing module
            from src.infrastructure.external.dataset_entrez import fetch_dataset_silent
            
            # Prepare parameters for the fetch function
            params = {
                "email": self.email,
                "db": db,
                "term": query,
                "n": retmax,
                "api_key": self.api_key,
            }
            
            # Add any additional parameters
            params.update(kwargs)
            
            logger.info(
                "Fetching dataset from NCBI: db=%s, query='%s', retmax=%d",
                db, query, retmax
            )
            
            # Fetch the dataset
            sequences, metadata = fetch_dataset_silent(params)
            
            if not sequences:
                raise DatasetNotFoundError(
                    f"No sequences found for query: {query} in database: {db}"
                )
            
            logger.info(
                "Successfully fetched %d sequences of length %d",
                len(sequences), len(sequences[0]) if sequences else 0
            )
            
            return sequences, metadata
            
        except Exception as e:
            if isinstance(e, (DatasetNotFoundError, InfrastructureError)):
                raise
            
            # Wrap other exceptions as infrastructure errors
            raise InfrastructureError(
                f"Failed to fetch dataset from NCBI: {str(e)}"
            ) from e

    def is_available(self) -> bool:
        """
        Check if Entrez service is available.
        
        Returns:
            bool: True if service is available
        """
        try:
            # Check if we have required dependencies and credentials
            if not self.email:
                return False
            
            # For basic availability check, just verify we have the required modules
            # The actual connectivity will be tested when fetching data
            return True
            
        except Exception as e:
            logger.warning("Entrez service availability check failed: %s", str(e))
            return False
