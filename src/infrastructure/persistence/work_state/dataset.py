from typing import Any, Optional
import json
from src.domain.dataset import Dataset


class DatasetMixin:
    def submit_dataset(self, id: str, dataset_obj: Dataset, meta: dict[str, any]) -> None:

        name = getattr(dataset_obj, "name", id)

        from src.infrastructure.persistence.work_state.utils.json_helpers import (
            serialize_to_json,
        )

        js = serialize_to_json(meta)

        self._execute(
            "INSERT OR IGNORE INTO datasets(id, name, meta_json) VALUES(?,?,?)",
            (id, name, js),
        )

        # Add sequences
        sequences = dataset_obj.get_sequences()
        self.add_sequences(id, sequences)

    def add_sequences(self, dataset_id: str, sequences: list[str]):
        """Add sequences for a dataset."""
        if not sequences:
            return
        rows = [(dataset_id, idx, seq) for idx, seq in enumerate(sequences)]
        self._executemany(
            "INSERT OR IGNORE INTO dataset_sequences(dataset_id, seq_index, sequence) VALUES(?,?,?)",
            rows,
        )

    def has_dataset(self, dataset_id: str) -> bool:
        """Check if dataset exists with sequences."""
        if not dataset_id:
            return False

        with self._lock:
            cursor = self._conn.cursor()
            cursor.execute(
                "SELECT COUNT(*) FROM dataset_sequences WHERE dataset_id=?",
                (dataset_id,),
            )
            count = cursor.fetchone()[0]
            return count > 0

    def get_dataset(self, dataset_id: str) -> Optional[Dataset]:
        """
        Retrieve a dataset by ID.

        Args:
            dataset_id: The unique identifier of the dataset

        Returns:
            Dataset object if found, None if not found or if dataset_id is empty

        Note:
            The returned Dataset object will have its sequences and metadata
            reconstructed from the database. The alphabet is extracted from
            the stored metadata if available.
        """
        if not dataset_id:
            return None

        with self._lock:
            cursor = self._conn.cursor()

            # Get dataset metadata
            cursor.execute(
                "SELECT id, name, meta_json FROM datasets WHERE id=?",
                (dataset_id,),
            )
            dataset_row = cursor.fetchone()

            if not dataset_row:
                return None

            ds_id, name, meta_json = dataset_row

            # Get sequences ordered by seq_index
            cursor.execute(
                "SELECT sequence FROM dataset_sequences WHERE dataset_id=? ORDER BY seq_index",
                (dataset_id,),
            )
            sequence_rows = cursor.fetchall()

            if not sequence_rows:
                return None

            sequences = [row[0] for row in sequence_rows]

            # Parse metadata
            metadata = {}
            if meta_json:
                try:
                    metadata = json.loads(meta_json)
                except json.JSONDecodeError:
                    metadata = {}

            # Extract alphabet from metadata if available
            alphabet = metadata.get("alphabet")

            # Create Dataset object
            dataset = Dataset(
                name=name or ds_id, sequences=sequences, alphabet=alphabet
            )

            return dataset
