"""Testes para normalização de status."""

import pytest

from src.domain.status import BaseStatus, normalize_status


class TestStatusNormalization:
    """Testes para helper de normalização de status."""

    def test_normalize_status_with_enum(self):
        """Testa normalização com instância de BaseStatus."""
        result = normalize_status(BaseStatus.RUNNING)
        assert result == "running"

        result = normalize_status(BaseStatus.COMPLETED)
        assert result == "completed"

        result = normalize_status(BaseStatus.FAILED)
        assert result == "failed"

    def test_normalize_status_with_string_lowercase(self):
        """Testa normalização com string já em lowercase."""
        result = normalize_status("running")
        assert result == "running"

        result = normalize_status("completed")
        assert result == "completed"

    def test_normalize_status_with_string_uppercase(self):
        """Testa normalização com string em uppercase."""
        result = normalize_status("RUNNING")
        assert result == "running"

        result = normalize_status("COMPLETED")
        assert result == "completed"

        result = normalize_status("FAILED")
        assert result == "failed"

    def test_normalize_status_with_string_mixedcase(self):
        """Testa normalização com string em case misto."""
        result = normalize_status("RunNing")
        assert result == "running"

        result = normalize_status("CompletED")
        assert result == "completed"

    def test_normalize_status_with_whitespace(self):
        """Testa normalização removendo whitespace."""
        result = normalize_status("  running  ")
        assert result == "running"

        result = normalize_status("\tCOMPLETED\n")
        assert result == "completed"

    def test_normalize_status_invalid_string(self):
        """Testa erro com string inválida."""
        # Generic assertion using substring for robustness
        with pytest.raises(ValueError) as exc:
            normalize_status("invalid_status")
        assert "Invalid status:" in str(exc.value)

        with pytest.raises(ValueError) as exc:
            normalize_status("UNKNOWN")
        assert "Invalid status:" in str(exc.value)

    def test_normalize_status_invalid_type(self):
        """Testa conversão de outros tipos para string."""
        # Números que correspondem a status inválidos devem falhar
        with pytest.raises(ValueError) as exc:
            normalize_status(123)
        assert "Invalid status:" in str(exc.value)

        with pytest.raises(ValueError) as exc:
            normalize_status(None)
        assert "Invalid status:" in str(exc.value)

    def test_all_enum_values_are_normalized_correctly(self):
        """Testa que todos os valores do enum são normalizados corretamente."""
        for status in BaseStatus:
            result = normalize_status(status)
            assert result == status.value
            assert isinstance(result, str)


# Testes legados mantidos para compatibilidade
def test_normalize_status_enum():
    assert normalize_status(BaseStatus.RUNNING) == "running"


def test_normalize_status_string_case():
    assert normalize_status("RUNNING") == "running"


def test_normalize_status_invalid():
    with pytest.raises(ValueError):
        normalize_status("not_a_status")
