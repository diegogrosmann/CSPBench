"""Testes para normalização de status."""

import pytest
from src.domain.status import BaseStatus, normalize_status
from src.infrastructure.persistence.work_state.utils.validation import validate_status


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
        with pytest.raises(ValueError, match="Status inválido"):
            normalize_status("invalid_status")

        with pytest.raises(ValueError, match="Status inválido"):
            normalize_status("UNKNOWN")

    def test_normalize_status_invalid_type(self):
        """Testa conversão de outros tipos para string."""
        # Números que correspondem a status inválidos devem falhar
        with pytest.raises(ValueError, match="Status inválido"):
            normalize_status(123)

        with pytest.raises(ValueError, match="Status inválido"):
            normalize_status(None)

    def test_all_enum_values_are_normalized_correctly(self):
        """Testa que todos os valores do enum são normalizados corretamente."""
        for status in BaseStatus:
            result = normalize_status(status)
            assert result == status.value
            assert isinstance(result, str)


class TestValidateStatus:
    """Testes para função validate_status que usa normalize_status."""

    def test_validate_status_with_enum(self):
        """Testa validação com enum - não deve levantar exceção."""
        validate_status(BaseStatus.RUNNING)
        validate_status(BaseStatus.COMPLETED)
        validate_status(BaseStatus.FAILED)

    def test_validate_status_with_valid_string(self):
        """Testa validação com string válida - não deve levantar exceção."""
        validate_status("running")
        validate_status("COMPLETED")
        validate_status("  failed  ")

    def test_validate_status_with_invalid_string(self):
        """Testa validação com string inválida - deve levantar exceção."""
        with pytest.raises(ValueError):
            validate_status("invalid")

        with pytest.raises(ValueError):
            validate_status("unknown")

        with pytest.raises(ValueError):
            validate_status("")

    def test_validate_status_integration_with_normalize(self):
        """Testa integração entre validate_status e normalize_status."""
        # Casos que devem passar
        test_cases = [
            BaseStatus.QUEUED,
            "running",
            "COMPLETED",
            "  paused  ",
            "CanceleD",
        ]

        for case in test_cases:
            # validate_status não deve levantar exceção
            validate_status(case)
            # normalize_status deve retornar string válida
            normalized = normalize_status(case)
            assert normalized in {s.value for s in BaseStatus}


# Testes legados mantidos para compatibilidade
def test_normalize_status_enum():
    assert normalize_status(BaseStatus.RUNNING) == "running"


def test_normalize_status_string_case():
    assert normalize_status("RUNNING") == "running"


def test_normalize_status_invalid():
    with pytest.raises(ValueError):
        normalize_status("not_a_status")


def test_validate_status_accepts_enum_and_string():
    # Não deve levantar
    validate_status(BaseStatus.QUEUED)
    validate_status("queued")
