import pytest
from unittest.mock import MagicMock, patch
httpx = pytest.importorskip("httpx")

from domestica.vendors.idt import IDTEvaluator
from domestica.vendors.thermofisher import ThermoFisherEvaluator


@patch("domestica.vendors.idt.httpx.Client")
def test_idt_evaluator_success(mock_client_cls):
    mock_client = MagicMock()
    mock_client_cls.return_value = mock_client

    # Mock token response and evaluation response
    token_response = MagicMock()
    token_response.status_code = 200
    token_response.json.return_value = {"access_token": "mock_idt_token", "expires_in": 3600}

    eval_response = MagicMock()
    eval_response.status_code = 200
    eval_response.json.return_value = [{"ComplexityScore": 5.0, "IsAcceptable": True}]

    mock_client.post.side_effect = [token_response, eval_response]

    evaluator = IDTEvaluator(product="eblocks")
    is_accepted, score = evaluator.evaluate("ATGCATGCATGC")

    assert is_accepted is True
    assert score == 5.0
    assert mock_client.post.call_count == 2


@patch("domestica.vendors.idt.httpx.Client")
def test_idt_evaluator_unauthorized_retry(mock_client_cls):
    mock_client = MagicMock()
    mock_client_cls.return_value = mock_client

    token_response_1 = MagicMock()
    token_response_1.status_code = 200
    token_response_1.json.return_value = {"access_token": "token_1", "expires_in": 3600}

    unauthorized_response = MagicMock()
    unauthorized_response.status_code = 401

    token_response_2 = MagicMock()
    token_response_2.status_code = 200
    token_response_2.json.return_value = {"access_token": "token_2", "expires_in": 3600}

    success_response = MagicMock()
    success_response.status_code = 200
    success_response.json.return_value = [{"ComplexityScore": 8.0, "IsAcceptable": True}]

    mock_client.post.side_effect = [
        token_response_1,
        unauthorized_response,
        token_response_2,
        success_response
    ]

    evaluator = IDTEvaluator(product="eblocks")
    is_accepted, score = evaluator.evaluate("ATGCATGCATGC")

    assert is_accepted is True
    assert score == 8.0


@patch("domestica.vendors.thermofisher.httpx.Client")
def test_thermofisher_evaluator_success(mock_client_cls):
    mock_client = MagicMock()
    mock_client_cls.return_value = mock_client

    token_response = MagicMock()
    token_response.status_code = 200
    token_response.json.return_value = {"access_token": "mock_tf_token", "expires_in": 3600}

    eval_response = MagicMock()
    eval_response.status_code = 200
    eval_response.json.return_value = {"content": {"complexity": "green"}}

    mock_client.post.side_effect = [token_response, eval_response]

    evaluator = ThermoFisherEvaluator(product="eblocks")
    is_accepted, score = evaluator.evaluate("ATGCATGCATGC")

    assert is_accepted is True
    assert score is None


@patch("domestica.vendors.thermofisher.httpx.Client")
def test_thermofisher_evaluator_red_complexity(mock_client_cls):
    mock_client = MagicMock()
    mock_client_cls.return_value = mock_client

    token_response = MagicMock()
    token_response.status_code = 200
    token_response.json.return_value = {"access_token": "mock_tf_token", "expires_in": 3600}

    eval_response = MagicMock()
    eval_response.status_code = 200
    eval_response.json.return_value = {"content": {"complexity": "red"}}

    mock_client.post.side_effect = [token_response, eval_response]

    evaluator = ThermoFisherEvaluator(product="eblocks")
    is_accepted, score = evaluator.evaluate("ATGCATGCATGC")

    assert is_accepted is False
    assert score is None