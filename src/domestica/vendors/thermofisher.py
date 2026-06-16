import time
import httpx
from typing import Tuple, Optional
from pydantic_settings import BaseSettings, SettingsConfigDict
from pydantic import Field, SecretStr
from domestica.vendors.base import ComplexityEvaluator, register_vendor


class ThermoFisherSettings(BaseSettings):
    client_id: str = Field(...)
    client_secret: SecretStr = Field(...)
    model_config = SettingsConfigDict(env_prefix="DOMESTICA_THERMOFISHER_", env_file=".env", extra="ignore")


@register_vendor("thermofisher")
class ThermoFisherEvaluator(ComplexityEvaluator):
    PRODUCT_MAPPING = {
        "dnastrings": "dnaStrings", "hqdnastrings": "hqDnaStrings",
        "eblocks": "dnaStrings", "gblocks": "dnaStrings", "genes": "dnaStrings"
    }

    def __init__(self, product: str):
        super().__init__(product)
        self.settings = ThermoFisherSettings()
        self.api_product_value = self.PRODUCT_MAPPING.get(product.lower(), "dnaStrings")
        self.http_client = httpx.Client(timeout=75.0)
        self._token = None
        self._token_expires_at = 0.0

    def _get_token(self, force_refresh: bool = False) -> str:
        if self._token and not force_refresh and time.time() < (self._token_expires_at - 60):
            return self._token

        response = self.http_client.post(
            "https://api.thermofisher.com/api/store/geneart/design-services/oauth2/token",
            json={
                "client_id": self.settings.client_id,
                "client_secret": self.settings.client_secret.get_secret_value(),
                "grant_type": "client_credentials"
            }
        )
        response.raise_for_status()
        body = response.json()
        self._token = body["access_token"]
        self._token_expires_at = time.time() + float(body.get("expires_in", 3600))
        return self._token

    def evaluate(self, sequence: str) -> Tuple[bool, Optional[float]]:
        cleaned_sequence = "".join(filter(str.isalpha, sequence.upper()))

        while True:
            response = self.http_client.post(
                "https://api.thermofisher.com/api/store/geneart/design-services/diagnostics/v1",
                json={"acgtSequence": cleaned_sequence, "product": self.api_product_value},
                headers={"Authorization": f"Bearer {self._get_token()}"},
                params={"waitSec": 60}
            )

            if response.status_code == 401:
                self._get_token(force_refresh=True)
                continue
            if response.status_code == 429:
                time.sleep(10)
                continue

            response.raise_for_status()
            complexity = response.json().get("content", {}).get("complexity", "red").lower()
            return complexity != "red", None