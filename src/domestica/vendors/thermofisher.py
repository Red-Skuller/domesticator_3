import time
import httpx
import logging
from typing import Tuple, Optional
from pydantic_settings import BaseSettings, SettingsConfigDict
from pydantic import Field, SecretStr
from domestica.vendors.base import ComplexityEvaluator, register_vendor

logger = logging.getLogger(__name__)


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

        logger.debug("Requesting new ThermoFisher OAuth2 access token credentials.")
        try:
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
            logger.debug("ThermoFisher OAuth2 token successfully established. Valid for %s seconds.", body.get("expires_in", 3600))
            return self._token
        except Exception:
            logger.exception("Critical exception thrown while requesting authentication token from ThermoFisher endpoint.")
            raise

    def evaluate(self, sequence: str) -> Tuple[bool, Optional[float]]:
        cleaned_sequence = "".join(filter(str.isalpha, sequence.upper()))
        url = "https://api.thermofisher.com/api/store/geneart/design-services/diagnostics/v1"

        while True:
            logger.debug("Dispatching post sequence validation request payload to ThermoFisher diagnostics API.")
            try:
                response = self.http_client.post(
                    url,
                    json={"acgtSequence": cleaned_sequence, "product": self.api_product_value},
                    headers={"Authorization": f"Bearer {self._get_token()}"},
                    params={"waitSec": 60}
                )

                if response.status_code == 401:
                    logger.warning("ThermoFisher API token invalidated (401 Unauthorized). Forcing token refresh routine sequence.")
                    self._get_token(force_refresh=True)
                    continue
                if response.status_code == 429:
                    logger.warning("ThermoFisher rate limit ceiling reached (429 Too Many Requests). Dormant cooling backoff for 10 seconds.")
                    time.sleep(10)
                    continue

                response.raise_for_status()
                complexity = response.json().get("content", {}).get("complexity", "red").lower()
                logger.debug("ThermoFisher evaluation structural assessment classification score returned: '%s'", complexity)
                return complexity != "red", None
            except Exception:
                logger.exception("Fatal processing breakdown when executing ThermoFisher API sequence screening requests.")
                raise