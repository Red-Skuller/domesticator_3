import random
import time
import httpx
import logging
from base64 import b64encode
from typing import Tuple, Optional
from pydantic_settings import BaseSettings, SettingsConfigDict
from pydantic import Field, SecretStr
from domestica.vendors.base import ComplexityEvaluator, register_vendor

logger = logging.getLogger(__name__)


class IDTSettings(BaseSettings):
    client_id: str = Field(...)
    client_secret: SecretStr = Field(...)
    username: str = Field(...)
    password: SecretStr = Field(...)
    model_config = SettingsConfigDict(env_prefix="DOMESTICA_IDT_", env_file=".env", extra="ignore")


@register_vendor("idt")
class IDTEvaluator(ComplexityEvaluator):
    PRODUCT_ENDPOINTS = {
        "eblocks": "https://www.idtdna.com/Restapi/v1/Complexities/ScreenEblockSequences",
        "gblocks": "https://www.idtdna.com/Restapi/v1/Complexities/ScreenGblockSequences",
        "genes": "https://www.idtdna.com/Restapi/v1/Complexities/ScreenGeneSequences"
    }

    def __init__(self, product: str):
        super().__init__(product)
        self.settings = IDTSettings()
        self.threshold = {"eblocks": 10.0, "gblocks": 10.0, "genes": 10.0}.get(product, 0.0)
        self.endpoint = self.PRODUCT_ENDPOINTS[product]
        self.http_client = httpx.Client(timeout=30.0)
        self._token = None
        self._token_expires_at = 0.0

    def _get_token(self, force_refresh: bool = False) -> str:
        if self._token and not force_refresh and time.time() < (self._token_expires_at - 60):
            return self._token

        logger.debug("Requesting new IDT IdentityServer verification token.")
        auth_str = b64encode(
            f"{self.settings.client_id}:{self.settings.client_secret.get_secret_value()}".encode()
        ).decode()

        try:
            response = self.http_client.post(
                "https://www.idtdna.com/Identityserver/connect/token",
                data={
                    "grant_type": "password", "scope": "test",
                    "username": self.settings.username,
                    "password": self.settings.password.get_secret_value()
                },
                headers={"Content-Type": "application/x-www-form-urlencoded", "Authorization": f"Basic {auth_str}"}
            )
            response.raise_for_status()
            body = response.json()
            self._token = body["access_token"]
            self._token_expires_at = time.time() + float(body.get("expires_in", 3600))
            logger.debug("IDT OAuth2 security connection established successfully. Valid for %s seconds.",
                         body.get("expires_in", 3600))
            return self._token
        except Exception:
            logger.exception("Critical exception thrown while requesting authentication token from IDT server.")
            raise

    def evaluate(self, sequence: str) -> Tuple[bool, Optional[float]]:
        payload = [{"Name": "Target", "Sequence": sequence}]
        for attempt in range(6):
            logger.debug("Dispatching sequence request to IDT endpoint. Attempt: %d/6. Sequence length: %d bp", attempt + 1, len(sequence))
            try:
                response = self.http_client.post(
                    self.endpoint, json=payload,
                    headers={"Authorization": f"Bearer {self._get_token()}"}
                )
                if response.status_code == 401:
                    logger.warning("IDT endpoint returned 401 Unauthorized status. Forcing token regeneration.")
                    self._get_token(force_refresh=True)
                    continue
                if response.status_code in (429, 500, 502, 503, 504):
                    backoff = random.uniform(0, min(60.0, 2.0 * (2 ** attempt)))
                    logger.warning("IDT server limits (Status: %d). Delaying exponential backoff: %.2f seconds.", response.status_code, backoff)
                    time.sleep(backoff)
                    continue

                logger.debug("IDT Response Status: %d", response.status_code)
                response.raise_for_status()
                res = response.json()[0]
                score = res.get("ComplexityScore")

                if score is not None:
                    score_val = float(score)
                    logger.debug("IDT returned complexity metric: %f (Threshold: <= %f)", score_val, self.threshold)
                    return score_val <= self.threshold, score_val

                is_acceptable = res.get("IsAcceptable", False)
                logger.debug("Complexity metrics missing. Falling back to explicit acceptance status: %s", is_acceptable)
                return is_acceptable, None
            except Exception:
                logger.exception("Exception encountered on evaluation cycle step %d of 6.", attempt + 1)
                if attempt == 5:
                    raise

        logger.error("All structural request retransmissions to IDT have failed.")
        return False, None