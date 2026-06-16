import random
import time
import httpx
from base64 import b64encode
from typing import Tuple, Optional
from pydantic_settings import BaseSettings, SettingsConfigDict
from pydantic import Field, SecretStr
from domestica.vendors.base import ComplexityEvaluator, register_vendor


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

        auth_str = b64encode(
            f"{self.settings.client_id}:{self.settings.client_secret.get_secret_value()}".encode()
        ).decode()

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
        return self._token

    def evaluate(self, sequence: str) -> Tuple[bool, Optional[float]]:
        for attempt in range(6):
            response = self.http_client.post(
                self.endpoint, json=[{"Name": "Target", "Sequence": sequence}],
                headers={"Authorization": f"Bearer {self._get_token()}"}
            )
            if response.status_code == 401:
                self._get_token(force_refresh=True)
                continue
            if response.status_code in (429, 500, 502, 503, 504):
                time.sleep(random.uniform(0, min(60.0, 2.0 * (2 ** attempt))))
                continue

            response.raise_for_status()
            res = response.json()[0]
            score = res.get("ComplexityScore")

            if score is not None:
                return float(score) <= self.threshold, float(score)
            return res.get("IsAcceptable", False), None

        return False, None