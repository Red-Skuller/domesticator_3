from typing import Any, Dict, Optional
from pydantic import BaseModel, Field, ConfigDict
from pydantic_settings import BaseSettings, SettingsConfigDict

class SequenceRecord(BaseModel):
    record_id: str
    protein_sequence: Optional[str] = None
    metadata: Dict[str, Any] = Field(default_factory=dict)

class ResultRow(BaseModel):
    record_id: str
    protein_sequence: Optional[str] = None
    naive_dna_sequence: Optional[str] = None
    optimized_sequence: Optional[str] = None
    vendor_score: Optional[float] = None
    accepted: bool = False
    status: str
    error_message: Optional[str] = None
    # Store the Biopython SeqRecord object; exclude from standard serializations
    optimized_record: Optional[Any] = Field(default=None, exclude=True)

    model_config = ConfigDict(arbitrary_types_allowed=True)

class PipelineConfig(BaseSettings):
    vendor_target: Optional[str] = None
    product: str = "eblocks"
    max_workers: int = Field(default=2, ge=1, le=64)
    min_length: int = Field(default=300, ge=0)

    model_config = SettingsConfigDict(
        env_prefix="DOMESTICA_",
        env_file=".env",
        extra="ignore"
    )