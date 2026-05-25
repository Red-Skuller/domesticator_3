from pydantic import BaseModel, Field
from typing import List, Optional
from pathlib import Path

class IDTConfig(BaseModel):
    username: str
    password: str
    client_id: str
    client_secret: str
    token_file_path: Optional[Path] = None

class PipelineConfig(BaseModel):
    input_path: Path
    output_path: Path
    params: List[str] = Field(default_factory=list)
    optimize: bool = False
    vector_path: Optional[Path] = None
    nstruct: int = 10
    skip_idt: bool = False
    ph: float = 7.4
    name_col: str = "Name"
    seq_col: str = "Sequence"
    idt_type: str = "gene"
    idt_credentials_dir: Path = Path("~/.idt_credentials")
    idt_threshold: float = 7.0
    n_tag: str = ""
    c_tag: str = ""
    out_cols: Optional[List[str]] = None

class ResultRow(BaseModel):
    name: str = Field(..., alias="Name")
    aa_seq: str = Field(..., alias="AA_seq")
    n_tag: str = Field("", alias="N_tag")
    c_tag: str = Field("", alias="C_tag")
    insert_name: Optional[str] = Field(None, alias="insert_name")
    aa_seq_final: Optional[str] = Field(None, alias="AA_Seq_Final")
    dna_seq: Optional[str] = Field(None, alias="DNA_seq")
    idt_score: Optional[float] = Field(None, alias="IDT_Score")
    dnachisel_score: Optional[float] = Field(None, alias="Dnachisel_Score")

    # Allow for extra metrics from AA_analysis
    model_config = {
        "extra": "allow",
        "populate_by_name": True
    }
