from __future__ import annotations

from typing import Dict

import pydantic


class LoggingFormatter(pydantic.BaseModel):
    file_format: str = pydantic.Field(
        alias="format"
    )  # Avoid shadowing Python `format`


class LoggingHandler(pydantic.BaseModel):
    file_class: str = pydantic.Field(
        alias="class"
    )  # Avoid shadowing Python `class`
    level: str
    formatter: str
    filename: str

    @pydantic.field_validator("level")
    @classmethod
    def validate_level(cls, v: str) -> str:
        if v not in ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]:
            raise ValueError(f"Invalid logging level: {v}")
        return v


class LoggingRoot(pydantic.BaseModel):
    level: str
    handlers: list[str]

    @pydantic.field_validator("level")
    @classmethod
    def validate_level(cls, v: str) -> str:
        if v not in ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]:
            raise ValueError(f"Invalid logging level: {v}")
        return v


class LoggingConfig(pydantic.BaseModel):
    version: int = 1
    disable_existing_loggers: bool = False
    formatters: Dict[str, LoggingFormatter]
    handlers: Dict[str, LoggingHandler]
    root: LoggingRoot


class CoralConfig(pydantic.BaseModel):
    logging: LoggingConfig
