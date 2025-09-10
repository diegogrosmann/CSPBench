"""
SQLAlchemy models for work state persistence.

Defines ORM models that replace the custom database driver implementation.
"""

import json
import time
from datetime import datetime
from typing import Any, Dict, List, Optional

from sqlalchemy import (
    CheckConstraint,
    Column,
    Float,
    ForeignKey,
    Integer,
    String,
    Text,
    UniqueConstraint,
    create_engine,
)
from sqlalchemy.orm import Session, relationship, sessionmaker, declarative_base
from sqlalchemy.types import TypeDecorator

Base = declarative_base()


class JSONType(TypeDecorator):
    """Custom type for storing JSON data in database."""
    
    impl = Text
    cache_ok = True

    def process_bind_param(self, value, dialect):
        """Convert Python dict/list to JSON string."""
        if value is not None:
            return json.dumps(value)
        return value

    def process_result_value(self, value, dialect):
        """Convert JSON string back to Python dict/list."""
        if value is not None:
            return json.loads(value)
        return value


class Work(Base):
    """Work table model."""
    
    __tablename__ = 'work'
    
    id = Column(String, primary_key=True, nullable=False)
    config_json = Column(JSONType, default=dict, nullable=False)
    status = Column(String, CheckConstraint("status IN ('queued', 'running', 'paused', 'canceled', 'completed', 'failed', 'error')"), nullable=False)
    created_at = Column(Float, nullable=False)
    updated_at = Column(Float)
    output_path = Column(String)
    error = Column(Text)
    extra_json = Column(JSONType, default=dict)
    
    # Relationships
    datasets = relationship("Dataset", back_populates="work", cascade="all, delete-orphan")
    combinations = relationship("Combination", back_populates="work", cascade="all, delete-orphan")
    events = relationship("Event", back_populates="work", cascade="all, delete-orphan")

    def to_dict(self) -> Dict[str, Any]:
        """Convert model to dictionary."""
        return {
            "id": self.id,
            "config_json": self.config_json or {},
            "status": self.status,
            "created_at": self.created_at,
            "updated_at": self.updated_at,
            "output_path": self.output_path,
            "error": self.error,
            "extra_json": self.extra_json or {},
        }


class Dataset(Base):
    """Dataset table model."""
    
    __tablename__ = 'datasets'
    
    id = Column(Integer, primary_key=True, autoincrement=True, nullable=False)
    dataset_id = Column(String, nullable=False)  # Original string ID
    work_id = Column(String, ForeignKey('work.id', ondelete='CASCADE'), nullable=False)
    name = Column(String, nullable=False)
    meta_json = Column(JSONType, default=dict)
    
    # Unique constraint for dataset_id within work
    __table_args__ = (
        UniqueConstraint('work_id', 'dataset_id', name='unique_dataset_per_work'),
    )
    
    # Relationships
    work = relationship("Work", back_populates="datasets")
    sequences = relationship("DatasetSequence", back_populates="dataset", cascade="all, delete-orphan")

    def to_dict(self) -> Dict[str, Any]:
        """Convert model to dictionary."""
        return {
            "id": self.id,
            "dataset_id": self.dataset_id,
            "work_id": self.work_id,
            "name": self.name,
            "meta": self.meta_json or {},
        }


class DatasetSequence(Base):
    """Dataset sequence table model."""
    
    __tablename__ = 'dataset_sequences'
    
    dataset_id = Column(Integer, ForeignKey('datasets.id', ondelete='CASCADE'), primary_key=True, nullable=False)
    seq_index = Column(Integer, primary_key=True, nullable=False)
    sequence = Column(Text, nullable=False)
    
    # Relationships
    dataset = relationship("Dataset", back_populates="sequences")

    def to_dict(self) -> Dict[str, Any]:
        """Convert model to dictionary."""
        return {
            "dataset_id": self.dataset_id,
            "seq_index": self.seq_index,
            "sequence": self.sequence,
        }


class Combination(Base):
    """Combination table model."""
    
    __tablename__ = 'combinations'
    
    id = Column(Integer, primary_key=True, nullable=False)
    work_id = Column(String, ForeignKey('work.id', ondelete='CASCADE'), nullable=False)
    task_id = Column(String, nullable=False)
    dataset_id = Column(String, nullable=False)  # Removed FK constraint for now
    preset_id = Column(String, nullable=False)
    algorithm_id = Column(String, nullable=False)
    mode = Column(String, nullable=False)
    status = Column(String, CheckConstraint("status IN ('queued', 'running', 'paused', 'canceled', 'completed', 'failed', 'error')"), nullable=False)
    total_sequences = Column(Integer, nullable=False)
    created_at = Column(Float, nullable=False)
    started_at = Column(Float)
    finished_at = Column(Float)
    
    # Constraints
    __table_args__ = (
        UniqueConstraint('work_id', 'task_id', 'dataset_id', 'preset_id', 'algorithm_id'),
    )
    
    # Relationships
    work = relationship("Work", back_populates="combinations")
    executions = relationship("Execution", back_populates="combination", cascade="all, delete-orphan")

    def to_dict(self) -> Dict[str, Any]:
        """Convert model to dictionary."""
        return {
            "id": self.id,
            "work_id": self.work_id,
            "task_id": self.task_id,
            "dataset_id": self.dataset_id,
            "preset_id": self.preset_id,
            "algorithm_id": self.algorithm_id,
            "mode": self.mode,
            "status": self.status,
            "total_sequences": self.total_sequences,
            "created_at": self.created_at,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
        }


class Execution(Base):
    """Execution table model."""
    
    __tablename__ = 'executions'
    
    id = Column(Integer, primary_key=True, nullable=False)
    unit_id = Column(String, nullable=False)
    combination_id = Column(Integer, ForeignKey('combinations.id', ondelete='CASCADE'), nullable=False)
    sequencia = Column(Integer, nullable=False)
    status = Column(String, CheckConstraint("status IN ('queued', 'running', 'paused', 'canceled', 'completed', 'failed', 'error')"), nullable=False)
    started_at = Column(Float)
    finished_at = Column(Float)
    params_json = Column(JSONType, default=dict)
    result_json = Column(JSONType, default=dict)
    objective = Column(Float)
    
    # Constraints
    __table_args__ = (
        UniqueConstraint('unit_id', 'combination_id', name='unique_unit_per_combination'),
    )
    
    # Relationships
    combination = relationship("Combination", back_populates="executions")
    progress = relationship("ExecutionProgress", back_populates="execution", cascade="all, delete-orphan")

    def to_dict(self) -> Dict[str, Any]:
        """Convert model to dictionary."""
        return {
            "id": self.id,
            "unit_id": self.unit_id,
            "combination_id": self.combination_id,
            "sequencia": self.sequencia,
            "status": self.status,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "params": self.params_json or {},
            "result": self.result_json or {},
            "objective": self.objective,
        }


class ExecutionProgress(Base):
    """Execution progress table model."""
    
    __tablename__ = 'execution_progress'
    
    id = Column(Integer, primary_key=True, nullable=False)
    execution_id = Column(Integer, ForeignKey('executions.id', ondelete='CASCADE'), nullable=False)
    progress = Column(Float, CheckConstraint("progress >= 0.0 AND progress <= 1.0"), nullable=False)
    message = Column(Text, nullable=False)
    timestamp = Column(Float, nullable=False)
    
    # Relationships
    execution = relationship("Execution", back_populates="progress")

    def to_dict(self) -> Dict[str, Any]:
        """Convert model to dictionary."""
        return {
            "id": self.id,
            "execution_id": self.execution_id,
            "progress": self.progress,
            "message": self.message,
            "timestamp": self.timestamp,
        }


class Event(Base):
    """Event table model."""
    
    __tablename__ = 'events'
    
    id = Column(Integer, primary_key=True, nullable=False)
    work_id = Column(String, ForeignKey('work.id', ondelete='CASCADE'), nullable=False)
    event_type = Column(String, CheckConstraint("event_type IN ('error', 'progress', 'warning')"), nullable=False)
    event_category = Column(String, CheckConstraint("event_category IN ('work', 'task', 'dataset', 'preset', 'combination', 'unit', 'algorithm', 'other')"), nullable=False)
    entity_data_json = Column(JSONType, default=dict, nullable=False)
    timestamp = Column(Float, nullable=False, default=time.time)
    
    # Relationships
    work = relationship("Work", back_populates="events")

    def to_dict(self) -> Dict[str, Any]:
        """Convert model to dictionary."""
        return {
            "id": self.id,
            "work_id": self.work_id,
            "event_type": self.event_type,
            "event_category": self.event_category,
            "entity_data": self.entity_data_json or {},
            "timestamp": self.timestamp,
        }
