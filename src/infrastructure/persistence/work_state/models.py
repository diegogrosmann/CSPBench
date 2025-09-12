"""
SQLAlchemy models for work state persistence.

This module defines ORM models that replace the custom database driver implementation,
providing type-safe database operations and automatic schema management.
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
    """
    Custom SQLAlchemy type for storing JSON data in database.
    
    This type automatically serializes Python dictionaries and lists to JSON strings
    for storage and deserializes them back when retrieving from the database.
    """
    
    impl = Text
    cache_ok = True

    def process_bind_param(self, value, dialect):
        """
        Convert Python dict/list to JSON string for database storage.
        
        Args:
            value: Python object to serialize
            dialect: SQLAlchemy dialect (unused)
            
        Returns:
            str: JSON string representation, or None if value is None
        """
        if value is not None:
            return json.dumps(value)
        return value

    def process_result_value(self, value, dialect):
        """
        Convert JSON string back to Python dict/list when retrieving from database.
        
        Args:
            value: JSON string from database
            dialect: SQLAlchemy dialect (unused)
            
        Returns:
            Any: Deserialized Python object, or None if value is None
        """
        if value is not None:
            return json.loads(value)
        return value


class Work(Base):
    """
    Work table model representing a benchmark work unit.
    
    A work represents a complete benchmark execution containing multiple datasets,
    algorithm combinations, and their results.
    
    Attributes:
        id: Unique identifier for the work
        config_json: JSON configuration for the work
        status: Current execution status
        created_at: Timestamp when work was created
        updated_at: Timestamp when work was last updated
        output_path: Path where results are stored
        error: Error message if work failed
        extra_json: Additional metadata
    """
    
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
        """
        Convert work model to dictionary representation.
        
        Returns:
            Dict[str, Any]: Dictionary containing all work attributes
        """
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
    """
    Dataset table model representing a collection of sequences.
    
    Datasets contain the input sequences used for CSP algorithm benchmarking.
    Each dataset belongs to a specific work and can contain multiple sequences.
    
    Attributes:
        id: Auto-incremented primary key
        dataset_id: Original string identifier
        work_id: Foreign key to the parent work
        name: Human-readable dataset name
        meta_json: Metadata about the dataset
    """
    
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
        """
        Convert dataset model to dictionary representation.
        
        Returns:
            Dict[str, Any]: Dictionary containing all dataset attributes
        """
        return {
            "id": self.id,
            "dataset_id": self.dataset_id,
            "work_id": self.work_id,
            "name": self.name,
            "meta": self.meta_json or {},
        }


class DatasetSequence(Base):
    """
    Dataset sequence table model representing individual sequences within a dataset.
    
    Each sequence is stored separately to enable efficient querying and processing
    of large datasets.
    
    Attributes:
        dataset_id: Foreign key to the parent dataset
        seq_index: Index of the sequence within the dataset
        sequence: The actual sequence string
    """
    
    __tablename__ = 'dataset_sequences'
    
    dataset_id = Column(Integer, ForeignKey('datasets.id', ondelete='CASCADE'), primary_key=True, nullable=False)
    seq_index = Column(Integer, primary_key=True, nullable=False)
    sequence = Column(Text, nullable=False)
    
    # Relationships
    dataset = relationship("Dataset", back_populates="sequences")

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert dataset sequence model to dictionary representation.
        
        Returns:
            Dict[str, Any]: Dictionary containing all sequence attributes
        """
        return {
            "dataset_id": self.dataset_id,
            "seq_index": self.seq_index,
            "sequence": self.sequence,
        }


class Combination(Base):
    """
    Combination table model representing an algorithm-dataset-preset combination.
    
    A combination defines a specific test case within a work, pairing a dataset
    with an algorithm configuration (preset) for execution.
    
    Attributes:
        id: Auto-incremented primary key
        work_id: Foreign key to the parent work
        task_id: Identifier for the task
        dataset_id: Identifier for the dataset
        preset_id: Identifier for the algorithm preset
        algorithm_id: Identifier for the algorithm
        mode: Execution mode
        status: Current execution status
        total_sequences: Total number of sequences to process
        created_at: Timestamp when combination was created
        started_at: Timestamp when execution started
        finished_at: Timestamp when execution finished
    """
    
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
        """
        Convert combination model to dictionary representation.
        
        Returns:
            Dict[str, Any]: Dictionary containing all combination attributes
        """
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
    """
    Execution table model representing a single algorithm execution.
    
    An execution represents the processing of one sequence by one algorithm
    configuration, containing the parameters used and results obtained.
    
    Attributes:
        id: Auto-incremented primary key
        unit_id: Unique identifier for the execution unit
        combination_id: Foreign key to the parent combination
        sequencia: Sequence index being processed
        status: Current execution status
        started_at: Timestamp when execution started
        finished_at: Timestamp when execution finished
        params_json: Algorithm parameters used
        result_json: Results obtained from execution
        objective: Objective function value
    """
    
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
        """
        Convert execution model to dictionary representation.
        
        Returns:
            Dict[str, Any]: Dictionary containing all execution attributes
        """
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
    """
    Execution progress table model for tracking algorithm execution progress.
    
    This table stores progress updates during algorithm execution, allowing
    real-time monitoring of long-running operations.
    
    Attributes:
        id: Auto-incremented primary key
        execution_id: Foreign key to the parent execution
        progress: Progress value between 0.0 and 1.0
        message: Progress message or description
        timestamp: Timestamp when progress was recorded
    """
    
    __tablename__ = 'execution_progress'
    
    id = Column(Integer, primary_key=True, nullable=False)
    execution_id = Column(Integer, ForeignKey('executions.id', ondelete='CASCADE'), nullable=False)
    progress = Column(Float, CheckConstraint("progress >= 0.0 AND progress <= 1.0"), nullable=False)
    message = Column(Text, nullable=False)
    timestamp = Column(Float, nullable=False)
    
    # Relationships
    execution = relationship("Execution", back_populates="progress")

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert execution progress model to dictionary representation.
        
        Returns:
            Dict[str, Any]: Dictionary containing all progress attributes
        """
        return {
            "id": self.id,
            "execution_id": self.execution_id,
            "progress": self.progress,
            "message": self.message,
            "timestamp": self.timestamp,
        }


class Event(Base):
    """
    Event table model for logging system events and errors.
    
    Events provide audit trail and debugging information for work execution,
    capturing errors, progress updates, and other significant occurrences.
    
    Attributes:
        id: Auto-incremented primary key
        work_id: Foreign key to the parent work
        event_type: Type of event (error, progress, warning)
        event_category: Category of the event source
        entity_data_json: JSON data about the event entity
        timestamp: Timestamp when event occurred
    """
    
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
        """
        Convert event model to dictionary representation.
        
        Returns:
            Dict[str, Any]: Dictionary containing all event attributes
        """
        return {
            "id": self.id,
            "work_id": self.work_id,
            "event_type": self.event_type,
            "event_category": self.event_category,
            "entity_data": self.entity_data_json or {},
            "timestamp": self.timestamp,
        }
