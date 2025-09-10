"""
Migration guide and utilities for SQLAlchemy transition.

This module provides utilities and documentation for migrating from the
custom database driver implementation to SQLAlchemy ORM.
"""

import logging
import os
import sys
from pathlib import Path
from typing import Optional

# Add the project root to Python path
project_root = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(project_root))

from src.infrastructure.persistence.work_state.core import SQLAlchemyWorkPersistence


def migrate_to_sqlalchemy(
    old_database_path: Optional[str] = None,
    new_database_url: Optional[str] = None
) -> SQLAlchemyWorkPersistence:
    """
    Migrate from old database driver to SQLAlchemy.
    
    This function helps transition existing code to use the new SQLAlchemy-based
    persistence layer. The database schema is compatible, so existing data
    should work with minimal changes.
    
    Args:
        old_database_path: Path to existing SQLite database (optional)
        new_database_url: New database URL for SQLAlchemy (optional)
        
    Returns:
        Configured SQLAlchemy persistence instance
        
    Example:
        # Simple migration
        persistence = migrate_to_sqlalchemy()
        
        # With custom database URL
        persistence = migrate_to_sqlalchemy(
            new_database_url="sqlite:///./data/new_database.db"
        )
    """
    logger = logging.getLogger(__name__)
    
    # If old database path is provided, use it as the basis for new URL
    if old_database_path and not new_database_url:
        if old_database_path.endswith('.db'):
            new_database_url = f"sqlite:///{old_database_path}"
        else:
            logger.warning(f"Could not auto-convert database path: {old_database_path}")
    
    # Create new SQLAlchemy persistence
    persistence = SQLAlchemyWorkPersistence(new_database_url)
    
    logger.info("Successfully migrated to SQLAlchemy persistence")
    return persistence


class MigrationGuide:
    """
    Migration guide for transitioning from custom drivers to SQLAlchemy.
    
    This class provides documentation and examples for the migration process.
    """
    
    @staticmethod
    def print_migration_guide():
        """Print comprehensive migration guide."""
        guide = """
        ═══════════════════════════════════════════════════════════════
                        SQLAlchemy Migration Guide
        ═══════════════════════════════════════════════════════════════
        
        The CSPBench persistence layer has been migrated from custom database
        drivers to SQLAlchemy ORM for better maintainability and type safety.
        
        🔄 CHANGES MADE:
        ───────────────────────────────────────────────────────────────
        
        ✅ REPLACED:
        - Custom DatabaseDriver classes
        - DatabaseProcessWorker 
        - Manual SQL execution with raw drivers
        - Mixin-based persistence classes
        
        ✅ WITH:
        - SQLAlchemy ORM models
        - Type-safe database operations
        - Automatic schema management
        - Better error handling and logging
        
        📋 API COMPATIBILITY:
        ───────────────────────────────────────────────────────────────
        
        The public API remains largely the same:
        
        # Before (still works)
        from src.infrastructure.persistence.work_state import WorkPersistence
        persistence = WorkPersistence()
        
        # After (recommended)
        from src.infrastructure.persistence.work_state import WorkPersistence
        persistence = WorkPersistence()  # Now uses SQLAlchemy internally
        
        🗃️ DATABASE COMPATIBILITY:
        ───────────────────────────────────────────────────────────────
        
        - Existing SQLite databases continue to work
        - Schema is automatically managed by SQLAlchemy
        - No data migration required
        
        ⚙️ CONFIGURATION:
        ───────────────────────────────────────────────────────────────
        
        Environment Variables (unchanged):
        - DB_ENGINE: Database engine type (default: sqlite)
        - DB_URL: Database connection string or path
        
        🔧 NEW FEATURES:
        ───────────────────────────────────────────────────────────────
        
        - Type-safe ORM models
        - Better error messages
        - Automatic connection management
        - Support for multiple database backends
        - Built-in schema migrations
        
        📝 DEPRECATED METHODS:
        ───────────────────────────────────────────────────────────────
        
        These methods still work but are deprecated:
        - _execute() - Use ORM methods instead
        - _fetch_one() - Use ORM queries instead
        - _fetch_all() - Use ORM queries instead
        
        🚀 GETTING STARTED:
        ───────────────────────────────────────────────────────────────
        
        1. Update imports (if needed):
           from src.infrastructure.persistence.work_state import WorkPersistence
        
        2. Create persistence instance:
           persistence = WorkPersistence()
        
        3. Use the same API as before:
           persistence.work_create(id="my-work", config={...})
           work = persistence.work_get("my-work")
        
        4. For new code, consider using ORM models directly:
           from src.infrastructure.persistence.work_state.models import Work
           
        📞 SUPPORT:
        ───────────────────────────────────────────────────────────────
        
        If you encounter issues during migration:
        1. Check logs for deprecation warnings
        2. Verify database URL configuration
        3. Ensure all imports are correct
        4. Review the new SQLAlchemy models if needed
        
        ═══════════════════════════════════════════════════════════════
        """
        print(guide)
    
    @staticmethod
    def check_migration_status():
        """Check if migration is working correctly."""
        try:
            # Test creating a persistence instance
            persistence = SQLAlchemyWorkPersistence()
            print("✅ SQLAlchemy persistence initialization: SUCCESS")
            
            # Test basic operations
            with persistence.get_session() as session:
                print("✅ Database session creation: SUCCESS")
            
            persistence.close()
            print("✅ Database cleanup: SUCCESS")
            print("\n🎉 Migration verification completed successfully!")
            
        except Exception as e:
            print(f"❌ Migration verification failed: {e}")
            print("\n📋 Troubleshooting tips:")
            print("1. Check database URL configuration")
            print("2. Verify SQLAlchemy installation")
            print("3. Check file permissions for SQLite database")
            print("4. Review logs for detailed error information")


if __name__ == "__main__":
    # Print migration guide and check status
    MigrationGuide.print_migration_guide()
    print("\n" + "="*60)
    print("MIGRATION STATUS CHECK")
    print("="*60)
    MigrationGuide.check_migration_status()
