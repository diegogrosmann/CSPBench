"""Database recovery utilities for handling SQLite locks and corruption."""

import logging
import sqlite3
import time
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)


class DatabaseRecovery:
    """Utilities for recovering from database lock and corruption issues."""

    @staticmethod
    def check_database_locks(db_path: Path) -> bool:
        """
        Check if database is locked or corrupted.
        
        Args:
            db_path: Path to the SQLite database
            
        Returns:
            True if database is accessible, False if locked/corrupted
        """
        if not db_path.exists():
            logger.warning(f"Database file does not exist: {db_path}")
            return False
            
        try:
            conn = sqlite3.connect(db_path, timeout=5.0)
            cursor = conn.cursor()
            cursor.execute("SELECT 1")
            cursor.fetchone()
            conn.close()
            return True
        except sqlite3.OperationalError as e:
            if "database is locked" in str(e).lower():
                logger.error(f"Database is locked: {db_path}")
                return False
            elif "database disk image is malformed" in str(e).lower():
                logger.error(f"Database is corrupted: {db_path}")
                return False
            else:
                logger.error(f"Database error: {e}")
                return False
        except Exception as e:
            logger.error(f"Unexpected error checking database: {e}")
            return False

    @staticmethod
    def force_unlock_database(db_path: Path) -> bool:
        """
        Force unlock a SQLite database by removing WAL and SHM files.
        
        WARNING: Only use this when you're sure no processes are using the database!
        
        Args:
            db_path: Path to the SQLite database
            
        Returns:
            True if unlock was attempted, False if database doesn't exist
        """
        if not db_path.exists():
            logger.warning(f"Database file does not exist: {db_path}")
            return False
            
        wal_file = db_path.with_suffix(db_path.suffix + "-wal")
        shm_file = db_path.with_suffix(db_path.suffix + "-shm")
        
        removed_files = []
        
        # Remove WAL file if exists
        if wal_file.exists():
            try:
                wal_file.unlink()
                removed_files.append("WAL")
                logger.info(f"Removed WAL file: {wal_file}")
            except Exception as e:
                logger.error(f"Failed to remove WAL file: {e}")
                
        # Remove SHM file if exists
        if shm_file.exists():
            try:
                shm_file.unlink()
                removed_files.append("SHM")
                logger.info(f"Removed SHM file: {shm_file}")
            except Exception as e:
                logger.error(f"Failed to remove SHM file: {e}")
        
        if removed_files:
            logger.info(f"Force unlock attempted by removing {', '.join(removed_files)} files")
            
            # Wait a moment and check if database is now accessible
            time.sleep(0.5)
            if DatabaseRecovery.check_database_locks(db_path):
                logger.info("Database is now accessible")
                return True
            else:
                logger.warning("Database is still not accessible after force unlock")
                return False
        else:
            logger.info("No WAL/SHM files found - database may not be locked by WAL mode")
            return True

    @staticmethod
    def backup_and_rebuild_database(db_path: Path, backup_suffix: str = ".backup") -> bool:
        """
        Create a backup and attempt to rebuild a corrupted database.
        
        Args:
            db_path: Path to the SQLite database
            backup_suffix: Suffix for backup file
            
        Returns:
            True if rebuild was successful, False otherwise
        """
        if not db_path.exists():
            logger.warning(f"Database file does not exist: {db_path}")
            return False
            
        backup_path = db_path.with_suffix(db_path.suffix + backup_suffix)
        
        try:
            # Create backup
            import shutil
            shutil.copy2(db_path, backup_path)
            logger.info(f"Created backup: {backup_path}")
            
            # Try to rebuild using VACUUM
            conn = sqlite3.connect(db_path, timeout=30.0)
            conn.execute("VACUUM")
            conn.close()
            
            logger.info("Database rebuilt successfully with VACUUM")
            return True
            
        except Exception as e:
            logger.error(f"Failed to rebuild database: {e}")
            
            # If VACUUM failed, try to dump and restore
            try:
                logger.info("Attempting dump and restore recovery...")
                temp_sql = db_path.with_suffix(".recovery.sql")
                
                # Dump to SQL
                with open(temp_sql, 'w') as f:
                    conn = sqlite3.connect(backup_path)
                    for line in conn.iterdump():
                        f.write('%s\n' % line)
                    conn.close()
                
                # Remove corrupted database
                db_path.unlink()
                
                # Restore from SQL
                conn = sqlite3.connect(db_path)
                with open(temp_sql, 'r') as f:
                    conn.executescript(f.read())
                conn.close()
                
                # Clean up
                temp_sql.unlink()
                
                logger.info("Database recovered using dump and restore")
                return True
                
            except Exception as e2:
                logger.error(f"Dump and restore also failed: {e2}")
                return False

    @staticmethod
    def get_database_status(db_path: Path) -> dict:
        """
        Get comprehensive status of a SQLite database.
        
        Args:
            db_path: Path to the SQLite database
            
        Returns:
            Dictionary with database status information
        """
        status = {
            "file_exists": db_path.exists() if db_path else False,
            "file_size": 0,
            "is_accessible": False,
            "wal_file_exists": False,
            "shm_file_exists": False,
            "wal_file_size": 0,
            "integrity_check": "unknown",
            "journal_mode": "unknown",
            "error": None
        }
        
        if not db_path or not db_path.exists():
            status["error"] = "Database file does not exist"
            return status
            
        try:
            # File size
            status["file_size"] = db_path.stat().st_size
            
            # Check WAL/SHM files
            wal_file = db_path.with_suffix(db_path.suffix + "-wal")
            shm_file = db_path.with_suffix(db_path.suffix + "-shm")
            
            status["wal_file_exists"] = wal_file.exists()
            status["shm_file_exists"] = shm_file.exists()
            
            if wal_file.exists():
                status["wal_file_size"] = wal_file.stat().st_size
            
            # Try to connect and get info
            conn = sqlite3.connect(db_path, timeout=5.0)
            cursor = conn.cursor()
            
            # Basic accessibility
            cursor.execute("SELECT 1")
            cursor.fetchone()
            status["is_accessible"] = True
            
            # Journal mode
            cursor.execute("PRAGMA journal_mode")
            status["journal_mode"] = cursor.fetchone()[0]
            
            # Integrity check (quick)
            cursor.execute("PRAGMA quick_check")
            result = cursor.fetchone()[0]
            status["integrity_check"] = result
            
            conn.close()
            
        except sqlite3.OperationalError as e:
            status["error"] = str(e)
            if "database is locked" in str(e).lower():
                status["is_accessible"] = False
            elif "database disk image is malformed" in str(e).lower():
                status["integrity_check"] = "corrupted"
        except Exception as e:
            status["error"] = str(e)
            
        return status

    @staticmethod
    def find_and_check_databases(base_path: Path) -> List[dict]:
        """
        Find all SQLite databases in a directory tree and check their status.
        
        Args:
            base_path: Base directory to search
            
        Returns:
            List of database status dictionaries
        """
        databases = []
        
        for db_file in base_path.rglob("*.db"):
            status = DatabaseRecovery.get_database_status(db_file)
            status["path"] = str(db_file)
            databases.append(status)
            
        return databases


def main():
    """CLI utility for database recovery operations."""
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(description="SQLite Database Recovery Utility")
    parser.add_argument("database", help="Path to SQLite database")
    parser.add_argument("--check", action="store_true", help="Check database status")
    parser.add_argument("--unlock", action="store_true", help="Force unlock database")
    parser.add_argument("--rebuild", action="store_true", help="Backup and rebuild database")
    parser.add_argument("--find-all", action="store_true", help="Find all databases in directory")
    
    args = parser.parse_args()
    
    db_path = Path(args.database)
    
    if args.find_all:
        databases = DatabaseRecovery.find_and_check_databases(db_path)
        for db_info in databases:
            print(f"\nDatabase: {db_info['path']}")
            print(f"  Accessible: {db_info['is_accessible']}")
            print(f"  Size: {db_info['file_size']} bytes")
            print(f"  Journal Mode: {db_info['journal_mode']}")
            print(f"  Integrity: {db_info['integrity_check']}")
            if db_info['error']:
                print(f"  Error: {db_info['error']}")
    
    elif args.check:
        status = DatabaseRecovery.get_database_status(db_path)
        print(f"Database Status for {db_path}:")
        for key, value in status.items():
            print(f"  {key}: {value}")
    
    elif args.unlock:
        print(f"Force unlocking database: {db_path}")
        success = DatabaseRecovery.force_unlock_database(db_path)
        if success:
            print("✅ Database unlock completed")
            sys.exit(0)
        else:
            print("❌ Database unlock failed")
            sys.exit(1)
    
    elif args.rebuild:
        print(f"Rebuilding database: {db_path}")
        success = DatabaseRecovery.backup_and_rebuild_database(db_path)
        if success:
            print("✅ Database rebuild completed")
            sys.exit(0)
        else:
            print("❌ Database rebuild failed")
            sys.exit(1)
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
