#!/usr/bin/env python3
"""
Script para limpar completamente o banco de dados do CSPBench (modo não-interativo).
Remove todas as tabelas e dados automaticamente, sem confirmação.
"""

import os
import sys
from pathlib import Path
from sqlalchemy import create_engine, text
from dotenv import load_dotenv

# Adicionar o diretório raiz ao PYTHONPATH
root_dir = Path(__file__).parent.parent
sys.path.insert(0, str(root_dir))

def clean_database_force():
    """Limpa completamente o banco de dados sem confirmação."""
    
    # Carregar variáveis de ambiente
    load_dotenv()
    
    db_url = os.getenv('DB_URL')
    if not db_url:
        print("❌ Erro: DB_URL não encontrada no arquivo .env")
        return False
    
    try:
        print(f"🔗 Conectando ao banco: {db_url.split('@')[1] if '@' in db_url else db_url}")
        
        # Criar engine de conexão
        engine = create_engine(db_url)
        
        with engine.connect() as connection:
            # Descobrir o tipo de banco
            is_sqlite = 'sqlite' in db_url
            is_postgres = 'postgresql' in db_url
            
            if is_sqlite:
                # Para SQLite
                result = connection.execute(
                    text("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%';")
                )
                tables = [row[0] for row in result.fetchall()]
                
            elif is_postgres:
                # Para PostgreSQL
                result = connection.execute(
                    text("SELECT tablename FROM pg_tables WHERE schemaname='public';")
                )
                tables = [row[0] for row in result.fetchall()]
                
            else:
                print("❌ Tipo de banco não suportado")
                return False
            
            if not tables:
                print("✅ Nenhuma tabela encontrada. Banco já está limpo!")
                return True
            
            print(f"🧹 Limpando {len(tables)} tabelas automaticamente...")
            
            # Remover tabelas
            if is_postgres:
                # Para PostgreSQL
                for table in tables:
                    try:
                        connection.execute(text(f'DROP TABLE IF EXISTS "{table}" CASCADE;'))
                        print(f"   ✓ {table}")
                    except Exception as e:
                        print(f"   ❌ Erro: {table} - {e}")
            else:
                # Para SQLite
                connection.execute(text('PRAGMA foreign_keys = OFF;'))
                for table in tables:
                    try:
                        connection.execute(text(f'DROP TABLE IF EXISTS "{table}";'))
                        print(f"   ✓ {table}")
                    except Exception as e:
                        print(f"   ❌ Erro: {table} - {e}")
                connection.execute(text('PRAGMA foreign_keys = ON;'))
            
            # Commit das mudanças
            connection.commit()
            print("✅ Banco de dados limpo com sucesso!")
            return True
            
    except Exception as e:
        print(f"❌ Erro ao limpar banco de dados: {e}")
        return False

if __name__ == "__main__":
    success = clean_database_force()
    sys.exit(0 if success else 1)
