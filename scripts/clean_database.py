#!/usr/bin/env python3
"""
Script para limpar completamente o banco de dados do CSPBench.
Remove todas as tabelas e dados, mantendo apenas a estrutura do banco.
"""

import os
import sys
from pathlib import Path
from sqlalchemy import create_engine, text, MetaData
from dotenv import load_dotenv

# Adicionar o diretório raiz ao PYTHONPATH
root_dir = Path(__file__).parent.parent
sys.path.insert(0, str(root_dir))

def clean_database():
    """Limpa completamente o banco de dados."""
    
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
            print("📋 Listando tabelas existentes...")
            
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
            
            print(f"📊 Encontradas {len(tables)} tabelas:")
            for table in tables:
                print(f"   - {table}")
            
            # Confirmar limpeza
            confirm = input("\n⚠️  Deseja realmente limpar TODAS as tabelas? (digite 'SIM' para confirmar): ")
            if confirm != 'SIM':
                print("❌ Operação cancelada pelo usuário")
                return False
            
            print("\n🧹 Limpando banco de dados...")
            
            # Remover tabelas
            if is_postgres:
                # Para PostgreSQL, desabilitar foreign key checks temporariamente
                for table in tables:
                    try:
                        connection.execute(text(f'DROP TABLE IF EXISTS "{table}" CASCADE;'))
                        print(f"   ✓ Tabela '{table}' removida")
                    except Exception as e:
                        print(f"   ❌ Erro ao remover tabela '{table}': {e}")
            else:
                # Para SQLite
                connection.execute(text('PRAGMA foreign_keys = OFF;'))
                for table in tables:
                    try:
                        connection.execute(text(f'DROP TABLE IF EXISTS "{table}";'))
                        print(f"   ✓ Tabela '{table}' removida")
                    except Exception as e:
                        print(f"   ❌ Erro ao remover tabela '{table}': {e}")
                connection.execute(text('PRAGMA foreign_keys = ON;'))
            
            # Commit das mudanças
            connection.commit()
            
            print("\n✅ Banco de dados limpo com sucesso!")
            print("💡 Execute as migrações para recriar as tabelas necessárias")
            
            return True
            
    except Exception as e:
        print(f"❌ Erro ao limpar banco de dados: {e}")
        return False

if __name__ == "__main__":
    print("🗃️  CSPBench - Limpeza de Banco de Dados")
    print("=" * 50)
    
    success = clean_database()
    
    if success:
        print("\n🎉 Limpeza concluída com sucesso!")
        sys.exit(0)
    else:
        print("\n💥 Falha na limpeza do banco de dados!")
        sys.exit(1)
