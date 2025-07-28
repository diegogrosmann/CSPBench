#!/usr/bin/env python3
"""
Test script to verify config_parser syntax without external dependencies
"""

import ast
import sys
from pathlib import Path

def check_syntax(file_path):
    """Check if a Python file has correct syntax"""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            source = f.read()
        
        # Parse the AST to check syntax
        ast.parse(source)
        print(f"✅ {file_path}: Syntax OK")
        return True
    except SyntaxError as e:
        print(f"❌ {file_path}: Syntax Error - {e}")
        return False
    except Exception as e:
        print(f"⚠️  {file_path}: Error - {e}")
        return False

if __name__ == "__main__":
    files_to_check = [
        "src/application/services/config_parser.py",
        "src/application/services/experiment_service.py"
    ]
    
    all_ok = True
    for file_path in files_to_check:
        if not check_syntax(file_path):
            all_ok = False
    
    if all_ok:
        print("\n🎉 All files have correct syntax!")
        sys.exit(0)
    else:
        print("\n💥 Some files have syntax errors!")
        sys.exit(1)
