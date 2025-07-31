#!/usr/bin/env python3
"""
Web Interface Launcher for CSPBench

This script provides a convenient way to start the web interface during development.
It follows CSPBench guidelines for environment configuration and security.
"""

import os
import sys
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

def main():
    """
    Launch the CSPBench web interface with proper configuration.
    
    Environment variables:
        WEB_HOST: Host to bind to (default: 0.0.0.0)
        WEB_PORT: Port to bind to (default: 8000)
        DEBUG: Enable debug mode (default: false)
        RELOAD: Enable auto-reload (default: true in debug)
    """
    try:
        import uvicorn
        from src.presentation.web.app import app
        
        # Configuration with environment variable fallbacks
        host = os.getenv("WEB_HOST", "0.0.0.0")
        port = int(os.getenv("WEB_PORT", "8000"))
        debug = os.getenv("DEBUG", "false").lower() == "true"
        reload = os.getenv("RELOAD", str(debug).lower()).lower() == "true"
        
        # Log configuration
        log_level = "debug" if debug else "info"
        
        print(f"🚀 Starting CSPBench Web Interface")
        print(f"📍 Host: {host}")
        print(f"🔌 Port: {port}")
        print(f"🐛 Debug: {debug}")
        print(f"🔄 Reload: {reload}")
        print(f"📊 Log Level: {log_level}")
        print(f"📂 Working Directory: {project_root}")
        print(f"🌐 URL: http://localhost:{port}")
        print("📊 Features:")
        print("   - Algorithm execution type selection")
        print("   - Modular component-based interface")
        print("   - Dataset upload, samples, and generation")
        print("   - Real-time execution monitoring")
        print("   - Results visualization and download")
        print("   - ZIP file export")
        print("=" * 50)
        
        # Start the web server
        uvicorn.run(
            app=app,
            host=host,
            port=port,
            reload=reload,
            log_level=log_level,
            access_log=True,
            use_colors=True,
            reload_dirs=[str(project_root / "src")] if reload else None
        )
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        print("💡 Make sure all dependencies are installed:")
        print("   pip install -r requirements.web.txt")
        sys.exit(1)
        
    except Exception as e:
        print(f"❌ Failed to start web interface: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
