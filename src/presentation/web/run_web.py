#!/usr/bin/env python3
"""
Web Interface Launcher for CSPBench Development and Production.

This script provides a convenient and configurable way to start the CSPBench web
interface with proper environment configuration, security settings, and development
features. It follows CSPBench guidelines for environment management and includes
comprehensive logging and error handling.
"""

import os
import sys
from pathlib import Path

# Add project root to Python path for development
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))


def main():
    """Launch the CSPBench web interface with comprehensive configuration.

    This function initializes and starts the web interface using uvicorn server
    with environment-based configuration and proper error handling for both
    development and production scenarios.

    Environment Variables:
        WEB_HOST (str): Host address to bind to (default: "0.0.0.0")
        PORT (str): Primary port configuration (default: "8000")
        WEB_PORT (str): Alternative port configuration (fallback)
        DEBUG (str): Enable debug mode with detailed logging (default: "false")
        RELOAD (str): Enable auto-reload for development (default: follows DEBUG)
        LOG_LEVEL (str): Uvicorn logging level (default: follows DEBUG)

    Features Enabled:
        - Algorithm execution type selection interface
        - Modular component-based UI architecture
        - Dataset upload, management, and generation tools
        - Real-time execution monitoring with WebSocket
        - Comprehensive results visualization and analysis
        - ZIP file export and bulk download capabilities
        - Batch configuration management and execution
        - Live progress tracking and status updates

    Raises:
        ImportError: If required dependencies (uvicorn, FastAPI) are not installed
        SystemExit: If configuration errors or startup failures occur

    Example:
        Run with default configuration::

            python src/presentation/web/run_web.py

        Run with custom configuration::

            WEB_HOST=localhost PORT=3000 DEBUG=true python run_web.py

        Production deployment::

            WEB_HOST=0.0.0.0 PORT=80 DEBUG=false python run_web.py

    Note:
        In production environments, consider using a proper WSGI/ASGI server
        like gunicorn with uvicorn workers for better performance and stability.
    """
    try:
        import uvicorn

        from src.presentation.web.app import app

        # Configuration with environment variable fallbacks and validation
        host = os.getenv("WEB_HOST", "0.0.0.0")
        port = int(os.getenv("PORT", os.getenv("WEB_PORT", "8000")))
        debug = os.getenv("DEBUG", "false").lower() == "true"
        reload = os.getenv("RELOAD", str(debug).lower()).lower() == "true"

        # Validate port range
        if not (1 <= port <= 65535):
            print(f"❌ Invalid port number: {port}. Must be between 1 and 65535.")
            sys.exit(1)

        # Log configuration with appropriate level
        log_level = "debug" if debug else "info"

        # Startup information display
        print("🚀 Starting CSPBench Web Interface")
        print("=" * 50)
        print(f"📍 Host: {host}")
        print(f"🔌 Port: {port}")
        print(f"🐛 Debug Mode: {debug}")
        print(f"🔄 Auto Reload: {reload}")
        print(f"📊 Log Level: {log_level}")
        print(f"📂 Working Directory: {project_root}")
        print(f"🌐 Access URL: http://localhost:{port}")
        print("\n🎯 Available Features:")
        print("   ├─ Algorithm execution type selection")
        print("   ├─ Modular component-based interface")
        print("   ├─ Dataset upload, samples, and generation")
        print("   ├─ Batch configuration management")
        print("   ├─ Real-time execution monitoring (WebSocket)")
        print("   ├─ Interactive progress visualization")
        print("   ├─ Results analysis and comparison")
        print("   └─ ZIP export and bulk download")
        print("=" * 50)
        print("📡 WebSocket Endpoints:")
        print(f"   ├─ General: ws://localhost:{port}/ws/{{client_id}}")
        print(f"   └─ Work Monitor: ws://localhost:{port}/ws/work/{{work_id}}")
        print("=" * 50)

        # Configure reload directories for development
        reload_dirs = [str(project_root / "src")] if reload else None

        # Start the uvicorn server with optimized configuration
        uvicorn.run(
            app=app,
            host=host,
            port=port,
            reload=reload,
            log_level=log_level,
            access_log=True,
            use_colors=True,
            reload_dirs=reload_dirs,
            # Additional production-ready configurations
            loop="auto",
            http="auto",
            ws_ping_interval=20,
            ws_ping_timeout=20,
        )

    except ImportError as e:
        print(f"❌ Import Error: {e}")
        print("\n💡 Required dependencies are missing. Please install them:")
        print("   pip install fastapi uvicorn[standard]")
        print("   # OR")
        print("   pip install -r requirements.web.txt")
        print("\n🔍 For development setup:")
        print("   pip install -e .")
        sys.exit(1)

    except ValueError as e:
        print(f"❌ Configuration Error: {e}")
        print("💡 Check your environment variables for correct format")
        sys.exit(1)

    except PermissionError as e:
        print(f"❌ Permission Error: {e}")
        print(
            "💡 Try using a different port (>1024) or run with appropriate permissions"
        )
        sys.exit(1)

    except Exception as e:
        print(f"❌ Startup Failed: {e}")
        print("💡 Check the logs above for detailed error information")
        print("🔧 Common solutions:")
        print("   - Ensure the port is not already in use")
        print("   - Check firewall settings")
        print("   - Verify all dependencies are installed")
        print("   - Check file system permissions")
        sys.exit(1)


if __name__ == "__main__":
    main()
