"""Start the alignment GUI: `python -m msldprep.gui [--port 8765]`."""

from __future__ import annotations

import argparse
import getpass
import os
import secrets
import socket


def _free_port(host: str) -> int:
    """A port nobody is listening on. The login node is shared, so a fixed
    default would often collide with another user's server."""
    with socket.socket(socket.AF_INET6 if ":" in host else socket.AF_INET) as probe:
        probe.bind((host, 0))
        return probe.getsockname()[1]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="python -m msldprep.gui", description=__doc__)
    parser.add_argument("--host", default="127.0.0.1", help="interface to bind (default: localhost only)")
    parser.add_argument("--port", type=int, default=0, help="port to listen on (default: any free port)")
    parser.add_argument("--token", default=None, help="access token (default: a random one)")
    parser.add_argument(
        "--reload", action="store_true",
        help="restart the server when a Python file of msldprep changes (for development; "
        "the reference, analogs and results are lost on each restart)",
    )
    args = parser.parse_args(argv)

    try:
        import uvicorn

        from .app import TOKEN_ENV, create_app
    except ImportError as exc:
        raise SystemExit(
            f"The GUI needs the web packages ({exc}). Install them with: pip install fastapi uvicorn"
        ) from exc

    port = args.port or _free_port(args.host)
    token = args.token or secrets.token_urlsafe(16)
    if args.host not in ("127.0.0.1", "localhost", "::1"):
        print(f"WARNING: binding {args.host}: anyone who can reach this port and has the token can use it.")
    print()
    print("Alignment GUI is starting. From your own computer, open a tunnel:")
    print(f"  ssh -L {port}:127.0.0.1:{port} {getpass.getuser()}@{socket.gethostname()}")
    print("then open this address in your browser:")
    print(f"  http://127.0.0.1:{port}/?token={token}")
    print("Nothing is saved on the server; use the download button in the page. Ctrl-C stops it.", flush=True)
    if args.reload:
        print("Reloading on code changes: after one, load the reference and analogs again.", flush=True)
    print(flush=True)
    if args.reload:
        # The reloaded server is a new process that imports the app again; it reads the
        # token from the environment, which (unlike the command line) other users cannot see.
        os.environ[TOKEN_ENV] = token
        package_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        uvicorn.run(
            "msldprep.gui.app:create_app_from_env", factory=True, reload=True, reload_dirs=[package_dir],
            host=args.host, port=port, log_level="warning", access_log=False,
        )
        return
    uvicorn.run(create_app(token), host=args.host, port=port, log_level="warning", access_log=False)


if __name__ == "__main__":
    main()
