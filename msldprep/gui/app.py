"""FastAPI shell over `service.AlignmentSession`.

One session per server, protected by a random access token: the page is
served only with `?token=...`, and every `/api/` call must carry it (header
`X-Token`, or `token=` in the query string for downloads). Requires the
`[web]` extra.
"""

from __future__ import annotations

import logging
import os
import secrets
import threading

from fastapi import Depends, FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from ..errors import MsldPrepError
from .service import AlignmentSession, SessionError

logger = logging.getLogger(__name__)

STATIC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")
MAX_TEXT = 20_000_000  # characters of one upload
TOKEN_ENV = "MSLDPREP_GUI_TOKEN"  # how `python -m msldprep.gui --reload` hands the token to the reloaded server


class ReferenceIn(BaseModel):
    name: str = Field(max_length=200)
    text: str = Field(max_length=MAX_TEXT)


class AnalogsIn(BaseModel):
    text: str = Field(max_length=MAX_TEXT)


class ParamsIn(BaseModel):
    complete_rings_only: bool | None = None
    ring_matches_ring_only: bool | None = None
    ignore_bond_orders: bool | None = None
    check_protonation: bool | None = None
    timeout: int | None = Field(default=None, ge=1, le=3600)
    seed: int | None = Field(default=None, ge=0, le=2**31 - 1)


class RealignIn(BaseModel):
    seed: int = Field(ge=0, le=2**31 - 1)


class DecisionIn(BaseModel):
    decision: str


def create_app(token: str, session: AlignmentSession | None = None) -> FastAPI:
    session = session or AlignmentSession()
    app = FastAPI(title="msld-prep alignment", docs_url=None, redoc_url=None, openapi_url=None)

    def check_token(request: Request) -> None:
        supplied = request.headers.get("x-token") or request.query_params.get("token") or ""
        if not secrets.compare_digest(supplied.encode(), token.encode()):
            raise HTTPException(status_code=403, detail="Missing or wrong access token.")

    auth = [Depends(check_token)]

    @app.exception_handler(MsldPrepError)
    async def _msldprep_error(_request: Request, exc: MsldPrepError) -> JSONResponse:
        status = 404 if isinstance(exc, SessionError) and "No analog named" in str(exc) else 422
        return JSONResponse(status_code=status, content={"detail": str(exc)})

    @app.get("/", dependencies=auth, include_in_schema=False)
    def index() -> FileResponse:
        return FileResponse(os.path.join(STATIC_DIR, "index.html"), headers={"Cache-Control": "no-store"})

    @app.get("/api/state", dependencies=auth)
    def state() -> dict:
        return session.summary()

    @app.put("/api/reference", dependencies=auth)
    def put_reference(body: ReferenceIn) -> dict:
        session.set_reference(body.name, body.text)
        return session.summary()

    @app.get("/api/reference", dependencies=auth)
    def get_reference() -> dict:
        return session.reference_view()

    @app.put("/api/analogs", dependencies=auth)
    def put_analogs(body: AnalogsIn) -> dict:
        session.set_analogs(body.text)
        return session.summary()

    @app.put("/api/params", dependencies=auth)
    def put_params(body: ParamsIn) -> dict:
        session.set_params(**body.model_dump())
        return session.summary()

    @app.post("/api/align", dependencies=auth, status_code=202)
    def align() -> dict:
        run = session.begin_alignment()
        threading.Thread(target=_run_logged, args=(run,), daemon=True).start()
        return session.summary()

    @app.get("/api/analogs/{name}", dependencies=auth)
    def get_analog(name: str) -> dict:
        return session.ligand_view(name)

    @app.post("/api/analogs/{name}/realign", dependencies=auth)
    def realign(name: str, body: RealignIn) -> dict:
        session.realign(name, body.seed)
        return session.ligand_view(name)

    @app.post("/api/analogs/{name}/decision", dependencies=auth)
    def decision(name: str, body: DecisionIn) -> dict:
        session.set_decision(name, body.decision)
        return session.summary()

    @app.get("/api/overlay", dependencies=auth)
    def overlay() -> list[dict]:
        return session.overlay()

    @app.get("/api/export", dependencies=auth)
    def export(
        format: str = "sdf", which: str = "aligned", reference: bool = True, per_file: bool = False, keep_names: bool = False
    ) -> Response:
        result = session.export(format, which, include_reference=reference, one_per_file=per_file, keep_names=keep_names)
        return Response(
            result.data,
            media_type=result.media_type,
            headers={"Content-Disposition": f'attachment; filename="{result.filename}"'},
        )

    @app.post("/api/reset", dependencies=auth)
    def reset() -> dict:
        session.reset()
        return session.summary()

    @app.middleware("http")
    async def _revalidate_static(request: Request, call_next):
        # The browser may keep the page's scripts but must ask whether they changed,
        # so an edited app.js or style.css shows on a plain refresh.
        response = await call_next(request)
        if request.url.path.startswith("/static/"):
            response.headers["Cache-Control"] = "no-cache"
        return response

    app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
    return app


def create_app_from_env() -> FastAPI:
    """`create_app` with the token from the environment: the factory a reloading
    server imports again after each change to the code."""
    token = os.environ.get(TOKEN_ENV)
    if not token:
        raise RuntimeError(f"{TOKEN_ENV} is not set.")
    return create_app(token)


def _run_logged(run) -> None:
    try:
        run()
    except Exception:  # the thread must not die silently; the job flag is cleared by `run`
        logger.exception("Alignment job failed")
