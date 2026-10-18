"""The committed Streamlit config binds the app to the loopback interface (audit C-16)."""

import tomllib
from pathlib import Path

CONFIG = Path(__file__).parent.parent / ".streamlit" / "config.toml"


def test_server_listens_on_localhost_only():
    cfg = tomllib.loads(CONFIG.read_text())
    assert cfg["server"]["address"] == "localhost"
    assert cfg["server"]["maxUploadSize"] == 200
    # headless is left to the command line: locally the browser should still open
    assert "headless" not in cfg["server"]
    assert cfg["theme"]["primaryColor"] == "#1F6FEB"  # theme kept
