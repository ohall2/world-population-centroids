"""Test av download_zenodo.py utan nätverk: avbruten hämtning återupptas
med HTTP Range, md5 kontrolleras och metadata skrivs."""
import hashlib, json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import download_zenodo as dz  # noqa: E402

CONTENT = bytes(range(256)) * 9000          # ca 2,3 MB
MD5 = hashlib.md5(CONTENT).hexdigest()


class Resp:
    def __init__(self, body, status=200, js=None):
        self.body, self.status_code, self._js = body, status, js
    def __enter__(self): return self
    def __exit__(self, *a): pass
    def raise_for_status(self): pass
    def json(self): return self._js
    def iter_content(self, n):
        for i in range(0, len(self.body), n):
            yield self.body[i:i + n]


def fake_get(fail_first: list):
    def get(url, headers=None, stream=False, timeout=None):
        if url.startswith("https://zenodo.org/api/records/") and "files" not in url:
            return Resp(b"", js={"doi": "10.5281/zenodo.1", "metadata": {"title": "T"},
                                 "files": [{"key": "CCNL_DMSP_1999_V1.tif", "size": len(CONTENT),
                                            "checksum": f"md5:{MD5}",
                                            "links": {"self": "https://zenodo.org/api/records/1/files/x/content"}},
                                           {"key": "README.txt", "size": 1, "checksum": "md5:x",
                                            "links": {"self": "u"}}]})
        start = int(headers["Range"][6:-1]) if headers else 0
        if fail_first:
            fail_first.pop()
            body = CONTENT[start:start + 1_000_000]
            class Broken(Resp):
                def iter_content(self, n):
                    yield from super().iter_content(n)
                    raise dz.requests.ConnectionError("avbruten")
            return Broken(body, 206 if start else 200)
        return Resp(CONTENT[start:], 206 if start else 200)
    return get


def test_resume_and_md5(tmp_path, monkeypatch):
    monkeypatch.setattr(dz.requests, "get", fake_get([True]))
    monkeypatch.setattr(dz.time, "sleep", lambda s: None)
    dz.main(["--product", "ccnl", "--out", str(tmp_path)])
    out = tmp_path / "CCNL_DMSP_1999_V1.tif"
    assert out.read_bytes() == CONTENT
    meta = json.loads(out.with_suffix(".json").read_text())
    assert meta["year"] == 1999 and meta["md5"] == MD5
    assert not list(tmp_path.glob("*.part")) and not (tmp_path / "README.txt").exists()
