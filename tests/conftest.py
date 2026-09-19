import os
import tempfile

import pytest

# 在任何 app 模块被收集导入前，指向临时数据库
_TMP = tempfile.mkdtemp(prefix="primer_test_")
os.environ.setdefault("PRIMER_DB_PATH", os.path.join(_TMP, "test.db"))


@pytest.fixture(scope="session")
def client():
    from app import db, main

    db.init_db()
    from fastapi.testclient import TestClient

    with TestClient(main.app) as c:
        yield c
