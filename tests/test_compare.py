"""批次风险对比接口测试：GET /api/batches/{base_id}/compare/{target_id}。"""

from __future__ import annotations


# base：3′阈值 10（极高），无 high；A-C、A-D、C-D 为 medium
BASE_PAYLOAD = {
    "primers": [
        {"name": "A", "sequence": "GGGGGAAAAA"},
        {"name": "B", "sequence": "ACGTACGTAC"},
        {"name": "C", "sequence": "TTTTGGGGGG"},
        {"name": "D", "sequence": "AAAACCCCGG"},
    ],
    "continuous_threshold": 4,
    "three_prime_threshold": 10,
}

# target：3′阈值 3；A 序列末尾 A->T（A-D 由 medium 降为 low），
# 删除 C、新增 E
TARGET_PAYLOAD = {
    "primers": [
        {"name": "A", "sequence": "GGGGGAAAAT"},
        {"name": "B", "sequence": "ACGTACGTAC"},
        {"name": "D", "sequence": "AAAACCCCGG"},
        {"name": "E", "sequence": "TTTTGGGGGG"},
    ],
    "continuous_threshold": 4,
    "three_prime_threshold": 3,
}


def _create(client, payload):
    resp = client.post("/api/batches", json=payload)
    assert resp.status_code == 201, resp.text
    return resp.json()["id"]


def test_compare_categories_and_summary(client):
    base_id = _create(client, BASE_PAYLOAD)
    target_id = _create(client, TARGET_PAYLOAD)
    resp = client.get(f"/api/batches/{base_id}/compare/{target_id}")
    assert resp.status_code == 200, resp.text
    data = resp.json()

    assert data["base_batch_id"] == base_id
    assert data["target_batch_id"] == target_id

    by_cat = {
        cat: [(d["pair_name"], d["base_risk_level"], d["target_risk_level"])
              for d in data[cat]]
        for cat in ("escalated", "deescalated", "introduced", "resolved")
    }

    # B、D 自二聚体 medium -> high
    assert by_cat["escalated"] == [
        ("B", "medium", "high"),
        ("D", "medium", "high"),
    ]
    # A-D medium -> low，且 A 序列变化（C-D 因 C 被删除属于 resolved）
    assert by_cat["deescalated"] == [("A|D", "medium", "low")]
    # 新增引物 E 带来的全部配对
    assert by_cat["introduced"] == [
        ("A|E", None, "medium"),
        ("B|E", None, "low"),
        ("D|E", None, "high"),
        ("E", None, "low"),
    ]
    # 删除引物 C 后消失的全部配对
    assert by_cat["resolved"] == [
        ("A|C", "medium", None),
        ("B|C", "low", None),
        ("C", "low", None),
        ("C|D", "medium", None),
    ]

    summary = data["summary"]
    assert summary == {
        "escalated": 2,
        "deescalated": 1,
        "introduced": 4,
        "resolved": 4,
        "unchanged": 3,  # A、A|B、B|D
        "base_pair_count": 10,
        "target_pair_count": 10,
    }

    for cat, details in (
        ("escalated", data["escalated"]),
        ("deescalated", data["deescalated"]),
        ("introduced", data["introduced"]),
        ("resolved", data["resolved"]),
    ):
        names = [d["pair_name"] for d in details]
        assert names == sorted(names), f"{cat} 未按配对名排序"
        for d in details:
            assert d["is_self_dimer"] == ("|" not in d["pair_name"])


def test_sequence_changed_flag(client):
    base_id = _create(client, BASE_PAYLOAD)
    target_id = _create(client, TARGET_PAYLOAD)
    data = client.get(
        f"/api/batches/{base_id}/compare/{target_id}"
    ).json()

    changed = {
        (cat, d["pair_name"])
        for cat in ("escalated", "deescalated", "introduced", "resolved")
        for d in data[cat]
        if d["sequence_changed"]
    }
    # 仅引物 A 的序列变化；涉及 A 的配对（A 自身风险仍为 low、A-B 仍为 low，
    # 二者归为 unchanged 不在明细中）为 deescalated 的 A|D、introduced 的 A|E，
    # 以及 resolved 的 A|C（C 被删除，但 A 在两批次间序列变化）
    assert changed == {
        ("deescalated", "A|D"),
        ("introduced", "A|E"),
        ("resolved", "A|C"),
    }

    # 显式核对 A|D 两档明细
    ad = next(d for d in data["deescalated"] if d["pair_name"] == "A|D")
    assert ad["sequence_changed"] is True
    assert ad["base_risk_level"] == "medium"
    assert ad["target_risk_level"] == "low"

    # 不涉及 A 的 resolved/introduced 不标记 sequence_changed
    assert all(
        d["sequence_changed"] is False
        for d in data["resolved"]
        if d["pair_name"] != "A|C"
    )
    assert all(
        d["sequence_changed"] is False
        for d in data["introduced"]
        if d["pair_name"] != "A|E"
    )
    assert all(
        d["sequence_changed"] is False for d in data["escalated"]
    )


def test_compare_hetero_pair_order_independent(client):
    # 同两个引物在两批次中提交顺序相反，应仍按同一配对键比较
    forward = {
        "primers": [
            {"name": "F", "sequence": "GGGGGAAAAA"},
            {"name": "R", "sequence": "GGGGGTTTTT"},
        ],
        "continuous_threshold": 4,
        "three_prime_threshold": 3,
    }
    reverse = {
        "primers": [
            {"name": "R", "sequence": "GGGGGTTTTT"},
            {"name": "F", "sequence": "GGGGGAAAAA"},
        ],
        "continuous_threshold": 4,
        "three_prime_threshold": 3,
    }
    fid = _create(client, forward)
    rid = _create(client, reverse)
    data = client.get(f"/api/batches/{fid}/compare/{rid}").json()

    # 共有配对为 F、R、F|R，全部 unchanged
    assert data["summary"]["unchanged"] == 3
    for cat in ("escalated", "deescalated", "introduced", "resolved"):
        assert data[cat] == []

    # 反向查询（交换 base/target）结论对称
    swapped = client.get(f"/api/batches/{rid}/compare/{fid}").json()
    assert swapped["summary"] == data["summary"]


def test_compare_same_batch_is_empty_diff(client):
    bid = _create(client, TARGET_PAYLOAD)
    data = client.get(f"/api/batches/{bid}/compare/{bid}").json()
    summary = data["summary"]
    assert summary["base_pair_count"] == 10
    assert summary["target_pair_count"] == 10
    assert summary["unchanged"] == 10
    assert summary["escalated"] == 0
    assert summary["deescalated"] == 0
    assert summary["introduced"] == 0
    assert summary["resolved"] == 0


def test_compare_404_when_either_batch_missing(client):
    target_id = _create(client, TARGET_PAYLOAD)

    missing_base = client.get(f"/api/batches/no-such-batch/compare/{target_id}")
    assert missing_base.status_code == 404

    bid = _create(client, BASE_PAYLOAD)
    missing_target = client.get(f"/api/batches/{bid}/compare/no-such-batch")
    assert missing_target.status_code == 404
