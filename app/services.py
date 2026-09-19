"""批次计算编排：枚举全部自二聚体与两两异二聚体并落库。"""

from __future__ import annotations

import json
import uuid

from . import db
from .dimer import PairResult, analyze_pair
from .pooling import build_pool_plan


def _pair_to_row(pair_index: int, r: PairResult) -> tuple:
    return (
        pair_index,
        r.name_top,
        r.name_bottom,
        1 if r.is_self_dimer else 0,
        r.alignments_checked,
        r.best_shift,
        r.overlap_len,
        r.max_complementary_run,
        r.top3_run,
        r.bottom3_run,
        r.terminal_run,
        json.dumps(list(r.involved_three_prime_ends), ensure_ascii=False),
        r.risk_level,
        json.dumps(list(r.risk_reasons), ensure_ascii=False),
        r.top_line,
        r.match_line,
        r.bottom_line,
    )


def run_batch(
    *,
    primers: list[tuple[str, str]],
    continuous_threshold: int,
    three_prime_threshold: int,
) -> str:
    """计算一个批次。

    primers 为 (name, normalized_sequence) 列表。
    枚举内容：N 个自二聚体 + N*(N-1)/2 个无序异二聚体。
    """
    batch_id = str(uuid.uuid4())
    created_at = db.utc_now_iso()

    results: list[PairResult] = []
    n = len(primers)

    # 自二聚体
    for name, seq in primers:
        results.append(
            analyze_pair(
                name, seq, name, seq,
                continuous_threshold, three_prime_threshold,
                is_self_dimer=True,
            )
        )
    # 两两异二聚体（无序对；top 固定为提交顺序中靠前者）
    for i in range(n):
        for j in range(i + 1, n):
            name_i, seq_i = primers[i]
            name_j, seq_j = primers[j]
            results.append(
                analyze_pair(
                    name_i, seq_i, name_j, seq_j,
                    continuous_threshold, three_prime_threshold,
                    is_self_dimer=False,
                )
            )

    risk_counts = {"high": 0, "medium": 0, "low": 0}
    for r in results:
        risk_counts[r.risk_level] += 1

    pair_rows = [
        _pair_to_row(idx, r) for idx, r in enumerate(results)
    ]
    primer_rows = [
        (pos, name, seq) for pos, (name, seq) in enumerate(primers)
    ]

    with db.get_conn() as conn:
        db.insert_batch(
            conn,
            batch_id=batch_id,
            created_at=created_at,
            continuous_threshold=continuous_threshold,
            three_prime_threshold=three_prime_threshold,
            primer_count=n,
            pair_count=len(results),
            risk_counts=risk_counts,
            primers=primer_rows,
            pair_rows=pair_rows,
        )
    return batch_id


def plan_batch_pools(
    *,
    batch_id: str,
    pool_count: int,
    pool_capacity: int,
    isolated_risk_levels: list[str],
) -> dict:
    """基于已保存批次的异二聚体结果构建冲突图并确定性分池。

    批次不存在返回 None（由路由层转 404）；容量总和不足抛出
    ValueError（由路由层转 422）。
    """
    with db.get_conn() as conn:
        batch = db.get_batch(conn, batch_id)
        if batch is None:
            return None
        primers = db.get_primers(conn, batch_id)
        pairs = db.get_pairs(conn, batch_id, risky_only=False)

    n = len(primers)
    if pool_count * pool_capacity < n:
        raise ValueError(
            f"池容量总和 {pool_count}×{pool_capacity}="
            f"{pool_count * pool_capacity} 不足以容纳 {n} 条引物"
        )

    # 批次内引物名称唯一（入库时已保证）
    position_by_name = {row["name"]: row["position"] for row in primers}
    names_by_position = {pos: name for name, pos in position_by_name.items()}

    # 仅取异二聚体（忽略自二聚体），且风险等级在隔离集合内
    levels = set(isolated_risk_levels)
    edges: list[tuple[int, int, str]] = []
    for row in pairs:
        if row["is_self_dimer"] or row["risk_level"] not in levels:
            continue
        edges.append(
            (
                position_by_name[row["name_top"]],
                position_by_name[row["name_bottom"]],
                row["risk_level"],
            )
        )

    plan = build_pool_plan(
        primer_count=n,
        edges=edges,
        pool_count=pool_count,
        pool_capacity=pool_capacity,
    )

    # 每个池的冲突边数
    conflict_count_by_pool = [0] * pool_count
    conflicts_out: list[dict] = []
    for pool_id, edge in plan.conflicts:
        conflict_count_by_pool[pool_id] += 1
        conflicts_out.append(
            {
                "pool": pool_id,
                "position_a": edge.position_a,
                "position_b": edge.position_b,
                "name_a": names_by_position[edge.position_a],
                "name_b": names_by_position[edge.position_b],
                "risk_level": edge.risk_level,
            }
        )

    pools_out = [
        {
            "pool": pool_id,
            "capacity": pool_capacity,
            "members": [names_by_position[pos] for pos in members],
            "member_positions": members,
            "conflict_count": conflict_count_by_pool[pool_id],
        }
        for pool_id, members in enumerate(plan.pool_members)
    ]

    return {
        "batch_id": batch_id,
        "pool_count": pool_count,
        "pool_capacity": pool_capacity,
        "isolated_risk_levels": isolated_risk_levels,
        "primer_count": n,
        "assigned_count": sum(len(m) for m in plan.pool_members),
        "total_conflicts": plan.total_conflicts,
        "pools": pools_out,
        "conflicts": conflicts_out,
    }
