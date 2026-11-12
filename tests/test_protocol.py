import copy
import json
import subprocess
import sys
from pathlib import Path

import pytest

from regjepa.grouping import stable_digest
from regjepa.protocol import file_sha256, freeze_protocol, load_split_plan, locked_input_sha256, prepare_protocol, verify_inputs, write_json
from regjepa.splits import build_split_plan, validate_split_plan

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def protocol():
    return prepare_protocol(ROOT, json.loads((ROOT / "configs/evaluation.json").read_text()))


def test_actual_cohort_groups_preserve_duplicate_source_and_named_series(protocol):
    plan, grouping, summary = protocol
    assert summary["molecules"] == summary["structures"] == 299
    assert summary["rows"] == 1264
    groups = grouping["compound_to_group"]
    assert len({groups[x] for x in ["comp178", "comp179", "comp180", "comp182"]}) == 1
    assert len({groups[x] for x in ["comp249", "comp250", "comp251"]}) == 1
    duplicate = next(g for g in grouping["groups"] if "comp13" in g["selected_members"])
    assert "comp598" in duplicate["all_source_members"]
    assert summary["maximum_between_group_tanimoto"] < .6
    validate_split_plan(plan["outer_folds"], sorted(groups), groups, plan["config"])


def test_unfrozen_proposal_cannot_be_used_for_an_experiment(protocol, tmp_path):
    plan, _, _ = protocol
    path = tmp_path / "proposal.json"
    write_json(path, plan)
    with pytest.raises(ValueError, match="explicitly frozen"):
        load_split_plan(ROOT, path)
    inspect = subprocess.run([sys.executable, "-B", str(ROOT / "scripts/run_linear_evaluation.py"), "--splits", str(path)],
                             capture_output=True, text=True)
    assert inspect.returncode == 0
    assert json.loads(inspect.stdout)["models_fitted"] == 0
    forbidden = subprocess.run([sys.executable, "-B", str(ROOT / "scripts/run_linear_evaluation.py"), "--splits", str(path), "--run"],
                               capture_output=True, text=True)
    assert forbidden.returncode != 0
    assert "explicitly frozen" in forbidden.stderr


def test_changed_locked_input_is_rejected_without_modifying_the_input(protocol):
    plan, _, _ = protocol
    bad = copy.deepcopy(plan)
    key = next(iter(bad["input_sha256"]))
    bad["input_sha256"][key] = "0" * 64
    with pytest.raises(ValueError, match="Locked input changed"):
        verify_inputs(ROOT, bad)


def test_changed_runtime_versions_require_a_new_proposal(protocol):
    plan, _, _ = protocol
    bad = copy.deepcopy(plan)
    bad["versions"]["numpy"] = "different"
    with pytest.raises(ValueError, match="Runtime versions differ"):
        verify_inputs(ROOT, bad)


def test_derived_text_input_locks_survive_git_line_endings_but_raw_bytes_do_not(tmp_path):
    derived = tmp_path / "data/processed/source.jsonl"
    raw = tmp_path / "data/raw/source.txt"
    derived.parent.mkdir(parents=True)
    raw.parent.mkdir(parents=True)
    derived.write_bytes(b'{"x":1}\r\n')
    raw.write_bytes(b'raw\r\n')
    before = locked_input_sha256(tmp_path, "data/processed/source.jsonl")
    raw_before = locked_input_sha256(tmp_path, "data/raw/source.txt")
    derived.write_bytes(b'{"x":1}\n')
    raw.write_bytes(b'raw\n')
    assert locked_input_sha256(tmp_path, "data/processed/source.jsonl") == before
    assert locked_input_sha256(tmp_path, "data/raw/source.txt") != raw_before


def test_split_edit_does_not_survive_identifier_check(protocol, tmp_path):
    plan, _, _ = protocol
    bad = copy.deepcopy(plan)
    bad["config"]["split_seed"] += 1
    path = tmp_path / "proposal.json"
    write_json(path, bad)
    with pytest.raises(ValueError, match="identifier"):
        load_split_plan(ROOT, path, require_frozen=False)


def test_explicit_freeze_reproduces_metadata_and_never_overwrites_different_contents(protocol, tmp_path):
    plan, grouping, _ = protocol
    proposal = tmp_path / "proposal.json"
    write_json(proposal, plan)
    write_json(tmp_path / "grouping.json", grouping)
    frozen_path = freeze_protocol(ROOT, proposal)
    frozen = load_split_plan(ROOT, frozen_path)
    assert frozen["status"] == "frozen"
    assert frozen["proposal_id"] == plan["proposal_id"]
    assert freeze_protocol(ROOT, proposal) == frozen_path
    with pytest.raises(FileExistsError):
        write_json(frozen_path, {"different": True}, exclusive=True)


def test_freeze_rejects_changed_grouping_even_if_fold_file_was_untouched(protocol, tmp_path):
    plan, grouping, _ = protocol
    path = tmp_path / "proposal.json"
    write_json(path, plan)
    bad = copy.deepcopy(grouping)
    bad["limitations"] = []
    write_json(tmp_path / "grouping.json", bad)
    with pytest.raises(ValueError, match="Grouping audit differs"):
        freeze_protocol(ROOT, path)


def bound_candidate(tmp_path, plan, grouping):
    proposal_path = tmp_path / "proposal.json"
    write_json(proposal_path, plan)
    write_json(tmp_path / "grouping.json", grouping)
    frozen = dict(plan, status="frozen", frozen_proposal_sha256=file_sha256(proposal_path))
    path = tmp_path / "frozen.json"
    write_json(path, frozen)
    return path


def refresh_identifier(plan):
    plan.pop("proposal_id", None)
    plan["proposal_id"] = stable_digest(plan)[:20]


def test_one_word_status_edit_cannot_bypass_the_freeze_contract(protocol, tmp_path):
    plan, grouping, _ = protocol
    write_json(tmp_path / "proposal.json", plan)
    write_json(tmp_path / "grouping.json", grouping)
    path = tmp_path / "frozen.json"
    write_json(path, dict(plan, status="frozen"))
    with pytest.raises(ValueError, match="not bound"):
        load_split_plan(ROOT, path)


def test_rehashed_dissolved_groups_fail_even_with_matching_forged_audit_and_binding(protocol, tmp_path):
    original, original_grouping, _ = protocol
    plan, grouping = copy.deepcopy(original), copy.deepcopy(original_grouping)
    mapping = {compound: "dissolved:" + compound for compound in plan["compound_to_group"]}
    plan["compound_to_group"] = mapping
    plan["outer_folds"] = build_split_plan(sorted(mapping), mapping, plan["config"])
    grouping["compound_to_group"] = mapping
    grouping["groups"] = [{"group_id": group, "selected_members": [compound], "all_source_members": [compound]}
                          for compound, group in mapping.items()]
    plan["grouping_digest"] = stable_digest(grouping)
    refresh_identifier(plan)
    path = bound_candidate(tmp_path, plan, grouping)
    with pytest.raises(ValueError, match="cannot be reproduced"):
        load_split_plan(ROOT, path)


@pytest.mark.parametrize("omission", ["all_training_sets", "label_budgets", "one_inner_fold"])
def test_rehashed_incomplete_training_schedule_is_rejected(protocol, tmp_path, omission):
    original, grouping, _ = protocol
    plan = copy.deepcopy(original)
    for fold in plan["outer_folds"]:
        if omission == "all_training_sets":
            del fold["training_sets"]
        elif omission == "label_budgets":
            fold["training_sets"] = [t for t in fold["training_sets"] if t["budget"] == "full"]
        else:
            fold["training_sets"][0]["inner_folds"].pop()
    refresh_identifier(plan)
    path = bound_candidate(tmp_path, plan, grouping)
    with pytest.raises(ValueError, match="Configured"):
        load_split_plan(ROOT, path)


def test_wrong_frozen_proposal_hash_is_rejected(protocol, tmp_path):
    plan, grouping, _ = protocol
    path = bound_candidate(tmp_path, plan, grouping)
    frozen = json.loads(path.read_text())
    frozen["frozen_proposal_sha256"] = "0" * 64
    write_json(path, frozen)
    with pytest.raises(ValueError, match="not bound"):
        load_split_plan(ROOT, path)


def test_bound_proposal_byte_edits_are_rejected_even_if_json_value_is_unchanged(protocol, tmp_path):
    plan, grouping, _ = protocol
    path = bound_candidate(tmp_path, plan, grouping)
    proposal = tmp_path / "proposal.json"
    proposal.write_bytes(proposal.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="not bound"):
        load_split_plan(ROOT, path)


def test_loading_a_frozen_artifact_requires_its_grouping_audit(protocol, tmp_path):
    plan, grouping, _ = protocol
    path = bound_candidate(tmp_path, plan, grouping)
    (tmp_path / "grouping.json").unlink()
    with pytest.raises(ValueError, match="missing its grouping audit"):
        load_split_plan(ROOT, path)
