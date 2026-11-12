import copy
import json
import subprocess
import sys
from pathlib import Path

import pytest

from regjepa.protocol import freeze_protocol, load_split_plan, locked_input_sha256, prepare_protocol, verify_inputs, write_json
from regjepa.splits import validate_split_plan

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
    validate_split_plan(plan["outer_folds"], sorted(groups), groups)


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
