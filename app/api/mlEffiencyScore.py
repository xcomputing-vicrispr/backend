import json
import os
import subprocess
import sys
import warnings

from rs3.seq import predict_seq

warnings.filterwarnings("ignore")

INVALID_SCORE = 999


def _valid_context(seq):
    seq = str(seq or "").upper()
    return len(seq) == 30 and set(seq) <= set("ACGT") and seq[25:27] == "GG"


def _metadata_value(metadata, *keys):
    if not isinstance(metadata, dict):
        return None
    for key in keys:
        value = metadata.get(key)
        if value not in (None, ""):
            return value
    return None


def _format_invalid_context(index, seq, pam, metadata):
    parts = [f"index={index}"]
    task_id = _metadata_value(metadata, "task_id", "idfile")
    stt = _metadata_value(metadata, "stt")
    pos = _metadata_value(metadata, "pos", "location")
    strand = _metadata_value(metadata, "strand")
    sgrna = _metadata_value(metadata, "sgRNA", "sequence")

    if task_id is not None:
        parts.append(f"task={task_id}")
    if stt is not None:
        parts.append(f"stt={stt}")
    if pos is not None:
        parts.append(f"pos={pos}")
    if strand is not None:
        parts.append(f"strand={strand}")
    if sgrna is not None:
        parts.append(f"sgRNA={sgrna}")

    parts.extend([f"seq[25:27]={pam!r}", f"mlseq={seq}"])
    return ", ".join(parts)


def _assert_ngg_contexts(seqlist, metadata_list=None):
    invalid = []
    for index, raw_seq in enumerate(seqlist):
        seq = str(raw_seq or "").upper()
        if len(seq) == 30 and set(seq) <= set("ACGT") and seq[25:27] != "GG":
            metadata = None
            if metadata_list and index < len(metadata_list):
                metadata = metadata_list[index]
            invalid.append((index, seq, seq[25:27], metadata))

    if invalid:
        examples = "; ".join(
            _format_invalid_context(index, seq, pam, metadata)
            for index, seq, pam, metadata in invalid[:5]
        )
        remaining = ""
        if len(invalid) > 5:
            remaining = f"; and {len(invalid) - 5} more invalid contexts"
        message = (
            "Critical Rule Set 2/3 input error: 30mer context without NGG "
            f"at seq[25:27]. {examples}{remaining}"
        )
        print(message, file=sys.stderr)
        raise RuntimeError(message)


def _merge_scores(seqlist, valid_scores):
    merged = []
    score_iter = iter(valid_scores)
    for seq in seqlist:
        if _valid_context(seq):
            merged.append(float(next(score_iter)))
        else:
            merged.append(INVALID_SCORE)
    return merged


def get_ml_score(seqlist, metadata_list=None):
    _assert_ngg_contexts(seqlist, metadata_list)
    valid_seqs = [str(seq).upper() for seq in seqlist if _valid_context(seq)]
    if not valid_seqs:
        return [INVALID_SCORE] * len(seqlist)

    seq_str = json.dumps(valid_seqs)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    test_file = os.path.abspath(os.path.join(current_dir, "worker", "getrs2.py"))
    test_dir = os.path.dirname(test_file)

    result = subprocess.Popen(
        [sys.executable, test_file, seq_str],
        cwd=test_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    stdout, stderr = result.communicate()
    stdout_text = stdout.decode(errors="replace")
    stderr_text = stderr.decode(errors="replace")

    print("=== ML Script STDOUT ===")
    print(valid_seqs)
    print(stdout_text)
    print("=== ML Script STDERR ===")
    print(stderr_text)

    if result.returncode != 0:
        raise RuntimeError(
            f"Rule Set 2 scorer failed with code {result.returncode}: {stderr_text.strip()}"
        )

    for line in reversed(stdout_text.splitlines()):
        line = line.strip()
        if line.startswith("[") and line.endswith("]"):
            values = json.loads(line)
            if len(values) != len(valid_seqs):
                raise RuntimeError(
                    f"Rule Set 2 returned {len(values)} scores for {len(valid_seqs)} valid sgRNAs"
                )
            return _merge_scores(seqlist, values)

    raise RuntimeError(
        f"Rule Set 2 scorer did not return a JSON array. STDOUT: {stdout_text.strip()} STDERR: {stderr_text.strip()}"
    )


def get_ml_score_azi3(seqlist, metadata_list=None):
    _assert_ngg_contexts(seqlist, metadata_list)
    valid_seqs = [str(seq).upper() for seq in seqlist if _valid_context(seq)]
    if not valid_seqs:
        return [INVALID_SCORE] * len(seqlist)

    res = predict_seq(valid_seqs, sequence_tracr="Hsu2013")
    values = [float(x) for x in res.tolist()]
    if len(values) != len(valid_seqs):
        raise RuntimeError(
            f"Rule Set 3 returned {len(values)} scores for {len(valid_seqs)} valid sgRNAs"
        )
    return _merge_scores(seqlist, values)
