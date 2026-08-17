#!/usr/bin/env python3
"""Stream CAPTURE-24 into wake-aligned 5-second ENMO functions.

The official archive stores one gzip-compressed CSV per participant inside a
ZIP archive.  Reading all raw 100 Hz rows into memory at once is unnecessary.
This helper processes one participant at a time in chunks, averages ENMO over
non-overlapping epochs, identifies the main sleep episode from the released
annotations, and circularly aligns the 24-hour curve at wake time.

The R pipeline owns all statistical analysis.  This helper only performs the
large, deterministic file-format conversion and writes simple binary matrices
that base R can read without an additional package.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import shutil
import sys
import time
import zipfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED_PARTICIPANTS = 151
EXPECTED_SAMPLE_RATE = 100
SECONDS_PER_DAY = 24 * 60 * 60


def archive_members(archive: Path) -> tuple[list[str], str, str]:
    """Return participant, metadata, and annotation-dictionary member names."""
    with zipfile.ZipFile(archive) as zf:
        names = zf.namelist()
    participant = sorted(
        name for name in names if re.search(r"(^|/)P\d{3}\.csv\.gz$", name)
    )
    metadata = next(
        (name for name in names if Path(name).name.lower() == "metadata.csv"), None
    )
    dictionary = next(
        (
            name
            for name in names
            if Path(name).name.lower() == "annotation-label-dictionary.csv"
        ),
        None,
    )
    if len(participant) != EXPECTED_PARTICIPANTS:
        raise RuntimeError(
            f"Expected {EXPECTED_PARTICIPANTS} participant files, found {len(participant)}."
        )
    if metadata is None or dictionary is None:
        raise RuntimeError("The archive is missing metadata.csv or the annotation dictionary.")
    return participant, metadata, dictionary


def read_small_csv(archive: Path, member: str) -> pd.DataFrame:
    with zipfile.ZipFile(archive) as zf, zf.open(member) as stream:
        return pd.read_csv(stream, dtype="string")


def sleep_annotation_keys(dictionary: pd.DataFrame) -> set[str]:
    """Find raw annotation values mapped to sleep by any released label system."""
    text = dictionary.fillna("").astype(str)
    sleep_rows = np.zeros(len(text), dtype=bool)
    for column in text.columns:
        sleep_rows |= text[column].str.contains(r"\bsleep", case=False, regex=True).to_numpy()
    if not sleep_rows.any():
        raise RuntimeError("No sleep category was found in the annotation dictionary.")
    annotation_column = next(
        (column for column in text.columns if column.lower() == "annotation"),
        text.columns[0],
    )
    return set(text.loc[sleep_rows, annotation_column].dropna().astype(str))


def close_short_false_gaps(values: np.ndarray, maximum_gap: int) -> np.ndarray:
    """Join sleep segments separated by at most maximum_gap non-sleep epochs."""
    result = values.astype(bool, copy=True)
    true_index = np.flatnonzero(result)
    if true_index.size < 2 or maximum_gap <= 0:
        return result
    gaps = np.diff(true_index) - 1
    for left, gap in zip(true_index[:-1], gaps):
        if 0 < gap <= maximum_gap:
            result[left + 1 : left + gap + 1] = True
    return result


def longest_true_run(values: np.ndarray) -> tuple[int, int]:
    """Return inclusive start and exclusive end indices of the longest True run."""
    padded = np.concatenate(([False], values.astype(bool), [False])).astype(np.int8)
    changes = np.diff(padded)
    starts = np.flatnonzero(changes == 1)
    ends = np.flatnonzero(changes == -1)
    if starts.size == 0:
        raise RuntimeError("No sleep episode was detected.")
    index = int(np.argmax(ends - starts))
    return int(starts[index]), int(ends[index])


def parse_epoch_times(values: list[str]) -> pd.DatetimeIndex:
    parsed = pd.to_datetime(pd.Series(values), errors="coerce", utc=True, format="mixed")
    if parsed.isna().any():
        parsed = pd.to_datetime(pd.Series(values), errors="coerce", utc=True)
    if parsed.isna().any():
        raise RuntimeError(f"Could not parse {int(parsed.isna().sum())} epoch timestamps.")
    return pd.DatetimeIndex(parsed)


def process_participant(
    archive_path: str,
    member: str,
    epoch_seconds: int,
    sleep_gap_minutes: float,
    chunk_rows: int,
    sleep_keys: set[str],
) -> dict:
    """Process one participant and return its aligned curve plus QC fields."""
    started = time.time()
    rows_per_epoch = EXPECTED_SAMPLE_RATE * epoch_seconds
    if chunk_rows % rows_per_epoch:
        chunk_rows = int(math.ceil(chunk_rows / rows_per_epoch) * rows_per_epoch)

    epoch_sum: list[np.ndarray] = []
    epoch_count: list[np.ndarray] = []
    sleep_count: list[np.ndarray] = []
    annotation_count: list[np.ndarray] = []
    epoch_times: list[str] = []
    raw_rows = 0

    with zipfile.ZipFile(archive_path) as zf, zf.open(member) as compressed:
        with gzip.GzipFile(fileobj=compressed) as stream:
            reader = pd.read_csv(
                stream,
                usecols=["time", "x", "y", "z", "annotation"],
                dtype={"x": "float32", "y": "float32", "z": "float32", "annotation": "string"},
                chunksize=chunk_rows,
            )
            carry = None
            for chunk in reader:
                if carry is not None:
                    chunk = pd.concat([carry, chunk], ignore_index=True)
                    carry = None
                complete_rows = (len(chunk) // rows_per_epoch) * rows_per_epoch
                if complete_rows < len(chunk):
                    carry = chunk.iloc[complete_rows:].copy()
                    chunk = chunk.iloc[:complete_rows]
                if chunk.empty:
                    continue

                groups = len(chunk) // rows_per_epoch
                xyz = chunk[["x", "y", "z"]].to_numpy(dtype=np.float32, copy=False)
                finite = np.isfinite(xyz).all(axis=1)
                norm = np.sqrt(np.square(xyz, dtype=np.float32).sum(axis=1))
                enmo_mg = np.maximum(norm - 1.0, 0.0) * 1000.0
                enmo_mg[~finite] = 0.0
                enmo_mg = enmo_mg.reshape(groups, rows_per_epoch)
                finite = finite.reshape(groups, rows_per_epoch)
                epoch_sum.append(enmo_mg.sum(axis=1, dtype=np.float64))
                epoch_count.append(finite.sum(axis=1, dtype=np.int32))

                annotations = chunk["annotation"]
                annotated = annotations.notna().to_numpy().reshape(groups, rows_per_epoch)
                annotation_text = annotations.fillna("").astype(str)
                is_sleep = annotation_text.isin(sleep_keys) | annotation_text.str.contains(
                    r"\bsleep", case=False, regex=True
                )
                sleep_count.append(is_sleep.to_numpy().reshape(groups, rows_per_epoch).sum(axis=1))
                annotation_count.append(annotated.sum(axis=1))
                epoch_times.extend(chunk["time"].iloc[::rows_per_epoch].astype(str).tolist())
                raw_rows += len(chunk)

            if carry is not None and len(carry):
                raw_rows += len(carry)

    sums = np.concatenate(epoch_sum)
    counts = np.concatenate(epoch_count)
    sleep_counts = np.concatenate(sleep_count)
    annotated_counts = np.concatenate(annotation_count)
    times = parse_epoch_times(epoch_times)
    valid_epoch = counts >= int(0.8 * rows_per_epoch)
    mean_enmo = np.divide(sums, counts, out=np.full_like(sums, np.nan), where=counts > 0)
    sleep_fraction = np.divide(
        sleep_counts,
        annotated_counts,
        out=np.zeros_like(sleep_counts, dtype=np.float64),
        where=annotated_counts > 0,
    )

    epoch_gaps = np.diff(times.asi8.astype(np.float64)) / 1e9
    sleep_epoch = sleep_fraction >= 0.5
    closed_sleep = close_short_false_gaps(
        sleep_epoch, int(round(sleep_gap_minutes * 60 / epoch_seconds))
    )
    # A true acquisition gap cannot be treated as a continuous sleep episode,
    # even if sleep labels happen to occur on both observed sides.
    acquisition_breaks = np.flatnonzero(epoch_gaps > sleep_gap_minutes * 60)
    closed_sleep[acquisition_breaks + 1] = False
    sleep_start, sleep_end = longest_true_run(closed_sleep)
    target_bins = SECONDS_PER_DAY // epoch_seconds
    wake_position = sleep_end
    wake_time = (
        times[wake_position]
        if wake_position < len(times)
        else times[-1] + pd.Timedelta(seconds=epoch_seconds)
    )
    # Select at most one actual 24-hour timestamp interval containing the main
    # sleep episode. Timestamp mapping preserves acquisition gaps instead of
    # collapsing them through row indices. A shorter record is never tiled or
    # repeated to manufacture a complete day.
    source_start = times[0]
    source_end = times[-1] + pd.Timedelta(seconds=epoch_seconds)
    day = pd.Timedelta(seconds=SECONDS_PER_DAY)
    if source_end - source_start >= day:
        latest_start = source_end - day
        desired_start = times[sleep_start]
        selected_start = max(source_start, min(desired_start, latest_start))
        selected_end = selected_start + day
        source_index = np.flatnonzero((times >= selected_start) & (times < selected_end))
    else:
        selected_start = source_start
        selected_end = source_end
        source_index = np.arange(len(times))
    relative_seconds = (
        times[source_index].asi8.astype(np.float64) - float(wake_time.value)
    ) / 1e9
    target = np.mod(np.rint(relative_seconds / epoch_seconds).astype(np.int64), target_bins)
    selected_valid = valid_epoch[source_index]
    selected_enmo = mean_enmo[source_index]
    curve_sum = np.bincount(
        target[selected_valid], weights=selected_enmo[selected_valid], minlength=target_bins
    )
    curve_count = np.bincount(target[selected_valid], minlength=target_bins)
    curve = np.divide(
        curve_sum,
        curve_count,
        out=np.full(target_bins, np.nan, dtype=np.float64),
        where=curve_count > 0,
    )
    curve = np.log1p(curve).astype(np.float32)
    coverage = (curve_count > 0).astype(np.uint8)

    elapsed = time.time() - started
    participant_id = Path(member).name.split(".")[0]
    irregular_epoch_gaps = int(np.sum(np.abs(epoch_gaps - epoch_seconds) > 0.25))
    source_duration_hours = (
        (times[-1] - times[0]).total_seconds() + epoch_seconds
    ) / 3600.0
    inferred_rate = raw_rows / max(source_duration_hours * 3600.0, 1.0)
    return {
        "participant_id": participant_id,
        "member": member,
        "curve": curve,
        "coverage": coverage,
        "raw_rows": raw_rows,
        "source_epochs": len(times),
        "source_duration_hours": source_duration_hours,
        "inferred_sample_rate_hz": inferred_rate,
        "maximum_epoch_gap_seconds": float(np.max(epoch_gaps)) if epoch_gaps.size else float("nan"),
        "irregular_epoch_gaps": irregular_epoch_gaps,
        "wake_time": wake_time.isoformat(),
        "main_sleep_start": times[sleep_start].isoformat(),
        "main_sleep_end": wake_time.isoformat(),
        "main_sleep_hours": (sleep_end - sleep_start) * epoch_seconds / 3600.0,
        "selected_window_start": selected_start.isoformat(),
        "selected_window_end": selected_end.isoformat(),
        "annotated_epoch_percent": 100.0 * np.mean(annotated_counts > 0),
        "observed_24h_bins": int(coverage.sum()),
        "missing_24h_minutes": float((target_bins - coverage.sum()) * epoch_seconds / 60.0),
        "mean_enmo_mg_observed": float(np.nanmean(np.expm1(curve))),
        "processing_seconds": elapsed,
    }


def write_outputs(results: list[dict], output_dir: Path, epoch_seconds: int, archive: Path) -> None:
    results.sort(key=lambda row: row["participant_id"])
    curves = np.vstack([row.pop("curve") for row in results]).astype("<f4", copy=False)
    coverage = np.vstack([row.pop("coverage") for row in results]).astype("u1", copy=False)
    curve_path = output_dir / "wake_aligned_log1p_enmo_mg_float32.bin"
    coverage_path = output_dir / "wake_aligned_coverage_uint8.bin"
    curves.tofile(curve_path)
    coverage.tofile(coverage_path)
    pd.DataFrame(results).to_csv(output_dir / "subject_qc.csv", index=False, quoting=csv.QUOTE_MINIMAL)
    manifest = {
        "dataset": "CAPTURE-24",
        "source_archive": str(archive.resolve()),
        "participants": int(curves.shape[0]),
        "grid_size": int(curves.shape[1]),
        "epoch_seconds": int(epoch_seconds),
        "curve_file": curve_path.name,
        "curve_dtype": "little-endian float32, C row-major",
        "coverage_file": coverage_path.name,
        "coverage_dtype": "uint8, C row-major",
        "transformation": "log1p(max(sqrt(x^2+y^2+z^2)-1,0)*1000) after 5-second averaging",
        "alignment": "circular 24-hour rotation to the end of the longest annotation-defined main sleep episode",
    }
    with open(output_dir / "preprocessing_manifest.json", "w", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--epoch-seconds", type=int, default=5)
    parser.add_argument("--sleep-gap-minutes", type=float, default=20.0)
    parser.add_argument("--chunk-rows", type=int, default=500_000)
    parser.add_argument("--n-jobs", type=int, default=4)
    args = parser.parse_args()

    archive = Path(args.archive)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    participant, metadata_member, dictionary_member = archive_members(archive)
    dictionary = read_small_csv(archive, dictionary_member)
    sleep_keys = sleep_annotation_keys(dictionary)

    with zipfile.ZipFile(archive) as zf:
        with zf.open(metadata_member) as source, open(output_dir / "metadata.csv", "wb") as target:
            shutil.copyfileobj(source, target)
        with zf.open(dictionary_member) as source, open(
            output_dir / "annotation-label-dictionary.csv", "wb"
        ) as target:
            shutil.copyfileobj(source, target)

    print(
        f"CAPTURE-24 preprocessing: {len(participant)} participants, "
        f"epoch={args.epoch_seconds}s, workers={args.n_jobs}",
        flush=True,
    )
    results: list[dict] = []
    with ProcessPoolExecutor(max_workers=max(1, args.n_jobs)) as executor:
        future_map = {
            executor.submit(
                process_participant,
                str(archive),
                member,
                args.epoch_seconds,
                args.sleep_gap_minutes,
                args.chunk_rows,
                sleep_keys,
            ): member
            for member in participant
        }
        for completed, future in enumerate(as_completed(future_map), start=1):
            member = future_map[future]
            try:
                result = future.result()
            except Exception as error:
                for pending in future_map:
                    pending.cancel()
                raise RuntimeError(f"Failed while processing {member}: {error}") from error
            results.append(result)
            print(
                f"CAPTURE-24 preprocess [{completed}/{len(participant)}] "
                f"{result['participant_id']} | source={result['source_duration_hours']:.2f}h | "
                f"missing={result['missing_24h_minutes']:.1f}min | "
                f"elapsed={result['processing_seconds']:.1f}s",
                flush=True,
            )

    write_outputs(results, output_dir, args.epoch_seconds, archive)
    print("CAPTURE-24 preprocessing complete.", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        raise
