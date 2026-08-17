"""Safe subprocess bridge to the rare-oncology 2.1 CLI."""

from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping


_GENE_PATTERN = re.compile(r"^[A-Z][A-Z0-9-]{0,30}$")
_ALLOWED_SOURCES = {"alphamissense", "gnomad", "panelapp"}


class RunnerError(RuntimeError):
    """Raised when 2.1 cannot be launched or returns a failure."""


def normalize_gene(value: str) -> str:
    """Role: Normalize gene into a stable internal representation."""
    gene = value.strip().upper()
    if not _GENE_PATTERN.fullmatch(gene):
        raise ValueError(
            "Gene symbols must start with a letter and contain only A-Z, 0-9 or '-'."
        )
    return gene


@dataclass(frozen=True)
class RunOptions:
    include_therapies: bool = True
    include_oncokb: bool = False
    evidence_sources: tuple[str, ...] = (
        "alphamissense",
        "gnomad",
        "panelapp",
    )
    max_clinvar_records: int | None = None
    max_therapy_records: int = 100
    include_string_network: bool = True
    string_score: int = 700
    string_partners: int = 20
    string_network_type: str = "functional"
    force_refresh: bool = False
    timeout_seconds: int = 3600


@dataclass(frozen=True)
class StringRunOptions:
    required_score: int = 700
    added_nodes: int = 20
    network_type: str = "functional"
    force_refresh: bool = False
    timeout_seconds: int = 300


@dataclass(frozen=True)
class StructureRunOptions:
    variants: tuple[str, ...]
    model_seeds: tuple[int, ...] = (1,)
    protein_sequence: str | None = None
    offline: bool = False
    timeout_seconds: int = 300


@dataclass(frozen=True)
class GeneratedStructureFile:
    platform: str
    protein_change: str | None
    filename: str
    content: bytes
    mime_type: str


@dataclass(frozen=True)
class StructureRunOutcome:
    gene: str
    command: tuple[str, ...]
    files: tuple[GeneratedStructureFile, ...]
    skipped: tuple[Mapping[str, Any], ...]
    stdout: str
    master_database: Path
    output_directory: Path

    @property
    def command_display(self) -> str:
        """Role: Return a shell-safe command preview for the dashboard."""
        return shlex.join(self.command)


@dataclass(frozen=True)
class RunOutcome:
    gene: str
    command: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str
    master_database: Path
    output_directory: Path

    @property
    def command_display(self) -> str:
        """Role: Return a shell-safe command preview for the dashboard."""
        return shlex.join(self.command)


class AnalysisRunner:
    def __init__(
        self,
        project_root: str | Path,
        master_database: str | Path,
        *,
        environment: Mapping[str, str] | None = None,
    ) -> None:
        """Role: Compute init for the runner workflow."""
        self.project_root = Path(project_root).expanduser().resolve()
        self.v21_root = self.project_root / "2.1"
        self.master_database = Path(master_database).expanduser().resolve()
        self.environment = dict(os.environ if environment is None else environment)

    def build_command(self, gene: str, options: RunOptions) -> tuple[str, ...]:
        """Role: Build command data for downstream use."""
        normalized_gene = normalize_gene(gene)
        if options.max_clinvar_records is not None and options.max_clinvar_records < 1:
            raise ValueError("max_clinvar_records must be positive or None")
        if options.max_therapy_records < 1:
            raise ValueError("max_therapy_records must be positive")
        if options.timeout_seconds < 1:
            raise ValueError("timeout_seconds must be positive")
        _validate_string_options(
            options.string_score,
            options.string_partners,
            options.string_network_type,
        )
        unknown = set(options.evidence_sources) - _ALLOWED_SOURCES
        if unknown:
            raise ValueError(f"Unsupported evidence sources: {', '.join(sorted(unknown))}")

        command = [*self._cli_prefix(), normalized_gene]
        command.extend(
            [
                "--program",
                "ALL",
                "--classification",
                "ALL",
                "--database",
                str(self.master_database),
                "--output-dir",
                str(self._output_directory(normalized_gene)),
                "--cache-max-age-hours",
                "8760",
            ]
        )
        for source in dict.fromkeys(options.evidence_sources):
            command.extend(["--source", source])
        if "gnomad" in options.evidence_sources:
            command.extend(["--gnomad-max-variants", "25"])
        if options.max_clinvar_records is not None:
            command.extend(
                ["--max-clinvar-records", str(options.max_clinvar_records)]
            )
        if options.include_therapies:
            if options.include_oncokb:
                if not self.environment.get("ONCOKB_TOKEN"):
                    raise RunnerError(
                        "ONCOKB_TOKEN is not configured in the dashboard environment."
                    )
                for source in ("civic", "oncokb", "clinicaltrials", "opentargets"):
                    command.extend(["--therapy-source", source])
            else:
                command.append("--therapies")
            command.extend(
                ["--max-therapy-records", str(options.max_therapy_records)]
            )
        if options.include_string_network:
            command.extend(
                [
                    "--string-network",
                    "--string-score",
                    str(options.string_score),
                    "--string-partners",
                    str(options.string_partners),
                    "--string-network-type",
                    options.string_network_type,
                ]
            )
        if options.force_refresh:
            command.append("--refresh")
        return tuple(command)

    def build_string_command(
        self, gene: str, options: StringRunOptions
    ) -> tuple[str, ...]:
        """Role: Build string command data for downstream use."""
        normalized_gene = normalize_gene(gene)
        _validate_string_options(
            options.required_score,
            options.added_nodes,
            options.network_type,
        )
        if options.timeout_seconds < 1:
            raise ValueError("timeout_seconds must be positive")
        command = [
            str(self._v21_python()),
            "-m",
            "rare_oncology.string_cli",
            normalized_gene,
            "--database",
            str(self.master_database),
            "--output-dir",
            str(self._output_directory(normalized_gene) / "string"),
            "--string-score",
            str(options.required_score),
            "--string-partners",
            str(options.added_nodes),
            "--string-network-type",
            options.network_type,
        ]
        if options.force_refresh:
            command.append("--refresh")
        return tuple(command)

    def run(self, gene: str, options: RunOptions = RunOptions()) -> RunOutcome:
        """Role: Execute the configured gene analysis and return its stored snapshot."""
        normalized_gene = normalize_gene(gene)
        command = self.build_command(normalized_gene, options)
        output_directory = self._output_directory(normalized_gene)
        self.master_database.parent.mkdir(parents=True, exist_ok=True)
        output_directory.mkdir(parents=True, exist_ok=True)
        try:
            child_environment = dict(self.environment)
            source_root = str(self.v21_root / "src")
            existing_pythonpath = child_environment.get("PYTHONPATH")
            child_environment["PYTHONPATH"] = (
                source_root
                if not existing_pythonpath
                else source_root + os.pathsep + existing_pythonpath
            )
            completed = subprocess.run(
                command,
                cwd=self.v21_root,
                env=child_environment,
                capture_output=True,
                text=True,
                timeout=options.timeout_seconds,
                check=False,
            )
        except FileNotFoundError as exc:
            raise RunnerError(f"rare-oncology 2.1 executable not found: {exc}") from exc
        except subprocess.TimeoutExpired as exc:
            raise RunnerError(
                f"Analysis timed out after {options.timeout_seconds} seconds."
            ) from exc
        outcome = RunOutcome(
            gene=normalized_gene,
            command=command,
            returncode=completed.returncode,
            stdout=completed.stdout,
            stderr=completed.stderr,
            master_database=self.master_database,
            output_directory=output_directory,
        )
        if completed.returncode != 0:
            message = completed.stderr.strip() or completed.stdout.strip()
            raise RunnerError(
                f"rare-oncology 2.1 failed with exit code {completed.returncode}: {message}"
            )
        return outcome

    def run_string_network(
        self,
        gene: str,
        options: StringRunOptions = StringRunOptions(),
    ) -> RunOutcome:
        """Role: Execute a STRING query and return the stored network snapshot."""
        normalized_gene = normalize_gene(gene)
        command = self.build_string_command(normalized_gene, options)
        output_directory = self._output_directory(normalized_gene) / "string"
        self.master_database.parent.mkdir(parents=True, exist_ok=True)
        output_directory.mkdir(parents=True, exist_ok=True)
        child_environment = dict(self.environment)
        source_root = str(self.v21_root / "src")
        existing_pythonpath = child_environment.get("PYTHONPATH")
        child_environment["PYTHONPATH"] = (
            source_root
            if not existing_pythonpath
            else source_root + os.pathsep + existing_pythonpath
        )
        try:
            completed = subprocess.run(
                command,
                cwd=self.v21_root,
                env=child_environment,
                capture_output=True,
                text=True,
                timeout=options.timeout_seconds,
                check=False,
            )
        except FileNotFoundError as exc:
            raise RunnerError(f"Python de rare-oncology 2.1 no encontrado: {exc}") from exc
        except subprocess.TimeoutExpired as exc:
            raise RunnerError(
                f"STRING timed out after {options.timeout_seconds} seconds."
            ) from exc
        outcome = RunOutcome(
            gene=normalized_gene,
            command=command,
            returncode=completed.returncode,
            stdout=completed.stdout,
            stderr=completed.stderr,
            master_database=self.master_database,
            output_directory=output_directory,
        )
        if completed.returncode != 0:
            message = completed.stderr.strip() or completed.stdout.strip()
            raise RunnerError(
                f"STRING 2.1 failed with exit code {completed.returncode}: {message}"
            )
        return outcome

    def build_structure_command(
        self,
        gene: str,
        options: StructureRunOptions,
        *,
        report_path: Path,
        sequence_path: Path | None,
    ) -> tuple[str, ...]:
        """Role: Build structure command data for downstream use."""
        normalized_gene = normalize_gene(gene)
        if not options.variants:
            raise ValueError("Select at least one protein variant")
        if any(not value.strip() or len(value) > 80 for value in options.variants):
            raise ValueError("Invalid protein variant")
        if not options.model_seeds or any(
            isinstance(seed, bool) or seed < 0 for seed in options.model_seeds
        ):
            raise ValueError("model_seeds must contain non-negative integers")
        if options.timeout_seconds < 1:
            raise ValueError("timeout_seconds must be positive")
        output = self._output_directory(normalized_gene) / "structures"
        command = [
            str(self._v21_python()),
            "-m",
            "rare_oncology.structure_cli",
            "--report-json",
            str(report_path),
            "--database",
            str(self.master_database),
            "--output-dir",
            str(output),
            "--max-variants",
            str(len(options.variants)),
        ]
        for variant in dict.fromkeys(options.variants):
            command.extend(["--variant", variant])
        for seed in dict.fromkeys(options.model_seeds):
            command.extend(["--model-seed", str(seed)])
        if sequence_path is not None:
            command.extend(["--protein-sequence-file", str(sequence_path)])
        if options.offline:
            command.append("--offline")
        return tuple(command)

    def run_structure_inputs(
        self,
        gene: str,
        analysis_result: Mapping[str, Any],
        options: StructureRunOptions,
    ) -> StructureRunOutcome:
        """Role: Generate structure-model inputs and return their artifacts."""
        normalized_gene = normalize_gene(gene)
        output_directory = self._output_directory(normalized_gene) / "structures"
        self.master_database.parent.mkdir(parents=True, exist_ok=True)
        output_directory.mkdir(parents=True, exist_ok=True)
        child_environment = self._child_environment()
        with tempfile.TemporaryDirectory(prefix="rare-oncology-structure-") as temporary:
            temporary_root = Path(temporary)
            report_path = temporary_root / "report.json"
            report_path.write_text(
                json.dumps(
                    dict(analysis_result),
                    ensure_ascii=False,
                    allow_nan=False,
                ),
                encoding="utf-8",
            )
            sequence_path: Path | None = None
            if options.protein_sequence is not None:
                sequence_path = temporary_root / "sequence.fasta"
                sequence_path.write_text(options.protein_sequence, encoding="utf-8")
            command = self.build_structure_command(
                normalized_gene,
                options,
                report_path=report_path,
                sequence_path=sequence_path,
            )
            try:
                completed = subprocess.run(
                    command,
                    cwd=self.v21_root,
                    env=child_environment,
                    capture_output=True,
                    text=True,
                    timeout=options.timeout_seconds,
                    check=False,
                )
            except FileNotFoundError as exc:
                raise RunnerError(
                    f"Python de rare-oncology 2.1 no encontrado: {exc}"
                ) from exc
            except subprocess.TimeoutExpired as exc:
                raise RunnerError(
                    f"Structure input generation timed out after "
                    f"{options.timeout_seconds} seconds."
                ) from exc
        if completed.returncode != 0:
            message = completed.stderr.strip() or completed.stdout.strip()
            raise RunnerError(
                f"Structure input generation failed with exit code "
                f"{completed.returncode}: {message}"
            )
        manifest_path = output_directory / f"{normalized_gene.lower()}_structure_inputs.json"
        manifest = _read_manifest(manifest_path)
        files: list[GeneratedStructureFile] = []
        for artifact in manifest.get("artifacts", []):
            if not isinstance(artifact, dict):
                continue
            artifact_path = Path(str(artifact.get("path", ""))).resolve()
            if (
                not artifact_path.is_file()
                or output_directory.resolve() not in artifact_path.parents
            ):
                raise RunnerError("Structure manifest references an invalid artifact path")
            platform = str(artifact.get("platform") or "unknown")
            files.append(
                GeneratedStructureFile(
                    platform=platform,
                    protein_change=(
                        str(artifact["protein_change"])
                        if artifact.get("protein_change")
                        else None
                    ),
                    filename=artifact_path.name,
                    content=artifact_path.read_bytes(),
                    mime_type=(
                        "application/json"
                        if artifact_path.suffix == ".json"
                        else "application/yaml"
                    ),
                )
            )
        files.append(
            GeneratedStructureFile(
                platform="manifest",
                protein_change=None,
                filename=manifest_path.name,
                content=_download_manifest(manifest, manifest_path.name),
                mime_type="application/json",
            )
        )
        skipped = tuple(
            item
            for item in manifest.get("skipped", [])
            if isinstance(item, dict)
        )
        return StructureRunOutcome(
            gene=normalized_gene,
            command=command,
            files=tuple(files),
            skipped=skipped,
            stdout=completed.stdout,
            master_database=self.master_database,
            output_directory=output_directory,
        )

    def _cli_prefix(self) -> tuple[str, ...]:
        """Role: Compute cli prefix for the runner workflow."""
        configured = self.environment.get("RARE_ONCOLOGY_CLI")
        if configured:
            return (str(Path(configured).expanduser()),)
        executable = self.v21_root / ".venv" / "bin" / "rare-oncology"
        if executable.is_file():
            return (str(executable),)
        package = self.v21_root / "src" / "rare_oncology"
        if package.is_dir():
            return (sys.executable, "-m", "rare_oncology")
        raise RunnerError(f"Could not locate rare-oncology 2.1 under {self.v21_root}")

    def _v21_python(self) -> Path:
        """Role: Compute v21 python for the runner workflow."""
        executable = self.v21_root / ".venv" / "bin" / "python"
        if executable.is_file():
            return executable
        if (self.v21_root / "src" / "rare_oncology").is_dir():
            return Path(sys.executable)
        raise RunnerError(f"Could not locate rare-oncology 2.1 under {self.v21_root}")

    def _child_environment(self) -> dict[str, str]:
        """Role: Compute child environment for the runner workflow."""
        child_environment = dict(self.environment)
        source_root = str(self.v21_root / "src")
        existing_pythonpath = child_environment.get("PYTHONPATH")
        child_environment["PYTHONPATH"] = (
            source_root
            if not existing_pythonpath
            else source_root + os.pathsep + existing_pythonpath
        )
        return child_environment

    def _output_directory(self, gene: str) -> Path:
        """Role: Compute output directory for the runner workflow."""
        return self.v21_root / "results" / "dashboard_runs" / gene


def _validate_string_options(score: int, partners: int, network_type: str) -> None:
    """Role: Validate string options constraints before execution."""
    if not 0 <= score <= 1000:
        raise ValueError("string_score must be between 0 and 1000")
    if not 0 <= partners <= 100:
        raise ValueError("string_partners must be between 0 and 100")
    if network_type not in {"functional", "physical"}:
        raise ValueError("string_network_type must be functional or physical")


def _read_manifest(path: Path) -> Mapping[str, Any]:
    """Role: Compute read manifest for the runner workflow."""
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RunnerError(f"Could not read structure manifest: {exc}") from exc
    if not isinstance(payload, dict):
        raise RunnerError("Structure manifest is not a JSON object")
    return payload


def _download_manifest(manifest: Mapping[str, Any], filename: str) -> bytes:
    """Remove host-local paths from the browser-facing manifest."""

    payload = dict(manifest)
    artifacts: list[dict[str, Any]] = []
    for value in manifest.get("artifacts", []):
        if not isinstance(value, dict):
            continue
        artifact = dict(value)
        if artifact.get("path"):
            artifact["path"] = Path(str(artifact["path"])).name
        artifacts.append(artifact)
    payload["artifacts"] = artifacts
    payload["manifest_path"] = filename
    return (
        json.dumps(
            payload,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
